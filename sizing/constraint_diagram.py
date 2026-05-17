"""Build an  constraint diagram with AeroSandbox.

This connects:
- the volume model, which determines `S_plan` from Kuechemann slenderness,
- the Raymer-style weight breakdown, which determines OEW/TOGW,
- the Zhang et al. constraint-analysis equation,
- and pyCycle engine-deck thrust availability.
"""

from dataclasses import dataclass, replace
from pathlib import Path
import csv
import sys

import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.tools.units as u

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from aero import (
    drag_build_up_coefficients,
    drag_geometry_from_planform_area,
    engine_deck_drag_point,
    leading_edge_sonic_mach,
    swept_wing_oswald_efficiency,
    swept_wing_tip_le_x_m,
)

try:
    from .aircraft import (
        Aircraft,
        DEFAULT_AIRCRAFT_JSON,
        FuelSystem,
        Payload,
        PropulsionSystem,
        build_geometric_asb_airplane as build_airplane,
        fuel_from_definition,
        load_aircraft_definition,
        interiors_from_definition,
        payload_from_definition,
        propulsion_from_definition,
        systems_from_definition,
    )
    from .constraint_equations import design_point_thrust_to_weight_from_wing_loading
    from .volume import aircraft_volume_breakdown
    from .weight import _weight_breakdown
except ImportError:
    from aircraft import (
        Aircraft,
        DEFAULT_AIRCRAFT_JSON,
        FuelSystem,
        Payload,
        PropulsionSystem,
        build_geometric_asb_airplane as build_airplane,
        fuel_from_definition,
        load_aircraft_definition,
        interiors_from_definition,
        payload_from_definition,
        propulsion_from_definition,
        systems_from_definition,
    )
    from constraint_equations import design_point_thrust_to_weight_from_wing_loading
    from volume import aircraft_volume_breakdown
    from weight import _weight_breakdown


plots_dir = Path("outputs/plots")


@dataclass(frozen=True)
class ConstraintDesignPoint:
    name: str
    mode: str
    mach: object
    altitude_m: object = None
    case: str = "generic"
    load_factor: object = 1.0
    beta: object = 1.0
    alpha: object = 1.0
    climb_rate_m_s: object = 0.0
    acceleration_m_s2: object = 0.0
    initial_mach: object = None
    final_mach: object = None
    initial_altitude_m: object = None
    final_altitude_m: object = None
    acceleration_time_s: object = 30.0
    drag_polar_k1: object = 0.05
    drag_polar_k2: object = 0.0
    cd0: object = 0.025
    cl_max: object = 1.8
    lift_coefficient: object = 0.5
    ground_roll_m: object = 900.0
    braking_roll_m: object = 800.0
    friction_coefficient: object = 0.03
    speed_ratio: object = 1.2
    climb_angle_deg: object = 3.0
    include_in_governing: bool = True


@dataclass(frozen=True)
class ConstraintDiagramConfig:
    aircraft_json: object = DEFAULT_AIRCRAFT_JSON
    engine_deck_csv: object = "propulsion/data/example_engine_deck.csv"
    kuechemann_tau: object = 0.0446
    void_volume_coefficient: object = 0.05
    thrust_scale: object = 1.0
    wing_form_factor: object = 1.15
    tail_form_factor: object = 1.15
    nacelle_form_factor: object = 1.30
    nacelle_length_to_diameter: object = 2.5
    main_wing_quarter_chord_sweep_rad: object = np.radians(60.0)
    airfoil_trailing_edge_angle_deg: object = 10.0
    design_cruise_mach: object = 3.0
    design_cruise_altitude_m: object = 70000.0 * u.foot
    drag_divergence_mach: object = 0.95
    drag_divergence_wave_cd: object = 0.002
    critical_mach_offset_from_mdd: object = 0.08
    supersonic_wave_drag_start_mach: object = 1.20
    drag_plot_cl_max: object = 1.8
    wing_loading_min_N_m2: object = 1000.0
    wing_loading_max_N_m2: object = 20000.0
    wing_loading_points: int = 250
    save_plot: object = plots_dir / "_constraint_diagram.png"
    show_plot: bool = False


def default_constraint_design_points():
    return (
        ConstraintDesignPoint(
            name="Case 1: constant-altitude/speed cruise",
            mode="ramjet",
            mach=3.0,
            altitude_m=60000.0 * u.foot,
            case="generic",
            cd0=0.032,
        ),
        ConstraintDesignPoint(
            name="Case 2: constant-speed climb",
            mode="fan_ab",
            mach=2.0,
            altitude_m=45000.0 * u.foot,
            case="generic",
            climb_rate_m_s=20.0,
            cd0=0.034,
        ),
        ConstraintDesignPoint(
            name="Case 3: constant-altitude/speed turn",
            mode="fan",
            mach=0.9,
            altitude_m=15000.0 * u.foot,
            case="generic",
            load_factor=2.5,
            cd0=0.028,
        ),
        ConstraintDesignPoint(
            name="Case 4a: transonic climb acceleration",
            mode="fan_ab",
            mach=1.2,
            case="climb_acceleration",
            initial_mach=0.7,
            final_mach=1.2,
            acceleration_time_s=30 * 60.0,
            initial_altitude_m=18000 * u.foot,
            final_altitude_m=25000 * u.foot,
            cd0=0.034,
        ),
        ConstraintDesignPoint(
            name="Case 4b: supersonic climb acceleration",
            mode="fan_ab",
            mach=2.2,
            case="climb_acceleration",
            initial_mach=1.7,
            final_mach=2.2,
            acceleration_time_s=30 * 60.0,
            initial_altitude_m=35000 * u.foot,
            final_altitude_m=45000 * u.foot,
            cd0=0.034,
        ),
        ConstraintDesignPoint(
            name="Case 4c: high-Mach climb acceleration",
            mode="ramjet",
            mach=3.1,
            case="climb_acceleration",
            initial_mach=2.2,
            final_mach=3.1,
            acceleration_time_s=30 * 60.0,
            initial_altitude_m=55000 * u.foot,
            final_altitude_m=70000 * u.foot,
            cd0=0.034,
        ),
        ConstraintDesignPoint(
            name="Case 5: takeoff ground roll ideal",
            mode="fan",
            mach=0.25,
            altitude_m=0.0,
            case="takeoff_ground_roll_ideal",
            cl_max=1.8,
            ground_roll_m=900.0,
            speed_ratio=1.2,
            cd0=0.040,
        ),
        ConstraintDesignPoint(
            name="Case 6: takeoff ground roll",
            mode="fan",
            mach=0.25,
            altitude_m=0.0,
            case="takeoff_ground_roll",
            cl_max=1.8,
            ground_roll_m=900.0,
            friction_coefficient=0.03,
            speed_ratio=1.2,
            cd0=0.040,
        ),
        ConstraintDesignPoint(
            name="Case 7: braking roll",
            mode="fan",
            mach=0.22,
            altitude_m=0.0,
            case="braking_roll",
            alpha=-1.0,
            cl_max=2.0,
            braking_roll_m=800.0,
            friction_coefficient=0.35,
            speed_ratio=1.3,
            cd0=0.050,
        ),
        ConstraintDesignPoint(
            name="Case 8: service ceiling",
            mode="ramjet",
            mach=3.0,
            altitude_m=80000.0 * u.foot,
            case="service_ceiling",
            lift_coefficient=0.5,
            climb_rate_m_s=0.508,
            cd0=0.045,
        ),
        ConstraintDesignPoint(
            name="Case 9: takeoff climb angle",
            mode="fan",
            mach=0.25,
            altitude_m=0.0,
            case="takeoff_climb_angle",
            cl_max=1.8,
            speed_ratio=1.2,
            climb_angle_deg=1.0,
            cd0=0.040,
        ),
    )


def config_propulsion(config, **overrides):
    return propulsion_from_definition(config.aircraft_json, **overrides)


def config_fuel(config, **overrides):
    return fuel_from_definition(config.aircraft_json, **overrides)


def config_payload(config, **overrides):
    return payload_from_definition(config.aircraft_json, **overrides)


def config_systems(config, **overrides):
    return systems_from_definition(config.aircraft_json, **overrides)


def config_interiors(config, **overrides):
    return interiors_from_definition(config.aircraft_json, **overrides)


def read_pycycle_engine_deck(engine_deck_csv):
    rows = []
    with Path(engine_deck_csv).open(newline="") as f:
        for row in csv.DictReader(f):
            if "altitude_ft" in row:
                altitude_m = float(row["altitude_ft"]) * 0.3048
                thrust_N = float(row["thrust_lbf"]) * 4.4482216152605
                fuel_flow_kg_s = float(row["fuel_flow_lbm_s"]) * 0.45359237
                electric_power_W = float(row.get("electric_power_hp") or 0.0) * 745.6998715822702
            else:
                altitude_m = float(row["altitude_m"])
                thrust_N = float(row["thrust_N"])
                fuel_flow_kg_s = float(row["fuel_flow_kg_s"])
                electric_power_W = float(row.get("electric_power_W") or 0.0)
            rows.append(
                {
                    "mode": row["mode"],
                    "mach": float(row["mach"]),
                    "altitude_m": altitude_m,
                    "throttle": float(row["throttle"]),
                    "thrust_N": thrust_N,
                    "fuel_flow_kg_s": fuel_flow_kg_s,
                    "electric_power_W": electric_power_W,
                }
            )
    if not rows:
        raise ValueError(f"No pyCycle engine-deck rows found in {engine_deck_csv}.")
    return rows


def pycycle_design_point_thrust_N(engine_deck_rows, design_point, thrust_scale=1.0):
    candidates = [row for row in engine_deck_rows if row["mode"] == design_point.mode]
    if not candidates:
        raise ValueError(f'Mode "{design_point.mode}" is missing from the pyCycle engine deck.')

    row = min(
        candidates,
        key=lambda item: (
            (item["mach"] - design_point.mach) ** 2
            + (
                (item["altitude_m"] - representative_altitude_m(design_point))
                / 10000.0
            ) ** 2
        ),
    )
    throttle = max(row["throttle"], 1e-9)
    return thrust_scale * row["thrust_N"] / throttle


def representative_altitude_m(design_point):
    if design_point.altitude_m is not None:
        return design_point.altitude_m
    if (
        design_point.initial_altitude_m is not None
        and design_point.final_altitude_m is not None
    ):
        return 0.5 * (design_point.initial_altitude_m + design_point.final_altitude_m)
    raise ValueError(
        f"{design_point.name} must define altitude_m or both "
        "initial_altitude_m and final_altitude_m."
    )


def design_point_flight_condition(design_point):
    atmosphere = asb.Atmosphere(altitude=representative_altitude_m(design_point))
    velocity_m_s = design_point.mach * atmosphere.speed_of_sound()
    dynamic_pressure_Pa = 0.5 * atmosphere.density() * velocity_m_s**2
    return velocity_m_s, dynamic_pressure_Pa


def propulsion_sizing_velocity_m_s(design_point):
    """Return representative velocity for converting thrust into propulsive power."""
    if (
        design_point.case in ("horizontal_acceleration", "climb_acceleration")
        and design_point.initial_mach is not None
        and design_point.final_mach is not None
    ):
        if (
            design_point.initial_altitude_m is not None
            and design_point.final_altitude_m is not None
        ):
            initial_atmosphere = asb.Atmosphere(altitude=design_point.initial_altitude_m)
            final_atmosphere = asb.Atmosphere(altitude=design_point.final_altitude_m)
            initial_velocity_m_s = (
                design_point.initial_mach * initial_atmosphere.speed_of_sound()
            )
            final_velocity_m_s = (
                design_point.final_mach * final_atmosphere.speed_of_sound()
            )
            return 0.5 * (initial_velocity_m_s + final_velocity_m_s)
        atmosphere = asb.Atmosphere(altitude=representative_altitude_m(design_point))
        return (
            0.5
            * (design_point.initial_mach + design_point.final_mach)
            * atmosphere.speed_of_sound()
        )
    return design_point_flight_condition(design_point)[0]


def propulsion_sizing_mach(design_point):
    """Return representative Mach number for electric propulsor sizing."""
    if (
        design_point.case in ("horizontal_acceleration", "climb_acceleration")
        and design_point.initial_mach is not None
        and design_point.final_mach is not None
    ):
        return 0.5 * (design_point.initial_mach + design_point.final_mach)
    return design_point.mach


def constraint_group(design_point):
    """Separate field constraints from in-flight propulsion mode constraints."""
    if design_point.case in {
        "takeoff_ground_roll_ideal",
        "takeoff_ground_roll",
        "takeoff_climb_angle",
        "braking_roll",
    }:
        return "field"
    return design_point.mode


def required_tw_for_design_point(
    design_point,
    wing_loading_N_m2,
    config=None,
    drag_geometry=None,
):
    """Return required sea-level static T/W for one literature constraint case."""
    atmosphere = asb.Atmosphere(altitude=representative_altitude_m(design_point))
    density_kg_m3 = atmosphere.density()
    gravity_m_s2 = 9.80665
    cd0 = design_point.cd0
    beta = design_point.beta
    alpha = design_point.alpha
    k1 = design_point.drag_polar_k1
    k2 = design_point.drag_polar_k2

    def apply_drag_build_up(altitude_m, velocity_m_s, lift_coefficient):
        if config is None or drag_geometry is None:
            return cd0, k1
        drag = drag_build_up_coefficients(
            config=config,
            drag_geometry=drag_geometry,
            altitude_m=altitude_m,
            velocity_m_s=velocity_m_s,
            lift_coefficient=lift_coefficient,
        )
        return drag["zero_lift_drag_coefficient"], drag["lift_dependent_k"]

    if design_point.case in ("horizontal_acceleration", "climb_acceleration"):
        initial_mach = design_point.initial_mach
        final_mach = design_point.final_mach
        if initial_mach is None or final_mach is None:
            raise ValueError(
                f"{design_point.name} must define initial_mach and final_mach."
            )
        if design_point.case == "climb_acceleration":
            has_altitude_segment = (
                design_point.initial_altitude_m is not None
                and design_point.final_altitude_m is not None
            )
            if has_altitude_segment:
                initial_altitude_m = design_point.initial_altitude_m
                final_altitude_m = design_point.final_altitude_m
                midpoint_altitude_m = 0.5 * (initial_altitude_m + final_altitude_m)
                atmosphere = asb.Atmosphere(altitude=midpoint_altitude_m)
                density_kg_m3 = atmosphere.density()
                initial_atmosphere = asb.Atmosphere(altitude=initial_altitude_m)
                final_atmosphere = asb.Atmosphere(altitude=final_altitude_m)
                initial_velocity_m_s = initial_mach * initial_atmosphere.speed_of_sound()
                final_velocity_m_s = final_mach * final_atmosphere.speed_of_sound()
                climb_rate_m_s = (
                    (final_altitude_m - initial_altitude_m)
                    / design_point.acceleration_time_s
                )
            else:
                speed_of_sound_m_s = atmosphere.speed_of_sound()
                initial_velocity_m_s = initial_mach * speed_of_sound_m_s
                final_velocity_m_s = final_mach * speed_of_sound_m_s
                climb_rate_m_s = design_point.climb_rate_m_s
        else:
            speed_of_sound_m_s = atmosphere.speed_of_sound()
            initial_velocity_m_s = initial_mach * speed_of_sound_m_s
            final_velocity_m_s = final_mach * speed_of_sound_m_s
            climb_rate_m_s = 0.0

        velocity_m_s = 0.5 * (initial_velocity_m_s + final_velocity_m_s)
        dynamic_pressure_Pa = 0.5 * density_kg_m3 * velocity_m_s**2
        cd0, k1 = apply_drag_build_up(
            representative_altitude_m(design_point),
            velocity_m_s,
            beta * wing_loading_N_m2 / dynamic_pressure_Pa,
        )
        acceleration_m_s2 = (
            (final_velocity_m_s - initial_velocity_m_s)
            / design_point.acceleration_time_s
        )
        return design_point_thrust_to_weight_from_wing_loading(
            wing_loading_N_m2=wing_loading_N_m2,
            dynamic_pressure_Pa=dynamic_pressure_Pa,
            velocity_m_s=velocity_m_s,
            installed_full_throttle_thrust_lapse=alpha,
            instantaneous_weight_fraction=beta,
            load_factor=1.0,
            drag_polar_k1=k1,
            drag_polar_k2=k2,
            zero_lift_drag_coefficient=cd0,
            climb_rate_m_s=climb_rate_m_s,
            acceleration_m_s2=acceleration_m_s2,
        )

    if design_point.case == "takeoff_ground_roll_ideal":
        cd0, k1 = apply_drag_build_up(
            representative_altitude_m(design_point),
            design_point.mach * atmosphere.speed_of_sound(),
            design_point.cl_max,
        )
        return (
            beta**2
            / alpha
            * design_point.speed_ratio**2
            / (
                design_point.ground_roll_m
                * density_kg_m3
                * gravity_m_s2
                * design_point.cl_max
            )
            * wing_loading_N_m2
        )

    if design_point.case == "takeoff_ground_roll":
        cd0, k1 = apply_drag_build_up(
            representative_altitude_m(design_point),
            design_point.mach * atmosphere.speed_of_sound(),
            design_point.cl_max / design_point.speed_ratio**2,
        )
        xi_to = cd0 - design_point.friction_coefficient * design_point.cl_max
        ideal = (
            beta**2
            / alpha
            * design_point.speed_ratio**2
            / (
                design_point.ground_roll_m
                * density_kg_m3
                * gravity_m_s2
                * design_point.cl_max
            )
            * wing_loading_N_m2
        )
        exponent = (
            design_point.ground_roll_m
            * density_kg_m3
            * gravity_m_s2
            * xi_to
            / (beta * wing_loading_N_m2)
        )
        denominator = 1.0 - np.exp(-exponent)
        full = (
            beta
            / alpha
            * (
                design_point.friction_coefficient
                + design_point.speed_ratio**2
                / design_point.cl_max
                * xi_to
                / denominator
            )
        )
        return np.where(np.abs(xi_to) < 1e-9, ideal, full)

    if design_point.case == "braking_roll":
        cd0, k1 = apply_drag_build_up(
            representative_altitude_m(design_point),
            design_point.mach * atmosphere.speed_of_sound(),
            design_point.cl_max / design_point.speed_ratio**2,
        )
        reverse_thrust_lapse = max(abs(alpha), 1e-9)
        ideal = (
            beta**2
            / reverse_thrust_lapse
            * design_point.speed_ratio**2
            / (
                design_point.braking_roll_m
                * density_kg_m3
                * gravity_m_s2
                * design_point.cl_max
            )
            * wing_loading_N_m2
        )
        xi_l = cd0 - design_point.friction_coefficient * design_point.cl_max
        exponent = (
            design_point.braking_roll_m
            * density_kg_m3
            * gravity_m_s2
            * xi_l
            / (beta * wing_loading_N_m2)
        )
        denominator = np.exp(exponent) - 1.0
        full = (
            beta
            / reverse_thrust_lapse
            * (
                design_point.speed_ratio**2
                / design_point.cl_max
                * xi_l
                / denominator
                - design_point.friction_coefficient
            )
        )
        return np.maximum(np.where(np.abs(xi_l) < 1e-9, ideal, full), 0.0)

    if design_point.case == "service_ceiling":
        lift_coefficient = design_point.lift_coefficient
        velocity_m_s = np.sqrt(
            2.0 * beta * wing_loading_N_m2 / (density_kg_m3 * lift_coefficient)
        )
        cd0, k1 = apply_drag_build_up(
            representative_altitude_m(design_point),
            velocity_m_s,
            lift_coefficient,
        )
        return (
            beta
            / alpha
            * (
                k1 * lift_coefficient
                + k2
                + cd0 / lift_coefficient
                + design_point.climb_rate_m_s / velocity_m_s
            )
        )

    if design_point.case == "takeoff_climb_angle":
        lift_coefficient = design_point.cl_max / design_point.speed_ratio**2
        cd0, k1 = apply_drag_build_up(
            representative_altitude_m(design_point),
            design_point.mach * atmosphere.speed_of_sound(),
            lift_coefficient,
        )
        return (
            beta
            / alpha
            * (
                k1 * lift_coefficient
                + k2
                + cd0 / lift_coefficient
                + np.sin(np.radians(design_point.climb_angle_deg))
            )
            * np.ones_like(wing_loading_N_m2)
        )

    velocity_m_s, dynamic_pressure_Pa = design_point_flight_condition(design_point)
    cd0, k1 = apply_drag_build_up(
        representative_altitude_m(design_point),
        velocity_m_s,
        design_point.load_factor * beta * wing_loading_N_m2 / dynamic_pressure_Pa,
    )
    return design_point_thrust_to_weight_from_wing_loading(
        wing_loading_N_m2=wing_loading_N_m2,
        dynamic_pressure_Pa=dynamic_pressure_Pa,
        velocity_m_s=velocity_m_s,
        installed_full_throttle_thrust_lapse=alpha,
        instantaneous_weight_fraction=beta,
        load_factor=design_point.load_factor,
        drag_polar_k1=k1,
        drag_polar_k2=k2,
        zero_lift_drag_coefficient=cd0,
        climb_rate_m_s=design_point.climb_rate_m_s,
        acceleration_m_s2=design_point.acceleration_m_s2,
    )


def weight_inputs_from_coupled_sizing(takeoff_mass_kg, planform_area_m2, config):
    fuselage_length_m = 4.0 * planform_area_m2**0.5
    propulsion = config_propulsion(
        config,
        mass_kg=0.0,
        engine_front_to_cockpit_length_m=0.35 * fuselage_length_m,
    )
    fuel = config_fuel(config)
    systems = config_systems(config)
    interiors = config_interiors(config)
    span_m = (planform_area_m2 * 3.0) ** 0.5
    root_chord_m = 2.0 * planform_area_m2 / (span_m * (1.0 + 0.25))
    tip_chord_m = 0.25 * root_chord_m
    vtail_area_m2 = 0.26 * planform_area_m2
    airplane = build_airplane(
        planform_area_m2=planform_area_m2,
        fuselage_length_m=fuselage_length_m,
        fuselage_height_m=0.12 * fuselage_length_m,
        fuselage_width_m=0.10 * fuselage_length_m,
        vtail_area_m2=vtail_area_m2,
        vtail_span_m=(vtail_area_m2 * 1.4) ** 0.5,
        vtail_dihedral_angle_deg=37.0,
        main_wing_tip_le_x_m=swept_wing_tip_le_x_m(
            span_m,
            root_chord_m,
            tip_chord_m,
            config.main_wing_quarter_chord_sweep_rad,
        ),
        sweep_25_rad=config.main_wing_quarter_chord_sweep_rad,
        vtail_le_x_m=0.55 * fuselage_length_m,
        tail_length_m=0.55 * fuselage_length_m,
        rudder_area_m2=0.025 * planform_area_m2,
        wing_mounted_control_area_m2=0.08 * planform_area_m2,
    )
    return Aircraft(
        airplane=airplane,
        mass_kg=takeoff_mass_kg,
        landing_mass_kg=0.85 * takeoff_mass_kg,
        propulsion=propulsion,
        fuel=fuel,
        payload=Payload(),
        systems=systems,
        interiors=interiors,
    ).to_weight_inputs()


def propulsion_system_sizing_breakdown(
    takeoff_mass_kg,
    wing_loading_N_m2,
    design_points,
    config,
    drag_geometry=None,
):
    """Return propulsor power split and component masses from a dry sizing point."""
    propulsion = config_propulsion(config)
    takeoff_weight_N = takeoff_mass_kg * 9.80665
    sizing_cases = []
    skipped_cases = []
    for design_point in design_points:
        if not design_point.include_in_governing:
            continue
        sizing_mach = propulsion_sizing_mach(design_point)
        if sizing_mach > propulsion.electric_propulsor_mach_limit:
            skipped_cases.append(
                {
                    "name": design_point.name,
                    "mach": sizing_mach,
                    "reason": "above electric propulsor Mach limit",
                }
            )
            continue
        required_thrust_to_weight = required_tw_for_design_point(
            design_point,
            wing_loading_N_m2,
            config=config,
            drag_geometry=drag_geometry,
        )
        required_thrust_N = required_thrust_to_weight * takeoff_weight_N
        velocity_m_s = propulsion_sizing_velocity_m_s(design_point)
        propulsive_power_W = (
            required_thrust_N * velocity_m_s / propulsion.propulsive_efficiency
        )
        sizing_cases.append(
            {
                "name": design_point.name,
                "mode": design_point.mode,
                "mach": sizing_mach,
                "required_thrust_to_weight": required_thrust_to_weight,
                "required_thrust_N": required_thrust_N,
                "velocity_m_s": velocity_m_s,
                "propulsive_power_W": propulsive_power_W,
            }
        )

    if not sizing_cases:
        raise ValueError(
            "No constraint design points remain for electric propulsor sizing. "
            "Increase electric_propulsor_mach_limit or add a low-Mach sizing case."
        )

    governing_case = max(sizing_cases, key=lambda item: item["propulsive_power_W"])
    governing_by_mode = {
        mode: max(
            (case for case in sizing_cases if case["mode"] == mode),
            key=lambda item: item["propulsive_power_W"],
        )
        for mode in sorted({case["mode"] for case in sizing_cases})
    }
    propulsive_power_W = governing_case["propulsive_power_W"]
    motor_power_each_W = propulsive_power_W / propulsion.number_propulsive_motors
    generator_electric_power_W = propulsive_power_W / propulsion.generator_efficiency
    generator_power_each_W = generator_electric_power_W / propulsion.number_generators
    turbine_shaft_power_W = (
        generator_electric_power_W / propulsion.turbine_mechanical_efficiency
    )
    motor_controller_mass_kg = (
        propulsion.number_propulsive_motors
        * motor_power_each_W
        / propulsion.motor_controller_power_density_W_kg
    )
    generator_mass_kg = (
        propulsion.number_generators
        * generator_power_each_W
        / propulsion.generator_power_density_W_kg
    )
    turbine_mass_kg = (
        propulsion.number_turbines
        * turbine_shaft_power_W
        / propulsion.turbine_power_density_W_kg
    )
    total_mass_kg = motor_controller_mass_kg + generator_mass_kg + turbine_mass_kg
    return {
        "governing_case": governing_case["name"],
        "governing_mode": governing_case["mode"],
        "governing_by_mode": governing_by_mode,
        "propulsive_power_W": propulsive_power_W,
        "motor_power_each_W": motor_power_each_W,
        "generator_power_each_W": generator_power_each_W,
        "turbine_shaft_power_W": turbine_shaft_power_W,
        "motor_controller_mass_kg": motor_controller_mass_kg,
        "generator_mass_kg": generator_mass_kg,
        "turbine_mass_kg": turbine_mass_kg,
        "total_mass_kg": total_mass_kg,
        "cases": sizing_cases,
        "skipped_cases": skipped_cases,
    }


def _solve_coupled_weight_volume_once(
    config,
    design_points,
    propulsion_mass_kg,
):
    design_points = design_points or default_constraint_design_points()
    opti = asb.Opti()
    log_planform_area_m2 = opti.variable(init_guess=np.log(80.0), scale=5.0)
    log_takeoff_mass_kg = opti.variable(init_guess=np.log(4000.0), scale=10.0)
    planform_area_m2 = np.exp(log_planform_area_m2)
    takeoff_mass_kg = np.exp(log_takeoff_mass_kg)
    opti.subject_to(log_planform_area_m2 >= np.log(1.0))
    opti.subject_to(log_takeoff_mass_kg >= np.log(100.0))

    fuselage_length_m = 4.0 * planform_area_m2**0.5
    span_m = (planform_area_m2 * 3.0) ** 0.5
    root_chord_m = 2.0 * planform_area_m2 / (span_m * (1.0 + 0.25))
    tip_chord_m = 0.25 * root_chord_m
    vtail_area_m2 = 0.26 * planform_area_m2
    airplane = build_airplane(
        planform_area_m2=planform_area_m2,
        fuselage_length_m=fuselage_length_m,
        fuselage_height_m=0.12 * fuselage_length_m,
        fuselage_width_m=0.10 * fuselage_length_m,
        vtail_area_m2=vtail_area_m2,
        vtail_span_m=(vtail_area_m2 * 1.4) ** 0.5,
        vtail_dihedral_angle_deg=37.0,
        main_wing_tip_le_x_m=swept_wing_tip_le_x_m(
            span_m,
            root_chord_m,
            tip_chord_m,
            config.main_wing_quarter_chord_sweep_rad,
        ),
        sweep_25_rad=config.main_wing_quarter_chord_sweep_rad,
        vtail_le_x_m=0.55 * fuselage_length_m,
        tail_length_m=0.55 * fuselage_length_m,
        rudder_area_m2=0.025 * planform_area_m2,
        wing_mounted_control_area_m2=0.08 * planform_area_m2,
    )
    propulsion = config_propulsion(
        config,
        mass_kg=propulsion_mass_kg,
        engine_front_to_cockpit_length_m=0.35 * fuselage_length_m,
    )
    fuel = config_fuel(config)
    payload = config_payload(config)
    systems = config_systems(config)
    interiors = config_interiors(config)
    aircraft = Aircraft(
        airplane=airplane,
        mass_kg=takeoff_mass_kg,
        landing_mass_kg=0.85 * takeoff_mass_kg,
        propulsion=propulsion,
        fuel=fuel,
        payload=payload,
        systems=systems,
        interiors=interiors,
    )
    volume_inputs = aircraft.to_volume_inputs()
    volume = aircraft_volume_breakdown(
        replace(
            volume_inputs,
            void_volume_coefficient=config.void_volume_coefficient,
        )
    )
    opti.subject_to(volume["kuechemann_slenderness_parameter"] == config.kuechemann_tau)

    weight = _weight_breakdown(aircraft.to_weight_inputs())
    opti.subject_to(takeoff_mass_kg == weight["total_aircraft_mass_kg"])
    opti.minimize(takeoff_mass_kg)

    try:
        sol = opti.solve()
    except RuntimeError:
        sol = opti.debug

    solved_planform_area_m2 = sol(planform_area_m2)
    solved_takeoff_mass_kg = sol(takeoff_mass_kg)
    solved_volume = {
        "total_aircraft_volume_m3": sol(volume["total_aircraft_volume_m3"]),
        "kuechemann_slenderness_parameter": sol(volume["kuechemann_slenderness_parameter"]),
    }
    solved_weight = {
        "operating_empty_without_engine_lb": sol(weight["operating_empty_without_engine_lb"]),
        "duality_weight_lb": sol(weight["duality_weight_lb"]),
        "tank_dry_weight_lb": sol(weight["tank_dry_weight_lb"]),
        "total_aircraft_weight_lb": sol(weight["total_aircraft_weight_lb"]),
        "total_aircraft_mass_kg": sol(weight["total_aircraft_mass_kg"]),
    }
    solved_wing_loading_N_m2 = solved_takeoff_mass_kg * 9.80665 / solved_planform_area_m2
    solved_drag_geometry = drag_geometry_from_planform_area(
        solved_planform_area_m2,
        config,
        config_propulsion(config).number_engines,
    )
    design_cruise_atmosphere = asb.Atmosphere(altitude=config.design_cruise_altitude_m)
    design_cruise_velocity_m_s = (
        config.design_cruise_mach * design_cruise_atmosphere.speed_of_sound()
    )
    design_cruise_dynamic_pressure_Pa = (
        0.5 * design_cruise_atmosphere.density() * design_cruise_velocity_m_s**2
    )
    solved_drag_geometry["design_lift_coefficient"] = (
        solved_wing_loading_N_m2 / design_cruise_dynamic_pressure_Pa
    )

    return {
        "planform_area_m2": solved_planform_area_m2,
        "takeoff_mass_kg": solved_takeoff_mass_kg,
        "takeoff_weight_N": solved_takeoff_mass_kg * 9.80665,
        "wing_loading_N_m2": solved_wing_loading_N_m2,
        "drag_geometry": solved_drag_geometry,
        "aircraft": Aircraft(
            airplane=airplane,
            mass_kg=solved_takeoff_mass_kg,
            landing_mass_kg=0.85 * solved_takeoff_mass_kg,
            propulsion=config_propulsion(
                config,
                mass_kg=propulsion_mass_kg,
                engine_front_to_cockpit_length_m=0.35 * sol(fuselage_length_m),
            ),
            fuel=fuel,
            payload=payload,
            systems=systems,
            interiors=interiors,
        ),
        "volume": solved_volume,
        "weight": solved_weight,
    }


def solve_coupled_weight_volume(config=ConstraintDiagramConfig(), design_points=None):
    design_points = design_points or default_constraint_design_points()
    dry = _solve_coupled_weight_volume_once(
        config=config,
        design_points=design_points,
        propulsion_mass_kg=0.0,
    )
    propulsion_sizing = propulsion_system_sizing_breakdown(
        takeoff_mass_kg=dry["takeoff_mass_kg"],
        wing_loading_N_m2=dry["wing_loading_N_m2"],
        design_points=design_points,
        config=config,
        drag_geometry=dry["drag_geometry"],
    )
    coupled = _solve_coupled_weight_volume_once(
        config=config,
        design_points=design_points,
        propulsion_mass_kg=propulsion_sizing["total_mass_kg"],
    )
    coupled["dry_without_propulsion"] = dry
    coupled["propulsion_sizing"] = propulsion_sizing
    return coupled


def build_constraint_diagram(config=ConstraintDiagramConfig(), design_points=None):
    design_points = design_points or default_constraint_design_points()
    coupled = solve_coupled_weight_volume(config, design_points=design_points)
    engine_deck_rows = read_pycycle_engine_deck(config.engine_deck_csv)

    wing_loading_N_m2 = np.linspace(
        config.wing_loading_min_N_m2,
        config.wing_loading_max_N_m2,
        config.wing_loading_points,
    )

    curves = {}
    for design_point in design_points:
        velocity_m_s, dynamic_pressure_Pa = design_point_flight_condition(design_point)
        drag_velocity_m_s = propulsion_sizing_velocity_m_s(design_point)
        required_tw = required_tw_for_design_point(
            design_point,
            wing_loading_N_m2,
            config=config,
            drag_geometry=coupled["drag_geometry"],
        )
        available_thrust_N = pycycle_design_point_thrust_N(
            engine_deck_rows,
            design_point,
            thrust_scale=config.thrust_scale,
        )
        available_tw = available_thrust_N / coupled["takeoff_weight_N"]
        required_at_coupled_wing_loading = required_tw_for_design_point(
            design_point,
            coupled["wing_loading_N_m2"],
            config=config,
            drag_geometry=coupled["drag_geometry"],
        )
        drag = drag_build_up_coefficients(
            config=config,
            drag_geometry=coupled["drag_geometry"],
            altitude_m=representative_altitude_m(design_point),
            velocity_m_s=drag_velocity_m_s,
            lift_coefficient=(
                coupled["wing_loading_N_m2"]
                / (
                    0.5
                    * asb.Atmosphere(
                        altitude=representative_altitude_m(design_point)
                    ).density()
                    * drag_velocity_m_s**2
                )
            ),
        )
        curves[design_point.name] = {
            "design_point": design_point,
            "velocity_m_s": velocity_m_s,
            "dynamic_pressure_Pa": dynamic_pressure_Pa,
            "parasite_cd0": drag["parasite_cd0"],
            "wave_cd0": drag["wave_cd0"],
            "zero_lift_drag_coefficient": drag["zero_lift_drag_coefficient"],
            "lift_dependent_k": drag["lift_dependent_k"],
            "required_thrust_to_weight": required_tw,
            "available_thrust_N": available_thrust_N,
            "available_thrust_to_weight": available_tw,
            "required_at_coupled_wing_loading": required_at_coupled_wing_loading,
            "margin_at_coupled_wing_loading": available_tw - required_at_coupled_wing_loading,
        }

    groups = sorted({constraint_group(point) for point in design_points if point.include_in_governing})
    governing_required_tw_by_group = {
        group: np.max(
            np.array(
                [
                    curves[point.name]["required_thrust_to_weight"]
                    for point in design_points
                    if point.include_in_governing and constraint_group(point) == group
                ]
            ),
            axis=0,
        )
        for group in groups
    }
    required_tw_at_coupled_by_group = {
        group: np.max(
            np.array(
                [
                    curves[point.name]["required_at_coupled_wing_loading"]
                    for point in design_points
                    if point.include_in_governing and constraint_group(point) == group
                ]
            )
        )
        for group in groups
    }
    required_stack = np.array(
        [
            curves[design_point.name]["required_thrust_to_weight"]
            for design_point in design_points
            if design_point.include_in_governing
        ]
    )
    governing_required_tw = np.max(required_stack, axis=0)
    coupled["required_thrust_to_weight_by_group"] = required_tw_at_coupled_by_group
    modes = sorted({point.mode for point in design_points if point.include_in_governing})
    required_tw_at_coupled_by_mode = {
        mode: np.max(
            np.array(
                [
                    curves[point.name]["required_at_coupled_wing_loading"]
                    for point in design_points
                    if point.include_in_governing and point.mode == mode
                ]
            )
        )
        for mode in modes
    }
    coupled["required_thrust_to_weight_by_mode"] = required_tw_at_coupled_by_mode
    coupled["required_thrust_to_weight"] = np.max(
        np.array(list(required_tw_at_coupled_by_group.values()))
    )

    return {
        "wing_loading_N_m2": wing_loading_N_m2,
        "governing_required_thrust_to_weight": governing_required_tw,
        "governing_required_thrust_to_weight_by_group": governing_required_tw_by_group,
        "curves": curves,
        "coupled": coupled,
    }


def engine_deck_aircraft_sizing(config=ConstraintDiagramConfig(), design_points=None):
    """Return the coupled aircraft sizing state used to build engine-deck requests."""
    design_points = design_points or default_constraint_design_points()
    result = build_constraint_diagram(config=config, design_points=design_points)
    coupled = result["coupled"]
    cruise_design_point = next(
        point
        for point in design_points
        if point.name == "Case 1: constant-altitude/speed cruise"
    )
    sizing_thrust_by_group_N = {
        mode: required_tw * coupled["takeoff_weight_N"]
        for mode, required_tw in coupled["required_thrust_to_weight_by_group"].items()
    }
    sizing_thrust_by_mode_N = {
        mode: required_tw * coupled["takeoff_weight_N"]
        for mode, required_tw in coupled["required_thrust_to_weight_by_mode"].items()
    }
    return {
        "constraint_result": result,
        "config": config,
        "design_points": design_points,
        "cruise_constraint_design_point": cruise_design_point,
        "cruise_constraint_mach": propulsion_sizing_mach(cruise_design_point),
        "cruise_constraint_altitude_ft": representative_altitude_m(cruise_design_point) / u.foot,
        "planform_area_m2": coupled["planform_area_m2"],
        "takeoff_mass_kg": coupled["takeoff_mass_kg"],
        "takeoff_weight_N": coupled["takeoff_weight_N"],
        "wing_loading_N_m2": coupled["wing_loading_N_m2"],
        "drag_geometry": coupled["drag_geometry"],
        "aircraft": coupled["aircraft"],
        "sizing_required_thrust_N": max(sizing_thrust_by_group_N.values()),
        "sizing_required_thrust_by_group_N": sizing_thrust_by_group_N,
        "sizing_required_thrust_by_mode_N": sizing_thrust_by_mode_N,
        "sizing_required_thrust_to_weight": coupled["required_thrust_to_weight"],
        "sizing_required_thrust_to_weight_by_group": coupled["required_thrust_to_weight_by_group"],
        "sizing_required_thrust_to_weight_by_mode": coupled["required_thrust_to_weight_by_mode"],
        "propulsion_sizing": coupled["propulsion_sizing"],
    }


def plot_constraint_diagram(result, save_plot=None, show_plot=False):
    import matplotlib.pyplot as plt

    wing_loading_N_m2 = result["wing_loading_N_m2"]
    wing_loading_lb_ft2 = wing_loading_N_m2 / (u.lbf / u.foot**2)
    coupled_wing_loading_lb_ft2 = result["coupled"]["wing_loading_N_m2"] / (
        u.lbf / u.foot**2
    )
    reference_wing_loadings_lb_ft2 = {
        "Concorde": 408000.0 / 3856.0,
        "SR-71 Blackbird": 84.0,
    }
    reference_thrust_to_weights = {
        "Concorde": 4.0 * 38050.0 / 408000.0,
        "SR-71 Blackbird": 2.0 * 32500.0 / 140000.0,
    }
    reference_colors = {
        "Concorde": "tab:purple",
        "SR-71 Blackbird": "tab:brown",
    }
    panels = (
        ("Fan mode constraints", ("fan",), "tab:blue"),
        ("Fan + afterburner mode constraints", ("fan_ab",), "tab:orange"),
        ("Ramjet mode constraints", ("ramjet",), "tab:red"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.6), sharex=True, sharey=True, constrained_layout=True)
    for ax, (title, groups, governing_color) in zip(axes, panels):
        for name, curve in result["curves"].items():
            design_point = curve["design_point"]
            if design_point.mode not in groups:
                continue
            ax.plot(
                wing_loading_lb_ft2,
                curve["required_thrust_to_weight"],
                linewidth=1.7,
                label=f"{design_point.name} ({design_point.mode}, M{float(design_point.mach):.2g})",
            )

        governing_curves = [
            curve["required_thrust_to_weight"]
            for curve in result["curves"].values()
            if curve["design_point"].include_in_governing and curve["design_point"].mode in groups
        ]
        if governing_curves:
            ax.plot(
                wing_loading_lb_ft2,
                np.max(np.array(governing_curves), axis=0),
                color=governing_color,
                linestyle="--",
                linewidth=2.4,
                label="Regime governing",
            )
        ax.axvline(
            coupled_wing_loading_lb_ft2,
            color="black",
            linestyle=":",
            linewidth=1.8,
            label="Coupled W/S",
        )
        for reference_name, reference_wing_loading_lb_ft2 in reference_wing_loadings_lb_ft2.items():
            ax.axvline(
                reference_wing_loading_lb_ft2,
                color=reference_colors[reference_name],
                linestyle="--",
                linewidth=1.0,
                alpha=0.6,
            )
        for reference_name, reference_thrust_to_weight in reference_thrust_to_weights.items():
            ax.axhline(
                reference_thrust_to_weight,
                color=reference_colors[reference_name],
                linestyle="--",
                linewidth=1.0,
                alpha=0.6,
            )
        ax.set_title(title)
        ax.set_xlabel("Wing loading W/S, lb/ft^2")
        ax.set_xlim(20.0, 150.0)
        ax.set_ylim(bottom=0.0)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", fontsize=7)

    axes[0].set_ylabel("Required installed thrust-to-weight at case condition")
    fig.suptitle("Constraint Requirements by Propulsion Mode")

    if save_plot is not None:
        save_plot = Path(save_plot)
        save_plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_plot, dpi=180)
    if show_plot:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_drag_build_up_vs_mach(
    result,
    config,
    altitude_m=60000.0 * u.foot,
    save_plot=plots_dir / "_drag_vs_mach.png",
    show_plot=False,
):
    """Plot the computed drag build-up versus Mach for the coupled geometry."""
    import matplotlib.pyplot as plt

    mach = np.linspace(0.1, 3.2, 250)
    atmosphere = asb.Atmosphere(altitude=altitude_m)
    velocity_m_s = mach * atmosphere.speed_of_sound()
    dynamic_pressure_Pa = 0.5 * atmosphere.density() * velocity_m_s**2
    lift_coefficient = result["coupled"]["wing_loading_N_m2"] / dynamic_pressure_Pa
    drag = drag_build_up_coefficients(
        config=config,
        drag_geometry=result["coupled"]["drag_geometry"],
        altitude_m=altitude_m,
        velocity_m_s=velocity_m_s,
    )
    drag_coefficient = (
        drag["zero_lift_drag_coefficient"]
        + drag["lift_dependent_k"] * lift_coefficient**2.0
    )
    valid_lift = lift_coefficient <= config.drag_plot_cl_max
    plotted_drag_coefficient = np.where(valid_lift, drag_coefficient, np.nan)
    plotted_cd0 = np.where(valid_lift, drag["zero_lift_drag_coefficient"], np.nan)
    plotted_parasite_cd0 = np.where(valid_lift, drag["parasite_cd0"], np.nan)
    plotted_wave_cd0 = np.where(valid_lift, drag["wave_cd0"], np.nan)

    fig, cd_axis = plt.subplots(figsize=(9.5, 5.5), constrained_layout=True)

    cd_axis.plot(
        mach,
        plotted_drag_coefficient,
        color="tab:purple",
        linewidth=2.4,
        label=f"CD total, CL <= {config.drag_plot_cl_max:g}",
    )
    cd_axis.plot(
        mach,
        plotted_cd0,
        color="black",
        linewidth=2.2,
        label="CD0 total",
    )
    cd_axis.plot(
        mach,
        plotted_parasite_cd0,
        color="tab:blue",
        linewidth=1.8,
        label="CD0 parasite",
    )
    cd_axis.plot(
        mach,
        plotted_wave_cd0,
        color="tab:red",
        linewidth=1.8,
        label="CD0 wave",
    )
    mcrit = config.drag_divergence_mach - config.critical_mach_offset_from_mdd
    markers = {
        "Mcrit": mcrit,
        "MDD": config.drag_divergence_mach,
        "M=1.2": config.supersonic_wave_drag_start_mach,
    }
    for label, marker_mach in markers.items():
        cd_axis.axvline(marker_mach, color="0.55", linestyle=":", linewidth=1.0)
        cd_axis.text(
            marker_mach,
            0.98,
            label,
            rotation=90,
            va="top",
            ha="right",
            color="0.35",
            transform=cd_axis.get_xaxis_transform(),
        )
    if np.any(valid_lift):
        first_valid_index = int(np.argmax(valid_lift))
        stall_limited_mach = mach[first_valid_index]
        cd_axis.axvline(
            stall_limited_mach,
            color="tab:purple",
            linestyle=":",
            linewidth=1.0,
        )
        cd_axis.text(
            stall_limited_mach,
            0.98,
            "CLmax",
            rotation=90,
            va="top",
            ha="left",
            color="tab:purple",
            transform=cd_axis.get_xaxis_transform(),
        )

    cd_axis.set_xlabel("Mach number")
    cd_axis.set_ylabel("Drag coefficient")
    cd_axis.set_title(f"Drag build-up vs Mach at {altitude_m / u.foot:.0f} ft")
    cd_axis.grid(True, alpha=0.3)
    if np.any(valid_lift):
        cd_axis.set_xlim(float(mach[int(np.argmax(valid_lift))]), float(mach[-1]))
    else:
        cd_axis.set_xlim(float(mach[0]), float(mach[-1]))
    cd_axis.set_ylim(bottom=0.0)

    lines = cd_axis.get_lines()
    cd_axis.legend(lines, [line.get_label() for line in lines], loc="upper left")

    if save_plot is not None:
        save_plot = Path(save_plot)
        save_plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_plot, dpi=180)
    if show_plot:
        plt.show()
    else:
        plt.close(fig)
    return fig, drag


def plot_cd0_vs_mach(
    result,
    config,
    altitude_m=60000.0 * u.foot,
    save_plot=plots_dir / "_cd0_vs_mach.png",
    show_plot=False,
):
    """Plot zero-lift drag components versus Mach for the coupled geometry."""
    import matplotlib.pyplot as plt

    mach = np.linspace(0.1, 3.2, 250)
    atmosphere = asb.Atmosphere(altitude=altitude_m)
    velocity_m_s = mach * atmosphere.speed_of_sound()
    drag = drag_build_up_coefficients(
        config=config,
        drag_geometry=result["coupled"]["drag_geometry"],
        altitude_m=altitude_m,
        velocity_m_s=velocity_m_s,
    )

    fig, ax = plt.subplots(figsize=(9.5, 5.5), constrained_layout=True)
    ax.plot(
        mach,
        drag["zero_lift_drag_coefficient"],
        color="black",
        linewidth=2.4,
        label="CD0 total",
    )
    ax.plot(
        mach,
        drag["parasite_cd0"],
        color="tab:blue",
        linewidth=1.8,
        label="CD0 parasite",
    )
    ax.plot(
        mach,
        drag["wave_cd0"],
        color="tab:red",
        linewidth=1.8,
        label="CD0 wave",
    )

    mcrit = config.drag_divergence_mach - config.critical_mach_offset_from_mdd
    markers = {
        "Mcrit": mcrit,
        "MDD": config.drag_divergence_mach,
        "M=1.2": config.supersonic_wave_drag_start_mach,
    }
    for label, marker_mach in markers.items():
        ax.axvline(marker_mach, color="0.55", linestyle=":", linewidth=1.0)
        ax.text(
            marker_mach,
            0.98,
            label,
            rotation=90,
            va="top",
            ha="right",
            color="0.35",
            transform=ax.get_xaxis_transform(),
        )

    ax.set_xlabel("Mach number")
    ax.set_ylabel("Zero-lift drag coefficient")
    ax.set_title(f"Zero-lift drag vs Mach at {altitude_m / u.foot:.0f} ft")
    ax.set_xlim(float(mach[0]), float(mach[-1]))
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left")

    if save_plot is not None:
        save_plot = Path(save_plot)
        save_plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_plot, dpi=180)
    if show_plot:
        plt.show()
    else:
        plt.close(fig)
    return fig, drag


def plot_level_flight_performance_vs_mach(
    result,
    config,
    altitude_m=60000.0 * u.foot,
    save_plot=plots_dir / "_level_flight_performance_vs_mach.png",
    show_plot=False,
):
    """Plot level-flight lift, drag, and L/D versus Mach for the coupled aircraft."""
    import matplotlib.pyplot as plt

    mach = np.linspace(0.1, 3.2, 250)
    atmosphere = asb.Atmosphere(altitude=altitude_m)
    velocity_m_s = mach * atmosphere.speed_of_sound()
    dynamic_pressure_Pa = 0.5 * atmosphere.density() * velocity_m_s**2
    wing_area_m2 = result["coupled"]["planform_area_m2"]
    weight_N = result["coupled"]["takeoff_weight_N"]
    lift_coefficient = weight_N / (dynamic_pressure_Pa * wing_area_m2)
    drag = drag_build_up_coefficients(
        config=config,
        drag_geometry=result["coupled"]["drag_geometry"],
        altitude_m=altitude_m,
        velocity_m_s=velocity_m_s,
        lift_coefficient=lift_coefficient,
    )
    drag_coefficient = (
        drag["zero_lift_drag_coefficient"]
        + drag["lift_dependent_k"] * lift_coefficient**2.0
    )
    drag_N = dynamic_pressure_Pa * wing_area_m2 * drag_coefficient
    lift_N = np.ones_like(mach) * weight_N
    lift_to_drag = lift_N / drag_N
    valid_lift = lift_coefficient <= config.drag_plot_cl_max

    plotted_drag_kN = np.where(valid_lift, drag_N / 1000.0, np.nan)
    plotted_lift_kN = np.where(valid_lift, lift_N / 1000.0, np.nan)
    plotted_lift_to_drag = np.where(valid_lift, lift_to_drag, np.nan)

    fig, force_axis = plt.subplots(figsize=(9.5, 5.5), constrained_layout=True)
    ld_axis = force_axis.twinx()

    force_axis.plot(
        mach,
        plotted_drag_kN,
        color="tab:red",
        linewidth=2.2,
        label="Drag",
    )
    force_axis.plot(
        mach,
        plotted_lift_kN,
        color="tab:blue",
        linewidth=1.8,
        label="Lift",
    )
    ld_axis.plot(
        mach,
        plotted_lift_to_drag,
        color="black",
        linewidth=2.0,
        linestyle="--",
        label="L/D",
    )

    if np.any(valid_lift):
        stall_limited_mach = mach[int(np.argmax(valid_lift))]
        force_axis.axvline(
            stall_limited_mach,
            color="0.45",
            linestyle=":",
            linewidth=1.0,
        )
        force_axis.text(
            stall_limited_mach,
            0.98,
            "CLmax",
            rotation=90,
            va="top",
            ha="left",
            color="0.35",
            transform=force_axis.get_xaxis_transform(),
        )
        force_axis.set_xlim(float(stall_limited_mach), float(mach[-1]))
    else:
        force_axis.set_xlim(float(mach[0]), float(mach[-1]))

    force_axis.set_xlabel("Mach number")
    force_axis.set_ylabel("Force (kN)")
    ld_axis.set_ylabel("L/D")
    force_axis.set_title(f"Level-flight performance vs Mach at {altitude_m / u.foot:.0f} ft")
    force_axis.grid(True, alpha=0.3)
    force_axis.set_ylim(bottom=0.0)
    ld_axis.set_ylim(bottom=0.0)

    lines = force_axis.get_lines() + ld_axis.get_lines()
    force_axis.legend(lines, [line.get_label() for line in lines], loc="upper left")

    if save_plot is not None:
        save_plot = Path(save_plot)
        save_plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_plot, dpi=180)
    if show_plot:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_altitude_sweep_lift_to_drag_vs_mach(
    result,
    config,
    altitude_ft_values=(30000.0, 45000.0, 60000.0, 70000.0, 80000.0),
    save_plot=plots_dir / "_altitude_sweep_lift_to_drag_vs_mach.png",
    show_plot=False,
):
    """Plot level-flight L/D versus Mach for several altitudes."""
    import matplotlib.pyplot as plt

    mach = np.linspace(0.1, 3.2, 250)
    wing_area_m2 = result["coupled"]["planform_area_m2"]
    weight_N = result["coupled"]["takeoff_weight_N"]

    fig, ax = plt.subplots(figsize=(9.5, 5.5), constrained_layout=True)

    for altitude_ft in altitude_ft_values:
        altitude_m = altitude_ft * u.foot
        atmosphere = asb.Atmosphere(altitude=altitude_m)
        velocity_m_s = mach * atmosphere.speed_of_sound()
        dynamic_pressure_Pa = 0.5 * atmosphere.density() * velocity_m_s**2
        lift_coefficient = weight_N / (dynamic_pressure_Pa * wing_area_m2)
        drag = drag_build_up_coefficients(
            config=config,
            drag_geometry=result["coupled"]["drag_geometry"],
            altitude_m=altitude_m,
            velocity_m_s=velocity_m_s,
            lift_coefficient=lift_coefficient,
        )
        drag_coefficient = (
            drag["zero_lift_drag_coefficient"]
            + drag["lift_dependent_k"] * lift_coefficient**2.0
        )
        lift_to_drag = lift_coefficient / drag_coefficient
        valid_lift = lift_coefficient <= config.drag_plot_cl_max
        ax.plot(
            mach,
            np.where(valid_lift, lift_to_drag, np.nan),
            linewidth=2.0,
            label=f"{altitude_ft / 1000.0:.0f} kft",
        )

    ax.set_xlabel("Mach number")
    ax.set_ylabel("L/D")
    ax.set_title("Level-flight L/D vs Mach by altitude")
    ax.set_xlim(float(mach[0]), float(mach[-1]))
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", title="Altitude")

    if save_plot is not None:
        save_plot = Path(save_plot)
        save_plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_plot, dpi=180)
    if show_plot:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_drag_terms_vs_mach(
    result,
    config,
    altitude_m=60000.0 * u.foot,
    save_plot=plots_dir / "_drag_terms_vs_mach.png",
    show_plot=False,
):
    """Plot CL, CD0, CDi, and total CD versus Mach for one altitude."""
    import matplotlib.pyplot as plt

    mach = np.linspace(0.1, 3.2, 250)
    atmosphere = asb.Atmosphere(altitude=altitude_m)
    velocity_m_s = mach * atmosphere.speed_of_sound()
    dynamic_pressure_Pa = 0.5 * atmosphere.density() * velocity_m_s**2
    wing_area_m2 = result["coupled"]["planform_area_m2"]
    weight_N = result["coupled"]["takeoff_weight_N"]
    lift_coefficient = weight_N / (dynamic_pressure_Pa * wing_area_m2)
    drag = drag_build_up_coefficients(
        config=config,
        drag_geometry=result["coupled"]["drag_geometry"],
        altitude_m=altitude_m,
        velocity_m_s=velocity_m_s,
        lift_coefficient=lift_coefficient,
    )
    induced_drag_coefficient = drag["lift_dependent_k"] * lift_coefficient**2.0
    total_drag_coefficient = (
        drag["zero_lift_drag_coefficient"] + induced_drag_coefficient
    )
    valid_lift = lift_coefficient <= config.drag_plot_cl_max

    fig, cd_axis = plt.subplots(figsize=(9.5, 5.5), constrained_layout=True)
    cl_axis = cd_axis.twinx()

    cd_axis.plot(
        mach,
        np.where(valid_lift, total_drag_coefficient, np.nan),
        color="black",
        linewidth=2.3,
        label="CD total",
    )
    cd_axis.plot(
        mach,
        np.where(valid_lift, drag["zero_lift_drag_coefficient"], np.nan),
        color="tab:blue",
        linewidth=1.8,
        label="CD0",
    )
    cd_axis.plot(
        mach,
        np.where(valid_lift, induced_drag_coefficient, np.nan),
        color="tab:red",
        linewidth=1.8,
        label="CDi",
    )
    cl_axis.plot(
        mach,
        np.where(valid_lift, lift_coefficient, np.nan),
        color="tab:green",
        linewidth=1.8,
        linestyle="--",
        label="CL",
    )

    mcrit = config.drag_divergence_mach - config.critical_mach_offset_from_mdd
    markers = {
        "Mcrit": mcrit,
        "MDD": config.drag_divergence_mach,
        "M=1.2": config.supersonic_wave_drag_start_mach,
    }
    for label, marker_mach in markers.items():
        cd_axis.axvline(marker_mach, color="0.55", linestyle=":", linewidth=1.0)
        cd_axis.text(
            marker_mach,
            0.98,
            label,
            rotation=90,
            va="top",
            ha="right",
            color="0.35",
            transform=cd_axis.get_xaxis_transform(),
        )

    if np.any(valid_lift):
        stall_limited_mach = mach[int(np.argmax(valid_lift))]
        cd_axis.axvline(
            stall_limited_mach,
            color="0.45",
            linestyle=":",
            linewidth=1.0,
        )
        cd_axis.text(
            stall_limited_mach,
            0.98,
            "CLmax",
            rotation=90,
            va="top",
            ha="left",
            color="0.35",
            transform=cd_axis.get_xaxis_transform(),
        )
        cd_axis.set_xlim(float(stall_limited_mach), float(mach[-1]))
    else:
        cd_axis.set_xlim(float(mach[0]), float(mach[-1]))

    cd_axis.set_xlabel("Mach number")
    cd_axis.set_ylabel("Drag coefficient")
    cl_axis.set_ylabel("Lift coefficient")
    cd_axis.set_title(f"Level-flight drag terms vs Mach at {altitude_m / u.foot:.0f} ft")
    cd_axis.set_ylim(bottom=0.0)
    cl_axis.set_ylim(bottom=0.0)
    cd_axis.grid(True, alpha=0.3)

    lines = cd_axis.get_lines() + cl_axis.get_lines()
    cd_axis.legend(lines, [line.get_label() for line in lines], loc="upper right")

    if save_plot is not None:
        save_plot = Path(save_plot)
        save_plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_plot, dpi=180)
    if show_plot:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_lift_dependent_k_vs_mach(
    result,
    config,
    altitude_m=60000.0 * u.foot,
    save_plot=plots_dir / "_k_vs_mach.png",
    show_plot=False,
):
    """Plot the lift-dependent drag factor K versus Mach."""
    import matplotlib.pyplot as plt

    mach = np.linspace(0.1, 3.2, 250)
    fig, ax = plt.subplots(figsize=(9.5, 5.5), constrained_layout=True)

    atmosphere = asb.Atmosphere(altitude=altitude_m)
    velocity_m_s = mach * atmosphere.speed_of_sound()
    dynamic_pressure_Pa = 0.5 * atmosphere.density() * velocity_m_s**2
    lift_coefficient = result["coupled"]["wing_loading_N_m2"] / dynamic_pressure_Pa
    drag = drag_build_up_coefficients(
        config=config,
        drag_geometry=result["coupled"]["drag_geometry"],
        altitude_m=altitude_m,
        velocity_m_s=velocity_m_s,
        lift_coefficient=lift_coefficient,
    )
    ax.plot(
        mach,
        drag["lift_dependent_k"],
        linewidth=2.0,
        label="K",
    )

    for label, marker_mach in {
        "M=1": 1.0,
        "M_LE": leading_edge_sonic_mach(result["coupled"]["drag_geometry"]),
    }.items():
        ax.axvline(marker_mach, color="0.55", linestyle=":", linewidth=1.0)
        ax.text(
            marker_mach,
            0.98,
            label,
            rotation=90,
            va="top",
            ha="right",
            color="0.35",
            transform=ax.get_xaxis_transform(),
        )

    ax.set_xlabel("Mach number")
    ax.set_ylabel("Lift-dependent drag factor K")
    ax.set_title(f"K vs Mach at {altitude_m / u.foot:.0f} ft")
    ax.set_xlim(float(mach[0]), float(mach[-1]))
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left")

    if save_plot is not None:
        save_plot = Path(save_plot)
        save_plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_plot, dpi=180)
    if show_plot:
        plt.show()
    else:
        plt.close(fig)
    return fig


def main():
    # Edit run options here.
    config = ConstraintDiagramConfig(
        engine_deck_csv="propulsion/data/example_engine_deck.csv",
        thrust_scale=1.0,
        save_plot=plots_dir / "_constraint_diagram.png",
        show_plot=False,
    )
    result = build_constraint_diagram(config)
    plot_constraint_diagram(result, save_plot=config.save_plot, show_plot=config.show_plot)
    drag_plot = plots_dir / "_drag_vs_mach.png"
    plot_drag_build_up_vs_mach(
        result,
        config,
        altitude_m=60000.0 * u.foot,
        save_plot=drag_plot,
        show_plot=config.show_plot,
    )
    cd0_plot = plots_dir / "_cd0_vs_mach.png"
    plot_cd0_vs_mach(
        result,
        config,
        altitude_m=60000.0 * u.foot,
        save_plot=cd0_plot,
        show_plot=config.show_plot,
    )
    performance_plot = plots_dir / "_level_flight_performance_vs_mach.png"
    plot_level_flight_performance_vs_mach(
        result,
        config,
        altitude_m=60000.0 * u.foot,
        save_plot=performance_plot,
        show_plot=config.show_plot,
    )
    altitude_sweep_plot = plots_dir / "_altitude_sweep_lift_to_drag_vs_mach.png"
    plot_altitude_sweep_lift_to_drag_vs_mach(
        result,
        config,
        save_plot=altitude_sweep_plot,
        show_plot=config.show_plot,
    )
    drag_terms_plot = plots_dir / "_drag_terms_vs_mach.png"
    plot_drag_terms_vs_mach(
        result,
        config,
        altitude_m=60000.0 * u.foot,
        save_plot=drag_terms_plot,
        show_plot=config.show_plot,
    )
    drag_terms_30k_plot = plots_dir / "_drag_terms_vs_mach_30kft.png"
    plot_drag_terms_vs_mach(
        result,
        config,
        altitude_m=30000.0 * u.foot,
        save_plot=drag_terms_30k_plot,
        show_plot=config.show_plot,
    )
    k_plot = plots_dir / "_k_vs_mach.png"
    plot_lift_dependent_k_vs_mach(
        result,
        config,
        altitude_m=60000.0 * u.foot,
        save_plot=k_plot,
        show_plot=config.show_plot,
    )
    k_30k_plot = plots_dir / "_k_vs_mach_30kft.png"
    plot_lift_dependent_k_vs_mach(
        result,
        config,
        altitude_m=30000.0 * u.foot,
        save_plot=k_30k_plot,
        show_plot=config.show_plot,
    )

    coupled = result["coupled"]
    print(" coupled constraint solution")
    print(f"S_plan: {coupled['planform_area_m2']:.3f} m^2")
    print(f"TO mass: {coupled['takeoff_mass_kg']:.3f} kg")
    print(f"OEW without engine: {coupled['weight']['operating_empty_without_engine_lb'] * u.lbm:.3f} kg")
    print(f"Duality mass: {coupled['weight']['duality_weight_lb'] * u.lbm:.3f} kg")
    print(f"Tank dry mass: {coupled['weight']['tank_dry_weight_lb'] * u.lbm:.3f} kg")
    print(f"W/S: {coupled['wing_loading_N_m2']:.3f} N/m^2")
    drag_geometry = coupled["drag_geometry"]
    oswald_efficiency = swept_wing_oswald_efficiency(
        drag_geometry["aspect_ratio"],
        drag_geometry["leading_edge_sweep_rad"],
    )
    print(f"AR: {drag_geometry['aspect_ratio']:.3f}")
    print(f"Leading-edge sweep: {np.degrees(drag_geometry['leading_edge_sweep_rad']):.3f} deg")
    print(f"Oswald efficiency e: {oswald_efficiency:.5f}")
    print(f"Governing required T/W: {coupled['required_thrust_to_weight']:.5f}")
    for mode, required_tw in coupled["required_thrust_to_weight_by_mode"].items():
        print(f"  {mode} governing T/W: {required_tw:.5f}")
    propulsion = coupled["propulsion_sizing"]
    print("Propulsion dry sizing:")
    print(f"  Governing case: {propulsion['governing_case']} ({propulsion['governing_mode']})")
    for mode, case in propulsion["governing_by_mode"].items():
        print(
            f"  {mode} electric sizing case: {case['name']}, "
            f"{case['propulsive_power_W'] / 1e6:.3f} MW"
        )
    print(f"  Propulsive power: {propulsion['propulsive_power_W'] / 1e6:.3f} MW")
    print(f"  Motor/controller mass: {propulsion['motor_controller_mass_kg']:.3f} kg")
    print(f"  Generator mass: {propulsion['generator_mass_kg']:.3f} kg")
    print(f"  Turbine mass: {propulsion['turbine_mass_kg']:.3f} kg")
    print(f"  Motor power each: {propulsion['motor_power_each_W'] / 1e6:.3f} MW")
    print(f"  Generator power each: {propulsion['generator_power_each_W'] / 1e6:.3f} MW")
    print(f"  Turbine shaft power: {propulsion['turbine_shaft_power_W'] / 1e6:.3f} MW")
    print(f"Plot: {Path(config.save_plot).resolve()}")
    print(f"Drag-vs-Mach plot: {Path(drag_plot).resolve()}")
    print(f"CD0-vs-Mach plot: {Path(cd0_plot).resolve()}")
    print(f"Level-flight performance plot: {Path(performance_plot).resolve()}")
    print(f"Altitude sweep L/D plot: {Path(altitude_sweep_plot).resolve()}")
    print(f"Drag-terms plot: {Path(drag_terms_plot).resolve()}")
    print(f"Drag-terms 30 kft plot: {Path(drag_terms_30k_plot).resolve()}")
    print(f"K-vs-Mach plot: {Path(k_plot).resolve()}")
    print(f"K-vs-Mach 30 kft plot: {Path(k_30k_plot).resolve()}")
    print("Drag build-up examples at 60,000 ft:")
    for mach_sample in (0.87, 0.95, 1.20, 2.00, 3.00):
        atmosphere = asb.Atmosphere(altitude=60000.0 * u.foot)
        velocity_m_s = mach_sample * atmosphere.speed_of_sound()
        dynamic_pressure_Pa = 0.5 * atmosphere.density() * velocity_m_s**2
        sample_drag = drag_build_up_coefficients(
            config=config,
            drag_geometry=coupled["drag_geometry"],
            altitude_m=60000.0 * u.foot,
            velocity_m_s=velocity_m_s,
            lift_coefficient=coupled["wing_loading_N_m2"] / dynamic_pressure_Pa,
        )
        lift_coefficient = coupled["wing_loading_N_m2"] / dynamic_pressure_Pa
        total_drag_coefficient = (
            sample_drag["zero_lift_drag_coefficient"]
            + sample_drag["lift_dependent_k"] * lift_coefficient**2.0
        )
        print(
            f"  M {mach_sample:.2f}: CD {total_drag_coefficient:.5f}, "
            f"CL {lift_coefficient:.5f}, "
            f"CD0 {sample_drag['zero_lift_drag_coefficient']:.5f}, "
            f"parasite {sample_drag['parasite_cd0']:.5f}, "
            f"wave {sample_drag['wave_cd0']:.5f}"
        )
    print("Design-point margins at coupled W/S:")
    for name, curve in result["curves"].items():
        print(
            f"  {name}: required {curve['required_at_coupled_wing_loading']:.5f}, "
            f"pyCycle available {curve['available_thrust_to_weight']:.5f}, "
            f"margin {curve['margin_at_coupled_wing_loading']:.5f}, "
            f"CD0 {curve['zero_lift_drag_coefficient']:.5f}, "
            f"CD0_parasite {curve['parasite_cd0']:.5f}, "
            f"CD0_wave {curve['wave_cd0']:.5f}"
        )


if __name__ == "__main__":
    main()
