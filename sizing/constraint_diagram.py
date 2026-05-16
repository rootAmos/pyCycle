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

from data.aero.interpolators import (
    cla_cla_theory_ratio as airfoil_cla_theory_ratio,
    leading_edge_suction_factor,
)

try:
    from .aircraft import (
        Aircraft,
        FuelSystem,
        Payload,
        PropulsionSystem,
        build_geometric_asb_airplane as build_airplane,
    )
    from .constraint_equations import design_point_thrust_to_weight_from_wing_loading
    from .volume import aircraft_volume_breakdown
    from .weight import _weight_breakdown
except ImportError:
    from aircraft import (
        Aircraft,
        FuelSystem,
        Payload,
        PropulsionSystem,
        build_geometric_asb_airplane as build_airplane,
    )
    from constraint_equations import design_point_thrust_to_weight_from_wing_loading
    from volume import aircraft_volume_breakdown
    from weight import _weight_breakdown


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
    engine_deck_csv: object = "coupled_mission/data/example_engine_deck.csv"
    fuel_mass_kg: object = 1200.0
    propulsion_volume_m3: object = 3.0
    payload_volume_m3: object = 5.0
    tank_dry_mass_kg: object = 350.0
    kuechemann_tau: object = 0.0446
    void_volume_coefficient: object = 0.05
    fuel_density_kg_m3: object = 422.0
    number_engines: object = 2.0
    number_propulsive_motors: object = 4.0
    number_generators: object = 4.0
    number_turbines: object = 1.0
    electric_propulsor_mach_limit: object = 1.2
    propulsive_efficiency: object = 0.75
    motor_controller_power_density_W_kg: object = 8000.0
    generator_power_density_W_kg: object = 8000.0
    turbine_power_density_W_kg: object = 6000.0
    generator_efficiency: object = 0.96
    turbine_mechanical_efficiency: object = 0.98
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
    save_plot: object = "_constraint_diagram.png"
    show_plot: bool = False


def default_constraint_design_points():
    return (
        ConstraintDesignPoint(
            name="Case 1: constant-altitude/speed cruise",
            mode="fan_ab",
            mach=3.0,
            altitude_m=60000.0 * u.foot,
            case="generic",
            cd0=0.032,
        ),
        ConstraintDesignPoint(
            name="Case 2: constant-speed climb",
            mode="fan_ab",
            mach=3.0,
            altitude_m=60000.0 * u.foot,
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
            mode="fan_ab",
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
            climb_angle_deg=3.0,
            cd0=0.040,
        ),
    )


def read_pycycle_engine_deck(engine_deck_csv):
    rows = []
    with Path(engine_deck_csv).open(newline="") as f:
        for row in csv.DictReader(f):
            rows.append(
                {
                    "mode": row["mode"],
                    "mach": float(row["mach"]),
                    "altitude_m": float(row["altitude_m"]),
                    "throttle": float(row["throttle"]),
                    "thrust_N": float(row["thrust_N"]),
                    "fuel_flow_kg_s": float(row["fuel_flow_kg_s"]),
                    "electric_power_W": float(row.get("electric_power_W") or 0.0),
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


def swept_wing_tip_le_x_m(span_m, root_chord_m, tip_chord_m, sweep_25_rad):
    """Return tip leading-edge x offset from quarter-chord sweep."""
    return (
        0.5 * span_m * np.tan(sweep_25_rad)
        + 0.25 * (root_chord_m - tip_chord_m)
    )


def drag_geometry_from_planform_area(planform_area_m2, config):
    """Return sizing geometry needed for parasite and wave drag build-up."""
    aspect_ratio = 3.0
    taper_ratio = 0.25
    root_thickness_to_chord = 0.06
    span_m = (planform_area_m2 * aspect_ratio) ** 0.5
    root_chord_m = 2.0 * planform_area_m2 / (span_m * (1.0 + taper_ratio))
    tip_chord_m = taper_ratio * root_chord_m
    main_wing_tip_le_x_m = swept_wing_tip_le_x_m(
        span_m,
        root_chord_m,
        tip_chord_m,
        config.main_wing_quarter_chord_sweep_rad,
    )
    leading_edge_sweep_rad = np.arctan(main_wing_tip_le_x_m / (0.5 * span_m))
    half_chord_sweep_rad = np.arctan(
        (
            main_wing_tip_le_x_m
            + 0.5 * tip_chord_m
            - 0.5 * root_chord_m
        )
        / (0.5 * span_m)
    )
    fuselage_length_m = 4.0 * planform_area_m2**0.5
    fuselage_height_m = 0.12 * fuselage_length_m
    fuselage_width_m = 0.10 * fuselage_length_m
    fuselage_radius_a_m = 0.5 * fuselage_width_m
    fuselage_radius_b_m = 0.5 * fuselage_height_m
    fuselage_perimeter_m = np.pi * (
        3.0 * (fuselage_radius_a_m + fuselage_radius_b_m)
        - (
            (3.0 * fuselage_radius_a_m + fuselage_radius_b_m)
            * (fuselage_radius_a_m + 3.0 * fuselage_radius_b_m)
        )
        ** 0.5
    )
    equivalent_fuselage_diameter_m = (fuselage_height_m * fuselage_width_m) ** 0.5
    engine_diameter_m = 2.0 * u.foot
    engine_length_m = config.nacelle_length_to_diameter * engine_diameter_m
    vtail_area_m2 = 0.26 * planform_area_m2

    return {
        "reference_area_m2": planform_area_m2,
        "aspect_ratio": aspect_ratio,
        "root_thickness_to_chord": root_thickness_to_chord,
        "tip_chord_m": tip_chord_m,
        "mean_aerodynamic_chord_m": 2.0
        / 3.0
        * root_chord_m
        * (1.0 + taper_ratio + taper_ratio**2)
        / (1.0 + taper_ratio),
        "leading_edge_sweep_rad": leading_edge_sweep_rad,
        "half_chord_sweep_rad": half_chord_sweep_rad,
        "wing_wetted_area_m2": 2.0
        * planform_area_m2
        * (1.0 + 0.25 * root_thickness_to_chord),
        "tail_wetted_area_m2": 2.0
        * vtail_area_m2
        * (1.0 + 0.25 * root_thickness_to_chord),
        "tail_mean_chord_m": vtail_area_m2 / (vtail_area_m2 * 1.4) ** 0.5,
        "fuselage_length_m": fuselage_length_m,
        "fuselage_wetted_area_m2": fuselage_perimeter_m * fuselage_length_m,
        "fuselage_fineness_ratio": fuselage_length_m / equivalent_fuselage_diameter_m,
        "max_cross_section_area_m2": 0.25
        * np.pi
        * fuselage_width_m
        * fuselage_height_m,
        "nacelle_wetted_area_m2": config.number_engines
        * np.pi
        * engine_diameter_m
        * engine_length_m,
        "nacelle_length_m": engine_length_m,
    }


def air_dynamic_viscosity_kg_m_s(temperature_K):
    """Sutherland-law dynamic viscosity for air."""
    reference_temperature_K = 273.15
    reference_viscosity_kg_m_s = 1.716e-5
    sutherland_temperature_K = 110.4
    return (
        reference_viscosity_kg_m_s
        * (temperature_K / reference_temperature_K) ** 1.5
        * (reference_temperature_K + sutherland_temperature_K)
        / (temperature_K + sutherland_temperature_K)
    )


def turbulent_skin_friction_coefficient(reynolds_number, mach):
    """Raymer-style turbulent flat-plate skin friction coefficient."""
    reynolds_number = np.maximum(reynolds_number, 1.0e5)
    return 0.455 / (
        np.log10(reynolds_number) ** 2.58 * (1.0 + 0.144 * mach**2) ** 0.65
    )


def swept_wing_oswald_efficiency(aspect_ratio, leading_edge_sweep_rad):
    """Raymer swept-wing Oswald efficiency correlation for Lambda_LE > 30 deg."""
    return (
        4.61
        * (1.0 - 0.045 * aspect_ratio**0.68)
        * np.cos(leading_edge_sweep_rad) ** 0.15
        - 3.1
    )


def smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x**2.0 * (3.0 - 2.0 * x)


def airfoil_theory_lift_curve_slope(config, drag_geometry):
    """Return theoretical 2D airfoil lift curve slope in 1/rad."""
    thickness_to_chord = drag_geometry["root_thickness_to_chord"]
    return (
        2.0 * np.pi
        + 4.7
        * thickness_to_chord
        * (1.0 + 0.00375 * config.airfoil_trailing_edge_angle_deg)
    )


def subsonic_finite_wing_lift_curve_slope(
    config,
    drag_geometry,
    mach,
    reynolds_number,
):
    """Return 3D subsonic CL_alpha using the airfoil-ratio data and finite-wing relation."""
    aspect_ratio = drag_geometry["aspect_ratio"]
    beta = np.sqrt(np.maximum(1.0 - mach**2.0, 1.0e-9))
    tan_half_te_ang = np.tan(np.radians(0.5 * config.airfoil_trailing_edge_angle_deg))
    clalpha_theory = airfoil_theory_lift_curve_slope(config, drag_geometry)
    clalpha_ratio = airfoil_cla_theory_ratio(
        tan_half_te_ang_deg=tan_half_te_ang,
        reynolds_number=reynolds_number,
    )
    airfoil_kappa = 1.05 * clalpha_ratio * clalpha_theory / (2.0 * np.pi)
    return (
        2.0
        * np.pi
        * aspect_ratio
        / (
            2.0
            + np.sqrt(
                aspect_ratio**2.0
                * beta**2.0
                / airfoil_kappa**2.0
                * (
                    1.0
                    + np.tan(drag_geometry["half_chord_sweep_rad"]) ** 2.0
                    / beta**2.0
                )
                + 4.0
            )
        )
    )


def supersonic_ackeret_lift_curve_slope(mach):
    """Return Ackeret 2D supersonic lift curve slope in 1/rad."""
    return 4.0 / np.sqrt(np.maximum(mach**2.0 - 1.0, 1.0e-9))


def blended_lift_curve_slope(config, drag_geometry, mach, reynolds_number):
    """Smoothly blend subsonic finite-wing CL_alpha to Ackeret CL_alpha."""
    subsonic_clalpha = subsonic_finite_wing_lift_curve_slope(
        config=config,
        drag_geometry=drag_geometry,
        mach=mach,
        reynolds_number=reynolds_number,
    )
    m_start = 1.0
    m_end = leading_edge_sonic_mach(drag_geometry)
    m_span = np.maximum(m_end - m_start, 1.0e-6)
    transonic_target_clalpha = supersonic_ackeret_lift_curve_slope(m_end)
    supersonic_clalpha = supersonic_ackeret_lift_curve_slope(np.maximum(mach, m_end))
    blend = smoothstep((mach - m_start) / m_span)
    transonic_clalpha = (
        (1.0 - blend) * subsonic_clalpha
        + blend * transonic_target_clalpha
    )
    return np.where(mach < m_end, transonic_clalpha, supersonic_clalpha)


def leading_edge_sonic_mach(drag_geometry):
    """Return Mach where the leading-edge normal component becomes sonic."""
    return 1.0 / np.maximum(np.cos(drag_geometry["leading_edge_sweep_rad"]), 1.0e-9)


def lift_dependent_drag_factor(
    config,
    drag_geometry,
    mach,
    reynolds_number,
    lift_coefficient,
):
    """Return K from leading-edge suction split between K100 and K0."""
    subsonic_clalpha = subsonic_finite_wing_lift_curve_slope(
        config=config,
        drag_geometry=drag_geometry,
        mach=mach,
        reynolds_number=reynolds_number,
    )
    m_start = 1.0
    m_end = leading_edge_sonic_mach(drag_geometry)
    m_span = np.maximum(m_end - m_start, 1.0e-6)
    supersonic_clalpha_at_transition = supersonic_ackeret_lift_curve_slope(m_end)
    supersonic_clalpha = supersonic_ackeret_lift_curve_slope(np.maximum(mach, m_end))
    design_cl = drag_geometry.get("design_lift_coefficient", 0.3)
    suction = leading_edge_suction_factor(
        cl=np.maximum(lift_coefficient, 0.0),
        cl_design=design_cl,
    )
    suction = np.clip(suction, 0.0, 1.0)
    aspect_ratio = drag_geometry["aspect_ratio"]
    k100 = 1.0 / (np.pi * aspect_ratio)
    subsonic_k = suction * k100 + (1.0 - suction) / subsonic_clalpha
    transition_supersonic_k = (
        suction * k100 + (1.0 - suction) / supersonic_clalpha_at_transition
    )
    supersonic_k = suction * k100 + (1.0 - suction) / supersonic_clalpha
    blend = smoothstep((mach - m_start) / m_span)
    transonic_k = (
        (1.0 - blend) * subsonic_k
        + blend * transition_supersonic_k
    )
    lift_dependent_k = np.where(mach < m_end, transonic_k, supersonic_k)
    clalpha = blended_lift_curve_slope(
        config=config,
        drag_geometry=drag_geometry,
        mach=mach,
        reynolds_number=reynolds_number,
    )
    return lift_dependent_k, clalpha, suction


def supersonic_wave_drag_coefficient(config, drag_geometry, mach):
    """Return the Ma >= 1.2 wave drag estimate from the supplied paper."""
    leading_edge_sweep_deg = np.degrees(drag_geometry["leading_edge_sweep_rad"])
    wave_factor = (
        1.0
        - 0.386
        * np.maximum(mach - config.supersonic_wave_drag_start_mach, 0.0) ** 0.57
        * (1.0 - np.pi * leading_edge_sweep_deg / 100.0) ** 2.0
    )
    return (
        1.5
        * np.maximum(wave_factor, 0.0)
        * 9.0
        * np.pi
        / 2.0
        * (drag_geometry["max_cross_section_area_m2"] / drag_geometry["fuselage_length_m"]) ** 2.0
        / drag_geometry["reference_area_m2"]
    )


def transonic_wave_drag_coefficient(config, drag_geometry, mach):
    """Bezier drag-rise estimate between Mcrit and the Ma 1.2 wave-drag model."""
    mdd = config.drag_divergence_mach
    mcrit = mdd - config.critical_mach_offset_from_mdd
    msup = config.supersonic_wave_drag_start_mach
    cdw_mdd = config.drag_divergence_wave_cd
    cdw_msup = supersonic_wave_drag_coefficient(config, drag_geometry, msup)

    t = np.clip((mach - mcrit) / (msup - mcrit), 0.0, 1.0)
    t_mdd = (mdd - mcrit) / (msup - mcrit)

    # Cubic Bezier ordinate. P0 is zero at Mcrit, P2 has the same CDw as P3
    # so the curve reaches the Ma 1.2 value with a flat tangent as in points B/A.
    p0 = 0.0
    p2 = cdw_msup
    p3 = cdw_msup
    p1_denominator = 3.0 * (1.0 - t_mdd) ** 2.0 * t_mdd
    p1_numerator = cdw_mdd - (
        3.0 * (1.0 - t_mdd) * t_mdd**2.0 * p2 + t_mdd**3.0 * p3
    )
    p1 = p1_numerator / p1_denominator
    wave_cd = (
        (1.0 - t) ** 3.0 * p0
        + 3.0 * (1.0 - t) ** 2.0 * t * p1
        + 3.0 * (1.0 - t) * t**2.0 * p2
        + t**3.0 * p3
    )
    return np.maximum(wave_cd, 0.0)


def drag_build_up_coefficients(
    config,
    drag_geometry,
    altitude_m,
    velocity_m_s,
    lift_coefficient=None,
):
    """Return condition-dependent CD0 and lift-dependent K."""
    atmosphere = asb.Atmosphere(altitude=altitude_m)
    density_kg_m3 = atmosphere.density()
    speed_of_sound_m_s = atmosphere.speed_of_sound()
    mach = velocity_m_s / speed_of_sound_m_s
    viscosity_kg_m_s = air_dynamic_viscosity_kg_m_s(atmosphere.temperature())

    reference_area_m2 = drag_geometry["reference_area_m2"]
    wing_re = density_kg_m3 * velocity_m_s * drag_geometry["mean_aerodynamic_chord_m"] / viscosity_kg_m_s
    tail_re = density_kg_m3 * velocity_m_s * drag_geometry["tail_mean_chord_m"] / viscosity_kg_m_s
    fuselage_re = density_kg_m3 * velocity_m_s * drag_geometry["fuselage_length_m"] / viscosity_kg_m_s
    nacelle_re = density_kg_m3 * velocity_m_s * drag_geometry["nacelle_length_m"] / viscosity_kg_m_s

    wing_cd0 = (
        turbulent_skin_friction_coefficient(wing_re, mach)
        * config.wing_form_factor
        * drag_geometry["wing_wetted_area_m2"]
        / reference_area_m2
    )
    tail_cd0 = (
        turbulent_skin_friction_coefficient(tail_re, mach)
        * config.tail_form_factor
        * drag_geometry["tail_wetted_area_m2"]
        / reference_area_m2
    )
    fuselage_form_factor = (
        1.0
        + 60.0 / drag_geometry["fuselage_fineness_ratio"] ** 3.0
        + drag_geometry["fuselage_fineness_ratio"] / 400.0
    )
    fuselage_cd0 = (
        turbulent_skin_friction_coefficient(fuselage_re, mach)
        * fuselage_form_factor
        * drag_geometry["fuselage_wetted_area_m2"]
        / reference_area_m2
    )
    nacelle_cd0 = (
        turbulent_skin_friction_coefficient(nacelle_re, mach)
        * config.nacelle_form_factor
        * drag_geometry["nacelle_wetted_area_m2"]
        / reference_area_m2
    )
    parasite_cd0 = wing_cd0 + tail_cd0 + fuselage_cd0 + nacelle_cd0

    transonic_wave_cd0 = transonic_wave_drag_coefficient(
        config,
        drag_geometry,
        mach,
    )
    supersonic_wave_cd0 = supersonic_wave_drag_coefficient(config, drag_geometry, mach)
    wave_cd0 = np.where(
        mach < config.drag_divergence_mach - config.critical_mach_offset_from_mdd,
        0.0,
        np.where(
            mach < config.supersonic_wave_drag_start_mach,
            transonic_wave_cd0,
            supersonic_wave_cd0,
        ),
    )

    if lift_coefficient is None:
        lift_coefficient = drag_geometry.get("design_lift_coefficient", 0.3)
    lift_dependent_k, lift_curve_slope, leading_edge_suction = lift_dependent_drag_factor(
        config=config,
        drag_geometry=drag_geometry,
        mach=mach,
        reynolds_number=wing_re,
        lift_coefficient=lift_coefficient,
    )
    oswald_efficiency = swept_wing_oswald_efficiency(
        drag_geometry["aspect_ratio"],
        drag_geometry["leading_edge_sweep_rad"],
    )

    return {
        "mach": mach,
        "parasite_cd0": parasite_cd0,
        "wave_cd0": wave_cd0,
        "zero_lift_drag_coefficient": parasite_cd0 + wave_cd0,
        "oswald_efficiency": oswald_efficiency,
        "lift_curve_slope": lift_curve_slope,
        "leading_edge_suction": leading_edge_suction,
        "lift_dependent_k": lift_dependent_k,
    }


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
        propulsion=PropulsionSystem(
            mass_kg=0.0,
            number_engines=config.number_engines,
            engine_front_to_cockpit_length_m=0.35 * fuselage_length_m,
        ),
        fuel=FuelSystem(
            mass_kg=config.fuel_mass_kg,
            volume_m3=config.fuel_mass_kg / config.fuel_density_kg_m3,
            fuel_density_kg_m3=config.fuel_density_kg_m3,
            tank_dry_mass_kg=config.tank_dry_mass_kg,
        ),
        payload=Payload(),
    ).to_weight_inputs()


def propulsion_system_sizing_breakdown(
    takeoff_mass_kg,
    wing_loading_N_m2,
    design_points,
    config,
    drag_geometry=None,
):
    """Return propulsor power split and component masses from a dry sizing point."""
    takeoff_weight_N = takeoff_mass_kg * 9.80665
    sizing_cases = []
    skipped_cases = []
    for design_point in design_points:
        if not design_point.include_in_governing:
            continue
        sizing_mach = propulsion_sizing_mach(design_point)
        if sizing_mach > config.electric_propulsor_mach_limit:
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
            required_thrust_N * velocity_m_s / config.propulsive_efficiency
        )
        sizing_cases.append(
            {
                "name": design_point.name,
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
    propulsive_power_W = governing_case["propulsive_power_W"]
    motor_power_each_W = propulsive_power_W / config.number_propulsive_motors
    generator_electric_power_W = propulsive_power_W / config.generator_efficiency
    generator_power_each_W = generator_electric_power_W / config.number_generators
    turbine_shaft_power_W = (
        generator_electric_power_W / config.turbine_mechanical_efficiency
    )
    motor_controller_mass_kg = (
        config.number_propulsive_motors
        * motor_power_each_W
        / config.motor_controller_power_density_W_kg
    )
    generator_mass_kg = (
        config.number_generators
        * generator_power_each_W
        / config.generator_power_density_W_kg
    )
    turbine_mass_kg = (
        config.number_turbines
        * turbine_shaft_power_W
        / config.turbine_power_density_W_kg
    )
    total_mass_kg = motor_controller_mass_kg + generator_mass_kg + turbine_mass_kg
    return {
        "governing_case": governing_case["name"],
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
    planform_area_m2 = opti.variable(init_guess=80.0, lower_bound=1.0, scale=100.0)
    takeoff_mass_kg = opti.variable(init_guess=4000.0, lower_bound=100.0, scale=5000.0)

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
    aircraft = Aircraft(
        airplane=airplane,
        mass_kg=takeoff_mass_kg,
        landing_mass_kg=0.85 * takeoff_mass_kg,
        propulsion=PropulsionSystem(
            mass_kg=propulsion_mass_kg,
            volume_m3=config.propulsion_volume_m3,
            number_engines=config.number_engines,
            engine_front_to_cockpit_length_m=0.35 * fuselage_length_m,
        ),
        fuel=FuelSystem(
            mass_kg=config.fuel_mass_kg,
            volume_m3=config.fuel_mass_kg / config.fuel_density_kg_m3,
            fuel_density_kg_m3=config.fuel_density_kg_m3,
            tank_dry_mass_kg=config.tank_dry_mass_kg,
        ),
        payload=Payload(volume_m3=config.payload_volume_m3),
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
        "custom_propulsion_weight_lb": sol(weight["custom_propulsion_weight_lb"]),
        "tank_dry_weight_lb": sol(weight["tank_dry_weight_lb"]),
        "total_aircraft_weight_lb": sol(weight["total_aircraft_weight_lb"]),
        "total_aircraft_mass_kg": sol(weight["total_aircraft_mass_kg"]),
    }
    solved_wing_loading_N_m2 = solved_takeoff_mass_kg * 9.80665 / solved_planform_area_m2
    solved_drag_geometry = drag_geometry_from_planform_area(
        solved_planform_area_m2,
        config,
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

    required_stack = np.array(
        [
            curves[design_point.name]["required_thrust_to_weight"]
            for design_point in design_points
            if design_point.include_in_governing
        ]
    )
    governing_required_tw = np.max(required_stack, axis=0)
    coupled["required_thrust_to_weight"] = np.max(
        np.array(
            [
                curves[design_point.name]["required_at_coupled_wing_loading"]
                for design_point in design_points
                if design_point.include_in_governing
            ]
        )
    )

    return {
        "wing_loading_N_m2": wing_loading_N_m2,
        "governing_required_thrust_to_weight": governing_required_tw,
        "curves": curves,
        "coupled": coupled,
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
    fig, ax = plt.subplots(figsize=(10.5, 6.5), constrained_layout=True)

    for name, curve in result["curves"].items():
        ax.plot(
            wing_loading_lb_ft2,
            curve["required_thrust_to_weight"],
            linewidth=1.8,
            label=f"{name} required",
        )
        ax.axhline(
            curve["available_thrust_to_weight"],
            linestyle="--",
            linewidth=1.0,
            alpha=0.5,
            label=f"{name} pyCycle available",
        )

    ax.plot(
        wing_loading_lb_ft2,
        result["governing_required_thrust_to_weight"],
        color="black",
        linewidth=2.4,
        label="Governing required",
    )
    ax.axvline(
        coupled_wing_loading_lb_ft2,
        color="black",
        linestyle=":",
        linewidth=2.0,
        label="Coupled weight/volume W/S",
    )
    ax.scatter(
        [coupled_wing_loading_lb_ft2],
        [result["coupled"]["required_thrust_to_weight"]],
        color="black",
        marker="D",
        zorder=5,
    )
    for reference_name, reference_wing_loading_lb_ft2 in reference_wing_loadings_lb_ft2.items():
        ax.axvline(
            reference_wing_loading_lb_ft2,
            color="0.35",
            linestyle="--",
            linewidth=1.2,
            alpha=0.8,
            label=f"{reference_name} W/S",
        )

    ax.set_xlabel("Wing loading W/S, lb/ft^2")
    ax.set_ylabel("Design thrust-to-weight T/W")
    ax.set_title(" constraint diagram")
    ax.set_xlim(0.0, 150.0)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0.0)
    ax.legend(loc="upper left", fontsize=8, ncols=2)

    if save_plot is not None:
        save_plot = Path(save_plot)
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
    save_plot="_drag_vs_mach.png",
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
    save_plot="_cd0_vs_mach.png",
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
    save_plot="_level_flight_performance_vs_mach.png",
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
    save_plot="_altitude_sweep_lift_to_drag_vs_mach.png",
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
    save_plot="_drag_terms_vs_mach.png",
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
    save_plot="_k_vs_mach.png",
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
        fig.savefig(save_plot, dpi=180)
    if show_plot:
        plt.show()
    else:
        plt.close(fig)
    return fig


def main():
    # Edit run options here.
    config = ConstraintDiagramConfig(
        engine_deck_csv="coupled_mission/data/example_engine_deck.csv",
        fuel_mass_kg=1200.0,
        propulsion_volume_m3=3.0,
        payload_volume_m3=5.0,
        tank_dry_mass_kg=350.0,
        thrust_scale=1.0,
        save_plot="_constraint_diagram.png",
        show_plot=False,
    )
    result = build_constraint_diagram(config)
    plot_constraint_diagram(result, save_plot=config.save_plot, show_plot=config.show_plot)
    drag_plot = "_drag_vs_mach.png"
    plot_drag_build_up_vs_mach(
        result,
        config,
        altitude_m=60000.0 * u.foot,
        save_plot=drag_plot,
        show_plot=config.show_plot,
    )
    cd0_plot = "_cd0_vs_mach.png"
    plot_cd0_vs_mach(
        result,
        config,
        altitude_m=60000.0 * u.foot,
        save_plot=cd0_plot,
        show_plot=config.show_plot,
    )
    performance_plot = "_level_flight_performance_vs_mach.png"
    plot_level_flight_performance_vs_mach(
        result,
        config,
        altitude_m=60000.0 * u.foot,
        save_plot=performance_plot,
        show_plot=config.show_plot,
    )
    altitude_sweep_plot = "_altitude_sweep_lift_to_drag_vs_mach.png"
    plot_altitude_sweep_lift_to_drag_vs_mach(
        result,
        config,
        save_plot=altitude_sweep_plot,
        show_plot=config.show_plot,
    )
    drag_terms_plot = "_drag_terms_vs_mach.png"
    plot_drag_terms_vs_mach(
        result,
        config,
        altitude_m=60000.0 * u.foot,
        save_plot=drag_terms_plot,
        show_plot=config.show_plot,
    )
    drag_terms_30k_plot = "_drag_terms_vs_mach_30kft.png"
    plot_drag_terms_vs_mach(
        result,
        config,
        altitude_m=30000.0 * u.foot,
        save_plot=drag_terms_30k_plot,
        show_plot=config.show_plot,
    )
    k_plot = "_k_vs_mach.png"
    plot_lift_dependent_k_vs_mach(
        result,
        config,
        altitude_m=60000.0 * u.foot,
        save_plot=k_plot,
        show_plot=config.show_plot,
    )
    k_30k_plot = "_k_vs_mach_30kft.png"
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
    print(f"Propulsion mass: {coupled['weight']['custom_propulsion_weight_lb'] * u.lbm:.3f} kg")
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
    propulsion = coupled["propulsion_sizing"]
    print("Propulsion dry sizing:")
    print(f"  Governing case: {propulsion['governing_case']}")
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
