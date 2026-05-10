"""Build an  constraint diagram with AeroSandbox.

This connects:
- the volume model, which determines `S_plan` from Kuechemann slenderness,
- the Raymer-style weight breakdown, which determines OEW/TOGW,
- the Zhang et al. constraint-analysis equation,
- and pyCycle engine-deck thrust availability.
"""

from dataclasses import dataclass
from pathlib import Path
import csv

import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.tools.units as u

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
    propulsion_mass_kg: object = 450.0
    propulsion_volume_m3: object = 3.0
    payload_volume_m3: object = 5.0
    kuechemann_tau: object = 0.0446
    fuel_density_kg_m3: object = 422.0
    number_engines: object = 2.0
    thrust_scale: object = 1.0
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
            name="Case 4: climb acceleration",
            mode="fan_ab",
            mach=1.2,
            case="climb_acceleration",
            initial_mach=1.0,
            final_mach=3.0,
            acceleration_time_s=30 * 60.0,
            initial_altitude_m=30000 * u.foot,
            final_altitude_m=45000 * u.foot,
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


def required_tw_for_design_point(design_point, wing_loading_N_m2):
    """Return required sea-level static T/W for one literature constraint case."""
    atmosphere = asb.Atmosphere(altitude=representative_altitude_m(design_point))
    density_kg_m3 = atmosphere.density()
    gravity_m_s2 = 9.80665
    cd0 = design_point.cd0
    beta = design_point.beta
    alpha = design_point.alpha
    k1 = design_point.drag_polar_k1
    k2 = design_point.drag_polar_k2

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
    vtail_area_m2 = 0.26 * planform_area_m2
    airplane = build_airplane(
        planform_area_m2=planform_area_m2,
        fuselage_length_m=fuselage_length_m,
        fuselage_height_m=0.12 * fuselage_length_m,
        fuselage_width_m=0.10 * fuselage_length_m,
        vtail_area_m2=vtail_area_m2,
        vtail_span_m=(vtail_area_m2 * 1.4) ** 0.5,
        vtail_dihedral_angle_deg=37.0,
        main_wing_tip_le_x_m=0.25 * root_chord_m,
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
            mass_kg=config.propulsion_mass_kg,
            number_engines=config.number_engines,
            engine_front_to_cockpit_length_m=0.35 * fuselage_length_m,
        ),
        fuel=FuelSystem(
            mass_kg=config.fuel_mass_kg,
            volume_m3=config.fuel_mass_kg / config.fuel_density_kg_m3,
            fuel_density_kg_m3=config.fuel_density_kg_m3,
        ),
        payload=Payload(),
    ).to_weight_inputs()


def solve_coupled_weight_volume(config=ConstraintDiagramConfig()):
    opti = asb.Opti()
    planform_area_m2 = opti.variable(init_guess=80.0, lower_bound=1.0, scale=100.0)
    takeoff_mass_kg = opti.variable(init_guess=4000.0, lower_bound=100.0, scale=5000.0)

    fuselage_length_m = 4.0 * planform_area_m2**0.5
    span_m = (planform_area_m2 * 3.0) ** 0.5
    root_chord_m = 2.0 * planform_area_m2 / (span_m * (1.0 + 0.25))
    vtail_area_m2 = 0.26 * planform_area_m2
    airplane = build_airplane(
        planform_area_m2=planform_area_m2,
        fuselage_length_m=fuselage_length_m,
        fuselage_height_m=0.12 * fuselage_length_m,
        fuselage_width_m=0.10 * fuselage_length_m,
        vtail_area_m2=vtail_area_m2,
        vtail_span_m=(vtail_area_m2 * 1.4) ** 0.5,
        vtail_dihedral_angle_deg=37.0,
        main_wing_tip_le_x_m=0.25 * root_chord_m,
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
            mass_kg=config.propulsion_mass_kg,
            volume_m3=config.propulsion_volume_m3,
            number_engines=config.number_engines,
            engine_front_to_cockpit_length_m=0.35 * fuselage_length_m,
        ),
        fuel=FuelSystem(
            mass_kg=config.fuel_mass_kg,
            volume_m3=config.fuel_mass_kg / config.fuel_density_kg_m3,
            fuel_density_kg_m3=config.fuel_density_kg_m3,
        ),
        payload=Payload(volume_m3=config.payload_volume_m3),
    )
    volume = aircraft_volume_breakdown(aircraft.to_volume_inputs())
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
        "total_aircraft_weight_lb": sol(weight["total_aircraft_weight_lb"]),
        "total_aircraft_mass_kg": sol(weight["total_aircraft_mass_kg"]),
    }
    solved_wing_loading_N_m2 = solved_takeoff_mass_kg * 9.80665 / solved_planform_area_m2

    return {
        "planform_area_m2": solved_planform_area_m2,
        "takeoff_mass_kg": solved_takeoff_mass_kg,
        "takeoff_weight_N": solved_takeoff_mass_kg * 9.80665,
        "wing_loading_N_m2": solved_wing_loading_N_m2,
        "volume": solved_volume,
        "weight": solved_weight,
    }


def build_constraint_diagram(config=ConstraintDiagramConfig(), design_points=None):
    design_points = design_points or default_constraint_design_points()
    coupled = solve_coupled_weight_volume(config)
    engine_deck_rows = read_pycycle_engine_deck(config.engine_deck_csv)

    wing_loading_N_m2 = np.linspace(
        config.wing_loading_min_N_m2,
        config.wing_loading_max_N_m2,
        config.wing_loading_points,
    )

    curves = {}
    for design_point in design_points:
        velocity_m_s, dynamic_pressure_Pa = design_point_flight_condition(design_point)
        required_tw = required_tw_for_design_point(design_point, wing_loading_N_m2)
        available_thrust_N = pycycle_design_point_thrust_N(
            engine_deck_rows,
            design_point,
            thrust_scale=config.thrust_scale,
        )
        available_tw = available_thrust_N / coupled["takeoff_weight_N"]
        required_at_coupled_wing_loading = required_tw_for_design_point(
            design_point,
            coupled["wing_loading_N_m2"],
        )
        curves[design_point.name] = {
            "design_point": design_point,
            "velocity_m_s": velocity_m_s,
            "dynamic_pressure_Pa": dynamic_pressure_Pa,
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
    fig, ax = plt.subplots(figsize=(10.5, 6.5), constrained_layout=True)

    for name, curve in result["curves"].items():
        ax.plot(
            wing_loading_N_m2,
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
        wing_loading_N_m2,
        result["governing_required_thrust_to_weight"],
        color="black",
        linewidth=2.4,
        label="Governing required",
    )
    ax.axvline(
        result["coupled"]["wing_loading_N_m2"],
        color="black",
        linestyle=":",
        linewidth=2.0,
        label="Coupled weight/volume W/S",
    )
    ax.scatter(
        [result["coupled"]["wing_loading_N_m2"]],
        [result["coupled"]["required_thrust_to_weight"]],
        color="black",
        marker="D",
        zorder=5,
    )
    ax.set_xlabel("Wing loading W/S, N/m^2")
    ax.set_ylabel("Design thrust-to-weight T/W")
    ax.set_title(" constraint diagram")
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


def main():
    # Edit run options here.
    config = ConstraintDiagramConfig(
        engine_deck_csv="coupled_mission/data/example_engine_deck.csv",
        fuel_mass_kg=1200.0,
        propulsion_mass_kg=450.0,
        propulsion_volume_m3=3.0,
        payload_volume_m3=5.0,
        thrust_scale=1.0,
        save_plot="_constraint_diagram.png",
        show_plot=False,
    )
    result = build_constraint_diagram(config)
    plot_constraint_diagram(result, save_plot=config.save_plot, show_plot=config.show_plot)

    coupled = result["coupled"]
    print(" coupled constraint solution")
    print(f"S_plan: {coupled['planform_area_m2']:.3f} m^2")
    print(f"TO mass: {coupled['takeoff_mass_kg']:.3f} kg")
    print(f"OEW without engine: {coupled['weight']['operating_empty_without_engine_lb'] * u.lbm:.3f} kg")
    print(f"W/S: {coupled['wing_loading_N_m2']:.3f} N/m^2")
    print(f"Governing required T/W: {coupled['required_thrust_to_weight']:.5f}")
    print(f"Plot: {Path(config.save_plot).resolve()}")
    print("Design-point margins at coupled W/S:")
    for name, curve in result["curves"].items():
        print(
            f"  {name}: required {curve['required_at_coupled_wing_loading']:.5f}, "
            f"pyCycle available {curve['available_thrust_to_weight']:.5f}, "
            f"margin {curve['margin_at_coupled_wing_loading']:.5f}"
        )


if __name__ == "__main__":
    main()
