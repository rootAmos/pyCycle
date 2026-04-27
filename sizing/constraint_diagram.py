"""Build an Astromechanic constraint diagram with AeroSandbox.

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
    from .engine_sizing import design_point_thrust_to_weight_from_wing_loading
    from .volume import AircraftVolumeInputs, aircraft_volume_breakdown
    from .weight import AstromechanicWeightInputs, astromechanic_weight_breakdown
except ImportError:
    from engine_sizing import design_point_thrust_to_weight_from_wing_loading
    from volume import AircraftVolumeInputs, aircraft_volume_breakdown
    from weight import AstromechanicWeightInputs, astromechanic_weight_breakdown


@dataclass(frozen=True)
class ConstraintDesignPoint:
    name: str
    mode: str
    mach: object
    altitude_m: object
    load_factor: object = 1.0
    beta: object = 1.0
    alpha: object = 1.0
    specific_excess_power_m_s: object = 0.0
    drag_polar_k1: object = 0.05
    drag_polar_k2: object = 0.0
    cd0: object = 0.025


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
    save_plot: object = "astromechanic_constraint_diagram.png"
    show_plot: bool = False


def default_constraint_design_points():
    return (
        ConstraintDesignPoint(
            name="Cruise subsonic",
            mode="fan",
            mach=0.8,
            altitude_m=30000.0 * u.foot,
            cd0=0.026,
        ),
        ConstraintDesignPoint(
            name="Cruise supersonic",
            mode="fan_ab",
            mach=1.6,
            altitude_m=45000.0 * u.foot,
            cd0=0.032,
        ),
        ConstraintDesignPoint(
            name="Cruise hypersonic",
            mode="ramjet",
            mach=5.0,
            altitude_m=80000.0 * u.foot,
            cd0=0.045,
        ),
        ConstraintDesignPoint(
            name="Climb 10k subsonic",
            mode="fan",
            mach=0.55,
            altitude_m=10000.0 * u.foot,
            specific_excess_power_m_s=10.0,
            cd0=0.028,
        ),
        ConstraintDesignPoint(
            name="Climb 30k supersonic",
            mode="fan_ab",
            mach=1.3,
            altitude_m=30000.0 * u.foot,
            specific_excess_power_m_s=20.0,
            cd0=0.034,
        ),
        ConstraintDesignPoint(
            name="Climb 50k hypersonic",
            mode="ramjet",
            mach=3.5,
            altitude_m=50000.0 * u.foot,
            specific_excess_power_m_s=30.0,
            cd0=0.043,
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
            + ((item["altitude_m"] - design_point.altitude_m) / 10000.0) ** 2
        ),
    )
    throttle = max(row["throttle"], 1e-9)
    return thrust_scale * row["thrust_N"] / throttle


def design_point_flight_condition(design_point):
    atmosphere = asb.Atmosphere(altitude=design_point.altitude_m)
    velocity = design_point.mach * atmosphere.speed_of_sound()
    dynamic_pressure = 0.5 * atmosphere.density() * velocity**2
    return velocity, dynamic_pressure


def weight_inputs_from_coupled_sizing(takeoff_mass_kg, planform_area_m2, config):
    aspect_ratio = 3.0
    taper_ratio = 0.25
    span_m = (planform_area_m2 * aspect_ratio) ** 0.5
    mean_chord_m = planform_area_m2 / span_m
    fuselage_length_m = 4.0 * planform_area_m2**0.5
    fuselage_depth_m = 0.12 * fuselage_length_m
    fuselage_width_m = 0.10 * fuselage_length_m
    horizontal_tail_area_m2 = 0.16 * planform_area_m2
    vertical_tail_area_m2 = 0.10 * planform_area_m2
    fuel_volume_m3 = config.fuel_mass_kg / config.fuel_density_kg_m3

    return AstromechanicWeightInputs(
        design_gross_weight_lb=takeoff_mass_kg / u.lbm,
        landing_design_gross_weight_lb=0.85 * takeoff_mass_kg / u.lbm,
        ultimate_load_factor=7.5,
        landing_ultimate_load_factor=4.5,
        mach=5.0,
        dynamic_pressure_lb_ft2=500.0,
        wing_area_ft2=planform_area_m2 / u.foot**2,
        aspect_ratio=aspect_ratio,
        taper_ratio=taper_ratio,
        sweep_25_rad=np.radians(60.0),
        root_thickness_to_chord=0.06,
        wing_mounted_control_area_ft2=0.08 * planform_area_m2 / u.foot**2,
        horizontal_tail_area_ft2=horizontal_tail_area_m2 / u.foot**2,
        horizontal_tail_span_ft=0.25 * span_m / u.foot,
        fuselage_width_at_htail_ft=fuselage_width_m / u.foot,
        vertical_tail_area_ft2=vertical_tail_area_m2 / u.foot**2,
        vertical_tail_aspect_ratio=1.4,
        vertical_tail_height_ft=(vertical_tail_area_m2 * 1.4) ** 0.5 / u.foot,
        horizontal_tail_height_ft=0.0,
        tail_length_ft=0.55 * fuselage_length_m / u.foot,
        rudder_area_ft2=0.25 * vertical_tail_area_m2 / u.foot**2,
        fuselage_structural_length_ft=fuselage_length_m / u.foot,
        fuselage_structural_depth_ft=fuselage_depth_m / u.foot,
        fuselage_structural_width_ft=fuselage_width_m / u.foot,
        main_gear_length_in=42.0,
        nose_gear_length_in=30.0,
        number_engines=config.number_engines,
        total_engine_thrust_lb=12000.0,
        thrust_per_engine_lb=12000.0 / config.number_engines,
        engine_diameter_ft=2.0,
        engine_front_to_cockpit_length_ft=0.35 * fuselage_length_m / u.foot,
        total_fuel_volume_gal=fuel_volume_m3 / u.gallon,
        number_mechanical_functions=1.0,
        number_generators=config.number_engines,
        fuel_weight_lb=config.fuel_mass_kg / u.lbm,
        custom_propulsion_weight_lb=config.propulsion_mass_kg / u.lbm,
    )


def solve_coupled_weight_volume(config=ConstraintDiagramConfig()):
    opti = asb.Opti()
    planform_area_m2 = opti.variable(init_guess=80.0, lower_bound=1.0, scale=100.0)
    takeoff_mass_kg = opti.variable(init_guess=4000.0, lower_bound=100.0, scale=5000.0)

    volume = aircraft_volume_breakdown(
        AircraftVolumeInputs(
            planform_area_m2=planform_area_m2,
            fuel_mass_kg=config.fuel_mass_kg,
            propulsion_volume_m3=config.propulsion_volume_m3,
            payload_volume_m3=config.payload_volume_m3,
            fuel_1_density_kg_m3=config.fuel_density_kg_m3,
        )
    )
    opti.subject_to(volume["kuechemann_slenderness_parameter"] == config.kuechemann_tau)

    weight = astromechanic_weight_breakdown(
        weight_inputs_from_coupled_sizing(takeoff_mass_kg, planform_area_m2, config)
    )
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
        velocity, dynamic_pressure = design_point_flight_condition(design_point)
        required_tw = design_point_thrust_to_weight_from_wing_loading(
            wing_loading=wing_loading_N_m2,
            dynamic_pressure=dynamic_pressure,
            velocity=velocity,
            installed_full_throttle_thrust_lapse=design_point.alpha,
            instantaneous_weight_fraction=design_point.beta,
            load_factor=design_point.load_factor,
            drag_polar_k1=design_point.drag_polar_k1,
            drag_polar_k2=design_point.drag_polar_k2,
            zero_lift_drag_coefficient=design_point.cd0,
            specific_excess_power=design_point.specific_excess_power_m_s,
        )
        available_thrust_N = pycycle_design_point_thrust_N(
            engine_deck_rows,
            design_point,
            thrust_scale=config.thrust_scale,
        )
        available_tw = available_thrust_N / coupled["takeoff_weight_N"]
        required_at_coupled_wing_loading = design_point_thrust_to_weight_from_wing_loading(
            wing_loading=coupled["wing_loading_N_m2"],
            dynamic_pressure=dynamic_pressure,
            velocity=velocity,
            installed_full_throttle_thrust_lapse=design_point.alpha,
            instantaneous_weight_fraction=design_point.beta,
            load_factor=design_point.load_factor,
            drag_polar_k1=design_point.drag_polar_k1,
            drag_polar_k2=design_point.drag_polar_k2,
            zero_lift_drag_coefficient=design_point.cd0,
            specific_excess_power=design_point.specific_excess_power_m_s,
        )
        curves[design_point.name] = {
            "design_point": design_point,
            "velocity_m_s": velocity,
            "dynamic_pressure_Pa": dynamic_pressure,
            "required_thrust_to_weight": required_tw,
            "available_thrust_N": available_thrust_N,
            "available_thrust_to_weight": available_tw,
            "required_at_coupled_wing_loading": required_at_coupled_wing_loading,
            "margin_at_coupled_wing_loading": available_tw - required_at_coupled_wing_loading,
        }

    required_stack = np.array(
        [curves[design_point.name]["required_thrust_to_weight"] for design_point in design_points]
    )
    governing_required_tw = np.max(required_stack, axis=0)
    coupled["required_thrust_to_weight"] = np.max(
        np.array(
            [
                curves[design_point.name]["required_at_coupled_wing_loading"]
                for design_point in design_points
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
        color="tab:red",
        linewidth=2.0,
        label="Coupled weight/volume W/S",
    )
    ax.scatter(
        [result["coupled"]["wing_loading_N_m2"]],
        [result["coupled"]["required_thrust_to_weight"]],
        color="tab:red",
        zorder=5,
    )
    ax.set_xlabel("Wing loading W/S, N/m^2")
    ax.set_ylabel("Design thrust-to-weight T/W")
    ax.set_title("Astromechanic constraint diagram")
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
        save_plot="astromechanic_constraint_diagram.png",
        show_plot=False,
    )
    result = build_constraint_diagram(config)
    plot_constraint_diagram(result, save_plot=config.save_plot, show_plot=config.show_plot)

    coupled = result["coupled"]
    print("Astromechanic coupled constraint solution")
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
