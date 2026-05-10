"""Close  fuel, tank, weight, volume, and mission sizing.

pyCycle and tank analyses are treated as outer-loop truth models. The sizing
loop uses pyCycle engine-deck samples for thrust/fuel-flow interpolation, runs
the 5-point mission, sizes the tank for the resulting fuel demand, and feeds
fuel/tank mass and volume back into the AeroSandbox weight/volume solve.
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
class MissionWaypoint:
    name: str
    altitude_m: object
    mach: object
    true_airspeed_m_s: object
    dynamic_pressure_Pa: object
    mass_fraction_guess: object
    mode: str


@dataclass(frozen=True)
class MissionClosureConfig:
    engine_deck_csv: object = "coupled_mission/data/example_engine_deck.csv"
    propellant: str = "LNG"
    fuel_density_kg_m3: object = 422.0
    initial_fuel_mass_kg: object = 1200.0
    initial_tank_dry_mass_kg: object = 350.0
    propulsion_mass_kg: object = 450.0
    propulsion_volume_m3: object = 3.0
    payload_volume_m3: object = 5.0
    kuechemann_tau: object = 0.0446
    reserve_fraction: object = 0.06
    tank_fill_fraction: object = 0.95
    tank_cyl_length_to_radius: object = 4.0
    number_engines: object = 2.0
    thrust_scale: object = 1.0
    max_iterations: int = 8
    convergence_tol: object = 1e-3
    segment_durations_s: tuple[object, object, object, object] = (
        20.0 * 60.0,
        35.0 * 60.0,
        3000.0 * u.nautical_mile / (2860.0 * u.knot),
        20.0 * 60.0,
    )
    run_tank_model: bool = False


def default_5_point_mission():
    return (
        MissionWaypoint("Takeoff", 0.0 * u.foot, 0.27, 180.0 * u.knot, 110.0 * u.psf, 1.000, "fan"),
        MissionWaypoint("Transonic accel", 40.0e3 * u.foot, 1.2, 688.0 * u.knot, 510.0 * u.psf, 0.940, "fan_ab"),
        MissionWaypoint("Begin cruise", 95.0e3 * u.foot, 5.0, 2860.0 * u.knot, 450.0 * u.psf, 0.850, "ramjet"),
        MissionWaypoint("End cruise", 95.0e3 * u.foot, 5.0, 2860.0 * u.knot, 450.0 * u.psf, 0.450, "ramjet"),
        MissionWaypoint("Landing", 0.0 * u.foot, 0.22, 145.0 * u.knot, 70.0 * u.psf, 0.420, "fan"),
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
                }
            )
    if not rows:
        raise ValueError(f"No pyCycle engine-deck rows found in {engine_deck_csv}.")
    return rows


def interpolate_pycycle_engine_deck(engine_deck_rows, mode, mach, altitude_m, thrust_scale=1.0):
    candidates = [row for row in engine_deck_rows if row["mode"] == mode]
    if not candidates:
        raise ValueError(f'Mode "{mode}" is missing from the pyCycle engine deck.')

    distances = np.array(
        [
            (row["mach"] - mach) ** 2 + ((row["altitude_m"] - altitude_m) / 10000.0) ** 2
            for row in candidates
        ]
    )
    weights = 1.0 / (distances + 1e-6)
    weights = weights / np.sum(weights)
    thrust_per_throttle = np.sum(
        weights
        * np.array(
            [
                thrust_scale * row["thrust_N"] / max(row["throttle"], 1e-9)
                for row in candidates
            ]
        )
    )
    fuel_flow_per_throttle = np.sum(
        weights
        * np.array(
            [
                row["fuel_flow_kg_s"] / max(row["throttle"], 1e-9)
                for row in candidates
            ]
        )
    )
    return {
        "thrust_per_throttle_N": thrust_per_throttle,
        "fuel_flow_per_throttle_kg_s": fuel_flow_per_throttle,
    }


def weight_inputs_from_sizing(
    takeoff_mass_kg,
    planform_area_m2,
    fuel_mass_kg,
    tank_dry_mass_kg,
    propulsion_mass_kg,
    number_engines,
    fuel_density_kg_m3,
):
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
            mass_kg=propulsion_mass_kg,
            number_engines=number_engines,
            engine_front_to_cockpit_length_m=0.35 * fuselage_length_m,
        ),
        fuel=FuelSystem(
            mass_kg=fuel_mass_kg,
            volume_m3=fuel_mass_kg / fuel_density_kg_m3,
            fuel_density_kg_m3=fuel_density_kg_m3,
            tank_dry_mass_kg=tank_dry_mass_kg,
        ),
        payload=Payload(),
    ).to_weight_inputs()


def solve_weight_volume_for_fuel_and_tank(fuel_mass_kg, tank_dry_mass_kg, tank_volume_m3, config):
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
            mass_kg=fuel_mass_kg,
            volume_m3=fuel_mass_kg / config.fuel_density_kg_m3,
            fuel_density_kg_m3=config.fuel_density_kg_m3,
            tank_dry_mass_kg=tank_dry_mass_kg,
        ),
        payload=Payload(volume_m3=config.payload_volume_m3),
    )
    volume = aircraft_volume_breakdown(aircraft.to_volume_inputs())
    opti.subject_to(volume["kuechemann_slenderness_parameter"] == config.kuechemann_tau)
    opti.subject_to(volume["fuel_volume_m3"] <= tank_volume_m3 * config.tank_fill_fraction)

    weight = _weight_breakdown(aircraft.to_weight_inputs())
    opti.subject_to(takeoff_mass_kg == weight["total_aircraft_mass_kg"])
    opti.minimize(takeoff_mass_kg)

    try:
        sol = opti.solve()
    except RuntimeError:
        sol = opti.debug

    return {
        "planform_area_m2": sol(planform_area_m2),
        "takeoff_mass_kg": sol(takeoff_mass_kg),
        "takeoff_weight_N": sol(takeoff_mass_kg) * 9.80665,
        "wing_loading_N_m2": sol(takeoff_mass_kg) * 9.80665 / sol(planform_area_m2),
        "total_volume_m3": sol(volume["total_aircraft_volume_m3"]),
        "oew_without_engine_kg": sol(weight["operating_empty_without_engine_lb"] * u.lbm),
    }


def size_tank_for_fuel(fuel_mass_kg, config):
    usable_volume_m3 = fuel_mass_kg / (config.fuel_density_kg_m3 * config.tank_fill_fraction)
    radius_m = (
        usable_volume_m3
        / (4.0 / 3.0 * np.pi + np.pi * config.tank_cyl_length_to_radius)
    ) ** (1.0 / 3.0)
    length_m = config.tank_cyl_length_to_radius * radius_m

    reference_volume_m3 = config.initial_fuel_mass_kg / (
        config.fuel_density_kg_m3 * config.tank_fill_fraction
    )
    tank_dry_mass_kg = config.initial_tank_dry_mass_kg * (
        usable_volume_m3 / max(reference_volume_m3, 1e-9)
    ) ** (2.0 / 3.0)

    if config.run_tank_model:
        from coupled_mission.tank_deck import TankCase, run_tank_case

        result = run_tank_case(
            TankCase(
                propellant=config.propellant,
                duration_h=sum(config.segment_durations_s) / 3600.0,
                radius_m=radius_m,
                length_m=length_m,
                m_dot_liq_out_kg_s=fuel_mass_kg / sum(config.segment_durations_s),
            )
        )
        tank_dry_mass_kg = result.tank_dry_mass_kg

    return {
        "radius_m": radius_m,
        "length_m": length_m,
        "usable_volume_m3": usable_volume_m3,
        "tank_dry_mass_kg": tank_dry_mass_kg,
    }


def run_5_point_mission_fuel(
    takeoff_mass_kg,
    planform_area_m2,
    engine_deck_rows,
    waypoints,
    config,
):
    takeoff_weight_N = takeoff_mass_kg * 9.80665
    current_weight_N = takeoff_weight_N
    fuel_burn_kg = 0.0
    segment_results = []

    for i in range(len(waypoints) - 1):
        start = waypoints[i]
        end = waypoints[i + 1]
        duration_s = config.segment_durations_s[i]
        mach = 0.5 * (start.mach + end.mach)
        altitude_m = 0.5 * (start.altitude_m + end.altitude_m)
        dynamic_pressure_Pa = 0.5 * (start.dynamic_pressure_Pa + end.dynamic_pressure_Pa)
        velocity_m_s = max(0.5 * (start.true_airspeed_m_s + end.true_airspeed_m_s), 1e-9)
        beta = current_weight_N / takeoff_weight_N
        wing_loading_N_m2 = takeoff_weight_N / planform_area_m2
        climb_rate_m_s = 0.0
        acceleration_m_s2 = 0.0
        if i in (0, 1):
            climb_rate_m_s = 15.0

        required_thrust_to_weight = design_point_thrust_to_weight_from_wing_loading(
            wing_loading_N_m2=wing_loading_N_m2,
            dynamic_pressure_Pa=dynamic_pressure_Pa,
            velocity_m_s=velocity_m_s,
            installed_full_throttle_thrust_lapse=1.0,
            instantaneous_weight_fraction=beta,
            load_factor=1.0,
            drag_polar_k1=0.05,
            drag_polar_k2=0.0,
            zero_lift_drag_coefficient=0.03 if mach < 2.0 else 0.045,
            climb_rate_m_s=climb_rate_m_s,
            acceleration_m_s2=acceleration_m_s2,
        )
        required_thrust_N = required_thrust_to_weight * takeoff_weight_N
        engine = interpolate_pycycle_engine_deck(
            engine_deck_rows,
            mode=end.mode,
            mach=mach,
            altitude_m=altitude_m,
            thrust_scale=config.thrust_scale,
        )
        throttle = required_thrust_N / max(engine["thrust_per_throttle_N"], 1e-9)
        fuel_flow_kg_s = engine["fuel_flow_per_throttle_kg_s"] * np.maximum(throttle, 0.0)
        segment_fuel_kg = fuel_flow_kg_s * duration_s
        current_weight_N = current_weight_N - segment_fuel_kg * 9.80665
        fuel_burn_kg = fuel_burn_kg + segment_fuel_kg
        segment_results.append(
            {
                "segment": f"{start.name} -> {end.name}",
                "mode": end.mode,
                "mach": mach,
                "altitude_m": altitude_m,
                "duration_s": duration_s,
                "required_thrust_N": required_thrust_N,
                "throttle": throttle,
                "fuel_flow_kg_s": fuel_flow_kg_s,
                "fuel_burn_kg": segment_fuel_kg,
            }
        )

    return {
        "fuel_burn_kg": fuel_burn_kg,
        "fuel_required_kg": fuel_burn_kg * (1.0 + config.reserve_fraction),
        "segment_results": segment_results,
    }


def close_mission_sizing(config=MissionClosureConfig(), waypoints=None):
    waypoints = waypoints or default_5_point_mission()
    engine_deck_rows = read_pycycle_engine_deck(config.engine_deck_csv)
    fuel_mass_kg = config.initial_fuel_mass_kg
    tank = size_tank_for_fuel(fuel_mass_kg, config)
    history = []

    for iteration in range(config.max_iterations):
        sizing = solve_weight_volume_for_fuel_and_tank(
            fuel_mass_kg=fuel_mass_kg,
            tank_dry_mass_kg=tank["tank_dry_mass_kg"],
            tank_volume_m3=tank["usable_volume_m3"],
            config=config,
        )
        mission = run_5_point_mission_fuel(
            takeoff_mass_kg=sizing["takeoff_mass_kg"],
            planform_area_m2=sizing["planform_area_m2"],
            engine_deck_rows=engine_deck_rows,
            waypoints=waypoints,
            config=config,
        )
        next_fuel_mass_kg = mission["fuel_required_kg"]
        next_tank = size_tank_for_fuel(next_fuel_mass_kg, config)
        fuel_change = abs(next_fuel_mass_kg - fuel_mass_kg) / max(fuel_mass_kg, 1e-9)
        tank_change = abs(next_tank["tank_dry_mass_kg"] - tank["tank_dry_mass_kg"]) / max(tank["tank_dry_mass_kg"], 1e-9)

        history.append(
            {
                "iteration": iteration,
                "sizing": sizing,
                "mission": mission,
                "tank": next_tank,
                "fuel_change": fuel_change,
                "tank_change": tank_change,
            }
        )

        fuel_mass_kg = next_fuel_mass_kg
        tank = next_tank
        if max(fuel_change, tank_change) < config.convergence_tol:
            break

    return history[-1] | {"history": history}


def main():
    # Edit run options here.
    config = MissionClosureConfig(
        engine_deck_csv="coupled_mission/data/example_engine_deck.csv",
        initial_fuel_mass_kg=1200.0,
        initial_tank_dry_mass_kg=350.0,
        propulsion_mass_kg=450.0,
        propulsion_volume_m3=3.0,
        payload_volume_m3=5.0,
        run_tank_model=True,
    )
    result = close_mission_sizing(config=config)
    sizing = result["sizing"]
    mission = result["mission"]
    tank = result["tank"]

    print(" mission closure")
    print(f"Iterations: {len(result['history'])}")
    print(f"TO mass: {sizing['takeoff_mass_kg']:.3f} kg")
    print(f"S_plan: {sizing['planform_area_m2']:.3f} m^2")
    print(f"W/S: {sizing['wing_loading_N_m2']:.3f} N/m^2")
    print(f"Fuel burn: {mission['fuel_burn_kg']:.3f} kg")
    print(f"Fuel required with reserve: {mission['fuel_required_kg']:.3f} kg")
    print(f"Tank dry mass: {tank['tank_dry_mass_kg']:.3f} kg")
    print(f"Tank radius: {tank['radius_m']:.3f} m")
    print(f"Tank cylindrical length: {tank['length_m']:.3f} m")
    print("Segments:")
    for segment in mission["segment_results"]:
        print(
            f"  {segment['segment']}: throttle {segment['throttle']:.3f}, "
            f"fuel {segment['fuel_burn_kg']:.3f} kg"
        )


if __name__ == "__main__":
    main()
