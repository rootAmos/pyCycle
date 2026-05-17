"""Generate a rubberized Duality engine deck from pyCycle baseline points."""

import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.duality_engine_deck import (
    preview_pycycle_engine_deck_setup,
    pycycle_engine_deck_conditions_from_operating_points,
    write_pycycle_engine_deck,
)
from aero import engine_deck_drag_point
from sizing.constraint_diagram import (
    ConstraintDiagramConfig,
    engine_deck_aircraft_sizing,
)
from sizing.electric_machines import build_duality_powertrain_deck
from scripts.scaled_turboshaft_deck import (
    ScaledTurboshaftDeck,
    expand_turboshaft_deck_altitude,
    load_turboshaft_csv,
    turboshaft_output_csv_for_target,
)


OUTPUT_CSV = Path("propulsion/data/duality_engine_deck.csv")
POWERTRAIN_OUTPUT_CSV = Path("propulsion/data/duality_powertrain_deck.csv")
POWERTRAIN_SIZING_SUMMARY_JSON = Path("propulsion/data/duality_powertrain_sizing_summary.json")
POWERTRAIN_SIZING_SUMMARY_CSV = Path("propulsion/data/duality_powertrain_sizing_summary.csv")
MOTOR_EFFICIENCY_MAP_CSV = Path("propulsion/data/duality_motor_efficiency_map.csv")
MOTOR_EFFICIENCY_MAP_PNG = Path("propulsion/data/duality_motor_efficiency_map.png")
MISSION_REF_JSON = Path("sizing/mission_ref.json")
CONSTRAINT_ENGINE_DECK_CSV = Path("propulsion/data/example_engine_deck.csv")
TURBOSHAFT_BASELINE_CSV = Path("propulsion/data/turbine/turboshaft_1120hp.csv")
TURBOSHAFT_INPUT_CSV = Path("propulsion/data/turbine/turboshaft_1120hp_expanded_100kft.csv")
TURBOSHAFT_EXPANDED_MAX_ALTITUDE_FT = 100000.0
SMOKE_TEST_ONLY = False
SMOKE_MAX_CASES = 18
INCLUDE_SEGMENT_MIDPOINTS = True
SIZE_ELECTRIC_MACHINES = True
SCALE_TURBOSHAFT_TO_GENERATOR_LOAD = True
SWEEP_ELECTRIC_POWER = True
ENGINE_POWER_SETTINGS = (0.35, 0.50, 0.70, 0.85, 1.00)


def condition_key(altitude_ft, mach):
    return (float(altitude_ft), float(mach))


def load_mission_ref(path):
    with Path(path).open("r", encoding="utf-8") as stream:
        return json.load(stream)


def add_operating_point(points, seen, altitude_ft, mach, source):
    mach = float(mach)
    if mach <= 0.0:
        return
    altitude_ft = float(altitude_ft)
    key = condition_key(altitude_ft, mach)
    if key in seen:
        return
    seen.add(key)
    points.append(
        {
            "altitude_ft": altitude_ft,
            "mach": mach,
            "source": source,
        }
    )


def mission_operating_points(mission):
    points = []
    seen = set()
    for segment in mission["segments"]:
        name = segment["name"]
        start_altitude_ft = float(segment["start_altitude_ft"])
        end_altitude_ft = float(segment["end_altitude_ft"])
        start_mach = float(segment["start_mach"])
        end_mach = float(segment["end_mach"])
        add_operating_point(points, seen, start_altitude_ft, start_mach, f"{name}:start")
        if INCLUDE_SEGMENT_MIDPOINTS:
            add_operating_point(
                points,
                seen,
                0.5 * (start_altitude_ft + end_altitude_ft),
                0.5 * (start_mach + end_mach),
                f"{name}:mid",
            )
        add_operating_point(points, seen, end_altitude_ft, end_mach, f"{name}:end")
    return points


def precompute_drag_points(aircraft_sizing, deck_conditions):
    drag_points = {}
    for row in deck_conditions:
        drag_points[condition_key(row["altitude_ft"], row["mach"])] = engine_deck_drag_point(
            aircraft_sizing=aircraft_sizing,
            mach=row["mach"],
            altitude_m=row["altitude_m"],
        )
    return drag_points


def propulsion_architecture(aircraft):
    propulsion = aircraft.propulsion
    number_duality_engines = int(propulsion.number_engines)
    number_motors = int(propulsion.number_propulsive_motors)
    number_generators = int(propulsion.number_generators)
    if number_duality_engines <= 0 or number_motors <= 0 or number_generators <= 0:
        raise ValueError("Aircraft propulsion architecture counts must be positive.")
    if number_motors % number_duality_engines:
        raise ValueError("number_propulsive_motors must divide evenly by number_engines.")
    if number_generators % number_motors:
        raise ValueError("number_generators must divide evenly by number_propulsive_motors.")
    return {
        "number_duality_engines": number_duality_engines,
        "motors_per_duality_engine": number_motors // number_duality_engines,
        "generators_per_motor": number_generators // number_motors,
        "number_turbines": int(propulsion.number_turbines),
        "generator_turbine_engine_face_mach": float(propulsion.generator_turbine_engine_face_mach),
        "generator_turbine_subsonic_pressure_recovery": float(propulsion.generator_turbine_subsonic_pressure_recovery),
        "generator_turbine_min_pressure_recovery": float(propulsion.generator_turbine_min_pressure_recovery),
        "generator_turbine_supersonic_recovery_coefficient": float(
            propulsion.generator_turbine_supersonic_recovery_coefficient
        ),
        "generator_turbine_supersonic_recovery_exponent": float(
            propulsion.generator_turbine_supersonic_recovery_exponent
        ),
    }


def turboshaft_mcp_throttle(input_csv):
    df, _ = load_turboshaft_csv(input_csv)
    return float(df["throttle"].max())


def check_turboshaft_design_point(input_csv, mach, altitude_ft, throttle):
    df, _ = load_turboshaft_csv(input_csv)
    mach_min, mach_max = float(df["mach"].min()), float(df["mach"].max())
    altitude_min_ft = float(df["altitude_ft"].min())
    altitude_max_ft = float(df["altitude_ft"].max())
    throttle_min, throttle_max = float(df["throttle"].min()), float(df["throttle"].max())
    if (
        mach_min <= mach <= mach_max
        and altitude_min_ft <= altitude_ft <= altitude_max_ft
        and throttle_min <= throttle <= throttle_max
    ):
        return
    raise ValueError(
        "Generator-turbine design point is outside the turboshaft deck bounds. "
        f"Requested M={mach:.3g}, h={altitude_ft:.0f} ft, throttle={throttle:.3g}; "
        f"deck supports M=[{mach_min}, {mach_max}], "
        f"h=[{altitude_min_ft}, {altitude_max_ft}] ft, "
        f"throttle=[{throttle_min}, {throttle_max}]."
    )


def ensure_turboshaft_input_csv():
    if TURBOSHAFT_INPUT_CSV.exists():
        return TURBOSHAFT_INPUT_CSV
    return expand_turboshaft_deck_altitude(
        TURBOSHAFT_BASELINE_CSV,
        TURBOSHAFT_INPUT_CSV,
        max_altitude_ft=TURBOSHAFT_EXPANDED_MAX_ALTITUDE_FT,
    )


def main():
    mission = load_mission_ref(MISSION_REF_JSON)
    operating_points = mission_operating_points(mission)
    deck_conditions = pycycle_engine_deck_conditions_from_operating_points(
        operating_points=operating_points,
    )
    constraint_config = ConstraintDiagramConfig(
        engine_deck_csv=CONSTRAINT_ENGINE_DECK_CSV,
    )
    aircraft_sizing = engine_deck_aircraft_sizing(config=constraint_config)
    drag_points_by_condition = precompute_drag_points(
        aircraft_sizing=aircraft_sizing,
        deck_conditions=deck_conditions,
    )

    if SMOKE_TEST_ONLY:
        preview = preview_pycycle_engine_deck_setup(
            altitudes_ft=None,
            mach_values=None,
            max_cases=SMOKE_MAX_CASES,
            drag_points_by_condition=drag_points_by_condition,
            operating_points=operating_points,
        )
        print("Duality engine-deck setup smoke test")
        print(f"Mission: {mission['mission_name']}")
        print(f"Candidate rows: {preview['total_candidate_rows']}")
        print(f"Candidate rows shown: {len(preview['rows'])}")
        print(f"Mode counts: {preview['mode_counts']}")
        print(
            "Aircraft sizing: "
            f"W/S={aircraft_sizing['wing_loading_N_m2']:.1f} N/m^2, "
            f"W={aircraft_sizing['takeoff_weight_N'] / 1000.0:.1f} kN, "
            f"T_req,max={aircraft_sizing['sizing_required_thrust_N'] / 1000.0:.1f} kN"
        )
        for mode, thrust_N in aircraft_sizing["sizing_required_thrust_by_mode_N"].items():
            print(f"  {mode} T_req={thrust_N / 1000.0:.1f} kN")
        for row in preview["rows"]:
            stall_note = ", stall-limited" if row["is_stall_limited"] else ""
            print(
                f"  M={row['mach']:.2f}, alt={row['altitude_ft']:.0f} ft -> "
                f"{row['point_name']} ({row['mode']}), "
                f"T_req={row['required_thrust_N'] / 1000.0:.1f} kN, "
                f"CL={row['lift_coefficient']:.3f}, "
                f"CD={row['total_drag_coefficient']:.4f}"
                f"{stall_note}"
            )
        return preview["rows"]

    rows = write_pycycle_engine_deck(
        output_csv=OUTPUT_CSV,
        operating_points=operating_points,
        power_settings=ENGINE_POWER_SETTINGS if SWEEP_ELECTRIC_POWER else (1.0,),
        drag_points_by_condition=drag_points_by_condition,
        sizing_required_thrust_N=aircraft_sizing["sizing_required_thrust_by_mode_N"].get(
            "ramjet",
            aircraft_sizing["sizing_required_thrust_N"],
        ),
    )
    print(f"Wrote {len(rows)} rubberized engine-deck rows to {OUTPUT_CSV.resolve()}")
    if SIZE_ELECTRIC_MACHINES:
        architecture = propulsion_architecture(aircraft_sizing["aircraft"])
        electric_summary = build_duality_powertrain_deck(
            OUTPUT_CSV,
            POWERTRAIN_OUTPUT_CSV,
            sizing_summary_json=POWERTRAIN_SIZING_SUMMARY_JSON,
            sizing_summary_csv=POWERTRAIN_SIZING_SUMMARY_CSV,
            motor_efficiency_map_csv=MOTOR_EFFICIENCY_MAP_CSV,
            motor_efficiency_map_png=MOTOR_EFFICIENCY_MAP_PNG,
            **architecture,
        )
        print("Electric machine sizing from Duality deck:")
        print(f"  Powertrain deck: {POWERTRAIN_OUTPUT_CSV.resolve()}")
        print(f"  Powertrain sizing summary: {POWERTRAIN_SIZING_SUMMARY_JSON.resolve()}")
        print(f"  Motor efficiency map data: {MOTOR_EFFICIENCY_MAP_CSV.resolve()}")
        print(f"  Motor efficiency map plot: {MOTOR_EFFICIENCY_MAP_PNG.resolve()}")
        print(f"  Motors: {electric_summary['number_motors']}")
        print(f"  Generators: {electric_summary['number_generators']}")
        print(f"  Turbines: {electric_summary['number_turbines']}")
        print(f"  Motor rated power: {electric_summary['unit_motor_rated_power_hp']:.1f} hp")
        print(f"  Generator shaft power: {electric_summary['unit_generator_shaft_power_hp']:.1f} hp")
        print(
            "  Max aircraft generator shaft power: "
            f"{electric_summary['aircraft_generator_shaft_power_hp']:.1f} hp"
        )
        if SCALE_TURBOSHAFT_TO_GENERATOR_LOAD:
            turboshaft_input_csv = ensure_turboshaft_input_csv()
            aircraft_design_mach = float(aircraft_sizing["cruise_constraint_mach"])
            design_mach = architecture["generator_turbine_engine_face_mach"]
            design_altitude_ft = float(aircraft_sizing["cruise_constraint_altitude_ft"])
            mcp_throttle = turboshaft_mcp_throttle(turboshaft_input_csv)
            target_mcp_kw = electric_summary["unit_turbine_shaft_power_kw"]
            output_csv = turboshaft_output_csv_for_target(target_mcp_kw)
            try:
                check_turboshaft_design_point(
                    turboshaft_input_csv,
                    design_mach,
                    design_altitude_ft,
                    mcp_throttle,
                )
                deck = ScaledTurboshaftDeck.from_csv(
                    turboshaft_input_csv,
                    target_mcp_kw=target_mcp_kw,
                    design_mach=design_mach,
                    design_altitude_ft=design_altitude_ft,
                    mcp_throttle=mcp_throttle,
                )
            except ValueError as exc:
                print(f"Skipped turboshaft fuel-flow post deck: {exc}")
            else:
                deck.write_scaled_csv(output_csv)
                electric_summary = build_duality_powertrain_deck(
                    OUTPUT_CSV,
                    POWERTRAIN_OUTPUT_CSV,
                    turboshaft_deck=deck,
                    sizing_summary_json=POWERTRAIN_SIZING_SUMMARY_JSON,
                    sizing_summary_csv=POWERTRAIN_SIZING_SUMMARY_CSV,
                    motor_efficiency_map_csv=MOTOR_EFFICIENCY_MAP_CSV,
                    motor_efficiency_map_png=MOTOR_EFFICIENCY_MAP_PNG,
                    **architecture,
                )
                print(
                    "Scaled turboshaft deck per turbine: "
                    f"{target_mcp_kw:.1f} kW at engine-face M={design_mach:.2f} "
                    f"(aircraft M={aircraft_design_mach:.2f}), "
                    f"h={design_altitude_ft:.0f} ft"
                )
                print(f"Scaled turboshaft deck: {output_csv.resolve()}")
                print(f"Updated powertrain deck with turboshaft fuel flow: {POWERTRAIN_OUTPUT_CSV.resolve()}")
                print(f"Powertrain sizing summary: {POWERTRAIN_SIZING_SUMMARY_JSON.resolve()}")
                print(f"Motor efficiency map data: {MOTOR_EFFICIENCY_MAP_CSV.resolve()}")
                print(f"Motor efficiency map plot: {MOTOR_EFFICIENCY_MAP_PNG.resolve()}")
                print(
                    "  Max aircraft generator-turbine inlet area: "
                    f"{electric_summary['aircraft_generator_turbine_inlet_area_in2']:.1f} in^2"
                )


if __name__ == "__main__":
    main()
