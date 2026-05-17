"""Sweep Duality fan modes over static and low-speed engine test conditions."""

from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[1]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from propulsion.duality_engine_deck import write_pycycle_engine_deck


output_csv = Path("propulsion/data/duality_static_engine_sweep.csv")
altitudes_ft = (0.0, 5000.0, 10000.0)
mach_values = (0.001, 0.03, 0.08, 0.15)
unit_motor_rated_power_W = 1.0e6
motors_per_duality_engine = 2
motor_peak_speed_rpm = 2500.0
motor_power_settings_W = (0.25e6, 0.50e6, 0.75e6, 1.00e6)
geometry_cases = (
    {"case_name": "baseline"},
    {
        "case_name": "inlet18_nozz12",
        "mode1_inlet_diameter_in": 18.0,
        "mode1_nozzle_throat_diameter_in": 12.0,
    },
    {
        "case_name": "inlet24_nozz14",
        "mode1_inlet_diameter_in": 24.0,
        "mode1_nozzle_throat_diameter_in": 14.0,
    },
    {
        "case_name": "inlet30_nozz16",
        "mode1_inlet_diameter_in": 30.0,
        "mode1_nozzle_throat_diameter_in": 16.0,
    },
)


def main():
    operating_points = [
        {"altitude_ft": altitude_ft, "mach": mach, "source": "static_engine_runner"}
        for altitude_ft in altitudes_ft
        for mach in mach_values
    ]
    rows = write_pycycle_engine_deck(
        output_csv=output_csv,
        operating_points=operating_points,
        shaft_power_settings_W=tuple(
            motor_power_W * motors_per_duality_engine
            for motor_power_W in motor_power_settings_W
        ),
        max_fan_shaft_power_W=unit_motor_rated_power_W * motors_per_duality_engine,
        max_fan_speed_rpm=motor_peak_speed_rpm,
        geometry_cases=geometry_cases,
    )
    failed = 0
    with output_csv.open("r", encoding="utf-8") as stream:
        failed = sum("run_model_failed" in line for line in stream)
    print(f"Wrote {output_csv.resolve()}")
    print(f"Converged rows returned: {len(rows)}")
    print(f"Failed rows written: {failed}")


if __name__ == "__main__":
    main()
