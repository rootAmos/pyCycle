"""Run the AeroSandbox coupled mission optimizer."""

import argparse
import importlib.util
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from coupled_mission.aerosandbox_mission import solve_aerosandbox_mission


def _load_airplane(module_path):
    if module_path is None:
        return None
    module_path = Path(module_path).resolve()
    spec = importlib.util.spec_from_file_location(module_path.stem, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "airplane"):
        raise AttributeError(f"{module_path} does not define an `airplane` object.")
    return module.airplane


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--airplane-module", help="Path to a Python file that defines `airplane`, like the AeroSandbox tutorial.")
    parser.add_argument("--engine-deck", default="propulsion/data/example_engine_deck.csv")
    parser.add_argument("--tank-deck", default=None)
    args = parser.parse_args()

    airplane = _load_airplane(args.airplane_module)
    model, sol = solve_aerosandbox_mission(
        engine_deck_csv=args.engine_deck,
        tank_deck_csv=args.tank_deck,
        airplane=airplane,
    )
    v = model.variables
    e = model.expressions
    print("AeroSandbox coupled mission solution")
    print(f"TOGW: {sol(v['mass_kg'][0]):.1f} kg")
    print(f"Wing area: {sol(v['wing_area_m2']):.2f} m^2")
    print(f"Aspect ratio: {sol(v['aspect_ratio']):.2f}")
    print(f"Tank radius: {sol(v['tank_radius_m']):.2f} m")
    print(f"Tank length: {sol(v['tank_length_m']):.2f} m")
    print(f"Initial fuel: {sol(v['fuel_initial_kg']):.1f} kg")
    print(f"Reserve fuel: {sol(e['fuel_remaining_kg']):.1f} kg")
    print(f"Final range: {sol(v['x_e'][-1]) / 1000:.1f} km")
    print(f"Final speed: {sol(v['speed'][-1]):.1f} m/s")
    print("Throttle:", [float(sol(x)) for x in v["throttle"]])


if __name__ == "__main__":
    main()
