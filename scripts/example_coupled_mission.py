"""Small analysis-only example using a CSV engine deck and tank initial fuel."""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from coupled_mission.engine_deck import EngineDeck
from coupled_mission.mission_sizing import MissionNode, evaluate_mission, fuel_burn_kg


def main():
    engine = EngineDeck.from_csv("data/propulsion/example_engine_deck.csv")
    nodes = [
        MissionNode(time_s=0.0, mach=0.3, altitude_m=0.0, required_thrust_N=4000.0, throttle=0.8, mode="fan"),
        MissionNode(time_s=600.0, mach=0.8, altitude_m=9000.0, required_thrust_N=3000.0, throttle=0.6, mode="fan_ab"),
        MissionNode(time_s=1800.0, mach=2.5, altitude_m=18000.0, required_thrust_N=5000.0, throttle=0.7, mode="ramjet"),
    ]
    history = evaluate_mission(nodes, engine, initial_fuel_mass_kg=500.0)
    print(f"Fuel burn: {fuel_burn_kg(history):.2f} kg")
    print(f"Final fuel: {history[-1]['fuel_mass_kg']:.2f} kg")
    print(f"Minimum thrust margin: {min(row['thrust_margin_N'] for row in history):.1f} N")


if __name__ == "__main__":
    main()
