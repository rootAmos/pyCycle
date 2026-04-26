"""Generate a small LNG/LH2 tank deck from the sibling HyTank checkout."""

from pathlib import Path
import argparse
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from coupled_mission import TankCase, write_tank_deck


def build_cases(propellant):
    cases = []
    for radius_m in [0.75, 1.0, 1.25]:
        for length_m in [1.0, 2.0, 3.0]:
            for m_dot in [0.0, 0.02, 0.05]:
                cases.append(
                    TankCase(
                        propellant=propellant,
                        num_nodes=31,
                        duration_h=1.0,
                        radius_m=radius_m,
                        length_m=length_m,
                        m_dot_liq_out_kg_s=m_dot,
                        T_env_K=300.0,
                    )
                )
    return cases


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--propellant", choices=["LNG", "LH2"], default="LNG")
    parser.add_argument("--out", default="coupled_mission/data/tank_deck.csv")
    args = parser.parse_args()

    output = Path(args.out)
    results = write_tank_deck(build_cases(args.propellant), output)
    print(f"Wrote {len(results)} tank cases -> {output}")


if __name__ == "__main__":
    main()
