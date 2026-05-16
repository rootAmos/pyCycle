"""Small smoke-test example for the digitized aero data interpolators."""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from sizing.aero_interpolators import (
    cla_cla_theory_ratio,
    leading_edge_suction_factor,
)


def scalar(value):
    return float(value[0])


def main():
    print("Aero data interpolation examples")
    print(
        "  leading_edge_suction_factor(cl=0.50, cl_design=0.30): "
        f"{scalar(leading_edge_suction_factor(cl=0.50, cl_design=0.30)):.5f}"
    )
    print(
        "  cla_cla_theory_ratio(tan_half_te_ang_deg=0.08, Re=1e7): "
        f"{scalar(cla_cla_theory_ratio(tan_half_te_ang_deg=0.08, reynolds_number=1e7)):.5f}"
    )


if __name__ == "__main__":
    main()
