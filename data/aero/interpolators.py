"""AeroSandbox-compatible interpolators for digitized aero reference data."""

from dataclasses import dataclass
from pathlib import Path
import csv
import re

import aerosandbox.numpy as np
import numpy as onp


AERO_DATA_ROOT = Path(__file__).resolve().parent


def _read_two_column_csv(path):
    rows = []
    with Path(path).open(newline="") as stream:
        for row in csv.reader(stream):
            if not row:
                continue
            rows.append((float(row[0]), float(row[1])))
    if not rows:
        raise ValueError(f"No data rows found in {path}.")
    data = onp.array(rows, dtype=float)
    order = onp.argsort(data[:, 0])
    return _deduplicate_sorted_xy(data[order, 0], data[order, 1])


def _deduplicate_sorted_xy(x, y):
    unique_x = onp.unique(x)
    if len(unique_x) == len(x):
        return x, y
    unique_y = onp.array([onp.mean(y[x == xi]) for xi in unique_x])
    return unique_x, unique_y


def _common_axis(curves):
    return onp.unique(onp.concatenate([curve[0] for curve in curves]))


def _structured_values_from_curves(curves, common_x):
    return onp.vstack(
        [
            onp.interp(common_x, curve_x, curve_y)
            for curve_x, curve_y in curves
        ]
    )


def _parse_decimal_filename(path):
    return float(Path(path).stem.replace("_", "."))


def _parse_reynolds_filename(path):
    match = re.fullmatch(r"re_10_to_the_([+-]?\d+(?:_\d+)?)", Path(path).stem)
    if match is None:
        raise ValueError(f"Could not parse Reynolds number from {path.name}.")
    exponent = float(match.group(1).replace("_", "."))
    return 10.0**exponent


def _interpn_xi(axis0_value, axis1_value):
    axis0_value, axis1_value = axis0_value + 0 * axis1_value, axis1_value + 0 * axis0_value
    return np.stack((axis0_value, axis1_value), axis=-1)


@dataclass(frozen=True)
class StructuredAeroInterpolator2D:
    """Regular-grid 2D interpolator backed by AeroSandbox's `np.interpn`."""

    axis0_values: object
    axis1_values: object
    values: object
    axis0_name: str
    axis1_name: str
    output_name: str
    method: str = "linear"

    def __call__(self, axis0_value, axis1_value):
        scalar_input = not hasattr(axis0_value, "__len__") and not hasattr(
            axis1_value,
            "__len__",
        )
        axis0_value = np.clip(
            axis0_value,
            float(self.axis0_values[0]),
            float(self.axis0_values[-1]),
        )
        axis1_value = np.clip(
            axis1_value,
            float(self.axis1_values[0]),
            float(self.axis1_values[-1]),
        )
        result = np.interpn(
            points=(self.axis0_values, self.axis1_values),
            values=self.values,
            xi=_interpn_xi(axis0_value, axis1_value),
            method=self.method,
        )
        if scalar_input:
            return result[0]
        return result


def load_leading_edge_suction_interpolator(
    data_dir=AERO_DATA_ROOT / "leading_edge_suction",
    method="linear",
):
    """Return `s = f(cl_design, cl)` from the leading-edge suction CSV curves."""
    csv_paths = sorted(Path(data_dir).glob("*.csv"), key=_parse_decimal_filename)
    if not csv_paths:
        raise FileNotFoundError(f"No CSV files found in {data_dir}.")
    design_cls = onp.array([_parse_decimal_filename(path) for path in csv_paths])
    curves = [_read_two_column_csv(path) for path in csv_paths]
    cl_values = _common_axis(curves)
    values = _structured_values_from_curves(curves, cl_values)
    return StructuredAeroInterpolator2D(
        axis0_values=design_cls,
        axis1_values=cl_values,
        values=values,
        axis0_name="cl_design",
        axis1_name="cl",
        output_name="leading_edge_suction_factor",
        method=method,
    )


def load_cla_theory_ratio_interpolator(
    data_dir=AERO_DATA_ROOT / "cla_theory_ratio",
    method="linear",
):
    """Return `cla/cla_theory = f(log10(Re), tan_half_te_ang_deg)` from CSV curves."""
    csv_paths = sorted(Path(data_dir).glob("*.csv"), key=_parse_reynolds_filename)
    if not csv_paths:
        raise FileNotFoundError(f"No CSV files found in {data_dir}.")
    reynolds_numbers = onp.array([_parse_reynolds_filename(path) for path in csv_paths])
    log10_reynolds = onp.log10(reynolds_numbers)
    curves = [_read_two_column_csv(path) for path in csv_paths]
    tan_half_te_ang_deg_values = _common_axis(curves)
    values = _structured_values_from_curves(curves, tan_half_te_ang_deg_values)
    return StructuredAeroInterpolator2D(
        axis0_values=log10_reynolds,
        axis1_values=tan_half_te_ang_deg_values,
        values=values,
        axis0_name="log10_reynolds_number",
        axis1_name="tan_half_te_ang_deg",
        output_name="cla_cla_theory_ratio",
        method=method,
    )


def leading_edge_suction_factor(cl, cl_design, method="linear"):
    """Evaluate the leading-edge suction factor table."""
    interpolator = load_leading_edge_suction_interpolator(method=method)
    return interpolator(cl_design, cl)


def cla_cla_theory_ratio(tan_half_te_ang_deg, reynolds_number, method="linear"):
    """Evaluate the airfoil lift-curve-slope ratio table."""
    interpolator = load_cla_theory_ratio_interpolator(method=method)
    return interpolator(np.log10(reynolds_number), tan_half_te_ang_deg)
