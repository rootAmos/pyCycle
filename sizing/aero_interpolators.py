"""Compatibility imports for aero data interpolators."""

from data.aero.interpolators import (
    StructuredAeroInterpolator2D,
    cla_cla_theory_ratio,
    leading_edge_suction_factor,
    load_cla_theory_ratio_interpolator,
    load_leading_edge_suction_interpolator,
)

__all__ = [
    "StructuredAeroInterpolator2D",
    "cla_cla_theory_ratio",
    "leading_edge_suction_factor",
    "load_cla_theory_ratio_interpolator",
    "load_leading_edge_suction_interpolator",
]
