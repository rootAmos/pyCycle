"""
scaled_turboshaft_deck.py

Case-B scaling of a baseline turboshaft engine deck to a target design MCP.

Method:
    1. Load a baseline empirical/simulated engine deck:
          Mach, altitude_ft, throttle -> shaft_power_hp, tailpipe_thrust_lbf, fuel_flow_lb_hr

    2. Choose the aircraft sizing condition:
          design_mach, design_altitude_ft, mcp_throttle

    3. Interpolate the baseline deck at that sizing condition:
          P_base_design = P_base(design_mach, design_altitude_ft, mcp_throttle)

    4. Define the scale factor:
          scale_factor = target_mcp_hp / P_base_design

    5. Scale power and fuel flow everywhere:
          P_scaled(M, h, throttle) = scale_factor * P_base(M, h, throttle)
          Wf_scaled(M, h, throttle) = scale_factor * Wf_base(M, h, throttle)

Rationale:
    This treats the baseline deck as an empirical off-design model of a geometrically
    similar engine family. The scale factor is set by the required Maximum Continuous
    Power at the aircraft sizing condition, not by sea-level-static power.

    The method preserves the baseline deck's normalized variation of power and SFC
    with Mach, altitude, and throttle.

Dependencies:
    pandas
    numpy
    scipy

Example:
    deck = ScaledTurboshaftDeck.from_csv(
        "turboshaft_1120hp.csv",
        target_mcp_kw=1500.0,
        design_mach=0.35,
        design_altitude_ft=15000.0,
        mcp_throttle=50.0,
    )

    result = deck.evaluate(mach=0.25, altitude_ft=10000.0, throttle=46.0)
    print(result)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Union

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator


HP_TO_KW = 0.745699872
KW_TO_HP = 1.0 / HP_TO_KW
LB_TO_KG = 0.45359237
FT_TO_M = 0.3048
M_TO_FT = 1.0 / FT_TO_M


COLUMN_NAMES = [
    "mach",
    "altitude_ft",
    "throttle",
    "shaft_power_hp",
    "tailpipe_thrust_lbf",
    "fuel_flow_lb_hr",
]


RUN_SCALER = "turboshaft"  # "turboshaft", "rubberized", or "expand_turboshaft"
USE_CONSTRAINT_DIAGRAM_DESIGN_POINT = False
CONSTRAINT_ENGINE_DECK_CSV = Path("data/propulsion/example_engine_deck.csv")
DESIGN_POINT_NAME = None

TURBOSHAFT_INPUT_CSV = Path("data/propulsion/turbine/turboshaft_1120hp.csv")
TURBOSHAFT_EXPANDED_CSV = Path("data/propulsion/turbine/turboshaft_1120hp_expanded_100kft.csv")
TURBOSHAFT_OUTPUT_DIR = Path("data/propulsion/turbine")
TURBOSHAFT_OUTPUT_STEM = "scaled_turboshaft_deck"
TURBOSHAFT_TARGET_MCP_KW = 1500.0
TURBOSHAFT_TARGET_MCP_HP = None
TURBOSHAFT_DESIGN_MACH = 0.35
TURBOSHAFT_DESIGN_ALTITUDE_FT = 15000.0
TURBOSHAFT_MCP_THROTTLE = 50.0
TURBOSHAFT_EXPANDED_MAX_ALTITUDE_FT = 100000.0
TURBOSHAFT_EXPANDED_ALTITUDE_STEP_FT = 5000.0

RUBBERIZED_INPUT_CSV = Path("data/propulsion/duality_engine_deck.csv")
RUBBERIZED_OUTPUT_CSV = Path("data/propulsion/scaled_duality_engine_deck.csv")
RUBBERIZED_TARGET_MCP_KW = 1500.0
RUBBERIZED_DESIGN_MACH = 0.35
RUBBERIZED_DESIGN_ALTITUDE_FT = 15000.0
RUBBERIZED_DESIGN_MODE = "fan"
RUBBERIZED_MCP_THROTTLE = 1.0
RUBBERIZED_POWER_COLUMN = "generator_shaft_power_hp"


def _read_metadata(csv_path: Union[str, Path]) -> Dict[str, float | str]:
    """
    Parse comment-line metadata of the form:
        # key: value
    from the top of the deck file.
    """
    metadata: Dict[str, float | str] = {}

    with open(csv_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                break

            clean = line.lstrip("#").strip()
            if ":" not in clean:
                continue

            key, value = clean.split(":", 1)
            key = key.strip().lower().replace(" ", "_")
            value = value.strip()

            try:
                metadata[key] = float(value)
            except ValueError:
                metadata[key] = value

    return metadata


def load_turboshaft_csv(csv_path: Union[str, Path]) -> tuple[pd.DataFrame, Dict[str, float | str]]:
    """
    Load the baseline engine deck.

    The uploaded deck has comment metadata lines beginning with '#',
    followed by one human-readable column-header line, followed by numeric CSV rows.
    """
    csv_path = Path(csv_path)
    metadata = _read_metadata(csv_path)

    rows = []
    with csv_path.open("r", encoding="utf-8", errors="replace", newline="") as stream:
        for row in pd.read_csv(
            stream,
            comment="#",
            header=None,
            skip_blank_lines=True,
            skipinitialspace=True,
        ).itertuples(index=False, name=None):
            if len(row) < len(COLUMN_NAMES):
                continue
            try:
                rows.append([float(value) for value in row[: len(COLUMN_NAMES)]])
            except (TypeError, ValueError):
                continue

    df = pd.DataFrame(rows, columns=COLUMN_NAMES)
    if df.empty:
        raise ValueError(f"No numeric turboshaft deck rows found in {csv_path}.")
    df = df.sort_values(["mach", "altitude_ft", "throttle"]).reset_index(drop=True)

    if (df["shaft_power_hp"] <= 0.0).any():
        raise ValueError("Deck contains non-positive shaft power values; SFC would be invalid.")

    if (df["fuel_flow_lb_hr"] <= 0.0).any():
        raise ValueError("Deck contains non-positive fuel-flow values.")

    return df, metadata


def isa_density_kg_m3(altitude_ft: float) -> float:
    altitude_m = float(altitude_ft) * FT_TO_M
    g0 = 9.80665
    gas_constant = 287.05287
    temperature = 288.15
    pressure = 101325.0
    base_altitude_m = 0.0
    layers = ((11000.0, -0.0065), (20000.0, 0.0), (32000.0, 0.0010))

    for top_altitude_m, lapse_rate in layers:
        next_altitude_m = min(altitude_m, top_altitude_m)
        delta_h = next_altitude_m - base_altitude_m
        if delta_h > 0.0:
            if abs(lapse_rate) < 1.0e-12:
                pressure *= np.exp(-g0 * delta_h / (gas_constant * temperature))
            else:
                next_temperature = temperature + lapse_rate * delta_h
                pressure *= (next_temperature / temperature) ** (-g0 / (lapse_rate * gas_constant))
                temperature = next_temperature
        base_altitude_m = next_altitude_m
        if altitude_m <= top_altitude_m:
            break

    return float(pressure / (gas_constant * temperature))


def expand_turboshaft_deck_altitude(
    input_csv: Union[str, Path],
    output_csv: Union[str, Path],
    max_altitude_ft: float = 100000.0,
    altitude_step_ft: float = 5000.0,
    lapse_exponent: float = 1.0,
) -> Path:
    df, metadata = load_turboshaft_csv(input_csv)
    top_altitude_ft = float(df["altitude_ft"].max())
    if max_altitude_ft <= top_altitude_ft:
        raise ValueError("max_altitude_ft must exceed the source deck altitude limit.")

    first_new_altitude_ft = (
        np.floor(top_altitude_ft / altitude_step_ft) * altitude_step_ft
        + altitude_step_ft
    )
    new_altitudes_ft = np.arange(
        first_new_altitude_ft,
        max_altitude_ft + 0.5 * altitude_step_ft,
        altitude_step_ft,
    )
    top_rows = df[df["altitude_ft"] == top_altitude_ft].copy()
    top_density = isa_density_kg_m3(top_altitude_ft)
    extended_rows = []
    for altitude_ft in new_altitudes_ft:
        row_block = top_rows.copy()
        lapse = (isa_density_kg_m3(float(altitude_ft)) / top_density) ** lapse_exponent
        row_block["altitude_ft"] = float(altitude_ft)
        for column in ("shaft_power_hp", "tailpipe_thrust_lbf", "fuel_flow_lb_hr"):
            row_block[column] = row_block[column] * lapse
        extended_rows.append(row_block)

    expanded = pd.concat([df, *extended_rows], ignore_index=True)
    expanded = expanded.sort_values(["mach", "altitude_ft", "throttle"]).reset_index(drop=True)

    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    metadata = dict(metadata)
    metadata.update(
        {
            "altitude_extension_method": "ISA density-lapse extrapolation from highest source altitude",
            "source_altitude_max_ft": top_altitude_ft,
            "expanded_altitude_max_ft": float(max_altitude_ft),
            "altitude_lapse_exponent": float(lapse_exponent),
        }
    )
    with output_csv.open("w", newline="") as stream:
        for key, value in metadata.items():
            stream.write(f"# {key}: {value}\n")
        expanded.to_csv(stream, index=False)
    return output_csv


def similarity_scale_summary(scale_factor: float) -> Dict[str, float | str]:
    """Return first-order same-cycle similarity scalings for interpreting the deck scale."""
    linear_scale = float(np.sqrt(scale_factor))
    return {
        "same_cycle_power_scale": float(scale_factor),
        "corrected_flow_scale": float(scale_factor),
        "approx_linear_scale": linear_scale,
        "approx_reynolds_scale": linear_scale,
        "pressure_ratio_scale": 1.0,
        "efficiency_scale": 1.0,
        "corrected_speed_scale": 1.0,
        "note": "Same-cycle deck scaling preserves normalized SFC; no Reynolds efficiency credit is applied.",
    }


def target_mcp_kw_from_inputs(target_mcp_kw, target_mcp_hp):
    if target_mcp_kw is not None:
        return float(target_mcp_kw)
    return float(target_mcp_hp) * HP_TO_KW


def mcp_filename_token(target_mcp_kw):
    target_mcp_kw = float(target_mcp_kw)
    if abs(target_mcp_kw - round(target_mcp_kw)) < 1.0e-9:
        return f"{int(round(target_mcp_kw))}kw"
    return f"{target_mcp_kw:.3f}".rstrip("0").rstrip(".").replace(".", "p") + "kw"


def turboshaft_output_csv_for_target(target_mcp_kw):
    return TURBOSHAFT_OUTPUT_DIR / f"{TURBOSHAFT_OUTPUT_STEM}_{mcp_filename_token(target_mcp_kw)}.csv"


@dataclass
class DeckBounds:
    mach_min: float
    mach_max: float
    altitude_min_ft: float
    altitude_max_ft: float
    throttle_min: float
    throttle_max: float

    def contains(self, mach: float, altitude_ft: float, throttle: float) -> bool:
        return (
            self.mach_min <= mach <= self.mach_max
            and self.altitude_min_ft <= altitude_ft <= self.altitude_max_ft
            and self.throttle_min <= throttle <= self.throttle_max
        )


class EngineDeckInterpolator:
    """
    Interpolates baseline deck values over Mach, altitude, and throttle.

    Uses RegularGridInterpolator if the input deck is a complete rectangular grid.
    Falls back to LinearNDInterpolator for scattered data.
    """

    def __init__(self, df: pd.DataFrame):
        self.df = df.copy()

        self.mach_grid = np.sort(self.df["mach"].unique())
        self.altitude_grid = np.sort(self.df["altitude_ft"].unique())
        self.throttle_grid = np.sort(self.df["throttle"].unique())

        self.bounds = DeckBounds(
            mach_min=float(self.mach_grid.min()),
            mach_max=float(self.mach_grid.max()),
            altitude_min_ft=float(self.altitude_grid.min()),
            altitude_max_ft=float(self.altitude_grid.max()),
            throttle_min=float(self.throttle_grid.min()),
            throttle_max=float(self.throttle_grid.max()),
        )

        self._use_regular_grid = self._is_complete_regular_grid()

        if self._use_regular_grid:
            self._build_regular_interpolators()
        else:
            self._build_scattered_interpolators()

    def _is_complete_regular_grid(self) -> bool:
        n_expected = (
            len(self.mach_grid)
            * len(self.altitude_grid)
            * len(self.throttle_grid)
        )
        n_actual = len(self.df.drop_duplicates(["mach", "altitude_ft", "throttle"]))
        return n_actual == n_expected

    def _build_regular_interpolators(self) -> None:
        indexed = self.df.set_index(["mach", "altitude_ft", "throttle"]).sort_index()

        shape = (len(self.mach_grid), len(self.altitude_grid), len(self.throttle_grid))

        def values_for(column: str) -> np.ndarray:
            arr = indexed[column].reindex(
                pd.MultiIndex.from_product(
                    [self.mach_grid, self.altitude_grid, self.throttle_grid],
                    names=["mach", "altitude_ft", "throttle"],
                )
            ).to_numpy()
            if np.isnan(arr).any():
                raise ValueError(f"Missing values detected while gridding {column}.")
            return arr.reshape(shape)

        self._power_interp = RegularGridInterpolator(
            (self.mach_grid, self.altitude_grid, self.throttle_grid),
            values_for("shaft_power_hp"),
            bounds_error=True,
        )
        self._fuel_interp = RegularGridInterpolator(
            (self.mach_grid, self.altitude_grid, self.throttle_grid),
            values_for("fuel_flow_lb_hr"),
            bounds_error=True,
        )
        self._tailpipe_interp = RegularGridInterpolator(
            (self.mach_grid, self.altitude_grid, self.throttle_grid),
            values_for("tailpipe_thrust_lbf"),
            bounds_error=True,
        )

    def _build_scattered_interpolators(self) -> None:
        points = self.df[["mach", "altitude_ft", "throttle"]].to_numpy()

        self._power_interp = LinearNDInterpolator(points, self.df["shaft_power_hp"].to_numpy())
        self._fuel_interp = LinearNDInterpolator(points, self.df["fuel_flow_lb_hr"].to_numpy())
        self._tailpipe_interp = LinearNDInterpolator(points, self.df["tailpipe_thrust_lbf"].to_numpy())

    def _eval_one(self, interpolator, mach: float, altitude_ft: float, throttle: float) -> float:
        if not self.bounds.contains(mach, altitude_ft, throttle):
            raise ValueError(
                "Requested point is outside the deck bounds: "
                f"M=[{self.bounds.mach_min}, {self.bounds.mach_max}], "
                f"h=[{self.bounds.altitude_min_ft}, {self.bounds.altitude_max_ft}] ft, "
                f"throttle=[{self.bounds.throttle_min}, {self.bounds.throttle_max}]"
            )

        x = np.array([[mach, altitude_ft, throttle]], dtype=float)
        val = interpolator(x)

        # RegularGridInterpolator returns shape (1,), LinearNDInterpolator may return array([nan])
        val = float(np.asarray(val).reshape(-1)[0])

        if not np.isfinite(val):
            raise ValueError("Interpolation returned NaN. Point may be outside scattered-data hull.")

        return val

    def evaluate_baseline(self, mach: float, altitude_ft: float, throttle: float) -> Dict[str, float]:
        shaft_power_hp = self._eval_one(self._power_interp, mach, altitude_ft, throttle)
        fuel_flow_lb_hr = self._eval_one(self._fuel_interp, mach, altitude_ft, throttle)
        tailpipe_thrust_lbf = self._eval_one(self._tailpipe_interp, mach, altitude_ft, throttle)

        return {
            "shaft_power_hp": shaft_power_hp,
            "shaft_power_kw": shaft_power_hp * HP_TO_KW,
            "fuel_flow_lb_hr": fuel_flow_lb_hr,
            "fuel_flow_kg_hr": fuel_flow_lb_hr * LB_TO_KG,
            "tailpipe_thrust_lbf": tailpipe_thrust_lbf,
            "sfc_lb_hp_hr": fuel_flow_lb_hr / shaft_power_hp,
            "sfc_lb_kw_hr": fuel_flow_lb_hr / (shaft_power_hp * HP_TO_KW),
        }


@dataclass
class ScaledTurboshaftDeck:
    """
    Case-B scaled turboshaft deck.

    Scaling factor:
        S = target_mcp_hp / P_base(design_mach, design_altitude_ft, mcp_throttle)

    Outputs:
        P_scaled  = S * P_base
        Wf_scaled = S * Wf_base

    The resulting SFC is preserved:
        SFC_scaled = Wf_scaled / P_scaled = Wf_base / P_base
    """

    interpolator: EngineDeckInterpolator
    target_mcp_hp: float
    design_mach: float
    design_altitude_ft: float
    mcp_throttle: float
    metadata: Optional[Dict[str, float | str]] = None

    def __post_init__(self) -> None:
        base_design = self.interpolator.evaluate_baseline(
            self.design_mach,
            self.design_altitude_ft,
            self.mcp_throttle,
        )
        self.base_design_power_hp = base_design["shaft_power_hp"]
        self.base_design_fuel_flow_lb_hr = base_design["fuel_flow_lb_hr"]

        if self.base_design_power_hp <= 0.0:
            raise ValueError("Baseline design-point power must be positive.")

        self.scale_factor = self.target_mcp_hp / self.base_design_power_hp
        self.implied_linear_scale = np.sqrt(self.scale_factor)

    @classmethod
    def from_csv(
        cls,
        csv_path: Union[str, Path],
        target_mcp_hp: Optional[float] = None,
        target_mcp_kw: Optional[float] = None,
        design_mach: float = 0.0,
        design_altitude_ft: float = 0.0,
        mcp_throttle: float = 50.0,
    ) -> "ScaledTurboshaftDeck":
        if (target_mcp_hp is None) == (target_mcp_kw is None):
            raise ValueError("Provide exactly one of target_mcp_hp or target_mcp_kw.")

        if target_mcp_kw is not None:
            target_mcp_hp = target_mcp_kw * KW_TO_HP

        assert target_mcp_hp is not None

        df, metadata = load_turboshaft_csv(csv_path)
        interpolator = EngineDeckInterpolator(df)

        return cls(
            interpolator=interpolator,
            target_mcp_hp=float(target_mcp_hp),
            design_mach=float(design_mach),
            design_altitude_ft=float(design_altitude_ft),
            mcp_throttle=float(mcp_throttle),
            metadata=metadata,
        )

    @property
    def target_mcp_kw(self) -> float:
        return self.target_mcp_hp * HP_TO_KW

    def evaluate(self, mach: float, altitude_ft: float, throttle: float) -> Dict[str, float]:
        base = self.interpolator.evaluate_baseline(mach, altitude_ft, throttle)

        shaft_power_hp = self.scale_factor * base["shaft_power_hp"]
        fuel_flow_lb_hr = self.scale_factor * base["fuel_flow_lb_hr"]
        tailpipe_thrust_lbf = self.scale_factor * base["tailpipe_thrust_lbf"]

        return {
            "mach": float(mach),
            "altitude_ft": float(altitude_ft),
            "throttle": float(throttle),

            "shaft_power_hp": shaft_power_hp,
            "shaft_power_kw": shaft_power_hp * HP_TO_KW,

            "fuel_flow_lb_hr": fuel_flow_lb_hr,
            "fuel_flow_kg_hr": fuel_flow_lb_hr * LB_TO_KG,

            "tailpipe_thrust_lbf": tailpipe_thrust_lbf,

            # SFC is preserved by same-factor scaling.
            "sfc_lb_hp_hr": fuel_flow_lb_hr / shaft_power_hp,
            "sfc_lb_kw_hr": fuel_flow_lb_hr / (shaft_power_hp * HP_TO_KW),

            # Useful for checking what the scaling did.
            "baseline_shaft_power_hp": base["shaft_power_hp"],
            "baseline_fuel_flow_lb_hr": base["fuel_flow_lb_hr"],
            "scale_factor": self.scale_factor,
        }

    def available_mcp(self, mach: float, altitude_ft: float) -> Dict[str, float]:
        """
        Convenience method: evaluate max continuous power at this flight condition.
        """
        return self.evaluate(mach, altitude_ft, self.mcp_throttle)

    def required_throttle_for_power(
        self,
        mach: float,
        altitude_ft: float,
        required_power_kw: float,
        throttle_min: Optional[float] = None,
        throttle_max: Optional[float] = None,
        n_grid: int = 200,
    ) -> Dict[str, float]:
        """
        Find the throttle required to produce required_power_kw at a given Mach/altitude.

        This uses a monotonic grid search + interpolation. It is robust enough for
        engine-deck work and avoids assuming the throttle-power curve is perfectly linear.
        """
        bounds = self.interpolator.bounds
        lo = bounds.throttle_min if throttle_min is None else throttle_min
        hi = bounds.throttle_max if throttle_max is None else throttle_max

        throttles = np.linspace(lo, hi, n_grid)
        powers_kw = np.array([
            self.evaluate(mach, altitude_ft, float(t))["shaft_power_kw"]
            for t in throttles
        ])

        # Make sure the requested power is within available range.
        p_min = float(np.nanmin(powers_kw))
        p_max = float(np.nanmax(powers_kw))

        if required_power_kw < p_min or required_power_kw > p_max:
            raise ValueError(
                f"Required power {required_power_kw:.3f} kW is outside available "
                f"range [{p_min:.3f}, {p_max:.3f}] kW at M={mach}, h={altitude_ft} ft."
            )

        # np.interp requires increasing x. If power is decreasing for some bad deck,
        # sort by power before interpolation.
        order = np.argsort(powers_kw)
        throttle_req = float(np.interp(required_power_kw, powers_kw[order], throttles[order]))

        result = self.evaluate(mach, altitude_ft, throttle_req)
        result["required_power_kw"] = required_power_kw
        result["required_throttle"] = throttle_req

        return result

    def summary(self) -> Dict[str, float]:
        summary = {
            "target_mcp_hp": self.target_mcp_hp,
            "target_mcp_kw": self.target_mcp_kw,
            "design_mach": self.design_mach,
            "design_altitude_ft": self.design_altitude_ft,
            "mcp_throttle": self.mcp_throttle,
            "base_design_power_hp": self.base_design_power_hp,
            "base_design_power_kw": self.base_design_power_hp * HP_TO_KW,
            "base_design_fuel_flow_lb_hr": self.base_design_fuel_flow_lb_hr,
            "scale_factor": self.scale_factor,
            "implied_linear_scale": self.implied_linear_scale,
        }
        summary.update(similarity_scale_summary(self.scale_factor))
        return summary

    def scaled_dataframe(self) -> pd.DataFrame:
        """Return the baseline turboshaft deck with scaled power/thrust/fuel columns."""
        df = self.interpolator.df.copy()
        for column in ("shaft_power_hp", "tailpipe_thrust_lbf", "fuel_flow_lb_hr"):
            df[column] = self.scale_factor * df[column]
        return df

    def write_scaled_csv(self, output_csv: Union[str, Path]) -> None:
        """Write a scaled turboshaft deck in the same CSV format expected by `from_csv`."""
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        metadata = dict(self.metadata or {})
        metadata.update(
            {
                "scaling_method": "Case B design-point MCP scaling",
                "target_mcp_kw": self.target_mcp_kw,
                "design_mach": self.design_mach,
                "design_altitude_ft": self.design_altitude_ft,
                "mcp_throttle": self.mcp_throttle,
                "base_design_power_kw": self.base_design_power_hp * HP_TO_KW,
                "scale_factor": self.scale_factor,
            }
        )
        with output_csv.open("w", newline="") as stream:
            for key, value in metadata.items():
                stream.write(f"# {key}: {value}\n")
            self.scaled_dataframe().to_csv(stream, index=False)


RUBBERIZED_REQUIRED_COLUMNS = {
    "mode",
    "mach",
    "altitude_ft",
    "throttle",
    "thrust_lbf",
    "fuel_flow_lbm_s",
}

RUBBERIZED_SCALE_COLUMNS = [
    "thrust_lbf",
    "fuel_flow_lbm_s",
    "electric_power_hp",
    "fan1_shaft_power_hp",
    "fan2_shaft_power_hp",
    "generator_shaft_power_hp",
    "inlet_area_in2",
    "nozzle_throat_area_in2",
]


def load_rubberized_engine_deck_csv(csv_path: Union[str, Path]) -> pd.DataFrame:
    """Load the repo-style imperial pyCycle/rubberized engine deck."""
    csv_path = Path(csv_path)
    df = pd.read_csv(csv_path)
    missing = RUBBERIZED_REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"{csv_path} is missing required columns: {sorted(missing)}")

    for col in df.columns:
        if col != "mode":
            df[col] = pd.to_numeric(df[col], errors="raise")

    return df.sort_values(["mode", "mach", "altitude_ft", "throttle"]).reset_index(drop=True)


def _interpolate_rubberized_column(
    df: pd.DataFrame,
    *,
    mode: str,
    mach: float,
    altitude_ft: float,
    throttle: float,
    column: str,
) -> float:
    """Interpolate a rubberized deck column for one mode, with nearest fallback."""
    if column not in df.columns:
        raise ValueError(f"Column {column!r} is not present in the engine deck.")

    mode_df = df[df["mode"] == mode]
    if mode_df.empty:
        raise ValueError(f"No rows found for mode={mode!r}.")

    points = mode_df[["mach", "altitude_ft", "throttle"]].to_numpy(dtype=float)
    values = mode_df[column].to_numpy(dtype=float)
    query = np.array([[mach, altitude_ft, throttle]], dtype=float)

    if len(mode_df) >= 4 and np.linalg.matrix_rank(points - points[0]) >= 3:
        interpolator = LinearNDInterpolator(points, values)
        value = float(np.asarray(interpolator(query)).reshape(-1)[0])
        if np.isfinite(value):
            return value

    distance = (
        (mode_df["mach"].to_numpy(dtype=float) - mach) ** 2
        + ((mode_df["altitude_ft"].to_numpy(dtype=float) - altitude_ft) / 3000.0) ** 2
        + (mode_df["throttle"].to_numpy(dtype=float) - throttle) ** 2
    )
    return float(values[int(np.argmin(distance))])


def rubberized_design_point_from_constraint_diagram(
    engine_deck_csv: Union[str, Path] = "data/propulsion/example_engine_deck.csv",
    design_point_name: Optional[str] = None,
) -> Dict[str, float | str]:
    """
    Return the current aircraft sizing condition and target MCP from constraint_diagram.py.

    The target MCP is the gas-turbine shaft power required to drive the generators
    at the current propulsion sizing point.
    """
    from sizing.constraint_diagram import (
        ConstraintDiagramConfig,
        build_constraint_diagram,
        default_constraint_design_points,
        propulsion_sizing_mach,
        representative_altitude_m,
    )

    config = ConstraintDiagramConfig(engine_deck_csv=str(engine_deck_csv))
    result = build_constraint_diagram(config)
    propulsion = result["coupled"]["propulsion_sizing"]
    governing_case = design_point_name or propulsion["governing_case"]
    design_points = {point.name: point for point in default_constraint_design_points()}
    if governing_case not in design_points:
        raise ValueError(
            f"Could not find design point {governing_case!r}. "
            f"Available points: {sorted(design_points)}"
        )

    design_point = design_points[governing_case]
    return {
        "target_mcp_kw": float(propulsion["turbine_shaft_power_W"]) / 1000.0,
        "design_mach": float(propulsion_sizing_mach(design_point)),
        "design_altitude_ft": float(representative_altitude_m(design_point)) * M_TO_FT,
        "design_mode": design_point.mode,
        "mcp_throttle": 1.0,
        "design_point_name": governing_case,
    }


def scale_rubberized_engine_deck(
    input_csv: Union[str, Path],
    output_csv: Union[str, Path],
    *,
    target_mcp_kw: float,
    design_mach: float,
    design_altitude_ft: float,
    design_mode: str,
    mcp_throttle: float = 1.0,
    power_column: str = "generator_shaft_power_hp",
    scale_columns: Optional[Iterable[str]] = None,
) -> Dict[str, float | str]:
    """
    Scale a repo-style imperial rubberized deck to hit target MCP at the sizing point.

    This applies Case-B scaling to the requested power column at the design point:
        S = target_mcp_hp / P_base(M_design, h_design, throttle_MCP)

    The same factor is then applied to the deck's power, fuel-flow, thrust, and
    flow-area columns. This intentionally treats the deck as a rubberized,
    geometrically similar engine/propulsor family for conceptual sizing.
    """
    input_csv = Path(input_csv)
    output_csv = Path(output_csv)
    df = load_rubberized_engine_deck_csv(input_csv)
    target_mcp_hp = float(target_mcp_kw) * KW_TO_HP

    base_design_power_hp = _interpolate_rubberized_column(
        df,
        mode=design_mode,
        mach=float(design_mach),
        altitude_ft=float(design_altitude_ft),
        throttle=float(mcp_throttle),
        column=power_column,
    )
    if base_design_power_hp <= 0.0:
        raise ValueError(
            f"Baseline {power_column} at the design point must be positive; "
            f"got {base_design_power_hp} hp."
        )

    scale_factor = target_mcp_hp / base_design_power_hp
    columns = list(scale_columns or RUBBERIZED_SCALE_COLUMNS)
    scaled = df.copy()
    for column in columns:
        if column in scaled.columns:
            scaled[column] = scaled[column] * scale_factor

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    scaled.to_csv(output_csv, index=False)

    return {
        "input_csv": str(input_csv),
        "output_csv": str(output_csv),
        "target_mcp_kw": float(target_mcp_kw),
        "design_mach": float(design_mach),
        "design_altitude_ft": float(design_altitude_ft),
        "design_mode": design_mode,
        "mcp_throttle": float(mcp_throttle),
        "power_column": power_column,
        "base_design_power_hp": float(base_design_power_hp),
        "base_design_power_kw": float(base_design_power_hp) * HP_TO_KW,
        "scale_factor": float(scale_factor),
        "implied_linear_scale": float(np.sqrt(scale_factor)),
        **similarity_scale_summary(scale_factor),
    }


def _print_summary(title, summary):
    print(title)
    for key, value in summary.items():
        print(f"  {key}: {value}")


def _constraint_diagram_design():
    design = rubberized_design_point_from_constraint_diagram(
        engine_deck_csv=CONSTRAINT_ENGINE_DECK_CSV,
        design_point_name=DESIGN_POINT_NAME,
    )
    print(f"Using constraint-diagram design point: {design['design_point_name']}")
    return design


def run_turboshaft_scaler():
    if USE_CONSTRAINT_DIAGRAM_DESIGN_POINT:
        design = _constraint_diagram_design()
        target_mcp_kw = design["target_mcp_kw"]
        target_mcp_hp = None
        design_mach = design["design_mach"]
        design_altitude_ft = design["design_altitude_ft"]
    else:
        if (TURBOSHAFT_TARGET_MCP_KW is None) == (TURBOSHAFT_TARGET_MCP_HP is None):
            raise ValueError("Set exactly one of TURBOSHAFT_TARGET_MCP_KW or TURBOSHAFT_TARGET_MCP_HP.")
        target_mcp_kw = TURBOSHAFT_TARGET_MCP_KW
        target_mcp_hp = TURBOSHAFT_TARGET_MCP_HP
        design_mach = TURBOSHAFT_DESIGN_MACH
        design_altitude_ft = TURBOSHAFT_DESIGN_ALTITUDE_FT

    output_csv = turboshaft_output_csv_for_target(
        target_mcp_kw_from_inputs(target_mcp_kw, target_mcp_hp)
    )
    deck = ScaledTurboshaftDeck.from_csv(
        TURBOSHAFT_INPUT_CSV,
        target_mcp_kw=target_mcp_kw,
        target_mcp_hp=target_mcp_hp,
        design_mach=design_mach,
        design_altitude_ft=design_altitude_ft,
        mcp_throttle=TURBOSHAFT_MCP_THROTTLE,
    )
    deck.write_scaled_csv(output_csv)
    summary = deck.summary()
    summary["output_csv"] = str(output_csv)
    _print_summary("Scaled turboshaft deck:", summary)


def run_turboshaft_expansion():
    output_csv = expand_turboshaft_deck_altitude(
        TURBOSHAFT_INPUT_CSV,
        TURBOSHAFT_EXPANDED_CSV,
        max_altitude_ft=TURBOSHAFT_EXPANDED_MAX_ALTITUDE_FT,
        altitude_step_ft=TURBOSHAFT_EXPANDED_ALTITUDE_STEP_FT,
    )
    df, _ = load_turboshaft_csv(output_csv)
    summary = {
        "output_csv": str(output_csv),
        "mach_min": float(df["mach"].min()),
        "mach_max": float(df["mach"].max()),
        "altitude_min_ft": float(df["altitude_ft"].min()),
        "altitude_max_ft": float(df["altitude_ft"].max()),
        "throttle_min": float(df["throttle"].min()),
        "throttle_max": float(df["throttle"].max()),
        "rows": int(len(df)),
    }
    _print_summary("Expanded turboshaft deck:", summary)


def run_rubberized_scaler():
    if USE_CONSTRAINT_DIAGRAM_DESIGN_POINT:
        design = _constraint_diagram_design()
        target_mcp_kw = design["target_mcp_kw"]
        design_mach = design["design_mach"]
        design_altitude_ft = design["design_altitude_ft"]
        design_mode = design["design_mode"]
        mcp_throttle = design["mcp_throttle"]
    else:
        target_mcp_kw = RUBBERIZED_TARGET_MCP_KW
        design_mach = RUBBERIZED_DESIGN_MACH
        design_altitude_ft = RUBBERIZED_DESIGN_ALTITUDE_FT
        design_mode = RUBBERIZED_DESIGN_MODE
        mcp_throttle = RUBBERIZED_MCP_THROTTLE

    summary = scale_rubberized_engine_deck(
        RUBBERIZED_INPUT_CSV,
        RUBBERIZED_OUTPUT_CSV,
        target_mcp_kw=target_mcp_kw,
        design_mach=design_mach,
        design_altitude_ft=design_altitude_ft,
        design_mode=design_mode,
        mcp_throttle=mcp_throttle,
        power_column=RUBBERIZED_POWER_COLUMN,
    )
    _print_summary("Scaled rubberized engine deck:", summary)


def main():
    if RUN_SCALER == "turboshaft":
        run_turboshaft_scaler()
    elif RUN_SCALER == "expand_turboshaft":
        run_turboshaft_expansion()
    elif RUN_SCALER == "rubberized":
        run_rubberized_scaler()
    else:
        raise ValueError(
            f"Unknown RUN_SCALER={RUN_SCALER!r}; use 'turboshaft', "
            "'expand_turboshaft', or 'rubberized'."
        )


if __name__ == "__main__":
    main()
