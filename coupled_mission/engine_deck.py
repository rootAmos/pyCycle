"""Engine deck interface for mission-level sizing models."""

from dataclasses import dataclass
import csv
from pathlib import Path


@dataclass
class EngineDeckRecord:
    """One engine operating point, expressed in SI units."""

    mode: str
    mach: float
    altitude_m: float
    throttle: float
    thrust_N: float
    fuel_flow_kg_s: float
    electric_power_W: float = 0.0
    fan1_shaft_power_W: float = 0.0
    fan2_shaft_power_W: float = 0.0
    generator_shaft_power_W: float = 0.0
    inlet_area_m2: float = 0.0
    nozzle_throat_area_m2: float = 0.0


class EngineDeck:
    """
    A small CSV-backed engine deck.

    This first implementation uses nearest-neighbor lookup. That is deliberate:
    it keeps the interface stable while we decide the exact pyCycle grid and
    then replace lookup with smooth interpolation or AeroSandbox surrogates.
    """

    si_required_columns = {
        "mode",
        "mach",
        "altitude_m",
        "throttle",
        "thrust_N",
        "fuel_flow_kg_s",
    }
    imperial_required_columns = {
        "mode",
        "mach",
        "altitude_ft",
        "throttle",
        "thrust_lbf",
        "fuel_flow_lbm_s",
    }

    required_columns = si_required_columns

    def __init__(self, records):
        self.records = list(records)
        if not self.records:
            raise ValueError("EngineDeck needs at least one record.")

    @staticmethod
    def _row_float(row, *names, default=0.0):
        for name in names:
            value = row.get(name)
            if value not in (None, ""):
                return float(value)
        return default

    @classmethod
    def from_csv(cls, csv_path):
        csv_path = Path(csv_path)
        with csv_path.open(newline="") as f:
            reader = csv.DictReader(f)
            fieldnames = set(reader.fieldnames or [])
            has_si = cls.si_required_columns <= fieldnames
            has_imperial = cls.imperial_required_columns <= fieldnames
            if has_si:
                units = "si"
                missing = set()
            elif has_imperial:
                units = "imperial"
                missing = set()
            else:
                missing = min(
                    cls.si_required_columns - fieldnames,
                    cls.imperial_required_columns - fieldnames,
                    key=len,
                )
            if missing:
                raise ValueError(f"{csv_path} is missing required columns: {sorted(missing)}")
            records = []
            for row in reader:
                if units == "imperial":
                    altitude_m = float(row["altitude_ft"]) * 0.3048
                    thrust_N = float(row["thrust_lbf"]) * 4.4482216152605
                    fuel_flow_kg_s = float(row["fuel_flow_lbm_s"]) * 0.45359237
                    electric_power_W = cls._row_float(row, "electric_power_hp") * 745.6998715822702
                    fan1_shaft_power_W = cls._row_float(row, "fan1_shaft_power_hp") * 745.6998715822702
                    fan2_shaft_power_W = cls._row_float(row, "fan2_shaft_power_hp") * 745.6998715822702
                    generator_shaft_power_W = cls._row_float(row, "generator_shaft_power_hp") * 745.6998715822702
                    inlet_area_m2 = cls._row_float(row, "inlet_area_in2") * 0.00064516
                    nozzle_throat_area_m2 = cls._row_float(row, "nozzle_throat_area_in2") * 0.00064516
                else:
                    altitude_m = float(row["altitude_m"])
                    thrust_N = float(row["thrust_N"])
                    fuel_flow_kg_s = float(row["fuel_flow_kg_s"])
                    electric_power_W = cls._row_float(row, "electric_power_W")
                    fan1_shaft_power_W = cls._row_float(row, "fan1_shaft_power_W")
                    fan2_shaft_power_W = cls._row_float(row, "fan2_shaft_power_W")
                    generator_shaft_power_W = cls._row_float(row, "generator_shaft_power_W")
                    inlet_area_m2 = cls._row_float(row, "inlet_area_m2")
                    nozzle_throat_area_m2 = cls._row_float(row, "nozzle_throat_area_m2")

                records.append(
                    EngineDeckRecord(
                        mode=row["mode"],
                        mach=float(row["mach"]),
                        altitude_m=altitude_m,
                        throttle=float(row["throttle"]),
                        thrust_N=thrust_N,
                        fuel_flow_kg_s=fuel_flow_kg_s,
                        electric_power_W=electric_power_W,
                        fan1_shaft_power_W=fan1_shaft_power_W,
                        fan2_shaft_power_W=fan2_shaft_power_W,
                        generator_shaft_power_W=generator_shaft_power_W,
                        inlet_area_m2=inlet_area_m2,
                        nozzle_throat_area_m2=nozzle_throat_area_m2,
                    )
                )
        return cls(records)

    def nearest(self, mach, altitude_m, throttle, mode=None):
        candidates = self.records if mode is None else [r for r in self.records if r.mode == mode]
        if not candidates:
            raise ValueError(f"No engine deck records found for mode={mode!r}.")

        def score(record):
            return (
                (record.mach - mach) ** 2
                + ((record.altitude_m - altitude_m) / 1000.0) ** 2
                + (record.throttle - throttle) ** 2
            )

        return min(candidates, key=score)

    def evaluate(self, mach, altitude_m, throttle, mode=None):
        """Return a dict matching the nearest deck point."""
        record = self.nearest(mach=mach, altitude_m=altitude_m, throttle=throttle, mode=mode)
        return {
            "mode": record.mode,
            "thrust_N": record.thrust_N,
            "fuel_flow_kg_s": record.fuel_flow_kg_s,
            "electric_power_W": record.electric_power_W,
            "fan1_shaft_power_W": record.fan1_shaft_power_W,
            "fan2_shaft_power_W": record.fan2_shaft_power_W,
            "generator_shaft_power_W": record.generator_shaft_power_W,
            "inlet_area_m2": record.inlet_area_m2,
            "nozzle_throat_area_m2": record.nozzle_throat_area_m2,
        }
