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
    inlet_area_m2: float = 0.0
    nozzle_throat_area_m2: float = 0.0


class EngineDeck:
    """
    A small CSV-backed engine deck.

    This first implementation uses nearest-neighbor lookup. That is deliberate:
    it keeps the interface stable while we decide the exact pyCycle grid and
    then replace lookup with smooth interpolation or AeroSandbox surrogates.
    """

    required_columns = {
        "mode",
        "mach",
        "altitude_m",
        "throttle",
        "thrust_N",
        "fuel_flow_kg_s",
    }

    def __init__(self, records):
        self.records = list(records)
        if not self.records:
            raise ValueError("EngineDeck needs at least one record.")

    @classmethod
    def from_csv(cls, csv_path):
        csv_path = Path(csv_path)
        with csv_path.open(newline="") as f:
            reader = csv.DictReader(f)
            missing = cls.required_columns - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{csv_path} is missing required columns: {sorted(missing)}")
            records = []
            for row in reader:
                records.append(
                    EngineDeckRecord(
                        mode=row["mode"],
                        mach=float(row["mach"]),
                        altitude_m=float(row["altitude_m"]),
                        throttle=float(row["throttle"]),
                        thrust_N=float(row["thrust_N"]),
                        fuel_flow_kg_s=float(row["fuel_flow_kg_s"]),
                        electric_power_W=float(row.get("electric_power_W") or 0.0),
                        inlet_area_m2=float(row.get("inlet_area_m2") or 0.0),
                        nozzle_throat_area_m2=float(row.get("nozzle_throat_area_m2") or 0.0),
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
            "inlet_area_m2": record.inlet_area_m2,
            "nozzle_throat_area_m2": record.nozzle_throat_area_m2,
        }
