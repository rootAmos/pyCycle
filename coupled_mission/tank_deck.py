"""HyTank/LNGTank execution helpers for mission-level tank decks."""

from dataclasses import dataclass, asdict
import csv
from pathlib import Path

import numpy as np

from .paths import ensure_hytank_on_path


@dataclass
class TankCase:
    """Inputs for one tank transient case."""

    propellant: str = "LNG"
    num_nodes: int = 51
    duration_h: float = 1.0
    radius_m: float = 1.0
    length_m: float = 1.0
    fill_level_init: float = 0.95
    ullage_T_init_K: float | None = None
    liquid_T_init_K: float | None = None
    ullage_P_init_Pa: float = 1.5e5
    T_env_K: float = 300.0
    N_layers: float = 20.0
    m_dot_liq_out_kg_s: float = 0.0
    m_dot_gas_out_kg_s: float = 0.0
    P_heater_W: float = 0.0
    environment_design_pressure_Pa: float = 1.5e5
    max_expected_operating_pressure_Pa: float = 1.0e6
    vacuum_gap_m: float = 0.05


@dataclass
class TankResult:
    """Selected outputs from one tank transient case."""

    propellant: str
    duration_h: float
    radius_m: float
    length_m: float
    m_dot_liq_out_kg_s: float
    m_dot_gas_out_kg_s: float
    tank_dry_mass_kg: float
    initial_total_mass_kg: float
    final_total_mass_kg: float
    initial_fuel_mass_kg: float
    final_fuel_mass_kg: float
    final_fill_level: float
    final_pressure_Pa: float
    final_liquid_temperature_K: float
    final_gas_temperature_K: float
    min_pressure_Pa: float
    max_pressure_Pa: float


def _tank_class(propellant):
    ensure_hytank_on_path()
    key = propellant.strip().upper()
    if key == "LH2":
        from hytank import LH2Tank

        return LH2Tank, 21.0, 20.0
    if key == "LNG":
        from lngtank import LNGTank

        return LNGTank, 120.0, 111.7
    raise ValueError(f'Unsupported propellant "{propellant}". Use "LH2" or "LNG".')


def run_tank_case(case):
    """
    Run one HyTank/LNGTank OpenMDAO case and return mission-relevant outputs.

    This uses constant extraction/heater/environment profiles across the case.
    For a full mission, generate one case per segment or extend this helper to
    accept vectors directly.
    """
    import openmdao.api as om

    tank_cls, default_ullage_T, default_liquid_T = _tank_class(case.propellant)
    ullage_T = default_ullage_T if case.ullage_T_init_K is None else case.ullage_T_init_K
    liquid_T = default_liquid_T if case.liquid_T_init_K is None else case.liquid_T_init_K

    nn = case.num_nodes
    p = om.Problem()
    p.model = tank_cls(
        num_nodes=nn,
        fill_level_init=case.fill_level_init,
        ullage_T_init=ullage_T,
        liquid_T_init=liquid_T,
        ullage_P_init=case.ullage_P_init_Pa,
    )
    p.model.nonlinear_solver = om.NewtonSolver(solve_subsystems=True, atol=1e-8, rtol=1e-8, iprint=-1)
    p.model.linear_solver = om.DirectSolver()
    p.setup()
    p.set_solver_print(level=-1)

    p.set_val("thermals.boil_off.integ.duration", case.duration_h, units="h")
    p.set_val("radius", case.radius_m, units="m")
    p.set_val("length", case.length_m, units="m")
    p.set_val("P_heater", np.full(nn, case.P_heater_W), units="W")
    p.set_val("m_dot_gas_out", np.full(nn, case.m_dot_gas_out_kg_s), units="kg/s")
    p.set_val("m_dot_liq_out", np.full(nn, case.m_dot_liq_out_kg_s), units="kg/s")
    p.set_val("T_env", np.full(nn, case.T_env_K), units="K")
    p.set_val("N_layers", case.N_layers)
    p.set_val("environment_design_pressure", case.environment_design_pressure_Pa, units="Pa")
    p.set_val("max_expected_operating_pressure", case.max_expected_operating_pressure_Pa, units="Pa")
    p.set_val("vacuum_gap", case.vacuum_gap_m, units="m")

    p.run_model()

    m_gas = p.get_val("m_gas", units="kg")
    m_liq = p.get_val("m_liq", units="kg")
    pressure = p.get_val("P", units="Pa")
    return TankResult(
        propellant=case.propellant.upper(),
        duration_h=case.duration_h,
        radius_m=case.radius_m,
        length_m=case.length_m,
        m_dot_liq_out_kg_s=case.m_dot_liq_out_kg_s,
        m_dot_gas_out_kg_s=case.m_dot_gas_out_kg_s,
        tank_dry_mass_kg=float(p.get_val("tank_weight", units="kg").item()),
        initial_total_mass_kg=float(p.get_val("total_weight", units="kg")[0]),
        final_total_mass_kg=float(p.get_val("total_weight", units="kg")[-1]),
        initial_fuel_mass_kg=float(m_gas[0] + m_liq[0]),
        final_fuel_mass_kg=float(m_gas[-1] + m_liq[-1]),
        final_fill_level=float(p.get_val("fill_level")[-1]),
        final_pressure_Pa=float(pressure[-1]),
        final_liquid_temperature_K=float(p.get_val("T_liq", units="K")[-1]),
        final_gas_temperature_K=float(p.get_val("T_gas", units="K")[-1]),
        min_pressure_Pa=float(np.min(pressure)),
        max_pressure_Pa=float(np.max(pressure)),
    )


def write_tank_deck(cases, output_csv):
    """Run cases and write one CSV row per tank case."""
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    results = [run_tank_case(case) for case in cases]
    fieldnames = list(asdict(results[0]).keys()) if results else list(TankResult.__dataclass_fields__)
    with output_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(asdict(result))
    return results
