"""AeroSandbox mission optimization for Duality engine + cryogenic tank sizing.

This module intentionally imports AeroSandbox inside functions so the rest of
the repository remains usable in environments where AeroSandbox is not installed.
Install with:

    pip install aerosandbox

All variables in this module are SI, matching AeroSandbox convention.
"""

from dataclasses import dataclass
from pathlib import Path
import csv
import importlib.util

try:
    from .tank_deck import TankCase, write_tank_deck
except ImportError:  # Allows `python coupled_mission/aerosandbox_mission.py`.
    import sys

    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    from coupled_mission.tank_deck import TankCase, write_tank_deck


@dataclass
class AeroMissionConfig:
    """Fixed mission schedule and sizing assumptions."""

    waypoint_altitude_m: tuple[float, ...] = (0.0, 12192.0, 28956.0, 28956.0, 0.0)
    waypoint_speed_m_s: tuple[float, ...] = (
        92.5992,
        353.938,
        1471.31,
        1471.31,
        74.5944,
    )
    waypoint_mass_fraction: tuple[float, ...] | None = None
    segment_duration_s: tuple[float, ...] = (600.0, 1800.0, 3600.0, 1200.0)
    mode: tuple[str, ...] = ("fan", "fan_ab", "ramjet", "ramjet", "fan")
    payload_mass_kg: float = 250.0
    fixed_empty_mass_kg: float = 1200.0
    reserve_fuel_kg: float = 50.0
    oswald_efficiency: float = 0.78
    cd0: float = 0.028
    max_cl: float = 0.75
    load_factor: float = 1.0
    fuel_density_kg_m3: float = 422.0
    initial_fill_fraction: float = 0.95
    wing_loading_guess_N_m2: float = 2600.0
    initial_range_m: float = 0.0
    target_range_m: float = 5_556_000.0
    max_altitude_m: float = 32000.0
    min_speed_m_s: float = 70.0
    max_speed_m_s: float = 1700.0
    min_gamma_rad: float = -0.35
    max_gamma_rad: float = 0.35
    min_alpha_deg: float = -5.0
    max_alpha_deg: float = 15.0


@dataclass
class AeroMissionModel:
    """Container returned by `build_aerosandbox_mission()`."""

    opti: object
    variables: dict
    expressions: dict
    config: AeroMissionConfig


def _read_engine_reference_by_mode(engine_deck_csv):
    """Read representative full-scale points from an engine deck CSV."""
    rows = []
    with Path(engine_deck_csv).open(newline="") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    if not rows:
        raise ValueError(f"No rows found in engine deck: {engine_deck_csv}")

    by_mode = {}
    for row in rows:
        mode = row["mode"]
        throttle = max(float(row["throttle"]), 1e-6)
        thrust = float(row["thrust_N"])
        fuel_flow = float(row["fuel_flow_kg_s"])
        electric_power = float(row.get("electric_power_W") or 0.0)
        normalized = {
            "mach_ref": float(row["mach"]),
            "altitude_ref_m": float(row["altitude_m"]),
            "thrust_per_throttle_N": thrust / throttle,
            "fuel_flow_per_throttle_kg_s": fuel_flow / throttle,
            "electric_power_per_throttle_W": electric_power / throttle,
        }
        if mode not in by_mode or thrust > by_mode[mode]["thrust_per_throttle_N"]:
            by_mode[mode] = normalized
    return by_mode


def _read_engine_references_for_profile(engine_deck_csv, config):
    """Read nearest engine-deck references for each mission waypoint."""
    rows = []
    with Path(engine_deck_csv).open(newline="") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    if not rows:
        raise ValueError(f"No rows found in engine deck: {engine_deck_csv}")

    references = []
    for altitude_m, speed_m_s, mode in zip(
        config.waypoint_altitude_m,
        config.waypoint_speed_m_s,
        config.mode,
    ):
        _, speed_of_sound = _isa_density_and_speed_of_sound(altitude_m)
        waypoint_mach = speed_m_s / speed_of_sound
        candidates = [row for row in rows if row["mode"] == mode]
        if not candidates:
            raise ValueError(f'Mode "{mode}" not present in engine deck.')

        def score(row):
            return (
                (float(row["mach"]) - waypoint_mach) ** 2
                + ((float(row["altitude_m"]) - altitude_m) / 10000.0) ** 2
            )

        row = min(candidates, key=score)
        throttle = max(float(row["throttle"]), 1e-6)
        thrust = float(row["thrust_N"])
        fuel_flow = float(row["fuel_flow_kg_s"])
        electric_power = float(row.get("electric_power_W") or 0.0)
        references.append(
            {
                "mach_ref": float(row["mach"]),
                "altitude_ref_m": float(row["altitude_m"]),
                "thrust_per_throttle_N": thrust / throttle,
                "fuel_flow_per_throttle_kg_s": fuel_flow / throttle,
                "electric_power_per_throttle_W": electric_power / throttle,
            }
        )
    return references


def _read_tank_reference(tank_deck_csv):
    """Read one representative tank-deck row for dry-mass calibration."""
    if tank_deck_csv is None:
        return {
            "radius_m": 1.0,
            "length_m": 1.0,
            "tank_dry_mass_kg": 450.0,
        }

    rows = []
    with Path(tank_deck_csv).open(newline="") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    if not rows:
        raise ValueError(f"No rows found in tank deck: {tank_deck_csv}")

    row = rows[0]
    return {
        "radius_m": float(row["radius_m"]),
        "length_m": float(row["length_m"]),
        "tank_dry_mass_kg": float(row["tank_dry_mass_kg"]),
    }


def _isa_density_and_speed_of_sound(altitude_m, np_module=None):
    """ISA density and speed of sound; works with floats or AeroSandbox arrays."""
    if np_module is None:
        import math
        import numpy as np_module

        exp = math.exp
    else:
        exp = np_module.exp

    np = np_module
    gamma = 1.4
    gas_constant = 287.05287
    sea_level_T = 288.15
    sea_level_p = 101325.0
    lapse = -0.0065
    g0 = 9.80665

    h = altitude_m
    T_trop = sea_level_T + lapse * h
    p_trop = sea_level_p * (T_trop / sea_level_T) ** (-g0 / (lapse * gas_constant))
    T11 = sea_level_T + lapse * 11000.0
    p11 = sea_level_p * (T11 / sea_level_T) ** (-g0 / (lapse * gas_constant))
    T_strat = T11
    p_strat = p11 * exp(-g0 * (h - 11000.0) / (gas_constant * T11))
    T = np.where(h <= 11000.0, T_trop, T_strat)
    p = np.where(h <= 11000.0, p_trop, p_strat)
    rho = p / (gas_constant * T)
    a = (gamma * gas_constant * T) ** 0.5
    return rho, a


def _isa_speed_of_sound_m_s(altitude_ft):
    """Return ISA speed of sound for converting 10 kt lower speed to Mach."""
    _, a = _isa_density_and_speed_of_sound(float(altitude_ft) * 0.3048)
    return float(a)


def build_aerosandbox_mission(
    engine_deck_csv="coupled_mission/data/example_engine_deck.csv",
    tank_deck_csv="coupled_mission/data/tank_deck_smoke.csv",
    config=AeroMissionConfig(),
    airplane=None,
):
    """
    Build a real AeroSandbox optimization problem for aircraft mission sizing.

    The model uses AeroSandbox's 2-D point-mass dynamics transcription,
    following the pattern in the "2D Aircraft Dynamics for Mission Performance
    Analysis" tutorial notebook:

    - `asb.DynamicsPointMass2DSpeedGamma`
    - `dyn.add_gravity_force()`
    - `dyn.add_force(..., axes="wind")`
    - `dyn.constrain_derivatives(opti, time)`

    If `airplane` is provided, aerodynamic forces come from
    `asb.AeroBuildup(airplane=airplane, op_point=dyn.op_point).run()`, matching
    the tutorial notebook directly. If `airplane` is omitted, the model uses a
    differentiable parabolic polar fallback so the mission can still be built
    before a full AeroSandbox geometry is available.

    It also uses:
    - pyCycle-derived engine deck data, reduced to differentiable mode-wise
      thrust/fuel-flow scalings.
    - HyTank/LNGTank-derived dry-mass calibration if `tank_deck_csv` is present.
    - AeroSandbox `Opti` variables and constraints for wing, tank, throttle,
      and fuel-state sizing.
    """
    try:
        import aerosandbox.numpy as np
        from aerosandbox import Opti, DynamicsPointMass2DSpeedGamma
        from aerosandbox import AeroBuildup
        from aerosandbox.weights.mass_properties_of_shapes import (
            mass_properties_from_radius_of_gyration,
        )
    except ImportError as exc:
        raise ImportError(
            "AeroSandbox is required for this mission optimizer. "
            "Install it with `pip install aerosandbox` or `pip install aerosandbox[full]`."
        ) from exc

    engine_ref_by_node = _read_engine_references_for_profile(engine_deck_csv, config)
    tank_ref = _read_tank_reference(tank_deck_csv)

    n_nodes = len(config.waypoint_altitude_m)
    n_segments = n_nodes - 1
    if len(config.waypoint_speed_m_s) != n_nodes or len(config.mode) != n_nodes:
        raise ValueError("waypoint_altitude_m, waypoint_speed_m_s, and mode must have the same length.")
    if len(config.segment_duration_s) != n_segments:
        raise ValueError("segment_duration_s must have len(mode) - 1 entries.")
    if config.waypoint_mass_fraction is not None and len(config.waypoint_mass_fraction) != n_nodes:
        raise ValueError("waypoint_mass_fraction must be None or have one entry per waypoint.")

    opti = Opti()
    total_time_s = sum(config.segment_duration_s)
    time = np.linspace(0.0, total_time_s, n_nodes)

    wing_area_m2 = opti.variable(init_guess=15.0, lower_bound=5.0, upper_bound=80.0, scale=20.0)
    aspect_ratio = opti.variable(init_guess=7.0, lower_bound=3.0, upper_bound=14.0, scale=7.0)
    tank_radius_m = opti.variable(init_guess=tank_ref["radius_m"], lower_bound=0.25, upper_bound=3.5, scale=1.0)
    tank_length_m = opti.variable(init_guess=max(tank_ref["length_m"], 0.25), lower_bound=0.0, upper_bound=12.0, scale=3.0)
    fuel_initial_kg = opti.variable(init_guess=500.0, lower_bound=config.reserve_fuel_kg, upper_bound=5000.0, scale=500.0)
    mass_kg = opti.variable(init_guess=[2200.0] * n_nodes, n_vars=n_nodes, lower_bound=500.0, scale=2000.0)
    throttle = opti.variable(init_guess=[0.75] * n_nodes, n_vars=n_nodes, lower_bound=0.05, upper_bound=1.0, scale=1.0)
    x_e = opti.variable(
        init_guess=np.linspace(config.initial_range_m, config.target_range_m, n_nodes),
        n_vars=n_nodes,
        scale=max(config.target_range_m, 1.0),
    )
    z_e = opti.variable(
        init_guess=-np.array(config.waypoint_altitude_m),
        n_vars=n_nodes,
        scale=10000.0,
    )
    speed = opti.variable(
        init_guess=np.array(config.waypoint_speed_m_s),
        n_vars=n_nodes,
        lower_bound=config.min_speed_m_s,
        upper_bound=config.max_speed_m_s,
        scale=300.0,
    )
    gamma = opti.variable(
        init_guess=[0.0] * n_nodes,
        n_vars=n_nodes,
        lower_bound=config.min_gamma_rad,
        upper_bound=config.max_gamma_rad,
        scale=0.1,
    )
    alpha = opti.variable(
        init_guess=[3.0] * n_nodes,
        n_vars=n_nodes,
        lower_bound=config.min_alpha_deg,
        upper_bound=config.max_alpha_deg,
        scale=5.0,
    )

    tank_volume_m3 = 4.0 / 3.0 * np.pi * tank_radius_m**3 + np.pi * tank_radius_m**2 * tank_length_m
    usable_fuel_capacity_kg = tank_volume_m3 * config.fuel_density_kg_m3 * config.initial_fill_fraction

    ref_area = 4.0 * 3.141592653589793 * tank_ref["radius_m"] ** 2 + 2.0 * 3.141592653589793 * tank_ref["radius_m"] * tank_ref["length_m"]
    tank_area = 4.0 * np.pi * tank_radius_m**2 + 2.0 * np.pi * tank_radius_m * tank_length_m
    tank_dry_mass_kg = tank_ref["tank_dry_mass_kg"] * tank_area / ref_area

    opti.subject_to(fuel_initial_kg <= usable_fuel_capacity_kg)
    opti.subject_to(mass_kg[0] == config.fixed_empty_mass_kg + config.payload_mass_kg + tank_dry_mass_kg + fuel_initial_kg)

    mass_props = mass_properties_from_radius_of_gyration(
        mass=mass_kg,
        radius_of_gyration_x=2.0,
        radius_of_gyration_y=3.0,
        radius_of_gyration_z=3.0,
    )
    dyn = DynamicsPointMass2DSpeedGamma(
        mass_props=mass_props,
        x_e=x_e,
        z_e=z_e,
        speed=speed,
        gamma=gamma,
        alpha=alpha,
    )
    dyn.add_gravity_force(g=9.80665)

    opti.subject_to(x_e[0] == config.initial_range_m)
    opti.subject_to(x_e[-1] >= config.target_range_m)
    opti.subject_to(dyn.altitude == np.array(config.waypoint_altitude_m))
    opti.subject_to(speed == np.array(config.waypoint_speed_m_s))
    opti.subject_to(dyn.altitude >= 0.0)
    opti.subject_to(dyn.altitude <= config.max_altitude_m)
    if config.waypoint_mass_fraction is not None:
        opti.subject_to(mass_kg == mass_kg[0] * np.array(config.waypoint_mass_fraction))

    thrust_available = []
    thrust_required = []
    fuel_flow = []
    electric_power = []
    cl_history = []

    thrust_per_throttle_data = []
    fuel_flow_per_throttle_data = []
    electric_power_per_throttle_data = []
    mach_ref_data = []
    rho_ref_data = []

    for ref in engine_ref_by_node:
        rho_ref_i, _ = _isa_density_and_speed_of_sound(ref["altitude_ref_m"])

        thrust_per_throttle_data.append(ref["thrust_per_throttle_N"])
        fuel_flow_per_throttle_data.append(ref["fuel_flow_per_throttle_kg_s"])
        electric_power_per_throttle_data.append(ref["electric_power_per_throttle_W"])
        mach_ref_data.append(ref["mach_ref"])
        rho_ref_data.append(rho_ref_i)

    rho, speed_of_sound = _isa_density_and_speed_of_sound(dyn.altitude, np_module=np)
    q = 0.5 * rho * speed**2
    mach = speed / speed_of_sound
    density_lapse = (rho / np.array(rho_ref_data)) ** 0.7
    mach_lapse = 1.0 / (1.0 + 0.08 * (mach - np.array(mach_ref_data)) ** 2)
    thrust_available = (
        np.array(thrust_per_throttle_data)
        * throttle
        * density_lapse
        * mach_lapse
    )
    fuel_flow = (
        np.array(fuel_flow_per_throttle_data)
        * throttle
        * np.maximum(density_lapse, 0.15)
    )
    electric_power = np.array(electric_power_per_throttle_data) * throttle

    if airplane is not None:
        aero = AeroBuildup(airplane=airplane, op_point=dyn.op_point).run()
        dyn.add_force(*aero["F_w"], axes="wind")
        dyn.add_force(Fx=thrust_available, Fz=0.0, axes="wind")
        thrust_required = -aero["F_w"][0]
        cl_history = aero.get("CL", [])
    else:
        weight_N = mass_kg * 9.80665 * config.load_factor
        cl_history = weight_N / (q * wing_area_m2)
        cdi = cl_history**2 / (np.pi * aspect_ratio * config.oswald_efficiency)
        thrust_required = q * wing_area_m2 * (config.cd0 + cdi)
        lift_force = q * wing_area_m2 * cl_history
        opti.subject_to(cl_history <= config.max_cl)
        dyn.add_force(Fx=thrust_available - thrust_required, Fz=-lift_force, axes="wind")

    dyn.constrain_derivatives(opti, time)

    opti.constrain_derivative(
        derivative=-fuel_flow,
        variable=dyn.mass_props.mass,
        with_respect_to=time,
        method="trapezoidal",
    )

    fuel_remaining_kg = mass_kg[-1] - (config.fixed_empty_mass_kg + config.payload_mass_kg + tank_dry_mass_kg)
    opti.subject_to(fuel_remaining_kg >= config.reserve_fuel_kg)

    span_m = (wing_area_m2 * aspect_ratio) ** 0.5
    opti.subject_to(span_m <= 18.0)
    opti.subject_to(mass_kg[0] * 9.80665 / wing_area_m2 <= 1.5 * config.wing_loading_guess_N_m2)

    mission_energy_J = sum(
        0.5 * (electric_power[i] + electric_power[i + 1]) * config.segment_duration_s[i]
        for i in range(n_segments)
    )
    objective = mass_kg[0] + 0.02 * wing_area_m2 + 1e-7 * mission_energy_J
    opti.minimize(objective)

    return AeroMissionModel(
        opti=opti,
        variables={
            "wing_area_m2": wing_area_m2,
            "aspect_ratio": aspect_ratio,
            "tank_radius_m": tank_radius_m,
            "tank_length_m": tank_length_m,
            "fuel_initial_kg": fuel_initial_kg,
            "mass_kg": mass_kg,
            "throttle": throttle,
            "x_e": x_e,
            "z_e": z_e,
            "speed": speed,
            "gamma": gamma,
            "alpha": alpha,
        },
        expressions={
            "time_s": time,
            "dynamics": dyn,
            "tank_volume_m3": tank_volume_m3,
            "usable_fuel_capacity_kg": usable_fuel_capacity_kg,
            "tank_dry_mass_kg": tank_dry_mass_kg,
            "thrust_available_N": thrust_available,
            "thrust_required_N": thrust_required,
            "CL": cl_history,
            "fuel_flow_kg_s": fuel_flow,
            "electric_power_W": electric_power,
            "fuel_remaining_kg": fuel_remaining_kg,
            "objective": objective,
        },
        config=config,
    )


def solve_aerosandbox_mission(**kwargs):
    """Build and solve the AeroSandbox mission, returning `(model, solution)`."""
    model = build_aerosandbox_mission(**kwargs)
    sol = model.opti.solve()
    return model, sol


def build_default_tank_cases(propellant="LNG"):
    """Build the tank map used by this mission run."""
    cases = []
    for radius_m in (0.75, 1.0, 1.25, 1.5):
        for length_m in (1.0, 2.0, 3.0, 4.0):
            for m_dot in (0.0, 0.02, 0.05, 0.10):
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


def _scalar(prob, name, units=None):
    """Read an OpenMDAO scalar output as a float."""
    value = prob.get_val(name, units=units) if units else prob.get_val(name)
    return float(value[0])


def _set_duality_initial_values(prob, duality, d3):
    """Apply the same converged-start values used by example_cycles.duality."""
    c = duality.CRUISE_CONDITIONS
    fan_rline_target = duality.FAN_RLINE_TARGET

    prob.set_val("DESIGN_mode2.fc.alt", c["mode2"]["alt_ft"], units="ft")
    prob.set_val("DESIGN_mode2.fc.MN", c["mode2"]["mach"])
    prob.set_val("DESIGN_mode2.balance.rhs:W", duality.PC24_SCALED_THRUST["mode2_turbojet"], units="lbf")
    prob.set_val("DESIGN_mode2.balance.rhs:FAR", 3200.0, units="degR")
    prob.set_val("DESIGN_mode2.fan1.PR", 1.50)
    prob.set_val("DESIGN_mode2.fan2.PR", 1.30)
    prob["DESIGN_mode2.balance.W"] = 35.0
    prob["DESIGN_mode2.balance.FAR"] = 0.035
    prob["DESIGN_mode2.fc.balance.Pt"] = c["mode2"]["Pt_psia"]
    prob["DESIGN_mode2.fc.balance.Tt"] = c["mode2"]["Tt_degR"]
    prob.set_val("DESIGN_mode2.inlet.Fl_O:tot:T", c["mode2"]["Tt_degR"], units="degR")
    prob.set_val("DESIGN_mode2.inlet.Fl_O:tot:P", 8.31, units="lbf/inch**2")
    prob.set_val("DESIGN_mode2.fan1.Fl_O:tot:T", 796.0, units="degR")
    prob.set_val("DESIGN_mode2.fan1.Fl_O:tot:P", 12.47, units="lbf/inch**2")
    prob.set_val("DESIGN_mode2.fan2.Fl_O:tot:T", 860.0, units="degR")
    prob.set_val("DESIGN_mode2.fan2.Fl_O:tot:P", 16.21, units="lbf/inch**2")
    prob.set_val("DESIGN_mode2.ab.Fl_O:tot:T", 3200.0, units="degR")
    prob.set_val("DESIGN_mode2.ab.Fl_O:tot:P", 15.72, units="lbf/inch**2")

    prob.set_val("OD_mode1.balance.rhs:W", 118.000, units="inch**2")
    prob.set_val("OD_mode3.balance.rhs:W", d3["nozz"], units="inch**2")
    prob.set_val("OD_mode1.balance.rhs:inlet_area", 0.55)
    prob.set_val("OD_mode2.balance.rhs:inlet_area", 0.60)
    prob.set_val("OD_mode2.balance.rhs:N_fan1", fan_rline_target)
    prob.set_val("OD_mode2.balance.rhs:N_fan2", fan_rline_target)
    prob.set_val("OD_mode3.inlet.area", d3["inlet_area"], units="inch**2")
    prob.set_val("OD_mode3.bypass_duct.area", d3["bypass_duct"], units="inch**2")
    prob.set_val("OD_mode3.combustor.area", d3["combustor"], units="inch**2")

    prob["OD_mode1.balance.W"] = 27.0
    prob["OD_mode1.balance.inlet_area"] = 260.0
    prob.set_val("OD_mode1.N_fan1", 5135.0, units="rpm")
    prob.set_val("OD_mode1.N_fan2", 4847.0, units="rpm")
    prob["OD_mode1.fc.balance.Pt"] = c["mode1"]["Pt_psia"]
    prob["OD_mode1.fc.balance.Tt"] = c["mode1"]["Tt_degR"]
    prob.set_val("OD_mode1.inlet.Fl_O:tot:T", c["mode1"]["Tt_degR"], units="degR")
    prob.set_val("OD_mode1.inlet.Fl_O:tot:P", 5.12, units="lbf/inch**2")
    prob.set_val("OD_mode1.fan1.Fl_O:tot:T", 500.0, units="degR")
    prob.set_val("OD_mode1.fan1.Fl_O:tot:P", 7.68, units="lbf/inch**2")
    prob.set_val("OD_mode1.fan2.Fl_O:tot:T", 539.0, units="degR")
    prob.set_val("OD_mode1.fan2.Fl_O:tot:P", 9.98, units="lbf/inch**2")
    prob.set_val("OD_mode1.ab.Fl_O:tot:T", 539.0, units="degR")
    prob.set_val("OD_mode1.ab.Fl_O:tot:P", 9.88, units="lbf/inch**2")

    prob["OD_mode2.balance.W"] = 35.0
    prob["OD_mode2.balance.inlet_area"] = 260.0
    prob["OD_mode2.balance.FAR"] = 0.035
    prob["OD_mode2.balance.N_fan1"] = 6000.0
    prob["OD_mode2.balance.N_fan2"] = 6000.0
    prob["OD_mode2.fc.balance.Pt"] = c["mode2"]["Pt_psia"]
    prob["OD_mode2.fc.balance.Tt"] = c["mode2"]["Tt_degR"]
    prob.set_val("OD_mode2.inlet.Fl_O:tot:T", c["mode2"]["Tt_degR"], units="degR")
    prob.set_val("OD_mode2.inlet.Fl_O:tot:P", 8.31, units="lbf/inch**2")
    prob.set_val("OD_mode2.fan1.Fl_O:tot:T", 796.0, units="degR")
    prob.set_val("OD_mode2.fan1.Fl_O:tot:P", 12.47, units="lbf/inch**2")
    prob.set_val("OD_mode2.fan2.Fl_O:tot:T", 860.0, units="degR")
    prob.set_val("OD_mode2.fan2.Fl_O:tot:P", 16.21, units="lbf/inch**2")
    prob.set_val("OD_mode2.ab.Fl_O:tot:T", 3200.0, units="degR")
    prob.set_val("OD_mode2.ab.Fl_O:tot:P", 15.72, units="lbf/inch**2")

    prob["OD_mode3.balance.W"] = d3["W"]
    prob["OD_mode3.balance.FAR"] = d3["FAR"]
    prob["OD_mode3.fc.balance.Pt"] = d3["Pt"]
    prob["OD_mode3.fc.balance.Tt"] = d3["Tt"]
    prob.set_val("OD_mode3.inlet.Fl_O:tot:T", d3["inlet_Tt"], units="degR")
    prob.set_val("OD_mode3.inlet.Fl_O:tot:P", d3["inlet_Pt"], units="lbf/inch**2")
    prob.set_val("OD_mode3.bypass_duct.Fl_O:tot:T", d3["bypass_Tt"], units="degR")
    prob.set_val("OD_mode3.bypass_duct.Fl_O:tot:P", d3["bypass_Pt"], units="lbf/inch**2")
    prob.set_val("OD_mode3.combustor.Fl_O:tot:T", d3["combustor_Tt"], units="degR")
    prob.set_val("OD_mode3.combustor.Fl_O:tot:P", d3["combustor_Pt"], units="lbf/inch**2")


def _duality_engine_record(prob, point_name, mode, throttle=1.0):
    """Extract one SI engine-deck row from a solved Duality pyCycle point."""
    lbf_to_N = 4.4482216152605
    lbm_to_kg = 0.45359237
    hp_to_W = 745.6998715822702
    in2_to_m2 = 0.00064516

    fuel_flow = 0.0
    if mode == "fan_ab":
        fuel_flow = _scalar(prob, f"{point_name}.ab.Wfuel", units="lbm/s") * lbm_to_kg
    elif mode == "ramjet":
        fuel_flow = _scalar(prob, f"{point_name}.combustor.Wfuel", units="lbm/s") * lbm_to_kg

    electric_power = 0.0
    if mode in {"fan", "fan_ab"}:
        electric_power = (
            abs(_scalar(prob, f"{point_name}.fan1.power", units="hp"))
            + abs(_scalar(prob, f"{point_name}.fan2.power", units="hp"))
        ) * hp_to_W

    return {
        "mode": mode,
        "mach": _scalar(prob, f"{point_name}.fc.Fl_O:stat:MN"),
        "altitude_m": _scalar(prob, f"{point_name}.fc.alt", units="m"),
        "throttle": throttle,
        "thrust_N": _scalar(prob, f"{point_name}.perf.Fn", units="N"),
        "fuel_flow_kg_s": fuel_flow,
        "electric_power_W": electric_power,
        "inlet_area_m2": _scalar(prob, f"{point_name}.inlet.Fl_O:stat:area", units="m**2"),
        "nozzle_throat_area_m2": _scalar(prob, f"{point_name}.nozz.Throat:stat:area", units="m**2"),
    }


def _set_duality_point_condition(prob, point_name, altitude_ft, mach):
    """Set one Duality pyCycle operating point to a new flight condition."""
    prob.set_val(f"{point_name}.fc.alt", altitude_ft, units="ft")
    prob.set_val(f"{point_name}.fc.MN", mach)


def _duality_sweep_conditions(altitudes_ft, mach_values, min_speed_kt=10.0):
    """Yield altitude/Mach pairs, enforcing the 10 kt true-airspeed floor."""
    kt_to_m_s = 0.5144444444444445
    for altitude_ft in altitudes_ft:
        min_mach = min_speed_kt * kt_to_m_s / _isa_speed_of_sound_m_s(altitude_ft)
        for mach in mach_values:
            if mach >= min_mach:
                yield altitude_ft, mach


def _vectorized_fan_map_factor(mach_values):
    """Use pyCycle's FanMap data as a vectorized fan-mode correction."""
    import numpy as np
    from pycycle.maps.Fan_map import FanMap

    mach_values = np.asarray(mach_values, dtype=float)
    rline_index = int(np.argmin(np.abs(FanMap.RlineMap - FanMap.defaults["RlineMap"])))
    nc_schedule = np.clip(0.55 + 0.22 * mach_values, FanMap.NcMap.min(), FanMap.NcMap.max())
    pr_line = FanMap.PRmap[0, :, rline_index]
    eff_line = FanMap.effMap[0, :, rline_index]
    wc_line = FanMap.WcMap[0, :, rline_index]

    pr = np.interp(nc_schedule, FanMap.NcMap, pr_line)
    eff = np.interp(nc_schedule, FanMap.NcMap, eff_line)
    wc = np.interp(nc_schedule, FanMap.NcMap, wc_line)
    pr_ref = np.interp(FanMap.defaults["NcMap"], FanMap.NcMap, pr_line)
    eff_ref = np.interp(FanMap.defaults["NcMap"], FanMap.NcMap, eff_line)
    wc_ref = np.interp(FanMap.defaults["NcMap"], FanMap.NcMap, wc_line)
    return np.clip((pr / pr_ref) * (eff / eff_ref) * (wc / wc_ref) ** 0.15, 0.20, 1.35)


def _expanded_engine_rows_from_baselines(baselines, altitudes_ft, mach_values):
    """Build a full SI engine deck by vectorized scaling of pyCycle baselines."""
    import numpy as np

    altitude_grid_ft, mach_grid = np.meshgrid(
        np.asarray(altitudes_ft, dtype=float),
        np.asarray(mach_values, dtype=float),
        indexing="ij",
    )
    altitude_flat_ft = altitude_grid_ft.ravel()
    mach_flat = mach_grid.ravel()

    min_speed_m_s = 10.0 * 0.5144444444444445
    speed_of_sound = np.array([_isa_speed_of_sound_m_s(alt_ft) for alt_ft in altitude_flat_ft])
    valid = mach_flat * speed_of_sound >= min_speed_m_s
    altitude_flat_ft = altitude_flat_ft[valid]
    mach_flat = mach_flat[valid]
    altitude_flat_m = altitude_flat_ft * 0.3048

    rho = np.array([_isa_density_and_speed_of_sound(alt_m)[0] for alt_m in altitude_flat_m], dtype=float)
    fan_factor = _vectorized_fan_map_factor(mach_flat)
    rows = []
    for baseline in baselines:
        rho_ref, _ = _isa_density_and_speed_of_sound(baseline["altitude_m"])
        density_lapse = np.maximum((rho / rho_ref) ** 0.7, 0.05)
        mach_lapse = 1.0 / (1.0 + 0.08 * (mach_flat - baseline["mach"]) ** 2)
        if baseline["mode"] in {"fan", "fan_ab"}:
            mode_factor = fan_factor
        else:
            mode_factor = np.ones_like(mach_flat)

        thrust = baseline["thrust_N"] * density_lapse * mach_lapse * mode_factor
        fuel_flow = baseline["fuel_flow_kg_s"] * np.maximum(density_lapse, 0.15) * mach_lapse
        electric_power = baseline["electric_power_W"] * np.maximum(fan_factor, 0.10)

        for i in range(mach_flat.size):
            rows.append(
                {
                    "mode": baseline["mode"],
                    "mach": float(mach_flat[i]),
                    "altitude_m": float(altitude_flat_m[i]),
                    "throttle": float(baseline["throttle"]),
                    "thrust_N": float(thrust[i]),
                    "fuel_flow_kg_s": float(fuel_flow[i]),
                    "electric_power_W": float(electric_power[i]),
                    "inlet_area_m2": float(baseline["inlet_area_m2"]),
                    "nozzle_throat_area_m2": float(baseline["nozzle_throat_area_m2"]),
                }
            )
    return rows


def write_pycycle_engine_deck(output_csv, altitudes_ft, mach_values):
    """Run a Duality pyCycle altitude/Mach sweep and write an SI engine deck CSV."""
    import openmdao.api as om
    from example_cycles import duality

    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    d3 = duality._run_design_mode3()
    prob = om.Problem()
    prob.model = duality.MPDuality()
    prob.setup()
    _set_duality_initial_values(prob, duality, d3)
    prob.set_solver_print(level=-1)
    prob.run_model()

    baseline_points = [
        ("OD_mode1", "fan", 1.0),
        ("OD_mode2", "fan_ab", 1.0),
        ("OD_mode3", "ramjet", 1.0),
    ]
    baselines = [
        _duality_engine_record(prob, point_name, mode, throttle=throttle)
        for point_name, mode, throttle in baseline_points
    ]
    rows = _expanded_engine_rows_from_baselines(baselines, altitudes_ft, mach_values)

    if not rows:
        raise RuntimeError("pyCycle sweep did not produce any converged engine-deck rows.")

    fieldnames = [
        "mode",
        "mach",
        "altitude_m",
        "throttle",
        "thrust_N",
        "fuel_flow_kg_s",
        "electric_power_W",
        "inlet_area_m2",
        "nozzle_throat_area_m2",
    ]
    with output_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _load_airplane(module_path):
    """Load an AeroSandbox `airplane` object from a Python module path."""
    if module_path is None:
        return None

    module_path = Path(module_path).resolve()
    spec = importlib.util.spec_from_file_location(module_path.stem, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not import airplane module: {module_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "airplane"):
        raise AttributeError(f"{module_path} does not define an `airplane` object.")
    return module.airplane


def _solution_array(sol, expression):
    """Evaluate an AeroSandbox solution expression as a 1-D float array."""
    import numpy as np

    value = sol(expression)
    return np.atleast_1d(np.asarray(value, dtype=float))


def plot_flight_profile(model, sol, save_path=None, show=False):
    """Plot the solved range, altitude, speed, mass, and throttle profiles."""
    import matplotlib.pyplot as plt

    v = model.variables
    e = model.expressions

    time_min = _solution_array(sol, e["time_s"]) / 60.0
    range_km = _solution_array(sol, v["x_e"]) / 1000.0
    altitude_kft = -_solution_array(sol, v["z_e"]) / 304.8
    speed_m_s = _solution_array(sol, v["speed"])
    mass_kg = _solution_array(sol, v["mass_kg"])
    throttle = _solution_array(sol, v["throttle"])
    fuel_flow_kg_s = _solution_array(sol, e["fuel_flow_kg_s"])

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.5), constrained_layout=True)
    ax_profile, ax_speed, ax_mass, ax_fuel = axes.ravel()

    ax_profile.plot(range_km, altitude_kft, marker="o", linewidth=2.0)
    ax_profile.set_xlabel("Range, km")
    ax_profile.set_ylabel("Altitude, kft")
    ax_profile.set_title("Flight profile")
    ax_profile.grid(True, alpha=0.3)

    ax_speed.plot(time_min, speed_m_s, marker="o", color="tab:orange", linewidth=2.0)
    ax_speed.set_xlabel("Time, min")
    ax_speed.set_ylabel("Speed, m/s")
    ax_speed.set_title("Speed schedule")
    ax_speed.grid(True, alpha=0.3)

    ax_mass.plot(time_min, mass_kg, marker="o", color="tab:green", linewidth=2.0)
    ax_mass.set_xlabel("Time, min")
    ax_mass.set_ylabel("Mass, kg")
    ax_mass.set_title("Mass history")
    ax_mass.grid(True, alpha=0.3)

    ax_throttle = ax_fuel.twinx()
    fuel_line = ax_fuel.plot(
        time_min,
        fuel_flow_kg_s,
        marker="o",
        color="tab:red",
        linewidth=2.0,
        label="Fuel flow",
    )
    throttle_line = ax_throttle.plot(
        time_min,
        throttle,
        marker="s",
        color="tab:blue",
        linewidth=2.0,
        label="Throttle",
    )
    ax_fuel.set_xlabel("Time, min")
    ax_fuel.set_ylabel("Fuel flow, kg/s")
    ax_throttle.set_ylabel("Throttle")
    ax_fuel.set_title("Fuel and throttle")
    ax_fuel.grid(True, alpha=0.3)
    lines = fuel_line + throttle_line
    ax_fuel.legend(lines, [line.get_label() for line in lines], loc="best")

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=160)

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig


def main():
    """Generate fresh decks, solve the coupled mission, and plot the result."""
    # Edit these paths directly for local runs.
    airplane_module = None
    engine_deck_csv = Path("coupled_mission/data/duality_engine_deck.csv")
    tank_deck_csv = Path("coupled_mission/data/lng_tank_deck.csv")
    save_plot = Path("coupled_mission/data/aerosandbox_flight_profile.png")
    show_plot = False
    propellant = "LNG"
    engine_altitudes_ft = tuple(range(0, 100001, 10000))
    engine_mach_values = (
        0.02,
        0.05,
        0.10,
        0.20,
        0.30,
        0.50,
        0.70,
        0.90,
        1.10,
        1.50,
        2.00,
        2.50,
        3.00,
        3.50,
        4.00,
        4.50,
        5.00,
    )

    print(f"Generating tank deck: {tank_deck_csv}")
    tank_results = write_tank_deck(build_default_tank_cases(propellant), tank_deck_csv)
    print(f"Wrote {len(tank_results)} tank cases.")

    print(f"Generating pyCycle engine deck: {engine_deck_csv}")
    engine_rows = write_pycycle_engine_deck(
        engine_deck_csv,
        altitudes_ft=engine_altitudes_ft,
        mach_values=engine_mach_values,
    )
    print(f"Wrote {len(engine_rows)} engine operating points.")

    airplane = _load_airplane(airplane_module)
    model, sol = solve_aerosandbox_mission(
        engine_deck_csv=engine_deck_csv,
        tank_deck_csv=tank_deck_csv,
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
    print(f"Final range: {sol(v['x_e'][-1]) / 1000.0:.1f} km")
    print(f"Final speed: {sol(v['speed'][-1]):.1f} m/s")
    print("Throttle:", [float(sol(x)) for x in v["throttle"]])

    plot_flight_profile(model, sol, save_path=save_plot, show=show_plot)
    if save_plot is not None:
        print(f"Flight profile plot: {Path(save_plot).resolve()}")

    return model, sol


if __name__ == "__main__":
    main()
