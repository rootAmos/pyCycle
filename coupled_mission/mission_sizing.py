"""Minimal AeroSandbox-facing mission sizing skeleton."""

from dataclasses import dataclass

import numpy as np


@dataclass
class MissionNode:
    """One mission collocation point."""

    time_s: float
    mach: float
    altitude_m: float
    required_thrust_N: float
    throttle: float
    mode: str


def evaluate_mission(nodes, engine_deck, initial_fuel_mass_kg):
    """
    Evaluate a fixed mission schedule against an engine deck.

    This is intentionally analysis-only. The next step is to replace fixed
    `nodes` with AeroSandbox `Opti` variables and constraints.
    """
    if len(nodes) < 2:
        raise ValueError("Mission evaluation needs at least two nodes.")

    fuel_mass = float(initial_fuel_mass_kg)
    history = []
    for i, node in enumerate(nodes):
        engine = engine_deck.evaluate(
            mach=node.mach,
            altitude_m=node.altitude_m,
            throttle=node.throttle,
            mode=node.mode,
        )
        thrust_margin_N = engine["thrust_N"] - node.required_thrust_N
        if i < len(nodes) - 1:
            dt = nodes[i + 1].time_s - node.time_s
            if dt < 0:
                raise ValueError("Mission node times must be monotonically increasing.")
            fuel_mass -= engine["fuel_flow_kg_s"] * dt
        history.append(
            {
                "time_s": node.time_s,
                "mode": node.mode,
                "mach": node.mach,
                "altitude_m": node.altitude_m,
                "throttle": node.throttle,
                "thrust_available_N": engine["thrust_N"],
                "required_thrust_N": node.required_thrust_N,
                "thrust_margin_N": thrust_margin_N,
                "fuel_flow_kg_s": engine["fuel_flow_kg_s"],
                "fuel_mass_kg": fuel_mass,
            }
        )
    return history


def fuel_burn_kg(history):
    """Return total fuel burned from an `evaluate_mission()` history."""
    masses = np.array([row["fuel_mass_kg"] for row in history], dtype=float)
    return float(masses[0] - masses[-1])
