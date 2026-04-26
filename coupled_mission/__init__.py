"""
Coupled aircraft mission scaffolding for pyCycle + HyTank + AeroSandbox.

The intended workflow is:

1. Run high-fidelity pyCycle and HyTank/LNGTank cases to generate decks.
2. Load those decks in an AeroSandbox mission/sizing model.
3. Keep unit conversion at this boundary so mission-level code stays SI.
"""

from .engine_deck import EngineDeck, EngineDeckRecord
from .tank_deck import TankCase, TankResult, run_tank_case, write_tank_deck
from .aerosandbox_mission import AeroMissionConfig, build_aerosandbox_mission

__all__ = [
    "AeroMissionConfig",
    "build_aerosandbox_mission",
    "EngineDeck",
    "EngineDeckRecord",
    "TankCase",
    "TankResult",
    "run_tank_case",
    "write_tank_deck",
]
