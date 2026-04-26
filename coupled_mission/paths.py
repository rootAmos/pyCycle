"""Repository path helpers for local pyCycle/HyTank development checkouts."""

from pathlib import Path
import sys


PYCYCLE_ROOT = Path(__file__).resolve().parents[1]
ASTROM_ROOT = PYCYCLE_ROOT.parent
HYTANK_ROOT = ASTROM_ROOT / "HyTank"


def ensure_hytank_on_path():
    """
    Make the sibling HyTank checkout importable without requiring installation.

    Returns
    -------
    pathlib.Path
        The expected HyTank repository root.
    """
    if HYTANK_ROOT.exists():
        hytank_str = str(HYTANK_ROOT)
        if hytank_str not in sys.path:
            sys.path.insert(0, hytank_str)
    return HYTANK_ROOT
