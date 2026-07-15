from __future__ import annotations

from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent.parent
PROJECT_DIR = SCRIPT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "output"
REFRACTION_SCRIPT = SCRIPT_DIR / "lvl_refraction.py"
