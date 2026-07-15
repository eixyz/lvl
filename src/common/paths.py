"""Central path resolution for the LVL project.

Single source of truth for the project root, data, and output directories.
Both `apps/` entry points and every `src/` module import from here so the
folder layout only needs to be defined once.
"""
from __future__ import annotations

from pathlib import Path

#: Root of the `lvl` project (parent of `apps/`, `data/`, `results/`, `src/`).
PROJECT_DIR: Path = Path(__file__).resolve().parents[2]

#: Raw SEG2 / geometry / field-report input data.
DATA_DIR: Path = PROJECT_DIR / "data"

#: Generated picks, exports, plots, and consolidated summaries.
OUTPUT_DIR: Path = PROJECT_DIR / "results"

#: Saved project sessions (project.json, per-project logs/results).
PROJECTS_DIR: Path = PROJECT_DIR / "projects"
