"""
Central path resolution for the LVL Seismic project.

This module is the single source of truth for the project directory
layout. No other module should construct paths manually.
"""

from __future__ import annotations

from pathlib import Path

# -----------------------------------------------------------------------------
# Project root
# -----------------------------------------------------------------------------

PROJECT_DIR = Path(__file__).resolve().parents[2]

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------

DATA_DIR = PROJECT_DIR / "data"

INPUT_DIR = DATA_DIR / "input"

RAW_DIR = INPUT_DIR / "raw"

GEOMETRY_DIR = INPUT_DIR / "geometry"

METADATA_DIR = INPUT_DIR / "metadata"

EXTERNAL_DIR = DATA_DIR / "external"

EXAMPLES_DIR = DATA_DIR / "examples"

IMPORT_PICKS_DIR = EXAMPLES_DIR / "picks"

# -----------------------------------------------------------------------------
# Processing
# -----------------------------------------------------------------------------

PROCESSING_DIR = DATA_DIR / "processing"

CACHE_DIR = PROCESSING_DIR / "cache"

PICKS_DIR = PROCESSING_DIR / "picks"

SESSIONS_DIR = PROCESSING_DIR / "sessions"

# -----------------------------------------------------------------------------
# Results
# -----------------------------------------------------------------------------

RESULTS_DIR = DATA_DIR / "results"

VELOCITY_RESULTS_DIR = RESULTS_DIR / "velocity"

PLOTS_DIR = RESULTS_DIR / "plots"

REPORTS_DIR = RESULTS_DIR / "reports"

# -----------------------------------------------------------------------------
# Projects
# -----------------------------------------------------------------------------

PROJECTS_DIR = PROJECT_DIR / "projects"

# =============================================================================
# Public API
# =============================================================================

__all__ = [
    "PROJECT_DIR",
    "DATA_DIR",
    "INPUT_DIR",
    "RAW_DIR",
    "GEOMETRY_DIR",
    "METADATA_DIR",
    "EXTERNAL_DIR",
    "EXAMPLES_DIR",
    "PROCESSING_DIR",
    "CACHE_DIR",
    "PICKS_DIR",
    "SESSIONS_DIR",
    "RESULTS_DIR",
    "VELOCITY_RESULTS_DIR",
    "PLOTS_DIR",
    "REPORTS_DIR",
    "PROJECTS_DIR",
]