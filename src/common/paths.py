"""
Central path resolution for the LVL Seismic project.

Two kinds of paths live here:

1. Application-level paths - fixed locations shipped with the installed
   copy of this repository (e.g. bundled default geometry templates, the
   old example dataset). These never change at runtime.

2. Project-level paths - the working folders of a single user *project*
   (one field survey). A project can live anywhere on disk; it is created
   or opened from the GUI (see `src.io.project_io`). Nothing about a
   project's location is hardcoded here - `ProjectPaths` simply resolves
   the standard sub-layout underneath whatever root folder the user
   picked, and this module keeps track of the currently "active" project
   so the rest of the app can resolve `processing/`, `results/`, etc.
   without threading a path through every function call.

No other module should construct project paths manually - always go
through `ProjectPaths` / `get_active_project()` / `require_active_project()`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

# -----------------------------------------------------------------------------
# Application root (the installed copy of this repository)
# -----------------------------------------------------------------------------

APP_ROOT = Path(__file__).resolve().parents[2]

# Default parent folder new projects are created under. The "New Project"
# dialog lets the user pick a different location; this is only the default.
PROJECTS_ROOT = APP_ROOT / "projects"

# -----------------------------------------------------------------------------
# Legacy / bundled application resources
# -----------------------------------------------------------------------------
# `data/` is the flat layout this app used before the project system
# existed. It is kept only for:
#   - bundled default receiver-geometry templates (geometry100.txt / 200.txt)
#   - the bundled example dataset referenced from the docs
# It is NOT where project output goes anymore - see ProjectPaths below.

LEGACY_DATA_DIR = APP_ROOT / "data"
LEGACY_INPUT_DIR = LEGACY_DATA_DIR / "input"

GEOM_TEMPLATES_DIR = LEGACY_INPUT_DIR / "geometry"   # geometry100.txt / geometry200.txt
LEGACY_METADATA_DIR = LEGACY_INPUT_DIR / "metadata"
LEGACY_RAW_DIR = LEGACY_INPUT_DIR / "raw"
LEGACY_EXTERNAL_DIR = LEGACY_DATA_DIR / "external"
LEGACY_EXAMPLES_DIR = LEGACY_DATA_DIR / "examples"
LEGACY_IMPORT_PICKS_DIR = LEGACY_EXAMPLES_DIR / "picks"

# -----------------------------------------------------------------------------
# Project name validation
# -----------------------------------------------------------------------------

_PROJECT_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 _.\-]{0,63}$")


def sanitize_project_name(name: str) -> str:
    """Validate/normalize a user-supplied project name.

    Raises ValueError with a human-readable message if the name is empty
    or contains characters that would be unsafe as a folder name.
    """
    cleaned = (name or "").strip()
    if not cleaned:
        raise ValueError("Project name must not be empty.")
    if not _PROJECT_NAME_RE.match(cleaned):
        raise ValueError(
            "Project name may only contain letters, numbers, spaces, "
            "'-', '_' and '.', must start with a letter or number, and "
            "be at most 64 characters long."
        )
    return cleaned


# -----------------------------------------------------------------------------
# Project-level layout
# -----------------------------------------------------------------------------

PROJECT_FILE_NAME = "project.json"


@dataclass(frozen=True)
class ProjectPaths:
    """Resolves every folder/file that belongs to a single project.

    Layout on disk::

        <root>/
            project.json
            processing/
                cache/
                picks/<profile>/picks.json
                sessions/<profile>/picks.session.json
                              .../layer_analysis.session.json
            results/
                plots/<profile>/...
                reports/<profile>/layer_analysis.json
                velocity/lvl_velocity_summary.xlsx

    Raw SEG2 data, geometry/coordinate files and field-report metadata are
    intentionally NOT part of this layout - they are large, usually live
    outside the project folder, and are only referenced by path (see
    `src.io.project_io.Project`).
    """

    root: Path

    def __post_init__(self) -> None:
        object.__setattr__(self, "root", Path(self.root).resolve())

    # -- top level ------------------------------------------------------
    @property
    def project_file(self) -> Path:
        return self.root / PROJECT_FILE_NAME

    @property
    def processing_dir(self) -> Path:
        return self.root / "processing"

    @property
    def results_dir(self) -> Path:
        return self.root / "results"

    # -- processing -------------------------------------------------------
    @property
    def cache_dir(self) -> Path:
        return self.processing_dir / "cache"

    @property
    def picks_dir(self) -> Path:
        return self.processing_dir / "picks"

    @property
    def sessions_dir(self) -> Path:
        return self.processing_dir / "sessions"

    # -- results ----------------------------------------------------------
    @property
    def plots_dir(self) -> Path:
        return self.results_dir / "plots"

    @property
    def reports_dir(self) -> Path:
        return self.results_dir / "reports"

    @property
    def velocity_dir(self) -> Path:
        return self.results_dir / "velocity"

    # -- per-profile helpers ------------------------------------------------
    def picks_json(self, profile: str) -> Path:
        return self.picks_dir / profile / "picks.json"

    def session_picks_json(self, profile: str) -> Path:
        return self.sessions_dir / profile / "picks.session.json"

    def layer_json(self, profile: str) -> Path:
        return self.reports_dir / profile / "layer_analysis.json"

    def layer_session_json(self, profile: str) -> Path:
        return self.sessions_dir / profile / "layer_analysis.session.json"

    def plots_dir_for(self, profile: str) -> Path:
        return self.plots_dir / profile

    def reports_dir_for(self, profile: str) -> Path:
        return self.reports_dir / profile

    def picks_dir_for(self, profile: str) -> Path:
        return self.picks_dir / profile

    # -- bookkeeping ------------------------------------------------------
    def all_dirs(self) -> list[Path]:
        return [
            self.processing_dir, self.cache_dir, self.picks_dir, self.sessions_dir,
            self.results_dir, self.plots_dir, self.reports_dir, self.velocity_dir,
        ]

    def ensure(self) -> "ProjectPaths":
        """Create every processing/results subfolder if missing. Idempotent."""
        for d in self.all_dirs():
            d.mkdir(parents=True, exist_ok=True)
        return self
    
    def ensure_dir(self, path: Path) -> Path:
        path.mkdir(parents=True, exist_ok=True)
        return path

    def relative(self, path: Path) -> str:
        """Best-effort path for display/logging, relative to the project root."""
        try:
            return str(Path(path).resolve().relative_to(self.root))
        except Exception:
            return str(path)


# -----------------------------------------------------------------------------
# Active project (the single "currently open" project for the running app)
# -----------------------------------------------------------------------------

_active_project: ProjectPaths | None = None


def set_active_project(root: Path | str) -> ProjectPaths:
    """Make `root` the active project, creating its subfolders if needed."""
    global _active_project
    _active_project = ProjectPaths(Path(root)).ensure()
    return _active_project


def get_active_project() -> ProjectPaths | None:
    """Return the active project's paths, or None if nothing is open."""
    return _active_project


def require_active_project() -> ProjectPaths:
    """Like `get_active_project`, but raises if no project is open.

    Used deep inside I/O helpers (picks, exports, ...) that need somewhere
    to read/write but don't want to receive a path on every call.
    """
    if _active_project is None:
        raise RuntimeError(
            "No active project. Create or open a project (File > New "
            "Project / Open Project) before running this operation."
        )
    return _active_project


def clear_active_project() -> None:
    global _active_project
    _active_project = None


__all__ = [
    "APP_ROOT",
    "PROJECTS_ROOT",
    "LEGACY_DATA_DIR",
    "GEOM_TEMPLATES_DIR",
    "LEGACY_METADATA_DIR",
    "LEGACY_RAW_DIR",
    "LEGACY_EXTERNAL_DIR",
    "LEGACY_EXAMPLES_DIR",
    "LEGACY_IMPORT_PICKS_DIR",
    "PROJECT_FILE_NAME",
    "ProjectPaths",
    "sanitize_project_name",
    "set_active_project",
    "get_active_project",
    "require_active_project",
    "clear_active_project",
]