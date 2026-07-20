"""
Project creation, discovery and persistence.

A "project" is a single field survey the user is working on. Its on-disk
footprint is a folder (anywhere the user likes, `PROJECTS_ROOT` by default)
containing a `project.json` manifest plus the `processing/` and `results/`
working folders described by `src.common.paths.ProjectPaths`.

Raw SEG2 data, geometry/coordinate files and field-report metadata are NOT
copied into the project - they usually live on a field laptop or network
share and can be large. The project only remembers *where* they are, so
the picker/exporters can find them again next time the project is opened.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from src.common.paths import (
    PROJECT_FILE_NAME,
    PROJECTS_ROOT,
    ProjectPaths,
    sanitize_project_name,
    set_active_project,
)


class ProjectError(Exception):
    """Base class for project-related errors."""


class ProjectExistsError(ProjectError):
    """Raised when creating a project that already exists."""


class ProjectNotFoundError(ProjectError):
    """Raised when opening a project that cannot be found or read."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _to_path(value) -> Path | None:
    return Path(value) if value else None


def _to_str(path: Path | None) -> str | None:
    return str(path) if path is not None else None


@dataclass
class Project:
    """In-memory representation of `project.json` plus its resolved paths."""

    name: str
    root: Path
    raw_folder: Path | None = None          # folder of SEG2 subfolders, one per profile
    geometry_folder: Path | None = None     # coordinate/geometry workbooks & shapefiles
    metadata_folder: Path | None = None     # field-report workbooks
    profiles: list[str] = field(default_factory=list)
    notes: str = ""
    created: str = field(default_factory=_now)
    modified: str = field(default_factory=_now)

    def __post_init__(self) -> None:
        self.root = Path(self.root).resolve()

    @property
    def paths(self) -> ProjectPaths:
        return ProjectPaths(self.root)

    # -- (de)serialization ------------------------------------------------
    def to_json_dict(self) -> dict:
        return {
            "name": self.name,
            "raw_folder": _to_str(self.raw_folder),
            "geometry_folder": _to_str(self.geometry_folder),
            "metadata_folder": _to_str(self.metadata_folder),
            "profiles": list(self.profiles),
            "notes": self.notes,
            "created": self.created,
            "modified": self.modified,
        }

    @classmethod
    def from_json_dict(cls, data: dict, root: Path) -> "Project":
        return cls(
            name=data.get("name") or Path(root).name,
            root=Path(root),
            raw_folder=_to_path(data.get("raw_folder")),
            geometry_folder=_to_path(data.get("geometry_folder")),
            metadata_folder=_to_path(data.get("metadata_folder")),
            profiles=list(data.get("profiles") or []),
            notes=data.get("notes", ""),
            created=data.get("created") or _now(),
            modified=data.get("modified") or _now(),
        )

    # -- convenience --------------------------------------------------------
    def register_profile(self, profile_name: str) -> None:
        if profile_name not in self.profiles:
            self.profiles.append(profile_name)
            self.profiles.sort()

    def activate(self) -> ProjectPaths:
        """Make this project the app-wide active project (see paths.py)."""
        return set_active_project(self.root)


# -----------------------------------------------------------------------------
# create / open / save / list
# -----------------------------------------------------------------------------

def create_project(
    name: str,
    parent_dir: Path | str = PROJECTS_ROOT,
    raw_folder: Path | str | None = None,
    geometry_folder: Path | str | None = None,
    metadata_folder: Path | str | None = None,
) -> Project:
    """Create a new project folder with the standard layout and manifest.

    Parameters
    ----------
    name : str
        Human-readable project name; also used as the folder name.
    parent_dir : Path
        Where the new project folder is created (default: PROJECTS_ROOT).
    raw_folder, geometry_folder, metadata_folder : Path, optional
        External folders the project should remember. If given, they must
        already exist (nothing is copied).

    Raises
    ------
    ProjectExistsError
        If `parent_dir/name` already contains a project.json.
    FileNotFoundError
        If one of the optional external folders does not exist.
    ValueError
        If `name` is empty or contains unsafe characters.
    """
    clean_name = sanitize_project_name(name)
    root = Path(parent_dir).resolve() / clean_name

    if (root / PROJECT_FILE_NAME).exists():
        raise ProjectExistsError(f"A project named '{clean_name}' already exists at {root}")

    for label, p in (
        ("raw_folder", raw_folder),
        ("geometry_folder", geometry_folder),
        ("metadata_folder", metadata_folder),
    ):
        if p is not None and not Path(p).exists():
            raise FileNotFoundError(f"{label} '{p}' does not exist.")

    project = Project(
        name=clean_name,
        root=root,
        raw_folder=Path(raw_folder).resolve() if raw_folder else None,
        geometry_folder=Path(geometry_folder).resolve() if geometry_folder else None,
        metadata_folder=Path(metadata_folder).resolve() if metadata_folder else None,
    )
    project.paths.ensure()
    save_project(project)
    return project


def open_project(path: Path | str) -> Project:
    """Open an existing project from its root folder (or its project.json)."""
    p = Path(path).resolve()
    project_file = p if p.name == PROJECT_FILE_NAME else p / PROJECT_FILE_NAME
    root = project_file.parent

    if not project_file.exists():
        raise ProjectNotFoundError(f"No '{PROJECT_FILE_NAME}' found at {root}")

    try:
        data = json.loads(project_file.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ProjectError(f"Could not read project file '{project_file}': {exc}") from exc

    project = Project.from_json_dict(data, root)
    project.paths.ensure()  # tolerate manually-deleted subfolders
    return project


def save_project(project: Project) -> Path:
    """Write project.json for `project`, updating the modified timestamp."""
    project.modified = _now()
    project.paths.root.mkdir(parents=True, exist_ok=True)
    project.paths.project_file.write_text(
        json.dumps(project.to_json_dict(), indent=2), encoding="utf-8"
    )
    return project.paths.project_file


def is_project_dir(path: Path | str) -> bool:
    """True if `path` is a project root (or a project.json file itself)."""
    p = Path(path)
    if p.is_dir():
        return (p / PROJECT_FILE_NAME).exists()
    return p.name == PROJECT_FILE_NAME and p.exists()


def list_projects(parent_dir: Path | str = PROJECTS_ROOT) -> list[Project]:
    """Discover every project folder directly under `parent_dir`.

    Unreadable project folders (corrupt project.json, etc.) are skipped
    rather than raising, so one bad project doesn't break project listing.
    """
    parent = Path(parent_dir)
    if not parent.exists():
        return []
    found: list[Project] = []
    for child in sorted(parent.iterdir()):
        if child.is_dir() and (child / PROJECT_FILE_NAME).exists():
            try:
                found.append(open_project(child))
            except ProjectError:
                continue
    return found


def delete_project(project: Project, *, remove_files: bool = False) -> None:
    """Forget a project. Only deletes files on disk if `remove_files=True`."""
    if remove_files:
        import shutil
        shutil.rmtree(project.root, ignore_errors=True)


__all__ = [
    "Project",
    "ProjectError",
    "ProjectExistsError",
    "ProjectNotFoundError",
    "create_project",
    "open_project",
    "save_project",
    "is_project_dir",
    "list_projects",
    "delete_project",
]