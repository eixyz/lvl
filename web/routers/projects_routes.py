from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel

from src.common.paths import PROJECTS_ROOT
from src.io import project_io as pio
from web.auth import get_current_user, grant_project_access, projects_for_user, user_role_for_project
from web.util import decode_project_id, encode_project_id

router = APIRouter(prefix="/api/projects", tags=["projects"])


class ProjectOut(BaseModel):
    id: str
    name: str
    role: str
    raw_folder: str | None
    geometry_folder: str | None
    metadata_folder: str | None
    profiles: list[str]
    created: str
    modified: str

    @classmethod
    def from_project(cls, project: pio.Project, role: str) -> "ProjectOut":
        return cls(
            id=encode_project_id(project.root),
            name=project.name,
            role=role,
            raw_folder=str(project.raw_folder) if project.raw_folder else None,
            geometry_folder=str(project.geometry_folder) if project.geometry_folder else None,
            metadata_folder=str(project.metadata_folder) if project.metadata_folder else None,
            profiles=project.profiles,
            created=project.created,
            modified=project.modified,
        )


class CreateProjectRequest(BaseModel):
    name: str
    raw_folder: str | None = None
    geometry_folder: str | None = None
    metadata_folder: str | None = None


def _open_authorized(project_id: str, user: dict) -> tuple[pio.Project, str]:
    """Decode a project_id, load it, and enforce the caller has access.

    Every route that touches a specific project should go through this,
    not `pio.open_project` directly - it's the only place access control
    is actually enforced.
    """
    try:
        root = decode_project_id(project_id)
    except Exception:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Invalid project id")

    role = user_role_for_project(user["id"], root)
    if role is None:
        # Same 404 whether the project doesn't exist or the user just can't
        # see it - don't leak which projects exist on disk to other users.
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Project not found")

    try:
        project = pio.open_project(root)
    except pio.ProjectError as exc:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(exc)) from exc

    return project, role


@router.get("", response_model=list[ProjectOut])
def list_projects(user: dict = Depends(get_current_user)):
    out = []
    for row in projects_for_user(user["id"]):
        try:
            project = pio.open_project(row["project_root"])
        except pio.ProjectError:
            continue  # project folder moved/deleted on disk; skip rather than 500
        out.append(ProjectOut.from_project(project, row["role"]))
    return out


@router.post("", response_model=ProjectOut, status_code=status.HTTP_201_CREATED)
def create_project_route(body: CreateProjectRequest, user: dict = Depends(get_current_user)):
    try:
        project = pio.create_project(
            body.name,
            parent_dir=PROJECTS_ROOT,
            raw_folder=Path(body.raw_folder) if body.raw_folder else None,
            geometry_folder=Path(body.geometry_folder) if body.geometry_folder else None,
            metadata_folder=Path(body.metadata_folder) if body.metadata_folder else None,
        )
    except pio.ProjectExistsError as exc:
        raise HTTPException(status_code=status.HTTP_409_CONFLICT, detail=str(exc)) from exc
    except (ValueError, FileNotFoundError) as exc:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(exc)) from exc

    grant_project_access(user["id"], project.root, project.name, role="owner")
    return ProjectOut.from_project(project, "owner")


@router.get("/{project_id}", response_model=ProjectOut)
def get_project(project_id: str, user: dict = Depends(get_current_user)):
    project, role = _open_authorized(project_id, user)
    return ProjectOut.from_project(project, role)


__all__ = ["router"]