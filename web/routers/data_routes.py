from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel

from src.common import paths as core_paths
from src.io.pick_reader import load_picks_json, load_session_picks_json, save_session_picks_json
from src.io.seg2_reader import read_seg2
from src.utils.geometry import discover_profile_folders
from web.auth import get_current_user, user_role_for_project
from web.db import role_at_least
from web.routers.projects_routes import _open_authorized

router = APIRouter(prefix="/api/projects/{project_id}", tags=["data"])


def _require_role(project_id: str, user: dict, minimum: str):
    project, role = _open_authorized(project_id, user)
    if not role_at_least(role, minimum):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=f"'{minimum}' access or higher required.",
        )
    return project, role


@router.get("/profiles")
def list_profiles(project_id: str, user: dict = Depends(get_current_user)):
    project, _role = _require_role(project_id, user, "viewer")
    if project.raw_folder is None:
        return {"profiles": list(project.profiles), "source": "registered"}
    found = discover_profile_folders(project.raw_folder)
    return {"profiles": found, "source": "raw_folder"}


@router.get("/profiles/{profile}/shots")
def list_shots(project_id: str, profile: str, user: dict = Depends(get_current_user)):
    project, _role = _require_role(project_id, user, "viewer")
    if project.raw_folder is None:
        raise HTTPException(status_code=400, detail="Project has no raw_folder set.")
    folder = project.raw_folder / profile
    if not folder.exists():
        raise HTTPException(status_code=404, detail=f"No such profile folder: {folder}")
    files = sorted(list(folder.glob("*.seg2")) + list(folder.glob("*.SEG2")))
    return {
        "profile": profile,
        "shots": [{"index": i, "filename": f.name} for i, f in enumerate(files)],
    }


class TraceGatherOut(BaseModel):
    profile: str
    shot_index: int
    filename: str
    dt_s: float
    n_traces: int
    n_samples: int
    delay_ms: float
    shot_pos_m: float | None
    ffid: int
    receiver_positions_m: list[float] | None
    traces: list[list[float]]
    picks: dict[str, float]        # trace index (1-based, str) -> pick time (ms)


@router.get("/profiles/{profile}/shots/{shot_index}", response_model=TraceGatherOut)
def get_shot_gather(project_id: str, profile: str, shot_index: int, user: dict = Depends(get_current_user)):
    """Full trace data for one shot gather, ready to plot in the browser.

    This is deliberately the same `read_seg2` used by the desktop app and
    CLI - one implementation of "read a SEG2 file", three front-ends.
    """
    project, _role = _require_role(project_id, user, "viewer")
    if project.raw_folder is None:
        raise HTTPException(status_code=400, detail="Project has no raw_folder set.")

    folder = project.raw_folder / profile
    files = sorted(list(folder.glob("*.seg2")) + list(folder.glob("*.SEG2")))
    if not (0 <= shot_index < len(files)):
        raise HTTPException(status_code=404, detail=f"No shot at index {shot_index} in profile '{profile}'.")
    path = files[shot_index]

    try:
        data, dt_s, n_traces, n_samples, shot_pos, ffid, delay_ms, recv_locs_m = read_seg2(path)
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Could not read {path.name}: {exc}") from exc

    # Existing picks for this shot, session first (in-progress work), then
    # finalized picks.json - same precedence the desktop app uses.
    project.activate()
    session_picks = load_session_picks_json(profile) or {}
    final_picks = load_picks_json(profile) or {}
    shot_key = str(ffid or (shot_index + 1))
    picks_for_shot = session_picks.get(shot_key) or final_picks.get(shot_key) or {}

    return TraceGatherOut(
        profile=profile,
        shot_index=shot_index,
        filename=path.name,
        dt_s=dt_s,
        n_traces=int(n_traces),
        n_samples=int(n_samples),
        delay_ms=float(delay_ms),
        shot_pos_m=(float(shot_pos) if shot_pos is not None else None),
        ffid=int(ffid),
        receiver_positions_m=(recv_locs_m.tolist() if recv_locs_m is not None else None),
        traces=data.tolist(),
        picks={str(k): float(v) for k, v in picks_for_shot.items()},
    )


class SavePicksRequest(BaseModel):
    picks: dict[str, float]   # trace index (1-based, str) -> pick time (ms)


@router.post("/profiles/{profile}/shots/{shot_index}/picks")
def save_shot_picks(
    project_id: str, profile: str, shot_index: int,
    body: SavePicksRequest, user: dict = Depends(get_current_user),
):
    project, _role = _require_role(project_id, user, "editor")
    if project.raw_folder is None:
        raise HTTPException(status_code=400, detail="Project has no raw_folder set.")

    folder = project.raw_folder / profile
    files = sorted(list(folder.glob("*.seg2")) + list(folder.glob("*.SEG2")))
    if not (0 <= shot_index < len(files)):
        raise HTTPException(status_code=404, detail=f"No shot at index {shot_index} in profile '{profile}'.")

    try:
        _, _, _, _, _, ffid, _, _ = read_seg2(files[shot_index])
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Could not read {files[shot_index].name}: {exc}") from exc

    project.activate()
    all_picks = load_session_picks_json(profile) or {}
    shot_key = str(ffid or (shot_index + 1))
    all_picks[shot_key] = dict(body.picks)
    save_session_picks_json(profile, all_picks)
    project.register_profile(profile)

    from src.io import project_io as pio
    pio.save_project(project)

    return {"ok": True, "shot": shot_key, "n_picks": len(body.picks)}


__all__ = ["router"]