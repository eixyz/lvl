"""Load/save first-break picks and layer-analysis results (JSON persistence)."""
from __future__ import annotations

import json
from pathlib import Path

from src.common.paths import require_active_project

def _project():
    return require_active_project()

def _picks_json_path(profile: str) -> Path:
    return _project().picks_json(profile)

def _session_picks_json_path(profile: str) -> Path:

    return _project().session_picks_json(profile)

def _layer_json_path(profile: str) -> Path:
    return _project().layer_json(profile)

def _layer_session_json_path(profile: str) -> Path:
    return _project().layer_session_json(profile)

def _coerce_layer_results(raw: dict) -> dict:
    """Normalize JSON-loaded layer results into int-shot keyed dict form."""
    out: dict = {}
    for sid, side_map in (raw or {}).items():
        try:
            shot_id = int(sid)
        except Exception:
            continue
        if not isinstance(side_map, dict):
            continue
        out[shot_id] = {str(side): dict(payload) for side, payload in side_map.items()
                        if isinstance(payload, dict)}
    return out

def load_layer_json(profile: str) -> dict:
    p = _layer_json_path(profile)
    if not p.exists():
        return {}
    try:
        with open(p, encoding="utf-8") as fh:
            raw = json.load(fh)
        return _coerce_layer_results(raw)
    except Exception as exc:
        print(f"  [WARN] Could not load layer_analysis.json: {exc}")
        return {}

def load_layer_session_json(profile: str) -> dict:
    p = _layer_session_json_path(profile)
    if not p.exists():
        return {}
    try:
        with open(p, encoding="utf-8") as fh:
            raw = json.load(fh)
        return _coerce_layer_results(raw)
    except Exception as exc:
        print(f"  [WARN] Could not load layer_analysis.session.json: {exc}")
        return {}

def save_layer_json(profile: str, layer_results: dict):
    """Persist final layer analysis results to layer_analysis.json."""
    p = _layer_json_path(profile)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump({str(k): v for k, v in (layer_results or {}).items()}, fh, indent=2)
    print(f"     Layer analysis saved -> {_project().relative(p)}")

def save_layer_session_json(profile: str, layer_results: dict):
    """Persist in-progress layer analysis to layer_analysis.session.json."""
    p = _layer_session_json_path(profile)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump({str(k): v for k, v in (layer_results or {}).items()}, fh, indent=2)
    print(f"     Layer session saved -> {_project().relative(p)}")

def clear_layer_session_json(profile: str):
    p = _layer_session_json_path(profile)
    try:
        if p.exists():
            p.unlink()
            print(f"     Layer session file cleared -> {_project().relative(p)}")
    except Exception as exc:
        print(f"  [WARN] Could not clear layer session file: {exc}")

def load_picks_json(profile: str) -> dict:
    """Load raw picks: {shot_id: {trace_idx: ms}}."""
    p = _picks_json_path(profile)
    if not p.exists():
        return {}
    try:
        with open(p, encoding="utf-8") as fh:
            raw = json.load(fh)
        return {int(k): {int(ti): float(tv) for ti, tv in v.items()}
                for k, v in raw.items()}
    except Exception as exc:
        print(f"  [WARN] Could not load picks.json: {exc}")
        return {}

def load_session_picks_json(profile: str) -> dict:
    """Load session picks: {shot_id: {trace_idx: ms}}."""
    p = _session_picks_json_path(profile)
    if not p.exists():
        return {}
    try:
        with open(p, encoding="utf-8") as fh:
            raw = json.load(fh)
        return {int(k): {int(ti): float(tv) for ti, tv in v.items()}
                for k, v in raw.items()}
    except Exception as exc:
        print(f"  [WARN] Could not load picks.session.json: {exc}")
        return {}

def save_picks_json(profile: str, all_picks: dict):
    """Persist finalized raw picks to picks.json."""
    p = _picks_json_path(profile)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump({str(k): {str(ti): tv for ti, tv in v.items()}
                   for k, v in all_picks.items()},
                  fh, indent=2)
    print(f"     Picks saved -> {_project().relative(p)}")

def save_session_picks_json(profile: str, all_picks: dict):
    """Persist in-progress raw picks to picks.session.json."""
    p = _session_picks_json_path(profile)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump({str(k): {str(ti): tv for ti, tv in v.items()}
                   for k, v in all_picks.items()},
                  fh, indent=2)
    print(f"     Session picks saved -> {_project().relative(p)}")

def clear_session_picks_json(profile: str):
    p = _session_picks_json_path(profile)
    try:
        if p.exists():
            p.unlink()
            print(f"     Session file cleared -> {_project().relative(p)}")
    except Exception as exc:
        print(f"  [WARN] Could not clear session file: {exc}")
