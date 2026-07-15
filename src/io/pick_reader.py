"""Load/save first-break picks and layer-analysis results (JSON persistence)."""
from __future__ import annotations

import json
from pathlib import Path

from src.common.paths import OUTPUT_DIR, PROJECT_DIR


def _picks_json_path(profile_name: str) -> Path:
    return OUTPUT_DIR / profile_name / "picks.json"

def _session_picks_json_path(profile_name: str) -> Path:
    return OUTPUT_DIR / profile_name / "picks.session.json"

def _layer_json_path(profile_name: str) -> Path:
    return OUTPUT_DIR / profile_name / "layer_analysis.json"

def _layer_session_json_path(profile_name: str) -> Path:
    return OUTPUT_DIR / profile_name / "layer_analysis.session.json"

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

def load_layer_json(profile_name: str) -> dict:
    p = _layer_json_path(profile_name)
    if not p.exists():
        return {}
    try:
        with open(p, encoding="utf-8") as fh:
            raw = json.load(fh)
        return _coerce_layer_results(raw)
    except Exception as exc:
        print(f"  [WARN] Could not load layer_analysis.json: {exc}")
        return {}

def load_layer_session_json(profile_name: str) -> dict:
    p = _layer_session_json_path(profile_name)
    if not p.exists():
        return {}
    try:
        with open(p, encoding="utf-8") as fh:
            raw = json.load(fh)
        return _coerce_layer_results(raw)
    except Exception as exc:
        print(f"  [WARN] Could not load layer_analysis.session.json: {exc}")
        return {}

def save_layer_json(profile_name: str, layer_results: dict):
    """Persist final layer analysis results to layer_analysis.json."""
    p = _layer_json_path(profile_name)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump({str(k): v for k, v in (layer_results or {}).items()}, fh, indent=2)
    print(f"     Layer analysis saved -> {p.relative_to(PROJECT_DIR)}")

def save_layer_session_json(profile_name: str, layer_results: dict):
    """Persist in-progress layer analysis to layer_analysis.session.json."""
    p = _layer_session_json_path(profile_name)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump({str(k): v for k, v in (layer_results or {}).items()}, fh, indent=2)
    print(f"     Layer session saved -> {p.relative_to(PROJECT_DIR)}")

def clear_layer_session_json(profile_name: str):
    p = _layer_session_json_path(profile_name)
    try:
        if p.exists():
            p.unlink()
            print(f"     Layer session file cleared -> {p.relative_to(PROJECT_DIR)}")
    except Exception as exc:
        print(f"  [WARN] Could not clear layer session file: {exc}")

def load_picks_json(profile_name: str) -> dict:
    """Load raw picks: {shot_id: {trace_idx: ms}}."""
    p = _picks_json_path(profile_name)
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

def load_session_picks_json(profile_name: str) -> dict:
    """Load session picks: {shot_id: {trace_idx: ms}}."""
    p = _session_picks_json_path(profile_name)
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

def save_picks_json(profile_name: str, all_picks: dict):
    """Persist finalized raw picks to picks.json."""
    p = _picks_json_path(profile_name)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump({str(k): {str(ti): tv for ti, tv in v.items()}
                   for k, v in all_picks.items()},
                  fh, indent=2)
    print(f"     Picks saved -> {p.relative_to(PROJECT_DIR)}")

def save_session_picks_json(profile_name: str, all_picks: dict):
    """Persist in-progress raw picks to picks.session.json."""
    p = _session_picks_json_path(profile_name)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump({str(k): {str(ti): tv for ti, tv in v.items()}
                   for k, v in all_picks.items()},
                  fh, indent=2)
    print(f"     Session picks saved -> {p.relative_to(PROJECT_DIR)}")

def clear_session_picks_json(profile_name: str):
    p = _session_picks_json_path(profile_name)
    try:
        if p.exists():
            p.unlink()
            print(f"     Session file cleared -> {p.relative_to(PROJECT_DIR)}")
    except Exception as exc:
        print(f"  [WARN] Could not clear session file: {exc}")
