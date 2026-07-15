"""Perpendicular-offset / inline-shift configuration resolution helpers."""
from __future__ import annotations


def parse_shot_value_map(text: str | None) -> dict:
    """
    Parse mappings like "1:0,2:3.5,3:0" into {1: 0.0, 2: 3.5, 3: 0.0}.
    """
    out: dict = {}
    if not text:
        return out
    for chunk in str(text).split(","):
        part = chunk.strip()
        if not part:
            continue
        if ":" not in part:
            raise ValueError(f"Invalid mapping '{part}'. Expected shot:value")
        k_s, v_s = part.split(":", 1)
        shot_id = int(k_s.strip())
        val = float(v_s.strip().replace(",", "."))
        out[shot_id] = val
    return out

def resolve_perp_by_shot(cfg: dict, shot_ids: list,
                         override_map: dict | None = None) -> dict:
    """
    Resolve perpendicular offsets per shot from profile config + CLI overrides.
    cfg['perp_m'] can be either scalar or dict {shot_id: perp_m}.
    """
    base = cfg.get("perp_m", 0.0)
    per_shot: dict = {}

    if isinstance(base, dict):
        for sid in shot_ids:
            per_shot[int(sid)] = float(base.get(int(sid), 0.0))
    else:
        base_val = float(base)
        for sid in shot_ids:
            per_shot[int(sid)] = base_val

    for sid, val in (override_map or {}).items():
        per_shot[int(sid)] = float(val)

    return per_shot
