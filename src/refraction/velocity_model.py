"""Layer-velocity averaging and travel-time prediction from fitted segments."""
from __future__ import annotations

import numpy as np

from src.refraction.travel_times import depth_2layer, depth_3layer


def compute_layer_averages(layer_results: dict, shot_pos_by_id: dict) -> dict:
    """
    Compute layer averages from off-end and center shots.

    Returns
    -------
    {
      "V0": {"off":..,"center":..,"avg":..},
      "V1": {"off":..,"center":..,"avg":..},
      "V2": {"off":..,"center":..,"avg":..},
      "ti1_ms": {"off":..,"center":..,"avg":..},
      "ti2_ms": {"off":..,"center":..,"avg":..},
      "h1_m": ..., "h2_m": ...,
      "off_shots": [...], "center_shots": [...]
    }
    """
    if not layer_results:
        return {}

    shot_ids = [sid for sid in shot_pos_by_id if sid in layer_results]
    if not shot_ids:
        return {}

    shot_ids_sorted = sorted(shot_ids, key=lambda sid: float(shot_pos_by_id.get(sid, sid)))
    if len(shot_ids_sorted) >= 2:
        off_shots = [shot_ids_sorted[0], shot_ids_sorted[-1]]
    else:
        off_shots = [shot_ids_sorted[0]]

    center_candidates = [sid for sid in shot_ids_sorted if sid not in off_shots]
    if not center_candidates:
        center_candidates = [sid for sid in shot_ids_sorted if sid in off_shots]

    def _collect(shots: list, key: str) -> list:
        vals: list = []
        for sid in shots:
            for side_res in (layer_results.get(sid, {}) or {}).values():
                v = float(side_res.get(key, 0.0))
                if v > 0.0:
                    vals.append(v)
        return vals

    def _mean(vals: list) -> float:
        return float(np.mean(vals)) if vals else 0.0

    def _avg_two(off_v: float, cen_v: float) -> float:
        picks = [v for v in (off_v, cen_v) if v > 0.0]
        return float(np.mean(picks)) if picks else 0.0

    out: dict = {"off_shots": off_shots, "center_shots": center_candidates}
    for key in ("V0_m_s", "V1_m_s", "V2_m_s", "ti1_ms", "ti2_ms"):
        off_v = _mean(_collect(off_shots, key))
        cen_v = _mean(_collect(center_candidates, key))
        avg_v = _avg_two(off_v, cen_v)
        short = key.replace("_m_s", "").replace("_ms", "_ms")
        out[short] = {"off": off_v, "center": cen_v, "avg": avg_v}

    V0 = out.get("V0", {}).get("avg", 0.0)
    V1 = out.get("V1", {}).get("avg", 0.0)
    V2 = out.get("V2", {}).get("avg", 0.0)
    ti1 = out.get("ti1_ms", {}).get("avg", 0.0)
    ti2 = out.get("ti2_ms", {}).get("avg", 0.0)

    h1 = depth_2layer(ti1, V0, V1) if (V0 > 0 and V1 > 0 and ti1 > 0) else None
    h2 = depth_3layer(ti2, V0, V1, V2, h1) if (h1 and V2 > 0 and ti2 > 0) else None
    out["h1_m"] = float(h1) if h1 else 0.0
    out["h2_m"] = float(h2) if h2 else 0.0
    return out

def _predict_time_from_fit(x_abs: float, fit_res: dict) -> float | None:
    """Predict t(ms) at corrected true offset x_abs from fitted layer windows."""
    windows = list((fit_res or {}).get("windows", []) or [])
    windows += [None] * (3 - len(windows))
    segs = list((fit_res or {}).get("segments", []) or [])
    segs += [{}] * (3 - len(segs))

    x = float(x_abs)
    for i in range(3):
        w = windows[i]
        seg = segs[i]
        if w is None:
            continue
        x0, x1 = float(w[0]), float(w[1])
        if not (x0 <= x <= x1):
            continue
        sl = float(seg.get("slope_ms_m", 0.0))
        ic = float(seg.get("intercept_ms", 0.0))
        if sl <= 0.0:
            return None
        return float(sl * x + ic)

    return None

def _choose_shot_fit(layer_by_side: dict) -> tuple[str | None, dict | None]:
    """
    Choose one representative side-fit per shot for summary export.

    Priority is:
    1) "ALL" side if present,
    2) side with the largest number of fitted points.
    """
    if not layer_by_side:
        return None, None
    if "ALL" in layer_by_side:
        return "ALL", layer_by_side.get("ALL")

    best_side = None
    best_res = None
    best_n = -1
    for side, res in layer_by_side.items():
        n = int((res or {}).get("n_points", 0))
        if n > best_n:
            best_n = n
            best_side = side
            best_res = res
    return best_side, best_res
