"""Corrected-offset data build and layer-window fitting/derived-metric helpers."""
from __future__ import annotations

from typing import Any

import numpy as np

from src.utils.math_utils import true_offset
from src.refraction.travel_times import fit_line, depth_2layer, depth_3layer
from src.refraction.velocity_model import _predict_time_from_fit, _choose_shot_fit
from src.picker.preprocessing import apply_bulk_static


def build_corrected_pick_data(shots_info: list, all_picks: dict,
                              recv_positions: Any,
                              perp_by_shot: dict | None = None,
                              inline_shift_by_shot: dict | None = None) -> dict:
    """
    Build per-shot corrected pick rows with true offset and corrected travel-time.

    Note
    ----
    Perpendicular/inline-shift corrections change the x-domain only. Pick times stay
    on their measured (bulk-corrected) basis, so fitting uses consistent (x_true, t).

    Returns
    -------
    {shot_id: [
        {
          trace_idx, trace_no, recv_pos_m, shot_pos_m,
          inline_signed_m, inline_abs_m, perp_m, true_off_m,
                    fb_raw_ms, fb_bulk_ms, fb_interp_inline_ms, side
        }, ...
    ]}
    """
    corrected: dict = {}
    for shot_id, shot_pos_m in shots_info:
        raw_picks = all_picks.get(shot_id, {})
        if not raw_picks:
            corrected[int(shot_id)] = []
            continue

        bulk_picks = apply_bulk_static(raw_picks)
        po = float((perp_by_shot or {}).get(int(shot_id), 0.0))
        x_shift = float((inline_shift_by_shot or {}).get(int(shot_id), 0.0))
        rows: list = []

        for trace_idx in sorted(raw_picks):
            if trace_idx >= len(recv_positions):
                continue
            recv_pos = float(recv_positions[trace_idx])
            inline_signed = recv_pos - float(shot_pos_m)
            inline_abs = abs(inline_signed)
            inline_corr = max(0.0, float(inline_abs + x_shift))
            true_off = float(true_offset(inline_corr, po))
            side = "L" if inline_signed < 0.0 else "R"

            rows.append({
                "trace_idx": int(trace_idx),
                "trace_no": int(trace_idx) + 1,
                "recv_pos_m": recv_pos,
                "shot_pos_m": float(shot_pos_m),
                "inline_signed_m": float(inline_signed),
                "inline_abs_m": float(inline_abs),
                "inline_shift_m": float(x_shift),
                "inline_corr_m": float(inline_corr),
                "perp_m": po,
                "true_off_m": true_off,
                "fb_raw_ms": float(raw_picks[trace_idx]),
                "fb_bulk_ms": float(bulk_picks.get(trace_idx, raw_picks[trace_idx])),
                "fb_interp_inline_ms": float(bulk_picks.get(trace_idx, raw_picks[trace_idx])),
                "fb_interp_geom_ms": float(bulk_picks.get(trace_idx, raw_picks[trace_idx])),
                "side": side,
            })

        corrected[int(shot_id)] = sorted(rows, key=lambda r: r["trace_idx"])

    return corrected

def _fit_layers_from_windows(x_vals: Any, t_vals: Any, windows: list) -> dict:
    """
    Fit up to three travel-time line segments from picked x-windows.

    Parameters
    ----------
    x_vals : Any
        Corrected true offsets (m).
    t_vals : Any
        Corrected first-break times (ms).
    windows : list
        Layer windows as [(x1, x2), (x3, x4), (x5, x6)] where each entry
        may be None if a layer was not selected.

    Returns
    -------
    dict
        Fitted segment metrics plus derived velocities/intercept times/depths.
    """
    x = np.asarray(x_vals, dtype=float)
    t = np.asarray(t_vals, dtype=float)

    # Handle cases where the shallowest visible layer is not picked (e.g. direct to L2):
    # shift non-empty windows to the left so fitting always starts at first available layer.
    wins_in = list((windows or []))[:3]
    wins_in += [None] * (3 - len(wins_in))
    normalized_windows = [w for w in wins_in if w is not None]
    normalized_windows += [None] * (3 - len(normalized_windows))

    segs: list = []
    for layer_idx in range(3):
        win = normalized_windows[layer_idx] if layer_idx < len(normalized_windows) else None
        if not win:
            segs.append({
                "layer": layer_idx + 1,
                "x0": 0.0,
                "x1": 0.0,
                "velocity_m_s": 0.0,
                "intercept_ms": 0.0,
                "slope_ms_m": 0.0,
                "r2": 0.0,
                "n": 0,
            })
            continue

        x0, x1 = float(win[0]), float(win[1])
        mask = (x >= x0) & (x <= x1)
        if int(mask.sum()) < 2:
            segs.append({
                "layer": layer_idx + 1,
                "x0": x0,
                "x1": x1,
                "velocity_m_s": 0.0,
                "intercept_ms": 0.0,
                "slope_ms_m": 0.0,
                "r2": 0.0,
                "n": int(mask.sum()),
            })
            continue

        slope, intercept, r2 = fit_line(x[mask], t[mask])
        vel = (1000.0 / slope) if slope > 1e-9 else 0.0
        segs.append({
            "layer": layer_idx + 1,
            "x0": x0,
            "x1": x1,
            "velocity_m_s": float(max(0.0, vel)),
            "intercept_ms": float(intercept),
            "slope_ms_m": float(slope),
            "r2": float(r2),
            "n": int(mask.sum()),
        })

    v0 = float(segs[0]["velocity_m_s"])
    v1 = float(segs[1]["velocity_m_s"])
    v2 = float(segs[2]["velocity_m_s"])
    ti1 = float(segs[1]["intercept_ms"]) if v1 > 0.0 else 0.0
    ti2 = float(segs[2]["intercept_ms"]) if v2 > 0.0 else 0.0

    h1 = depth_2layer(ti1, v0, v1) if (v0 > 0.0 and v1 > 0.0 and ti1 > 0.0) else None
    h2 = depth_3layer(ti2, v0, v1, v2, h1) if (h1 and v2 > 0.0 and ti2 > 0.0) else None

    return {
        "segments": segs,
        "V0_m_s": v0,
        "V1_m_s": v1,
        "V2_m_s": v2,
        "ti1_ms": ti1,
        "ti2_ms": ti2,
        "h1_m": float(h1) if h1 else 0.0,
        "h2_m": float(h2) if h2 else 0.0,
    }

def _compute_fit_rms(x_vals: Any, t_vals: Any, fit_res: dict) -> float:
    """Compute RMS(ms) between observed and fitted times for given x-values."""
    x = np.asarray(x_vals, dtype=float)
    t = np.asarray(t_vals, dtype=float)
    if x.size == 0 or t.size == 0:
        return 0.0
    pred: list = []
    obs: list = []
    for xv, tv in zip(x, t):
        tp = _predict_time_from_fit(float(xv), fit_res)
        if tp is None:
            continue
        pred.append(float(tp))
        obs.append(float(tv))
    if not pred:
        return 0.0
    pa = np.asarray(pred, dtype=float)
    oa = np.asarray(obs, dtype=float)
    return float(np.sqrt(np.mean((oa - pa) ** 2)))

def _recompute_fit_derived(fit_res: dict) -> dict:
    """Recompute V/ti/depth summary fields after manual segment edits."""
    segs = list((fit_res or {}).get("segments", []) or [])
    segs += [{}, {}, {}]
    segs = segs[:3]

    v0 = float(segs[0].get("velocity_m_s", 0.0) or 0.0)
    v1 = float(segs[1].get("velocity_m_s", 0.0) or 0.0)
    v2 = float(segs[2].get("velocity_m_s", 0.0) or 0.0)
    ti1 = float(segs[1].get("intercept_ms", 0.0) or 0.0) if v1 > 0.0 else 0.0
    ti2 = float(segs[2].get("intercept_ms", 0.0) or 0.0) if v2 > 0.0 else 0.0

    h1 = depth_2layer(ti1, v0, v1) if (v0 > 0.0 and v1 > 0.0 and ti1 > 0.0) else None
    h2 = depth_3layer(ti2, v0, v1, v2, h1) if (h1 and v2 > 0.0 and ti2 > 0.0) else None

    fit_res["segments"] = segs
    fit_res["V0_m_s"] = v0
    fit_res["V1_m_s"] = v1
    fit_res["V2_m_s"] = v2
    fit_res["ti1_ms"] = ti1
    fit_res["ti2_ms"] = ti2
    fit_res["h1_m"] = float(h1) if h1 else 0.0
    fit_res["h2_m"] = float(h2) if h2 else 0.0
    return fit_res

def _drop_fit_layer(fit_res: dict, layer_idx: int) -> dict:
    """Disable one fitted layer (0-based index) and recompute summaries."""
    segs = list((fit_res or {}).get("segments", []) or [])
    while len(segs) < 3:
        segs.append({})
    i = int(layer_idx)
    if 0 <= i < 3:
        seg = dict(segs[i])
        seg["velocity_m_s"] = 0.0
        seg["intercept_ms"] = 0.0
        seg["slope_ms_m"] = 0.0
        seg["r2"] = 0.0
        segs[i] = seg
        wins = list((fit_res or {}).get("windows", []) or [])
        while len(wins) < 3:
            wins.append(None)
        wins[i] = None
        fit_res["windows"] = wins[:3]
        fit_res["segments"] = segs[:3]
    return _recompute_fit_derived(fit_res)

def build_analysis_from_layers(corrected_by_shot: dict,
                               layer_results: dict) -> dict:
    """
    Build per-shot analysis payload for Excel/T-X export.

    Output shape per shot:
      {
        "segments": [...],
        "depths": {...},
        "rms_ms": float,
        "n_fit": int,
        "fit_side": "L"|"R"|"ALL"
      }

    Parameters
    ----------
    corrected_by_shot : dict
        Corrected pick rows produced by build_corrected_pick_data().
    layer_results : dict
        Interactive fit results per shot/side.

    Returns
    -------
    dict
        Compact per-shot analysis payload used by Excel and QC exports.
    """
    out: dict = {}
    for shot_id, rows in (corrected_by_shot or {}).items():
        sides = (layer_results or {}).get(shot_id, {}) or {}
        fit_side, fit_res = _choose_shot_fit(sides)
        if not fit_res:
            continue

        if fit_side == "ALL":
            rows_fit = list(rows)
        else:
            rows_fit = [r for r in rows if r.get("side") == fit_side]
            if not rows_fit:
                rows_fit = list(rows)

        obs: list = []
        pred: list = []
        for r in rows_fit:
            x_abs = float(r.get("true_off_m", 0.0))
            t_obs = float(r.get("fb_interp_inline_ms", 0.0))
            t_pred = _predict_time_from_fit(x_abs, fit_res)
            if t_pred is None:
                continue
            obs.append(t_obs)
            pred.append(float(t_pred))

        if obs:
            oa = np.asarray(obs, dtype=float)
            pa = np.asarray(pred, dtype=float)
            rms = float(np.sqrt(np.mean((oa - pa) ** 2)))
        else:
            rms = 0.0

        segs_raw = list((fit_res or {}).get("segments", []) or [])
        segs: list = []
        for i, seg in enumerate(segs_raw):
            s = dict(seg)
            w = None
            wins = list((fit_res or {}).get("windows", []) or [])
            if i < len(wins):
                w = wins[i]
            if w is not None:
                s["x_start"] = float(w[0])
                s["x_end"] = float(w[1])
            else:
                s["x_start"] = float(seg.get("x0", 0.0))
                s["x_end"] = float(seg.get("x1", 0.0))
            segs.append(s)

        out[int(shot_id)] = {
            "segments": segs,
            "depths": {
                "ti1_ms": float((fit_res or {}).get("ti1_ms", 0.0)),
                "ti2_ms": float((fit_res or {}).get("ti2_ms", 0.0)),
                "h1_m": float((fit_res or {}).get("h1_m", 0.0)),
                "h2_m": float((fit_res or {}).get("h2_m", 0.0)),
            },
            "rms_ms": rms,
            "n_fit": int(len(obs)),
            "fit_side": fit_side or "",
        }

    return out
