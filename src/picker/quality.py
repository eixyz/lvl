"""
Picker Engine V2 - physical consistency validation.

`confidence.py` scores how trustworthy a single pick looks on its own
trace. This module checks the *gather* as a whole against basic physics:
first-arrival travel time should increase smoothly with offset, implied
apparent velocity should stay within a plausible range for near-surface
refraction, and neighbouring picks shouldn't jump around unrealistically.
None of these checks require re-touching the waveform - they only look
at the chosen pick samples plus geometry.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from src.common.settings import AUTO_VMAX_M_S, AUTO_VMIN_M_S


def validate_pick_window(
    sample: int | None,
    win_start_s: float,
    win_end_s: float,
    dt_s: float,
) -> bool:
    """True if `sample` falls inside [win_start_s, win_end_s] (inclusive)."""
    if sample is None:
        return False
    t = float(sample) * dt_s
    return win_start_s <= t <= win_end_s


def detect_outlier_pick(
    picks: list[int | None],
    window: int = 5,
    max_jump_samples: float = 50.0,
) -> list[bool]:
    """Non-destructive query version of `optimizer.remove_outliers`.

    Returns one bool per trace: True if that pick deviates from its local
    neighbourhood median by more than `max_jump_samples`.
    """
    n = len(picks)
    flags = [False] * n
    vals = np.array([np.nan if s is None else float(s) for s in picks])
    half = max(1, int(window)) // 2

    for i in range(n):
        if np.isnan(vals[i]):
            continue
        lo, hi = max(0, i - half), min(n, i + half + 1)
        neighbours = np.delete(vals[lo:hi], i - lo)
        neighbours = neighbours[~np.isnan(neighbours)]
        if neighbours.size == 0:
            continue
        med = float(np.median(neighbours))
        flags[i] = abs(vals[i] - med) > max_jump_samples
    return flags


def detect_crossing_picks(
    picks: list[int | None],
    offsets_m: list[float],
    tolerance_samples: float = 0.0,
) -> list[bool]:
    """Flag picks that break the expected offset -> travel-time ordering.

    For a single refracting medium, travel time should be non-decreasing
    with |offset|. This sorts traces by |offset| and flags any pick whose
    travel time is smaller than the running maximum seen so far at lower
    offsets (beyond `tolerance_samples`), which usually indicates a
    mis-pick rather than real geology.
    """
    n = len(picks)
    flags = [False] * n
    if n == 0 or len(offsets_m) != n:
        return flags

    order = sorted(range(n), key=lambda i: abs(offsets_m[i]))
    running_max = -np.inf
    for i in order:
        s = picks[i]
        if s is None:
            continue
        if s < running_max - tolerance_samples:
            flags[i] = True
        else:
            running_max = max(running_max, s)
    return flags


def enforce_velocity_limits(
    picks: list[int | None],
    offsets_m: list[float],
    dt_s: float,
    vmin_m_s: float = AUTO_VMIN_M_S,
    vmax_m_s: float = AUTO_VMAX_M_S,
    min_offset_m: float = 1.0,
) -> list[bool]:
    """Flag picks whose implied apparent velocity (|offset| / time) falls
    outside [vmin_m_s, vmax_m_s].

    Reuses the same bounds as the legacy `AUTO_USE_VELOCITY_GATE` picking
    gate in `src.common.settings`, so v2 stays consistent with v1 unless
    explicitly overridden. Traces within `min_offset_m` of the source are
    skipped (never flagged) - apparent velocity is undefined/meaningless
    right at zero offset, so this avoids a spurious flag on every
    gather's near-offset trace.
    """
    n = len(picks)
    flags = [False] * n
    if len(offsets_m) != n:
        return flags

    for i in range(n):
        s = picks[i]
        if s is None or abs(offsets_m[i]) < min_offset_m:
            continue
        t = float(s) * dt_s
        if t <= 1e-9:
            continue
        v = abs(offsets_m[i]) / t
        flags[i] = not (vmin_m_s <= v <= vmax_m_s)
    return flags


def detect_unrealistic_jump(
    picks: list[int | None],
    max_jump_samples: float = 60.0,
) -> list[bool]:
    """Flag picks that differ from their *immediate* neighbour by more
    than `max_jump_samples` (cheaper, more local than `detect_outlier_pick`).
    """
    n = len(picks)
    flags = [False] * n
    for i in range(n):
        if picks[i] is None:
            continue
        neighbours = [picks[j] for j in (i - 1, i + 1) if 0 <= j < n and picks[j] is not None]
        if not neighbours:
            continue
        if all(abs(picks[i] - nb) > max_jump_samples for nb in neighbours):
            flags[i] = True
    return flags


def check_pick_consistency(
    picks: list[int | None],
    offsets_m: list[float],
    dt_s: float,
    vmin_m_s: float = AUTO_VMIN_M_S,
    vmax_m_s: float = AUTO_VMAX_M_S,
    max_jump_samples: float = 50.0,
) -> dict[str, list[bool]]:
    """Run every gather-level check and return the flags side by side."""
    return {
        "outlier": detect_outlier_pick(picks, max_jump_samples=max_jump_samples),
        "crossing": detect_crossing_picks(picks, offsets_m),
        "velocity": enforce_velocity_limits(picks, offsets_m, dt_s, vmin_m_s, vmax_m_s),
        "unrealistic_jump": detect_unrealistic_jump(picks, max_jump_samples=max_jump_samples * 1.2),
    }


def quality_report(
    picks: list[int | None],
    offsets_m: list[float],
    dt_s: float,
    vmin_m_s: float = AUTO_VMIN_M_S,
    vmax_m_s: float = AUTO_VMAX_M_S,
    max_jump_samples: float = 50.0,
) -> list[dict[str, Any]]:
    """Per-trace consolidated report: `{sample, valid, flags: [...]}`.

    `valid` is False only when the pick failed a hard physical check
    (crossing order or outside the velocity gate); `outlier` /
    `unrealistic_jump` are softer signals surfaced as flags but don't by
    themselves invalidate a pick, since edge traces or real near-surface
    velocity contrasts can legitimately look unusual.
    """
    checks = check_pick_consistency(picks, offsets_m, dt_s, vmin_m_s, vmax_m_s, max_jump_samples)
    n = len(picks)
    report: list[dict[str, Any]] = []
    for i in range(n):
        flags = [name.upper() for name, arr in checks.items() if arr[i]]
        if picks[i] is None:
            flags.append("NO_PICK")
        valid = picks[i] is not None and not checks["crossing"][i] and not checks["velocity"][i]
        report.append({"index": i, "sample": picks[i], "valid": valid, "flags": flags})
    return report


__all__ = [
    "validate_pick_window", "detect_outlier_pick", "detect_crossing_picks",
    "enforce_velocity_limits", "detect_unrealistic_jump",
    "check_pick_consistency", "quality_report",
]