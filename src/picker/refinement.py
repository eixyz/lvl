"""Local refinement of candidate pick times (zero-crossing snap)."""
from __future__ import annotations

from typing import Any

import numpy as np


def _zero_crossing_before(x: np.ndarray, peak_idx: int, dt_s: float, i0: int) -> float:
    """Zero-crossing time immediately before `peak_idx`, searching backward to `i0`."""
    if x[peak_idx] == 0.0:
        return peak_idx * dt_s
    for right in range(peak_idx, i0, -1):
        left = right - 1
        y0, y1 = float(x[left]), float(x[right])
        if y0 == 0.0:
            return left * dt_s
        if y1 == 0.0:
            return right * dt_s
        if y0 * y1 < 0.0:
            frac = -y0 / (y1 - y0)
            return (left + frac) * dt_s
    return peak_idx * dt_s


def _zero_crossing_after(x: np.ndarray, peak_idx: int, dt_s: float, i1: int) -> float:
    """Zero-crossing time immediately after `peak_idx`, searching forward to `i1`."""
    if x[peak_idx] == 0.0:
        return peak_idx * dt_s
    for left in range(peak_idx, i1):
        right = left + 1
        y0, y1 = float(x[left]), float(x[right])
        if y0 == 0.0:
            return left * dt_s
        if y1 == 0.0:
            return right * dt_s
        if y0 * y1 < 0.0:
            frac = -y0 / (y1 - y0)
            return (left + frac) * dt_s
    return peak_idx * dt_s


def _zero_crossing_from_extremum_samples(
    samples: Any,
    dt_s: float,
    win_start_s: float,
    win_end_s: float,
    search_direction: str = "backward",
    use_abs_peak: bool = True,
) -> float | None:
    """Return zero-crossing time (s) near local extremum in a window.

    Finds the single largest-amplitude sample anywhere in the window and
    snaps to the nearest zero-crossing. Appropriate for REFINEMENT of an
    already-approximate pick in a narrow, already-localized window (e.g.
    snapping a manual click) - NOT appropriate as a standalone first-break
    detector over a wide velocity-gated search window, where a later,
    louder reverberation cycle will usually beat the true (often smaller,
    emergent) first onset. Use
    `first_significant_extremum_zero_crossing` for that instead.
    """
    x = np.asarray(samples, dtype=float)
    if x.size < 2:
        return None

    i0 = max(0, int(np.floor(win_start_s / dt_s)))
    i1 = min(x.size - 1, int(np.ceil(win_end_s / dt_s)))
    if i1 <= i0:
        return None

    seg = x[i0:i1 + 1]
    if seg.size == 0:
        return None

    if use_abs_peak:
        rel_idx = int(np.argmax(np.abs(seg)))
    else:
        rel_idx = int(np.argmax(seg))
    peak_idx = i0 + rel_idx

    direction = str(search_direction).strip().lower()
    if direction not in {"backward", "forward"}:
        direction = "backward"

    if direction == "backward":
        return _zero_crossing_before(x, peak_idx, dt_s, i0)
    return _zero_crossing_after(x, peak_idx, dt_s, i1)


def first_significant_extremum_zero_crossing(
    samples: Any,
    dt_s: float,
    win_start_s: float,
    win_end_s: float,
    threshold_ratio: float = 0.25,
    noise_win_start_s: float | None = None,
    noise_win_end_s: float | None = None,
    search_direction: str = "backward",
) -> float | None:
    """First-break-appropriate alternative to
    `_zero_crossing_from_extremum_samples`: finds the EARLIEST local
    extremum in the window whose absolute amplitude exceeds
    `threshold_ratio` of the window's own peak amplitude (or, if a noise
    window is given, `threshold_ratio` x noise std above the noise
    floor) - then snaps to the zero-crossing at that onset, not at
    whichever sample happens to be loudest overall.

    This is what makes it suitable as a standalone picker over a wide
    velocity-gated window: a genuine first break is, by definition, the
    FIRST energy to arrive - it doesn't have to be the loudest part of
    the trace, and in practice often isn't (later reverberation cycles
    routinely have larger amplitude than the initial onset).
    """
    x = np.asarray(samples, dtype=float)
    if x.size < 2:
        return None

    i0 = max(0, int(np.floor(win_start_s / dt_s)))
    i1 = min(x.size - 1, int(np.ceil(win_end_s / dt_s)))
    if i1 <= i0:
        return None

    seg = x[i0:i1 + 1]
    if seg.size < 2:
        return None

    if noise_win_start_s is not None and noise_win_end_s is not None:
        n0 = max(0, int(np.floor(noise_win_start_s / dt_s)))
        n1 = min(x.size - 1, int(np.ceil(noise_win_end_s / dt_s)))
        noise_seg = x[n0:n1 + 1] if n1 > n0 else seg
        floor = float(np.std(noise_seg)) if noise_seg.size >= 2 else float(np.std(seg))
        threshold = max(threshold_ratio * float(np.max(np.abs(seg))), 4.0 * floor)
    else:
        threshold = threshold_ratio * float(np.max(np.abs(seg)))

    if threshold <= 1e-20:
        return None

    above = np.where(np.abs(seg) >= threshold)[0]
    if above.size == 0:
        return None
    onset_rel = int(above[0])
    peak_idx = i0 + onset_rel

    direction = str(search_direction).strip().lower()
    if direction not in {"backward", "forward"}:
        direction = "backward"

    if direction == "backward":
        return _zero_crossing_before(x, peak_idx, dt_s, i0)
    return _zero_crossing_after(x, peak_idx, dt_s, i1)