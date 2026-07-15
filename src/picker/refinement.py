"""Local refinement of candidate pick times (zero-crossing snap)."""
from __future__ import annotations

from typing import Any

import numpy as np


def _zero_crossing_from_extremum_samples(
    samples: Any,
    dt_s: float,
    win_start_s: float,
    win_end_s: float,
    search_direction: str = "backward",
    use_abs_peak: bool = True,
) -> float | None:
    """Return zero-crossing time (s) near local extremum in a window."""
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
        if x[peak_idx] == 0.0:
            return peak_idx * dt_s
        for right in range(peak_idx, i0, -1):
            left = right - 1
            y0 = float(x[left])
            y1 = float(x[right])
            if y0 == 0.0:
                return left * dt_s
            if y1 == 0.0:
                return right * dt_s
            if y0 * y1 < 0.0:
                frac = -y0 / (y1 - y0)
                return (left + frac) * dt_s
    else:
        if x[peak_idx] == 0.0:
            return peak_idx * dt_s
        for left in range(peak_idx, i1):
            right = left + 1
            y0 = float(x[left])
            y1 = float(x[right])
            if y0 == 0.0:
                return left * dt_s
            if y1 == 0.0:
                return right * dt_s
            if y0 * y1 < 0.0:
                frac = -y0 / (y1 - y0)
                return (left + frac) * dt_s

    return peak_idx * dt_s
