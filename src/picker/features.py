"""Feature computations for first-break picking (Hilbert envelope, onset/peak).

These functions only compute feature arrays / candidate times; they do not
perform the final pick placement (see `src.picker.refinement`).
"""
from __future__ import annotations

from typing import Any

import numpy as np
from obspy import Trace
from scipy.signal import hilbert as _scipy_hilbert

from src.common.settings import HILBERT_ONSET_PCT


def hilbert_envelope_pick(
    trace: Trace,
    win_start_s: float,
    win_end_s: float,
    onset_pct: float | None = HILBERT_ONSET_PCT,
    noise_window_s: tuple[float, float] | None = None,
) -> tuple[Any, float | None, float | None]:
    """
    Compute Hilbert envelope and pick arrival in a time window.

    Returns
    -------
    envelope : ndarray
        Instantaneous amplitude envelope (same sample length as input trace).
    peak_time_s : float | None
        Time of envelope maximum in the search window (seconds from trace start).
    onset_time_s : float | None
        Optional onset time where envelope first exceeds threshold before peak.
    """
    x = np.asarray(trace.data, dtype=float)
    if x.size < 2:
        return np.asarray([], dtype=np.float32), None, None

    dt_s = float(trace.stats.delta)
    env = np.abs(_scipy_hilbert(x))

    i0 = max(0, int(np.floor(float(win_start_s) / dt_s)))
    i1 = min(x.size - 1, int(np.ceil(float(win_end_s) / dt_s)))
    if i1 <= i0:
        return env.astype(np.float32), None, None

    seg = env[i0:i1 + 1]
    if seg.size == 0:
        return env.astype(np.float32), None, None

    peak_rel = int(np.argmax(seg))
    peak_idx = i0 + peak_rel
    peak_time_s = peak_idx * dt_s

    onset_time_s: float | None = None
    if onset_pct is not None:
        pct = float(max(0.0, min(1.0, onset_pct)))
        if noise_window_s is not None:
            n0 = max(0, int(np.floor(float(noise_window_s[0]) / dt_s)))
            n1 = min(x.size - 1, int(np.ceil(float(noise_window_s[1]) / dt_s)))
            if n1 > n0:
                noise_floor = float(np.median(env[n0:n1 + 1]))
            else:
                noise_floor = float(np.median(env[: max(1, i0)])) if i0 > 0 else 0.0
        else:
            noise_floor = float(np.median(env[: max(1, i0)])) if i0 > 0 else 0.0

        peak_val = float(env[peak_idx])
        thr = noise_floor + pct * max(0.0, peak_val - noise_floor)
        for j in range(i0 + 1, peak_idx + 1):
            if env[j] >= thr and env[j - 1] < thr:
                y0 = float(env[j - 1])
                y1 = float(env[j])
                frac = 0.0 if abs(y1 - y0) < 1e-30 else (thr - y0) / (y1 - y0)
                onset_time_s = ((j - 1) + frac) * dt_s
                break
        if onset_time_s is None:
            onset_time_s = float(i0 * dt_s)

    return env.astype(np.float32), float(peak_time_s), onset_time_s
