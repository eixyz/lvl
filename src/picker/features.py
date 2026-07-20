"""Feature computations for first-break picking (Hilbert envelope, onset/peak).

These functions only compute feature arrays / candidate times; they do not
perform the final pick placement (see `src.picker.refinement`).

Feature extraction for the LVL probabilistic first-break picker.

All returned features are normalized to [0,1].

Pipeline

Trace
    ↓
Hilbert
STA/LTA
AIC
Gradient
Energy
    ↓
Normalization
    ↓
FeatureSet
"""

from __future__ import annotations

from typing import Any

import numpy as np
from obspy import Trace
from obspy.signal.trigger import classic_sta_lta

from scipy.ndimage import gaussian_filter1d
from scipy.signal import hilbert

from .models import FeatureSet

from src.common.settings import HILBERT_ONSET_PCT


EPS = 1e-12


# -------------------------------------------------------------
# Normalization
# -------------------------------------------------------------

def normalize_feature(
    x: np.ndarray,
    invert: bool = False,
    lower_percentile: float = 2.0,
    upper_percentile: float = 98.0,
) -> np.ndarray:
    """
    Robust normalization to [0,1].

    Uses percentiles instead of min/max to reduce sensitivity
    to spikes.
    """

    x = np.asarray(x, dtype=np.float64)

    if x.size == 0:
        return np.zeros(0, dtype=np.float32)

    lo = np.percentile(x, lower_percentile)
    hi = np.percentile(x, upper_percentile)

    if hi <= lo:
        out = np.zeros_like(x)

    else:
        out = (x - lo) / (hi - lo)

    out = np.clip(out, 0.0, 1.0)

    if invert:
        out = 1.0 - out

    return out.astype(np.float32)


# -------------------------------------------------------------
# Hilbert envelope
# -------------------------------------------------------------

def compute_hilbert_feature(
    samples: np.ndarray,
) -> np.ndarray:

    env = np.abs(hilbert(samples))

    env = gaussian_filter1d(env, sigma=2)

    return normalize_feature(env)


# -------------------------------------------------------------
# STA/LTA
# -------------------------------------------------------------

def compute_stalta_feature(
    samples: np.ndarray,
    fs: float,
    sta_ms: float = 2.0,
    lta_ms: float = 20.0,
) -> np.ndarray:

    sta = max(2, int(sta_ms * fs / 1000))
    lta = max(sta + 1, int(lta_ms * fs / 1000))

    ratio = classic_sta_lta(samples.astype(np.float64), sta, lta)

    ratio[np.isnan(ratio)] = 0

    return normalize_feature(ratio)


# -------------------------------------------------------------
# Fast AIC
# -------------------------------------------------------------

def compute_aic_feature(
    samples: np.ndarray,
) -> np.ndarray:
    """
    Fast O(N) Akaike Information Criterion.
    """

    x = np.asarray(samples, dtype=np.float64)

    n = len(x)

    if n < 10:
        return np.zeros(n, dtype=np.float32)

    cs = np.cumsum(x)
    cs2 = np.cumsum(x * x)

    aic = np.full(n, np.nan)

    for k in range(2, n - 2):

        left_n = k
        right_n = n - k

        left_mean = cs[k - 1] / left_n

        right_mean = (cs[-1] - cs[k - 1]) / right_n

        left_var = (
            cs2[k - 1]
            - 2 * left_mean * cs[k - 1]
            + left_n * left_mean**2
        ) / left_n

        right_var = (
            cs2[-1]
            - cs2[k - 1]
            - 2 * right_mean * (cs[-1] - cs[k - 1])
            + right_n * right_mean**2
        ) / right_n

        left_var = max(left_var, EPS)
        right_var = max(right_var, EPS)

        aic[k] = (
            left_n * np.log(left_var)
            + right_n * np.log(right_var)
        )

    aic[np.isnan(aic)] = np.nanmax(aic)

    aic = gaussian_filter1d(aic, sigma=2)

    return normalize_feature(aic, invert=True)


# -------------------------------------------------------------
# Gradient
# -------------------------------------------------------------

def compute_gradient_feature(
    samples: np.ndarray,
) -> np.ndarray:

    grad = np.abs(np.gradient(samples))

    grad = gaussian_filter1d(grad, sigma=1)

    return normalize_feature(grad)


# -------------------------------------------------------------
# Energy
# -------------------------------------------------------------

def compute_energy_feature(
    samples: np.ndarray,
    window: int = 12,
) -> np.ndarray:

    x2 = samples * samples

    kernel = np.ones(window) / window

    energy = np.convolve(
        x2,
        kernel,
        mode="same",
    )

    return normalize_feature(energy)


# -------------------------------------------------------------
# Main extraction
# -------------------------------------------------------------

def extract_features(
    trace: Trace,
) -> FeatureSet:
    """
    Compute every feature required by the picker.

    Returns
    -------
    FeatureSet
    """

    samples = np.asarray(trace.data, dtype=np.float64)

    fs = trace.stats.sampling_rate

    hilb = compute_hilbert_feature(samples)

    stalta = compute_stalta_feature(samples, fs)

    aic = compute_aic_feature(samples)

    gradient = compute_gradient_feature(samples)

    energy = compute_energy_feature(samples)

    return FeatureSet(
        hilbert=hilb,
        stalta=stalta,
        aic=aic,
        gradient=gradient,
        energy=energy,
        coherence=None,
    )

#################################

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
    env = np.abs(hilbert(x))

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