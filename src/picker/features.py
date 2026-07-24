"""Feature computations for first-break picking (Hilbert envelope, onset/peak).

These functions only compute feature arrays / candidate times; they do not
perform the final pick placement (see `src.picker.refinement`).

The functions below `compute_aic` / `compute_gradient` are the start of the
"Picker Engine V2" feature layer: every `compute_*_feature` function returns
a per-sample curve normalized to [0, 1] where **higher means more likely to
be the arrival**, so they can be fused directly by `likelihood.py` without
each caller worrying about sign conventions or scale.
"""
from __future__ import annotations

from typing import Any

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from obspy import Trace
from scipy.signal import hilbert as _scipy_hilbert

from src.common.settings import HILBERT_ONSET_PCT
from src.picker.models import FeatureSet


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

def compute_aic(samples: Any) -> np.ndarray:
    """
    Compute AIC (Akaike Information Criterion) function for a 1D signal.

    Classic Maeda/Akaike-picker formula:
        AIC[k] = k * log(var(x[:k])) + (n-k) * log(var(x[k:]))

    Computed via cumulative sums so each expanding-window variance is
    O(1) instead of recomputing `np.var` from scratch at every k (which
    made this O(n^2) - the dominant cost in the whole picking pipeline
    on real, thousands-of-samples-long traces).

    Parameters
    ----------
    samples : array-like
        Input 1D signal samples.

    Returns
    -------
    aic : ndarray
        AIC values for each sample in the input signal (0 at the edges,
        matching the original loop's range(1, n-1)).
    """
    x = np.asarray(samples, dtype=np.float64)
    n = x.size
    if n < 2:
        return np.asarray([], dtype=np.float32)
    if n < 3:
        return np.zeros(n, dtype=np.float32)

    cs1 = np.concatenate(([0.0], np.cumsum(x)))
    cs2 = np.concatenate(([0.0], np.cumsum(x * x)))

    k = np.arange(1, n - 1)
    sum1, sum1_sq = cs1[k], cs2[k]
    mean1 = sum1 / k
    var1 = np.maximum(sum1_sq / k - mean1 ** 2, 0.0)

    n_k = n - k
    sum2 = cs1[n] - cs1[k]
    sum2_sq = cs2[n] - cs2[k]
    mean2 = sum2 / n_k
    var2 = np.maximum(sum2_sq / n_k - mean2 ** 2, 0.0)

    aic = np.zeros(n, dtype=np.float64)
    aic[k] = k * np.log(var1 + 1e-10) + n_k * np.log(var2 + 1e-10)
    return aic.astype(np.float32)

def compute_gradient(samples: Any) -> np.ndarray:
    """
    Compute the gradient (first derivative) of a 1D signal.

    Parameters
    ----------
    samples : array-like
        Input 1D signal samples.

    Returns
    -------
    gradient : ndarray
        Gradient values for each sample in the input signal.
    """
    x = np.asarray(samples, dtype=float)
    if x.size < 2:
        return np.asarray([], dtype=np.float32)

    grad = np.gradient(x)
    return grad.astype(np.float32)


# -----------------------------------------------------------------------------
# Normalization
# -----------------------------------------------------------------------------

def normalize_feature(values: Any, invert: bool = False) -> np.ndarray:
    """Min-max normalize a 1D array to [0, 1].

    Parameters
    ----------
    values : array-like
        Raw feature curve.
    invert : bool, optional
        If True, high raw values become LOW output (used for features like
        AIC, where the arrival sits at a *minimum*, not a maximum).

    A constant (or empty/too-short) input returns an all-zero array of the
    same length rather than raising, so a degenerate trace just contributes
    nothing to the fused likelihood instead of crashing the pipeline.
    """
    x = np.asarray(values, dtype=np.float64)
    if x.size == 0:
        return np.asarray([], dtype=np.float32)

    lo, hi = float(np.min(x)), float(np.max(x))
    if hi - lo < 1e-20:
        return np.zeros_like(x, dtype=np.float32)

    out = (x - lo) / (hi - lo)
    if invert:
        out = 1.0 - out
    return out.astype(np.float32)


def _rolling_windows(x: np.ndarray, win: int) -> np.ndarray:
    """Return a (n, win) view of overlapping windows, edge-padded by reflection."""
    win = max(3, int(win))
    if win % 2 == 0:
        win += 1
    pad = win // 2
    xp = np.pad(x, (pad, pad), mode="reflect")
    return sliding_window_view(xp, win)


# -----------------------------------------------------------------------------
# Per-feature curves (each returns a normalized [0, 1] array, same length as
# the input trace)
# -----------------------------------------------------------------------------

def compute_hilbert_feature(trace: Any, dt_s: float | None = None) -> np.ndarray:
    """Normalized Hilbert-envelope amplitude for the whole trace."""
    x = np.asarray(trace, dtype=float) if not hasattr(trace, "data") else np.asarray(trace.data, dtype=float)
    if x.size < 2:
        return np.asarray([], dtype=np.float32)
    env = np.abs(_scipy_hilbert(x))
    return normalize_feature(env)


def compute_stalta_feature(samples: Any, dt_s: float, sta_ms: float, lta_ms: float) -> np.ndarray:
    """Normalized short-term/long-term-average ratio curve.

    Same characteristic function as the legacy `_auto_pick_stalta` trigger
    in lvl_refraction.py, but returned as a full per-sample curve rather
    than a single trigger sample.
    """
    x = np.asarray(samples, dtype=float)
    n = x.size
    if n < 4:
        return np.asarray([], dtype=np.float32)

    n_sta = max(1, int(sta_ms / 1000.0 / dt_s))
    n_lta = max(3, int(lta_ms / 1000.0 / dt_s))

    energy = x.astype(np.float64) ** 2
    sta_kernel = np.ones(n_sta) / n_sta
    lta_kernel = np.ones(n_lta) / n_lta
    sta = np.convolve(energy, sta_kernel, mode="same")
    lta = np.convolve(energy, lta_kernel, mode="same")
    lta = np.where(lta > 1e-20, lta, 1e-20)

    ratio = sta / lta
    return normalize_feature(ratio)


def compute_aic_feature(samples: Any) -> np.ndarray:
    """Normalized AIC feature (the AIC *minimum* marks the likely arrival)."""
    raw = compute_aic(samples)
    if raw.size == 0:
        return raw
    return normalize_feature(raw, invert=True)


def compute_gradient_feature(samples: Any) -> np.ndarray:
    """Normalized absolute gradient (sharp waveform changes score highest)."""
    raw = compute_gradient(samples)
    if raw.size == 0:
        return raw
    return normalize_feature(np.abs(raw))


def compute_energy_feature(samples: Any, window: int = 9) -> np.ndarray:
    """Normalized short-window energy (mean squared amplitude)."""
    x = np.asarray(samples, dtype=float)
    if x.size < 2:
        return np.asarray([], dtype=np.float32)
    windows = _rolling_windows(x, window)
    energy = np.mean(windows ** 2, axis=1)
    return normalize_feature(energy)


def compute_snr_feature(samples: Any, window: int = 9, noise_samples: int | None = None) -> np.ndarray:
    """Normalized local SNR: short-window RMS divided by a fixed noise floor.

    `noise_samples` should be the number of leading samples that are
    "before the shot" (pure noise). Defaults to the first 10% of the
    trace if not given.
    """
    x = np.asarray(samples, dtype=float)
    n = x.size
    if n < 2:
        return np.asarray([], dtype=np.float32)

    if noise_samples is None or noise_samples < 2:
        noise_samples = max(2, n // 10)
    noise_floor = float(np.std(x[:noise_samples])) if noise_samples < n else float(np.std(x))
    noise_floor = max(noise_floor, 1e-20)

    windows = _rolling_windows(x, window)
    local_rms = np.sqrt(np.mean(windows ** 2, axis=1))
    snr = local_rms / noise_floor
    return normalize_feature(snr)


def compute_kurtosis_feature(samples: Any, window: int = 21) -> np.ndarray:
    """Normalized rolling-window excess kurtosis.

    Kurtosis rises sharply when a window straddles the impulsive onset of
    a seismic arrival, which is a classic onset-detection statistic
    (Baillard et al., 2014-style kurtosis picker).
    """
    x = np.asarray(samples, dtype=float)
    if x.size < 5:
        return np.asarray([], dtype=np.float32)
    windows = _rolling_windows(x, window)
    mean = np.mean(windows, axis=1)
    std = np.std(windows, axis=1)
    std = np.where(std > 1e-20, std, 1e-20)
    centered = windows - mean[:, None]
    m4 = np.mean(centered ** 4, axis=1)
    kurt = m4 / (std ** 4) - 3.0  # excess kurtosis (0 for a Gaussian window)
    return normalize_feature(kurt)


def compute_skewness_feature(samples: Any, window: int = 21) -> np.ndarray:
    """Normalized rolling-window skewness."""
    x = np.asarray(samples, dtype=float)
    if x.size < 5:
        return np.asarray([], dtype=np.float32)
    windows = _rolling_windows(x, window)
    mean = np.mean(windows, axis=1)
    std = np.std(windows, axis=1)
    std = np.where(std > 1e-20, std, 1e-20)
    centered = windows - mean[:, None]
    m3 = np.mean(centered ** 3, axis=1)
    skew = m3 / (std ** 3)
    return normalize_feature(np.abs(skew))


def compute_teager_energy(samples: Any) -> np.ndarray:
    """Teager-Kaiser energy operator: psi[n] = x[n]^2 - x[n-1]*x[n+1].

    Sensitive to both amplitude and instantaneous frequency changes,
    making it a good complement to plain energy for onset detection.
    """
    x = np.asarray(samples, dtype=float)
    n = x.size
    if n < 3:
        return np.asarray([], dtype=np.float32)
    psi = np.zeros(n, dtype=np.float64)
    psi[1:-1] = x[1:-1] ** 2 - x[:-2] * x[2:]
    psi[0] = psi[1]
    psi[-1] = psi[-2]
    return normalize_feature(psi)


def compute_instantaneous_frequency(samples: Any, dt_s: float) -> np.ndarray:
    """Normalized instantaneous frequency from the analytic signal's phase.

    A jump in instantaneous frequency often accompanies the first arrival
    as the trace transitions from low-frequency noise to the higher-
    frequency P-wave onset.
    """
    x = np.asarray(samples, dtype=float)
    n = x.size
    if n < 3:
        return np.asarray([], dtype=np.float32)
    analytic = _scipy_hilbert(x)
    phase = np.unwrap(np.angle(analytic))
    inst_freq = np.diff(phase) / (2.0 * np.pi * dt_s)
    inst_freq = np.concatenate([[inst_freq[0]], inst_freq])  # pad back to length n
    return normalize_feature(np.abs(inst_freq))


# -----------------------------------------------------------------------------
# Aggregator
# -----------------------------------------------------------------------------

def extract_features(
    samples: Any,
    dt_s: float,
    sta_ms: float = 3.0,
    lta_ms: float = 20.0,
    energy_window: int = 9,
    stat_window: int = 21,
    noise_samples: int | None = None,
) -> FeatureSet:
    """Compute every normalized feature curve for one trace.

    Any individual feature that fails to compute (e.g. trace too short)
    is left as `None` on the returned `FeatureSet` rather than raising,
    so one bad trace doesn't abort picking a whole profile.
    """
    x = np.asarray(samples, dtype=float)
    fs = FeatureSet()

    def _safe(fn, *args, **kwargs):
        try:
            result = fn(*args, **kwargs)
            return result if result.size == x.size and result.size > 0 else None
        except Exception:
            return None

    fs.hilbert = _safe(compute_hilbert_feature, x)
    fs.stalta = _safe(compute_stalta_feature, x, dt_s, sta_ms, lta_ms)
    fs.aic = _safe(compute_aic_feature, x)
    fs.gradient = _safe(compute_gradient_feature, x)
    fs.energy = _safe(compute_energy_feature, x, energy_window)
    fs.snr = _safe(compute_snr_feature, x, energy_window, noise_samples)
    fs.kurtosis = _safe(compute_kurtosis_feature, x, stat_window)
    fs.skewness = _safe(compute_skewness_feature, x, stat_window)
    # `coherence` is intentionally left None here - it needs neighbouring
    # traces and is filled in by src.picker.coherence, not per-trace.
    return fs


__all__ = [
    "hilbert_envelope_pick",
    "compute_aic", "compute_gradient",
    "normalize_feature",
    "compute_hilbert_feature", "compute_stalta_feature",
    "compute_aic_feature", "compute_gradient_feature",
    "compute_energy_feature", "compute_snr_feature",
    "compute_kurtosis_feature", "compute_skewness_feature",
    "compute_teager_energy", "compute_instantaneous_frequency",
    "extract_features",
]