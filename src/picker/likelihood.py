"""
Picker Engine V2 - likelihood fusion.

Turns a per-trace `FeatureSet` (see `src.picker.features`) into a single
arrival-likelihood curve by weighted fusion, then finds candidate arrival
samples on that curve. This module does not know about neighbouring
traces (see `src.picker.coherence` for that) or about optimizing a path
across a whole shot gather (see `src.picker.optimizer`) - it only answers
"for this one trace, which samples look like the arrival, and how much".
"""
from __future__ import annotations

from typing import Any

import numpy as np
from scipy.signal import find_peaks as _scipy_find_peaks

from src.picker.models import FeatureSet
from src.picker.settings import PickerSettings


# -----------------------------------------------------------------------------
# Fusion primitives
# -----------------------------------------------------------------------------

def weighted_sum(feature_set: FeatureSet, weights: dict[str, float]) -> np.ndarray:
    """Combine the available features into one weighted-sum curve.

    Features that are `None` (not computed for this trace) are skipped
    entirely, and the remaining weights are renormalized to sum to 1 -
    so a trace missing e.g. `coherence` isn't unfairly penalized relative
    to one that has it.
    """
    available = {name: arr for name, arr in feature_set.as_dict().items() if arr is not None}
    if not available:
        return np.asarray([], dtype=np.float32)

    n = len(next(iter(available.values())))
    total_weight = sum(max(0.0, float(weights.get(name, 0.0))) for name in available)
    if total_weight <= 1e-20:
        # No usable weights - fall back to an unweighted average so the
        # curve is still meaningful rather than all-zero.
        total_weight = float(len(available))
        weights = {name: 1.0 for name in available}

    out = np.zeros(n, dtype=np.float64)
    for name, arr in available.items():
        w = max(0.0, float(weights.get(name, 0.0))) / total_weight
        if w > 0.0:
            out += w * np.asarray(arr, dtype=np.float64)
    return out.astype(np.float32)


def adaptive_weights(feature_set: FeatureSet, base_weights: dict[str, float]) -> dict[str, float]:
    """Boost/penalize each feature's base weight by how "peaky" it is on
    this specific trace.

    A feature with one sharp, well-defined maximum (high peak-to-median
    ratio) is more trustworthy for this trace than one that's flat or
    noisy, so it gets weighted up; a nearly-flat feature gets weighted
    down rather than dragging the fused curve toward noise.
    """
    out: dict[str, float] = {}
    for name, arr in feature_set.as_dict().items():
        if arr is None or name not in base_weights:
            continue
        x = np.asarray(arr, dtype=np.float64)
        if x.size == 0:
            continue
        peak = float(np.max(x))
        med = float(np.median(x))
        spread = max(peak - med, 1e-6)
        # Reliability in (0, ~2]: sharp/prominent peak -> higher; flat -> lower.
        reliability = float(np.clip(spread / (med + 1e-6), 0.2, 2.0))
        out[name] = max(0.0, float(base_weights[name])) * reliability
    return out


def normalize_probability(values: Any) -> np.ndarray:
    """Rescale a fused curve to [0, 1].

    This is a *relative* likelihood curve (its peak marks the most likely
    arrival sample), not a probability mass function - it does not sum to 1.
    """
    x = np.asarray(values, dtype=np.float64)
    if x.size == 0:
        return np.asarray([], dtype=np.float32)
    x = np.clip(x, 0.0, None)
    lo, hi = float(np.min(x)), float(np.max(x))
    if hi - lo < 1e-20:
        return np.zeros_like(x, dtype=np.float32)
    return ((x - lo) / (hi - lo)).astype(np.float32)


def smooth_probability(values: Any, window: int = 5) -> np.ndarray:
    """Light smoothing (Hann-weighted moving average) to suppress single-
    sample spikes before peak-picking. `window` is in samples; even
    values are bumped up to the next odd number.
    """
    x = np.asarray(values, dtype=np.float64)
    if x.size == 0:
        return np.asarray([], dtype=np.float32)
    win = max(1, int(window))
    if win % 2 == 0:
        win += 1
    if win <= 1 or x.size < win:
        return x.astype(np.float32)

    kernel = np.hanning(win + 2)[1:-1]  # drop the zero end-taps
    kernel = kernel / kernel.sum()
    pad = win // 2
    xp = np.pad(x, (pad, pad), mode="edge")
    smoothed = np.convolve(xp, kernel, mode="valid")
    return smoothed.astype(np.float32)


# -----------------------------------------------------------------------------
# High-level entry points
# -----------------------------------------------------------------------------

def compute_arrival_likelihood(
    feature_set: FeatureSet,
    settings: PickerSettings | None = None,
    adaptive: bool = True,
    smooth_window: int = 5,
) -> np.ndarray:
    """Fuse a `FeatureSet` into one normalized, smoothed likelihood curve."""
    settings = settings or PickerSettings()
    base_weights = settings.feature_weights()
    weights = adaptive_weights(feature_set, base_weights) if adaptive else base_weights

    raw = weighted_sum(feature_set, weights)
    if raw.size == 0:
        return raw
    prob = normalize_probability(raw)
    prob = smooth_probability(prob, smooth_window)
    return normalize_probability(prob)  # re-normalize after smoothing flattens the peak slightly


def compute_candidate_peaks(
    probability: Any,
    min_prominence: float = 0.05,
    min_distance_samples: int = 5,
    max_candidates: int = 5,
) -> list[tuple[int, float]]:
    """Find local maxima on a likelihood curve, sorted best-first.

    Returns a list of (sample_index, score) tuples, highest score first,
    capped at `max_candidates`.
    """
    x = np.asarray(probability, dtype=np.float64)
    if x.size == 0:
        return []

    peaks, props = _scipy_find_peaks(
        x, prominence=max(0.0, float(min_prominence)),
        distance=max(1, int(min_distance_samples)),
    )
    if peaks.size == 0:
        # Fall back to the global maximum so callers always get at least
        # one candidate for a non-empty curve.
        idx = int(np.argmax(x))
        return [(idx, float(x[idx]))]

    scored = sorted(((int(p), float(x[p])) for p in peaks), key=lambda t: t[1], reverse=True)
    return scored[:max(1, int(max_candidates))]


def compute_arrival_probability(
    samples: Any,
    dt_s: float,
    settings: PickerSettings | None = None,
    sta_ms: float | None = None,
    lta_ms: float | None = None,
    adaptive: bool = True,
    smooth_window: int = 5,
) -> tuple[np.ndarray, FeatureSet]:
    """Full per-trace entry point: extract features, then fuse them.

    Returns (probability_curve, feature_set) - the feature set is handed
    back too because `src.picker.coherence` and `confidence.py` reuse the
    individual per-feature curves, not just the fused result.
    """
    # Local import avoids a hard import-time dependency between features.py
    # and likelihood.py in either direction.
    from src.picker.features import extract_features

    settings = settings or PickerSettings()
    feature_set = extract_features(
        samples, dt_s,
        sta_ms=sta_ms if sta_ms is not None else settings.sta_window,
        lta_ms=lta_ms if lta_ms is not None else settings.lta_window,
    )
    probability = compute_arrival_likelihood(
        feature_set, settings, adaptive=adaptive, smooth_window=smooth_window
    )
    return probability, feature_set


__all__ = [
    "weighted_sum", "adaptive_weights", "normalize_probability", "smooth_probability",
    "compute_arrival_likelihood", "compute_candidate_peaks", "compute_arrival_probability",
]