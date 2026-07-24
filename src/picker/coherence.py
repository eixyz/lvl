"""
Picker Engine V2 - neighbour coherence.

Everything in `src.picker.features` / `src.picker.likelihood` works on one
trace at a time. This module is the first stage that looks *across*
traces in a shot gather: real first breaks line up smoothly from
receiver to receiver (a moveout curve), while noise picks don't. A
sample whose local waveform shape resembles its neighbours' at the same
two-way time is more likely to be a genuine, aligned arrival.

Design principle (see picker v2 spec): coherence *improves* the fused
likelihood curve, it never replaces it - a sample with strong single-
trace evidence but poor coherence (e.g. a receiver near a geometry
discontinuity) should still be pickable, just not preferentially boosted.
"""
from __future__ import annotations

from typing import Any

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from src.picker.likelihood import normalize_probability


# -----------------------------------------------------------------------------
# Building blocks
# -----------------------------------------------------------------------------

def cross_correlation(trace_a: Any, trace_b: Any, max_lag: int = 20) -> tuple[int, float]:
    """Normalized cross-correlation between two traces.

    Returns (best_lag, best_score) where `best_lag` is the sample shift
    of `trace_b` relative to `trace_a` (positive = trace_b lags behind)
    that maximizes the normalized correlation coefficient, and
    `best_score` is that coefficient in [-1, 1].
    """
    a = np.asarray(trace_a, dtype=np.float64)
    b = np.asarray(trace_b, dtype=np.float64)
    n = min(a.size, b.size)
    if n < 2:
        return 0, 0.0
    a, b = a[:n], b[:n]

    a = a - a.mean()
    b = b - b.mean()
    denom = (np.linalg.norm(a) * np.linalg.norm(b))
    if denom < 1e-20:
        return 0, 0.0

    max_lag = max(0, min(int(max_lag), n - 1))
    best_lag, best_score = 0, -1.0
    for lag in range(-max_lag, max_lag + 1):
        if lag >= 0:
            seg_a, seg_b = a[lag:], b[:n - lag]
        else:
            seg_a, seg_b = a[:n + lag], b[-lag:]
        if seg_a.size < 2:
            continue
        num = float(np.dot(seg_a, seg_b))
        d = float(np.linalg.norm(seg_a) * np.linalg.norm(seg_b))
        score = num / d if d > 1e-20 else 0.0
        if score > best_score:
            best_score, best_lag = score, lag
    return best_lag, float(best_score)


def local_similarity(window_a: Any, window_b: Any) -> float:
    """Normalized correlation coefficient between two equal-length windows.

    Simpler/cheaper than `cross_correlation` - no lag search, used inside
    the sliding-window coherence computation where windows are assumed
    already roughly aligned.
    """
    a = np.asarray(window_a, dtype=np.float64)
    b = np.asarray(window_b, dtype=np.float64)
    n = min(a.size, b.size)
    if n < 2:
        return 0.0
    a, b = a[:n] - a[:n].mean(), b[:n] - b[:n].mean()
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom < 1e-20:
        return 0.0
    return float(np.dot(a, b) / denom)


def semblance(windows: Any) -> float:
    """Classic multi-trace semblance for a stack of equal-length windows.

    `windows` has shape (n_traces_local, window_len) - short, same-length
    excerpts from several neighbouring traces at the same sample position.
    Returns a value in [0, 1]: 1.0 means the traces are perfectly aligned
    and in phase at this window, 0.0 means they're incoherent.
    """
    w = np.asarray(windows, dtype=np.float64)
    if w.ndim != 2 or w.shape[0] < 2 or w.shape[1] < 1:
        return 0.0
    n_traces = w.shape[0]
    stack = np.sum(w, axis=0)
    num = float(np.sum(stack ** 2))
    den = float(n_traces * np.sum(w ** 2))
    if den < 1e-20:
        return 0.0
    return float(np.clip(num / den, 0.0, 1.0))


def stack_neighbours(
    data: Any,
    center_idx: int,
    radius: int = 2,
    shift_samples: dict[int, int] | None = None,
) -> np.ndarray:
    """Average the traces around `center_idx` (inclusive of it) within
    `radius`, optionally applying a per-neighbour integer sample shift
    first (e.g. a moveout correction from `cross_correlation`).

    `shift_samples` maps trace_index -> shift (positive = shift trace
    forward in time before stacking).
    """
    d = np.asarray(data, dtype=np.float64)
    n_traces, n_samples = d.shape
    lo = max(0, center_idx - radius)
    hi = min(n_traces, center_idx + radius + 1)

    acc = np.zeros(n_samples, dtype=np.float64)
    count = 0
    for i in range(lo, hi):
        tr = d[i]
        shift = 0 if shift_samples is None else int(shift_samples.get(i, 0))
        if shift:
            tr = np.roll(tr, shift)
        acc += tr
        count += 1
    if count == 0:
        return np.zeros(n_samples, dtype=np.float32)
    return (acc / count).astype(np.float32)


# -----------------------------------------------------------------------------
# Per-trace, per-sample coherence curve
# -----------------------------------------------------------------------------

def compute_trace_coherence(
    data: Any,
    trace_idx: int,
    radius: int = 2,
    window: int = 21,
    align: bool = True,
    max_shift: int = 30,
) -> np.ndarray:
    """Per-sample semblance curve for one trace against its neighbours.

    For every sample, takes a short window centred on it from
    `trace_idx` and from each neighbour within `radius`, and computes
    their semblance.

    If `align=True` (default), each neighbour is first shifted by a
    constant lag - found once via whole-trace `cross_correlation` against
    `trace_idx` - before computing semblance. This matters a lot on real
    data: real first-break moveout across even 2-3 receiver spacings is
    often several samples to tens of samples, and comparing UNALIGNED
    windows means comparing genuinely different parts of the waveform
    across traces, which actively degrades the fused pick rather than
    reinforcing it (this was the likely cause of Picker V2 underperforming
    a plain single-feature picker on real, noisier data - see picker
    module docs). The alignment here is a single constant shift per
    neighbour (not time-varying), which is an approximation, but is far
    better than none.

    Returns a curve normalized to [0, 1], same length as the trace.
    """
    d = np.asarray(data, dtype=np.float64)
    n_traces, n_samples = d.shape
    if n_traces < 2 or n_samples < 3:
        return np.zeros(n_samples, dtype=np.float32)

    lo = max(0, trace_idx - radius)
    hi = min(n_traces, trace_idx + radius + 1)
    neighbours = [i for i in range(lo, hi)]
    if len(neighbours) < 2:
        return np.zeros(n_samples, dtype=np.float32)

    win = max(3, int(window))
    if win % 2 == 0:
        win += 1
    half = win // 2

    block = d[neighbours].copy()  # (n_neighbours, n_samples)

    if align:
        ref = d[trace_idx]
        for row, nb_idx in enumerate(neighbours):
            if nb_idx == trace_idx:
                continue
            lag, score = cross_correlation(ref, d[nb_idx], max_lag=max_shift)
            # Only trust the alignment if it's a genuinely good match;
            # a weak/noisy correlation is more likely to misalign than help.
            if score > 0.3 and lag != 0:
                block[row] = np.roll(block[row], lag)

    padded = np.pad(block, ((0, 0), (half, half)), mode="reflect")

    # Vectorized semblance for every sample at once, instead of a Python
    # loop calling semblance() per sample (that loop is what made this
    # hang on real traces with thousands of samples - O(n_samples) Python-
    # level calls per trace, times every trace in the gather).
    windows_all = sliding_window_view(padded, win, axis=1)  # (n_neighbours, n_samples, win)
    stack = np.sum(windows_all, axis=0)                     # (n_samples, win)
    num = np.sum(stack ** 2, axis=1)                         # (n_samples,)
    den = float(len(neighbours)) * np.sum(windows_all ** 2, axis=(0, 2))  # (n_samples,)
    den = np.where(den > 1e-20, den, 1e-20)
    coh = np.clip(num / den, 0.0, 1.0)

    return normalize_probability(coh)


# -----------------------------------------------------------------------------
# Combining coherence with the fused likelihood curve
# -----------------------------------------------------------------------------

def coherence_weight(coherence: Any, sharpen: float = 1.0) -> np.ndarray:
    """Turn a raw [0, 1] coherence curve into a reinforcement multiplier.

    `sharpen` > 1 makes only strongly-coherent samples count (power-law
    sharpening); `sharpen` = 1 leaves the curve unchanged.
    """
    c = np.asarray(coherence, dtype=np.float64)
    if c.size == 0:
        return np.asarray([], dtype=np.float32)
    c = np.clip(c, 0.0, 1.0)
    if sharpen != 1.0:
        c = c ** max(0.1, float(sharpen))
    return c.astype(np.float32)


def update_probability(probability: Any, coherence: Any, weight: float = 0.75) -> np.ndarray:
    """Reinforce a per-trace likelihood curve with a coherence curve.

    Multiplicative, one-sided boost: `updated = probability * (1 + weight
    * coherence)`, renormalized to [0, 1]. Coherence can only raise a
    sample's relative standing, never suppress one that already had
    strong single-trace evidence (see module docstring).
    """
    p = np.asarray(probability, dtype=np.float64)
    c = np.asarray(coherence, dtype=np.float64)
    if p.size == 0:
        return np.asarray([], dtype=np.float32)
    if c.size != p.size:
        # Coherence unavailable/mismatched - probability passes through
        # unchanged rather than raising, so callers can always compute
        # this even for a gather too small to have neighbours.
        return normalize_probability(p)

    weight = float(np.clip(weight, 0.0, 1.0))
    boosted = p * (1.0 + weight * np.clip(c, 0.0, 1.0))
    return normalize_probability(boosted)


def compute_coherence_probability(
    data: Any,
    probabilities: list[np.ndarray],
    radius: int = 2,
    window: int = 21,
    weight: float = 0.75,
    sharpen: float = 1.0,
    align: bool = True,
    max_shift: int = 30,
) -> list[np.ndarray]:
    """Full per-gather entry point.

    Parameters
    ----------
    data : Any
        Trace matrix (n_traces, n_samples) - the same preprocessed data
        the per-trace likelihood curves in `probabilities` were computed
        from.
    probabilities : list[np.ndarray]
        One fused likelihood curve per trace (from
        `src.picker.likelihood.compute_arrival_probability`), same order
        as `data`'s rows.
    align : bool, optional
        Shift-align neighbours before computing semblance (see
        `compute_trace_coherence`) - matters a lot on real data with
        real moveout across the receiver spread.

    Returns
    -------
    list[np.ndarray]
        One coherence-reinforced likelihood curve per trace, same order.
    """
    d = np.asarray(data, dtype=np.float64)
    n_traces = d.shape[0]
    if len(probabilities) != n_traces:
        raise ValueError(
            f"probabilities has {len(probabilities)} entries but data has {n_traces} traces."
        )

    updated: list[np.ndarray] = []
    for i in range(n_traces):
        coh = compute_trace_coherence(d, i, radius=radius, window=window, align=align, max_shift=max_shift)
        coh = coherence_weight(coh, sharpen=sharpen)
        updated.append(update_probability(probabilities[i], coh, weight=weight))
    return updated


__all__ = [
    "cross_correlation", "local_similarity", "semblance", "stack_neighbours",
    "compute_trace_coherence", "coherence_weight", "update_probability",
    "compute_coherence_probability",
]