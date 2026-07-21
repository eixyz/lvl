"""
Picker Engine V2 - path optimization across a shot gather.

`src.picker.likelihood` scores one trace at a time; `src.picker.coherence`
lets neighbours reinforce (but never override) each other. This module is
where a genuinely bad trace (dead channel, high noise, wrong polarity...)
gets pulled back into line: real first breaks form a smooth moveout curve
across offset, so the best pick per trace is chosen jointly, as the path
through all traces' candidate peaks that best balances "match the local
evidence" against "stay close to your neighbours" (Viterbi / dynamic
programming over per-trace candidate lists).
"""
from __future__ import annotations

import math
from typing import Any, Callable

import numpy as np

from src.picker.settings import PickerSettings

Candidate = tuple[int | None, float]


# -----------------------------------------------------------------------------
# Candidate generation
# -----------------------------------------------------------------------------

def find_local_maxima(
    probabilities: list[Any],
    min_prominence: float = 0.05,
    min_distance_samples: int = 5,
    max_candidates: int = 5,
) -> list[list[Candidate]]:
    """Batch version of `likelihood.compute_candidate_peaks` over a whole gather.

    Returns one candidate list per trace; a trace with an empty/invalid
    probability curve gets `[(None, 0.0)]` (a single "no pick" placeholder)
    rather than an empty list, so downstream DP code doesn't need a
    separate empty-trace special case.
    """
    from src.picker.likelihood import compute_candidate_peaks

    out: list[list[Candidate]] = []
    for p in probabilities:
        arr = np.asarray(p, dtype=np.float64)
        if arr.size == 0:
            out.append([(None, 0.0)])
            continue
        peaks = compute_candidate_peaks(
            arr, min_prominence=min_prominence,
            min_distance_samples=min_distance_samples,
            max_candidates=max_candidates,
        )
        out.append(peaks if peaks else [(None, 0.0)])
    return out


# -----------------------------------------------------------------------------
# Generic DP core + Viterbi-flavoured wrapper
# -----------------------------------------------------------------------------

def dynamic_programming(
    candidates: list[list[Candidate]],
    transition_cost_fn: Callable[[int | None, int | None], float],
    emission_cost_fn: Callable[[int | None, float], float] | None = None,
) -> tuple[list[int | None], float]:
    """Generic min-cost path through per-trace candidate lists.

    `transition_cost_fn(sample_a, sample_b)` scores moving from a chosen
    candidate on one trace to a candidate on the next trace.
    `emission_cost_fn(sample, score)` scores choosing that candidate at
    all (defaults to `-log(score)`, i.e. plain Viterbi emission cost).

    Returns (path, total_cost) where `path[i]` is the chosen sample for
    trace i (or None if that trace has no usable candidate).
    """
    if emission_cost_fn is None:
        emission_cost_fn = lambda s, score: (-math.log(max(float(score), 1e-6)) if s is not None else 5.0)

    n = len(candidates)
    if n == 0:
        return [], 0.0

    dp: list[list[float]] = []
    bp: list[list[int]] = []

    first = candidates[0] or [(None, 0.0)]
    dp.append([emission_cost_fn(s, sc) for s, sc in first])
    bp.append([-1] * len(first))

    for i in range(1, n):
        prev = candidates[i - 1] or [(None, 0.0)]
        cur = candidates[i] or [(None, 0.0)]
        dp_i: list[float] = []
        bp_i: list[int] = []
        for (s_b, sc_b) in cur:
            emit = emission_cost_fn(s_b, sc_b)
            best_cost, best_j = None, 0
            for j, (s_a, _sc_a) in enumerate(prev):
                trans = transition_cost_fn(s_a, s_b)
                cost = dp[i - 1][j] + trans + emit
                if best_cost is None or cost < best_cost:
                    best_cost, best_j = cost, j
            dp_i.append(best_cost if best_cost is not None else emit)
            bp_i.append(best_j)
        dp.append(dp_i)
        bp.append(bp_i)

    last_idx = int(np.argmin(dp[-1]))
    total_cost = float(dp[-1][last_idx])

    path_idx = [0] * n
    path_idx[-1] = last_idx
    for i in range(n - 1, 0, -1):
        path_idx[i - 1] = bp[i][path_idx[i]]

    path: list[int | None] = []
    for i in range(n):
        ci = candidates[i] or [(None, 0.0)]
        s, _sc = ci[path_idx[i]]
        path.append(s)
    return path, total_cost


def viterbi_path(
    candidates: list[list[Candidate]],
    smoothness_penalty: float = 0.1,
    jump_penalty: float = 0.2,
    jump_threshold: float = 40.0,
) -> tuple[list[int | None], float]:
    """Picking-specific application of `dynamic_programming`.

    Transition cost is a Huber-style penalty on the sample jump between
    consecutive traces: a small constant cost per sample
    (`smoothness_penalty`), plus an extra linear penalty (`jump_penalty`)
    once the jump exceeds `jump_threshold` samples. A gap (either trace
    has no candidate) costs a small fixed amount rather than being
    forbidden, so a genuinely unpickable trace doesn't break the path.
    """
    def transition(s_a: int | None, s_b: int | None) -> float:
        if s_a is None or s_b is None:
            return 1.0
        d = abs(int(s_b) - int(s_a))
        cost = smoothness_penalty * d
        if d > jump_threshold:
            cost += jump_penalty * (d - jump_threshold)
        return cost

    return dynamic_programming(candidates, transition)


# -----------------------------------------------------------------------------
# Path post-processing
# -----------------------------------------------------------------------------

def smooth_path(path_samples: list[int | None], window: int = 3) -> list[int | None]:
    """Median-filter the chosen pick samples across traces (None-aware)."""
    n = len(path_samples)
    if n == 0:
        return []
    vals = np.array([np.nan if s is None else float(s) for s in path_samples])
    win = max(1, int(window))
    if win % 2 == 0:
        win += 1
    half = win // 2

    out = vals.copy()
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        seg = vals[lo:hi]
        seg = seg[~np.isnan(seg)]
        if seg.size:
            out[i] = float(np.median(seg))
    return [None if np.isnan(v) else int(round(v)) for v in out]


def remove_outliers(
    path_samples: list[int | None],
    max_jump_samples: float = 50.0,
    window: int = 5,
) -> list[int | None]:
    """Null out picks that deviate too far from their local neighbourhood
    median (a lightweight physical-plausibility check on the moveout curve).
    """
    n = len(path_samples)
    if n == 0:
        return []
    vals = np.array([np.nan if s is None else float(s) for s in path_samples])
    out = vals.copy()
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
        if abs(vals[i] - med) > max_jump_samples:
            out[i] = np.nan
    return [None if np.isnan(v) else int(round(v)) for v in out]


def interpolate_missing(path_samples: list[int | None]) -> list[int | None]:
    """Linearly interpolate gaps (None); edges are held constant rather
    than extrapolated. Returns all-None unchanged if nothing is valid.
    """
    n = len(path_samples)
    if n == 0:
        return []
    vals = np.array([np.nan if s is None else float(s) for s in path_samples])
    valid = ~np.isnan(vals)
    if not valid.any():
        return [None] * n
    idx = np.arange(n)
    interp = np.interp(idx, idx[valid], vals[valid])  # edge-holds outside the valid range
    return [int(round(v)) for v in interp]


# -----------------------------------------------------------------------------
# Public entry points
# -----------------------------------------------------------------------------

def optimize_pick(
    probability: Any,
    min_prominence: float = 0.05,
    min_distance_samples: int = 5,
) -> Candidate:
    """Single-trace convenience: best candidate on one likelihood curve,
    with no cross-trace context (used e.g. for an isolated/preview pick).
    """
    from src.picker.likelihood import compute_candidate_peaks

    arr = np.asarray(probability, dtype=np.float64)
    if arr.size == 0:
        return None, 0.0
    peaks = compute_candidate_peaks(
        arr, min_prominence=min_prominence,
        min_distance_samples=min_distance_samples, max_candidates=1,
    )
    return peaks[0] if peaks else (None, 0.0)


def optimize_profile(
    probabilities: list[Any],
    settings: PickerSettings | None = None,
    min_prominence: float = 0.05,
    min_distance_samples: int = 5,
    max_candidates: int = 5,
    jump_threshold: float = 40.0,
    smoothing_window: int = 3,
    outlier_window: int = 5,
    max_jump_samples: float | None = None,
) -> tuple[list[int | None], float]:
    """Full per-gather entry point: candidates -> Viterbi path -> smooth ->
    drop outliers -> interpolate gaps.

    Returns (picks, path_cost) where `picks[i]` is the chosen sample index
    for trace i (or None if it truly couldn't be picked/interpolated).
    """
    settings = settings or PickerSettings()

    candidates = find_local_maxima(
        probabilities, min_prominence=min_prominence,
        min_distance_samples=min_distance_samples, max_candidates=max_candidates,
    )
    path, cost = viterbi_path(
        candidates,
        smoothness_penalty=settings.smoothness_penalty,
        jump_penalty=settings.jump_penalty,
        jump_threshold=jump_threshold,
    )
    path = smooth_path(path, window=smoothing_window)
    path = remove_outliers(
        path, max_jump_samples=(max_jump_samples if max_jump_samples is not None else jump_threshold * 2),
        window=outlier_window,
    )
    path = interpolate_missing(path)
    return path, cost


__all__ = [
    "find_local_maxima", "dynamic_programming", "viterbi_path",
    "smooth_path", "remove_outliers", "interpolate_missing",
    "optimize_pick", "optimize_profile",
]