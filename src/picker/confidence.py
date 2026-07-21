"""
Picker Engine V2 - per-pick confidence estimation.

Everything upstream (features -> likelihood -> coherence -> optimizer)
produces one chosen sample per trace. This module answers a different
question: *how much should you trust that sample*? It looks at how sharp
and unambiguous the winning peak was, how well independent features
agreed with each other, and how well-supported it was by neighbouring
traces - then reduces all of that to one confidence number plus a small
set of human-readable quality flags.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import warnings
from scipy.signal import peak_prominences as _scipy_peak_prominences, peak_widths as _scipy_peak_widths

from src.picker.models import FeatureSet
from src.picker.settings import PickerSettings

# Flag names (kept as plain strings, per the picker v2 spec)
HIGH_CONFIDENCE = "HIGH_CONFIDENCE"
LOW_SNR = "LOW_SNR"
LOW_COHERENCE = "LOW_COHERENCE"
MULTIPLE_PEAKS = "MULTIPLE_PEAKS"
WEAK_AIC = "WEAK_AIC"
LOW_STALTA = "LOW_STALTA"
MANUAL_REVIEW = "MANUAL_REVIEW"

# Default thresholds. All operate on already-normalized [0, 1] curves, so
# these are dimensionless and shouldn't need per-survey retuning often.
_HIGH_CONFIDENCE_THRESHOLD = 0.7
_LOW_SNR_THRESHOLD = 0.3
_LOW_COHERENCE_THRESHOLD = 0.3
_MULTIPLE_PEAKS_RATIO = 0.8
_WEAK_AIC_THRESHOLD = 0.3
_LOW_STALTA_THRESHOLD = 0.3


def _sample_at(curve: Any, sample: int | None) -> float:
    """Safe curve[sample] lookup: out-of-range or missing -> 0.0."""
    if curve is None or sample is None:
        return 0.0
    arr = np.asarray(curve)
    if arr.size == 0 or not (0 <= sample < arr.size):
        return 0.0
    return float(arr[sample])


# -----------------------------------------------------------------------------
# Individual sub-scores
# -----------------------------------------------------------------------------

def peak_prominence(probability: Any, sample: int | None) -> float:
    """Scipy prominence of the curve at `sample`, 0.0 if not a local peak."""
    if sample is None:
        return 0.0
    arr = np.asarray(probability, dtype=np.float64)
    if arr.size == 0 or not (0 <= sample < arr.size):
        return 0.0
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            prom, _, _ = _scipy_peak_prominences(arr, [int(sample)])
        return float(prom[0])
    except Exception:
        return 0.0


def peak_width(probability: Any, sample: int | None, rel_height: float = 0.5) -> float:
    """Width (in samples) of the peak at `sample`, at `rel_height` of its prominence."""
    if sample is None:
        return 0.0
    arr = np.asarray(probability, dtype=np.float64)
    if arr.size == 0 or not (0 <= sample < arr.size):
        return 0.0
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            widths, *_ = _scipy_peak_widths(arr, [int(sample)], rel_height=rel_height)
        return float(widths[0])
    except Exception:
        return 0.0


def second_peak_ratio(candidates: list[tuple[int | None, float]], chosen_sample: int | None) -> float:
    """Score of the best *other* candidate divided by the chosen one's score.

    Close to 1.0 means the pick was ambiguous (a near-tied second peak
    elsewhere on the trace); close to 0.0 means it clearly dominated.
    """
    if not candidates or chosen_sample is None:
        return 0.0
    chosen_score = next((sc for s, sc in candidates if s == chosen_sample), None)
    if chosen_score is None or chosen_score <= 1e-9:
        return 0.0
    others = [sc for s, sc in candidates if s != chosen_sample]
    if not others:
        return 0.0
    return float(np.clip(max(others) / chosen_score, 0.0, 1.0))


def feature_agreement(feature_set: FeatureSet, sample: int | None, tolerance_samples: int = 10) -> float:
    """Fraction of the individually-computed features whose own maximum
    falls within `tolerance_samples` of the chosen sample.

    Independent agreement is one of the strongest confidence signals: a
    pick that AIC, STA/LTA and the Hilbert envelope all "agree" on is far
    more trustworthy than one only one feature liked.
    """
    if sample is None:
        return 0.0
    available = [arr for arr in feature_set.as_dict().values() if arr is not None]
    if not available:
        return 0.0
    agree = 0
    for arr in available:
        a = np.asarray(arr)
        if a.size == 0:
            continue
        if abs(int(np.argmax(a)) - int(sample)) <= tolerance_samples:
            agree += 1
    return agree / len(available)


def coherence_score(coherence: Any, sample: int | None) -> float:
    """Coherence curve value at `sample` (already normalized [0, 1])."""
    return _sample_at(coherence, sample)


def snr_score(feature_set: FeatureSet, sample: int | None) -> float:
    """Normalized SNR feature value at `sample`."""
    return _sample_at(feature_set.snr, sample)


# -----------------------------------------------------------------------------
# Fused confidence + flags
# -----------------------------------------------------------------------------

def estimate_confidence(
    feature_set: FeatureSet,
    probability: Any,
    sample: int | None,
    candidates: list[tuple[int | None, float]] | None = None,
    coherence: Any = None,
    agreement_tolerance_samples: int = 10,
) -> float:
    """Combine peak sharpness, feature agreement, coherence and SNR into
    one confidence value in [0, 1].
    """
    if sample is None:
        return 0.0

    arr = np.asarray(probability, dtype=np.float64)
    prom = peak_prominence(arr, sample)
    prom_norm = float(np.clip(prom, 0.0, 1.0))  # curve is already in [0,1], so prominence is too

    sp_ratio = second_peak_ratio(candidates or [], sample)
    ambiguity_score = 1.0 - sp_ratio  # unambiguous -> high score

    agreement = feature_agreement(feature_set, sample, agreement_tolerance_samples)
    coh_score = coherence_score(coherence, sample) if coherence is not None else 0.5
    snr = snr_score(feature_set, sample)

    # Equal-weighted blend of five independent signals. Simple by design
    # (principle: thresholds/weights should be easy to retune per-survey
    # without touching the rest of the engine).
    parts = [prom_norm, ambiguity_score, agreement, coh_score, snr]
    return float(np.clip(np.mean(parts), 0.0, 1.0))


def quality_flags(
    feature_set: FeatureSet,
    sample: int | None,
    confidence: float,
    candidates: list[tuple[int | None, float]] | None = None,
    coherence: Any = None,
    settings: PickerSettings | None = None,
) -> list[str]:
    """Human-readable quality flags for one pick, per the picker v2 spec."""
    settings = settings or PickerSettings()
    flags: list[str] = []

    if sample is None:
        flags.append(MANUAL_REVIEW)
        return flags

    if confidence >= _HIGH_CONFIDENCE_THRESHOLD:
        flags.append(HIGH_CONFIDENCE)

    if snr_score(feature_set, sample) < _LOW_SNR_THRESHOLD:
        flags.append(LOW_SNR)

    if coherence is not None and coherence_score(coherence, sample) < _LOW_COHERENCE_THRESHOLD:
        flags.append(LOW_COHERENCE)

    if second_peak_ratio(candidates or [], sample) > _MULTIPLE_PEAKS_RATIO:
        flags.append(MULTIPLE_PEAKS)

    if feature_set.aic is not None and _sample_at(feature_set.aic, sample) < _WEAK_AIC_THRESHOLD:
        flags.append(WEAK_AIC)

    if feature_set.stalta is not None and _sample_at(feature_set.stalta, sample) < _LOW_STALTA_THRESHOLD:
        flags.append(LOW_STALTA)

    if confidence < settings.minimum_confidence:
        flags.append(MANUAL_REVIEW)

    return flags


__all__ = [
    "HIGH_CONFIDENCE", "LOW_SNR", "LOW_COHERENCE", "MULTIPLE_PEAKS",
    "WEAK_AIC", "LOW_STALTA", "MANUAL_REVIEW",
    "peak_prominence", "peak_width", "second_peak_ratio", "feature_agreement",
    "coherence_score", "snr_score", "estimate_confidence", "quality_flags",
]