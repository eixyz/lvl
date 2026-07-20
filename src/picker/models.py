from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np


@dataclass(slots=True)
class FeatureSet:
    """
    Normalized feature vectors.
    Every feature has exactly the same length as the trace.
    Values are normalized to [0,1].
    """

    hilbert: np.ndarray
    stalta: np.ndarray
    aic: np.ndarray
    gradient: np.ndarray
    energy: np.ndarray
    coherence: np.ndarray | None = None


@dataclass(slots=True)
class PickResult:
    """
    Final picker output.
    """

    sample: int
    time: float

    confidence: float

    likelihood: np.ndarray

    features: FeatureSet

    method: str = "Likelihood"

    metadata: dict = field(default_factory=dict)

@dataclass(slots=True)
class PickerSettings:

    hilbert_weight: float = 1.0
    stalta_weight: float = 1.0
    aic_weight: float = 1.5
    gradient_weight: float = 0.8
    energy_weight: float = 0.8
    coherence_weight: float = 2.0

    smoothness_penalty: float = 0.15

    coherence_radius: int = 2

    min_confidence: float = 0.35