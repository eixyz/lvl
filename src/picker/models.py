"""
Picker Engine V2 - data structures.

Pure data containers used throughout the picking pipeline
(preprocessing -> features -> likelihood -> coherence -> optimizer ->
confidence -> quality -> picker). No algorithm code lives here.

See docs/picker_v2/ (or the design note it was generated from) for the
overall pipeline description.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class FeatureSet:
    """Per-sample, normalized (0-1) feature curves for one trace.

    Every array has the same length as the input trace. Features that
    aren't computed for a given call are left as None rather than being
    filled with zeros, so downstream code can tell "not computed" apart
    from "computed and equal to zero".
    """

    hilbert: Any = None
    stalta: Any = None
    aic: Any = None
    gradient: Any = None
    energy: Any = None
    snr: Any = None
    kurtosis: Any = None
    skewness: Any = None
    coherence: Any = None

    def as_dict(self) -> dict:
        return {
            "hilbert": self.hilbert,
            "stalta": self.stalta,
            "aic": self.aic,
            "gradient": self.gradient,
            "energy": self.energy,
            "snr": self.snr,
            "kurtosis": self.kurtosis,
            "skewness": self.skewness,
            "coherence": self.coherence,
        }

    def available(self) -> list[str]:
        """Names of the features that were actually computed (not None)."""
        return [name for name, value in self.as_dict().items() if value is not None]


@dataclass
class PickResult:
    """Outcome of picking a single trace."""

    sample: int | None = None
    time: float | None = None                 # seconds, from trace start
    confidence: float | None = None            # 0-1
    likelihood: float | None = None            # arrival probability at `sample`
    arrival_probability: Any = None            # full per-sample probability vector

    candidate_samples: list[int] = field(default_factory=list)
    candidate_scores: list[float] = field(default_factory=list)

    quality_flags: list[str] = field(default_factory=list)
    features: FeatureSet | None = None

    method: str | None = None                  # e.g. "stalta", "hilbert_env", "v2"
    metadata: dict = field(default_factory=dict)

    @property
    def is_valid(self) -> bool:
        return self.sample is not None and self.time is not None

    def add_flag(self, flag: str) -> None:
        if flag not in self.quality_flags:
            self.quality_flags.append(flag)


__all__ = ["FeatureSet", "PickResult"]