"""
Picker Engine V2 - tunable settings.

`PickerSettings` bundles the feature weights and processing parameters the
v2 engine needs. Defaults are seeded from the existing legacy picking
constants in `src.common.settings` where a direct equivalent already
exists (STA/LTA windows, Hilbert onset threshold, ...), so the new engine
starts from behaviour close to the current picker rather than from
arbitrary numbers. New parameters that have no legacy equivalent
(feature weights, coherence radius, smoothness/jump penalties, minimum
confidence) get sensible defaults that should be tuned against real data.
"""
from __future__ import annotations

from dataclasses import dataclass

from src.common.settings import (
    STA_MS, LTA_MS, HILBERT_ONSET_PCT,
    AUTO_USE_VELOCITY_GATE, AUTO_VMIN_M_S, AUTO_VMAX_M_S, AUTO_GATE_PAD_MS,
)


@dataclass
class PickerSettings:
    # -- Feature weights (relative importance when fusing features into a
    #    single arrival likelihood in likelihood.py). These are normalized
    #    internally, so only relative magnitudes matter. ------------------
    hilbert_weight: float = 1.0
    stalta_weight: float = 1.0
    aic_weight: float = 1.0
    gradient_weight: float = 0.5
    energy_weight: float = 0.5
    snr_weight: float = 0.5
    kurtosis_weight: float = 0.25
    skewness_weight: float = 0.25
    coherence_weight: float = 0.75

    # -- Processing parameters -----------------------------------------
    sta_window: float = STA_MS            # ms
    lta_window: float = LTA_MS            # ms
    hilbert_onset_pct: float = HILBERT_ONSET_PCT
    use_velocity_gate: bool = AUTO_USE_VELOCITY_GATE
    vmin_m_s: float = AUTO_VMIN_M_S
    vmax_m_s: float = AUTO_VMAX_M_S
    gate_pad_ms: float = AUTO_GATE_PAD_MS

    coherence_radius: int = 2             # number of neighbouring traces each side
    coherence_align: bool = True          # shift-align neighbours before semblance (real moveout)
    coherence_max_shift: int = 30         # samples; cap on the alignment search
    smoothness_penalty: float = 0.1       # dynamic-programming path smoothness cost
    jump_penalty: float = 0.2             # cost for large sample-to-sample pick jumps
    minimum_confidence: float = 0.3       # picks below this are flagged MANUAL_REVIEW

    def feature_weights(self) -> dict[str, float]:
        """Feature name -> weight, matching `FeatureSet` field names."""
        return {
            "hilbert": self.hilbert_weight,
            "stalta": self.stalta_weight,
            "aic": self.aic_weight,
            "gradient": self.gradient_weight,
            "energy": self.energy_weight,
            "snr": self.snr_weight,
            "kurtosis": self.kurtosis_weight,
            "skewness": self.skewness_weight,
            "coherence": self.coherence_weight,
        }


__all__ = ["PickerSettings"]