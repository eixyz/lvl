"""
Picker Engine V2 - main entry point.

`FirstBreakPicker` is the *only* public interface to the v2 pipeline
(preprocessing -> features -> likelihood -> coherence -> optimizer ->
confidence -> quality). GUI and CLI code should call this class, not the
individual pipeline modules directly - that's what keeps the algorithm
itself swappable and improvable without touching either front-end
(design principle #2 in the picker v2 spec).

This intentionally does NOT touch matplotlib, PyQt, tkinter, or any GUI
state - it takes plain arrays in and returns `PickResult`s out, so it can
be called the same way from the desktop app, the CLI, and eventually a
web backend.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from src.picker import coherence, confidence, likelihood, optimizer, preprocessing, quality
from src.picker.features import extract_features
from src.picker.models import FeatureSet, PickResult
from src.picker.settings import PickerSettings


class FirstBreakPicker:
    """Modular, probabilistic first-break picker.

    Examples
    --------
    >>> picker = FirstBreakPicker()
    >>> result = picker.pick_trace(samples, dt_s)
    >>> results = picker.pick_profile(shot_gather, dt_s, offsets_m=offsets)
    """

    def __init__(self, settings: PickerSettings | None = None):
        self.settings = settings or PickerSettings()

    # -- pipeline stages (exposed individually for testing/introspection) --

    def preprocess(self, samples: Any, dt_s: float, **kwargs) -> np.ndarray:
        return preprocessing.preprocess_trace(samples, dt_s, **kwargs)

    def extract_features(self, samples: Any, dt_s: float) -> FeatureSet:
        return extract_features(
            samples, dt_s,
            sta_ms=self.settings.sta_window,
            lta_ms=self.settings.lta_window,
        )

    def compute_probability(self, feature_set: FeatureSet, adaptive: bool = True) -> np.ndarray:
        return likelihood.compute_arrival_likelihood(feature_set, self.settings, adaptive=adaptive)

    def compute_coherence(self, data: Any, trace_idx: int, radius: int | None = None,
                          window: int = 21) -> np.ndarray:
        radius = self.settings.coherence_radius if radius is None else radius
        raw = coherence.compute_trace_coherence(data, trace_idx, radius=radius, window=window)
        return coherence.coherence_weight(raw)

    def optimize(self, probabilities: list[np.ndarray], **kwargs) -> tuple[list[int | None], float]:
        return optimizer.optimize_profile(probabilities, settings=self.settings, **kwargs)

    def estimate_confidence(self, feature_set: FeatureSet, probability: Any, sample: int | None,
                            candidates=None, coherence_curve=None) -> float:
        return confidence.estimate_confidence(feature_set, probability, sample, candidates, coherence_curve)

    def quality_control(self, picks: list[int | None], offsets_m: list[float], dt_s: float, **kwargs):
        return quality.quality_report(picks, offsets_m, dt_s, **kwargs)

    # -- public API -------------------------------------------------------

    def pick_trace(self, samples: Any, dt_s: float, do_preprocess: bool = True,
                   adaptive: bool = True) -> PickResult:
        """Pick a single, isolated trace - no cross-trace coherence or
        gather-wide path optimization (use `pick_profile` for that).
        """
        x = self.preprocess(samples, dt_s) if do_preprocess else np.asarray(samples, dtype=np.float32)
        fset = self.extract_features(x, dt_s)
        prob = self.compute_probability(fset, adaptive=adaptive)

        if prob.size == 0:
            return PickResult(method="v2_single_trace", features=fset, quality_flags=["NO_PICK"])

        candidates = likelihood.compute_candidate_peaks(prob, max_candidates=5)
        sample, score = candidates[0] if candidates else (None, 0.0)
        conf_value = self.estimate_confidence(fset, prob, sample, candidates)
        flags = confidence.quality_flags(fset, sample, conf_value, candidates, None, self.settings)

        return PickResult(
            sample=sample,
            time=(sample * dt_s if sample is not None else None),
            confidence=conf_value,
            likelihood=score,
            arrival_probability=prob,
            candidate_samples=[s for s, _ in candidates if s is not None],
            candidate_scores=[sc for _, sc in candidates],
            quality_flags=flags,
            features=fset,
            method="v2_single_trace",
        )

    def pick_profile(
        self,
        data: Any,
        dt_s: float,
        offsets_m: list[float] | None = None,
        do_preprocess: bool = True,
        adaptive: bool = True,
        use_coherence: bool = True,
    ) -> list[PickResult]:
        """Pick every trace in one shot gather, reinforced by cross-trace
        coherence and a gather-wide Viterbi path (see `src.picker.optimizer`).

        If `offsets_m` is given (one value per trace), gather-level
        physical-consistency flags (velocity gate, crossing picks, ...)
        are attached to each `PickResult` too.
        """
        d = np.asarray(data, dtype=np.float64)
        if d.ndim != 2:
            raise ValueError(f"pick_profile expects a (n_traces, n_samples) array, got shape {d.shape}.")
        n_traces = d.shape[0]

        if do_preprocess:
            d = np.array([self.preprocess(d[i], dt_s) for i in range(n_traces)])

        feature_sets: list[FeatureSet] = []
        probabilities: list[np.ndarray] = []
        for i in range(n_traces):
            fset = self.extract_features(d[i], dt_s)
            feature_sets.append(fset)
            probabilities.append(self.compute_probability(fset, adaptive=adaptive))

        coherences: list[Any] = [None] * n_traces
        if use_coherence and n_traces > 1:
            coherences = [self.compute_coherence(d, i) for i in range(n_traces)]
            probabilities = [
                coherence.update_probability(probabilities[i], coherences[i], weight=self.settings.coherence_weight)
                for i in range(n_traces)
            ]

        # Candidates are recomputed with the same defaults `optimize()` uses
        # internally, so they describe the exact curve the chosen path came
        # from (kept as a separate call only so confidence scoring below can
        # see the full candidate list, not just the winning sample).
        candidates = optimizer.find_local_maxima(probabilities)
        picks, path_cost = self.optimize(probabilities)

        results: list[PickResult] = []
        for i in range(n_traces):
            sample = picks[i]
            score = next((sc for s, sc in candidates[i] if s == sample), 0.0)
            coh_curve = coherences[i]
            conf_value = self.estimate_confidence(
                feature_sets[i], probabilities[i], sample, candidates[i], coh_curve
            )
            flags = confidence.quality_flags(
                feature_sets[i], sample, conf_value, candidates[i], coh_curve, self.settings
            )
            results.append(PickResult(
                sample=sample,
                time=(sample * dt_s if sample is not None else None),
                confidence=conf_value,
                likelihood=score,
                arrival_probability=probabilities[i],
                candidate_samples=[s for s, _ in candidates[i] if s is not None],
                candidate_scores=[sc for _, sc in candidates[i]],
                quality_flags=list(flags),
                features=feature_sets[i],
                method="v2_profile",
                metadata={"path_cost": path_cost, "trace_index": i},
            ))

        if offsets_m is not None and len(offsets_m) == n_traces:
            report = self.quality_control(picks, offsets_m, dt_s)
            for i, r in enumerate(report):
                for f in r["flags"]:
                    results[i].add_flag(f)
                results[i].metadata["gather_valid"] = r["valid"]

        return results

    def pick_dataset(
        self,
        gathers: dict[str, Any],
        dt_s_map: dict[str, float],
        offsets_map: dict[str, list[float]] | None = None,
        **kwargs,
    ) -> dict[str, list[PickResult]]:
        """Pick every shot gather in a dataset (e.g. all shots on a profile).

        `gathers` maps an arbitrary key (shot id, filename, ...) to a
        (n_traces, n_samples) array; `dt_s_map` gives the sample interval
        for each key (SEG2 files on the same profile are usually all the
        same, but this doesn't assume that).
        """
        out: dict[str, list[PickResult]] = {}
        for key, data in gathers.items():
            if key not in dt_s_map:
                raise ValueError(f"No dt_s given for gather '{key}'.")
            offsets_m = (offsets_map or {}).get(key)
            out[key] = self.pick_profile(data, dt_s_map[key], offsets_m=offsets_m, **kwargs)
        return out


__all__ = ["FirstBreakPicker"]