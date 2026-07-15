"""Trace preprocessing: bandpass filtering, gain/AGC, and bulk static shift."""
from __future__ import annotations

from typing import Any

import numpy as np
from scipy.signal import butter as _scipy_butter, sosfiltfilt as _scipy_sosfiltfilt

from src.common.settings import (
    BP_F1, BP_F2, BP_F3, BP_F4, BP_FFT_PAD, BP_REAPPLY, BUTTER_ORDER,
    BULK_SHIFT_MS, GAIN_MODE, AGC_WINDOW_MS, AGC_STAT,
)


def _ormsby_response(freqs: Any,
                     f1: float, f2: float, f3: float, f4: float) -> Any:
    """Trapezoidal Ormsby amplitude response (linear ramps)."""
    H = np.zeros(len(freqs), dtype=float)
    m_lo = (freqs > f1) & (freqs <= f2)
    m_pb = (freqs > f2) & (freqs <= f3)
    m_hi = (freqs > f3) & (freqs < f4)
    H[m_lo] = (freqs[m_lo] - f1) / (f2 - f1)
    H[m_pb] = 1.0
    H[m_hi] = (f4 - freqs[m_hi]) / (f4 - f3)
    return H

def ormsby(trace: Any, dt_s: float,
           f1: float = BP_F1, f2: float = BP_F2,
           f3: float = BP_F3, f4: float = BP_F4,
           pad_pct: float = BP_FFT_PAD,
           reapply: bool = BP_REAPPLY) -> Any:
    """
    Zero-phase Ormsby bandpass filter  (ProMAX-compatible).

    Algorithm
    ---------
    1. Zero-pad to next power-of-2  (pad_pct=0.25 -> 25 % extra = ProMAX default)
    2. FFT -> multiply by trapezoidal amplitude response (zero-phase: real filter)
    3. IFFT -> trim to original length
    4. Optionally repeat once  (reapply=True squares the amplitude response)
    """
    n     = len(trace)
    n_pad = int(n * (1.0 + pad_pct))
    n_fft = 1 << (n_pad - 1).bit_length()          # next power of 2
    spec  = np.fft.rfft(trace.astype(np.float64), n=n_fft)
    freqs = np.fft.rfftfreq(n_fft, d=dt_s)
    H     = _ormsby_response(freqs, f1, f2, f3, f4)
    if reapply:
        H = H ** 2
    filtered = np.fft.irfft(spec * H, n=n_fft)
    return filtered[:n].astype(np.float32)

def apply_ormsby_all(data: Any, dt_s: float) -> Any:
    """Apply Ormsby bandpass to every trace in data (n_traces x n_samples)."""
    out = np.empty_like(data)
    for i in range(data.shape[0]):
        out[i] = ormsby(data[i], dt_s)
    return out

def apply_ormsby_all_params(data: Any, dt_s: float,
                            f1: float, f2: float, f3: float, f4: float) -> Any:
    """Apply Ormsby with explicit frequency parameters to all traces."""
    out = np.empty_like(data)
    for i in range(data.shape[0]):
        out[i] = ormsby(data[i], dt_s, f1=f1, f2=f2, f3=f3, f4=f4)
    return out

def butterworth_bandpass(trace: Any, dt_s: float,
                         low_hz: float, high_hz: float,
                         order: int = BUTTER_ORDER) -> Any:
    """Zero-phase Butterworth bandpass using scipy SOS filter."""
    x = np.asarray(trace, dtype=np.float64)
    if x.size < 8:
        return x.astype(np.float32)

    fs = 1.0 / max(float(dt_s), 1e-12)
    nyq = 0.5 * fs
    lo = max(0.001, float(low_hz))
    hi = min(float(high_hz), nyq * 0.999)
    if not (0.0 < lo < hi < nyq):
        return x.astype(np.float32)

    sos = _scipy_butter(int(max(1, order)), [lo, hi], btype="bandpass", fs=fs, output="sos")
    try:
        y = _scipy_sosfiltfilt(sos, x)
    except Exception:
        y = x
    return np.asarray(y, dtype=np.float32)

def apply_butterworth_all_params(data: Any, dt_s: float,
                                 low_hz: float, high_hz: float,
                                 order: int = BUTTER_ORDER) -> Any:
    """Apply Butterworth bandpass to every trace in data."""
    out = np.empty_like(data)
    for i in range(data.shape[0]):
        out[i] = butterworth_bandpass(data[i], dt_s, low_hz=low_hz, high_hz=high_hz, order=order)
    return out

def butterworth_highpass(trace: Any, dt_s: float, low_hz: float, order: int = BUTTER_ORDER) -> Any:
    """Zero-phase Butterworth high-pass (low-cut)."""
    x = np.asarray(trace, dtype=np.float64)
    if x.size < 8:
        return np.asarray(trace, dtype=np.float32)
    fs = 1.0 / max(float(dt_s), 1e-12)
    nyq = 0.5 * fs
    lo = max(0.001, min(float(low_hz), nyq * 0.999))
    if not (0.0 < lo < nyq):
        return np.asarray(trace, dtype=np.float32)
    sos = _scipy_butter(int(max(1, order)), lo, btype="highpass", fs=fs, output="sos")
    try:
        y = _scipy_sosfiltfilt(sos, x)
    except Exception:
        y = x
    return np.asarray(y, dtype=np.float32)

def butterworth_lowpass(trace: Any, dt_s: float, high_hz: float, order: int = BUTTER_ORDER) -> Any:
    """Zero-phase Butterworth low-pass (high-cut)."""
    x = np.asarray(trace, dtype=np.float64)
    if x.size < 8:
        return np.asarray(trace, dtype=np.float32)
    fs = 1.0 / max(float(dt_s), 1e-12)
    nyq = 0.5 * fs
    hi = max(0.001, min(float(high_hz), nyq * 0.999))
    if not (0.0 < hi < nyq):
        return np.asarray(trace, dtype=np.float32)
    sos = _scipy_butter(int(max(1, order)), hi, btype="lowpass", fs=fs, output="sos")
    try:
        y = _scipy_sosfiltfilt(sos, x)
    except Exception:
        y = x
    return np.asarray(y, dtype=np.float32)

def apply_cutpass_all_params(data: Any, dt_s: float,
                             low_cut_hz: float, high_cut_hz: float,
                             order: int = BUTTER_ORDER) -> Any:
    """Apply low-cut then high-cut (high-pass + low-pass) to all traces."""
    out = np.empty_like(data)
    for i in range(data.shape[0]):
        y = butterworth_highpass(data[i], dt_s, low_hz=low_cut_hz, order=order)
        y = butterworth_lowpass(y, dt_s, high_hz=high_cut_hz, order=order)
        out[i] = y
    return out

def apply_bulk_static(picks_raw: dict) -> dict:
    """
    Apply BULK_SHIFT_MS to a picks dict {trace_idx: time_ms}.
    Call at analysis / export time; JSON always stores RAW picks.
    """
    if abs(BULK_SHIFT_MS) < 1e-9:
        return picks_raw
    return {k: max(0.0, float(v) + BULK_SHIFT_MS) for k, v in picks_raw.items()}

def apply_gain(data: Any, dt_s: float,
               mode: str = GAIN_MODE,
               window_ms: float = AGC_WINDOW_MS,
               stat: str = AGC_STAT) -> Any:
    """
    Apply display-time amplitude normalization per trace.

    Parameters
    ----------
    data : Any
        Trace matrix with shape (n_traces, n_samples).
    dt_s : float
        Sample interval in seconds.
    mode : str, optional
        Gain mode: "none", "norm", or "agc".
    window_ms : float, optional
        AGC moving window length in milliseconds.
    stat : str, optional
        AGC statistic: "rms" (default) or "mean".

    Returns
    -------
    Any
        Gain-adjusted trace matrix as float32.
    """
    mode_l = str(mode).lower()
    out = np.asarray(data, dtype=np.float32).copy()

    if mode_l == "none":
        return out

    if mode_l == "norm":
        for i in range(out.shape[0]):
            mx = float(np.max(np.abs(out[i])))
            if mx > 1e-20:
                out[i] /= mx
        return out

    if mode_l == "agc":
        win = max(3, int((window_ms / 1000.0) / dt_s))
        if win % 2 == 0:
            win += 1
        ker = np.ones(win, dtype=np.float64) / float(win)
        for i in range(out.shape[0]):
            tr = out[i].astype(np.float64)
            if str(stat).lower() == "mean":
                env = np.convolve(np.abs(tr), ker, mode="same")
            else:
                env = np.sqrt(np.convolve(tr * tr, ker, mode="same"))
            env = np.where(env > 1e-12, env, 1e-12)
            out[i] = (tr / env).astype(np.float32)
        return out

    return out
