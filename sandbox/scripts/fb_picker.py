"""
fb_picker.py -- Clean first-break picker for LVL profiles
=========================================================

Purpose
-------
This module contains only first-break picking functionality:
- SEG2 reading
- Geometry assignment
- Optional Ormsby filtering
- Interactive picking UI
- Session and final JSON persistence

It intentionally excludes refraction-analysis and export logic.
"""

from __future__ import annotations

import argparse
import importlib
import json
import subprocess
import sys
import tkinter as tk
from pathlib import Path
from tkinter import filedialog
from typing import Any


def _ensure(pkg: str, mod: str = "") -> Any:
    """Import a module, installing it into the active environment if missing."""
    module_name = mod or pkg
    try:
        return importlib.import_module(module_name)
    except ImportError:
        print(f"  Installing '{pkg}' ...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "-q"])
        return importlib.import_module(module_name)


np = _ensure("numpy")
_ensure("obspy")
_ensure("matplotlib")
_ensure("scipy")
pd = _ensure("pandas")
_ensure("openpyxl")
_ensure("xlrd")

from obspy import Trace, read as _read_obspy  # type: ignore[import-untyped]
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.widgets import Button, Slider
from scipy.signal import butter as _scipy_butter, hilbert as _scipy_hilbert, sosfiltfilt as _scipy_sosfiltfilt


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CWD = Path(__file__).parent
DATA_DIR = CWD.parent / "data"
OUTPUT_DIR = CWD.parent / "output"

GEOM_FILES: dict = {
    100: DATA_DIR / "geometry100.txt",
    200: DATA_DIR / "geometry200.txt",
}

PROFILES: dict = {
    "120": {"geom": 100, "line_no": 1, "perp_m": 0.0, "shots": "auto"},
    "214_A": {"geom": 100, "line_no": 2, "perp_m": 0.0, "shots": "auto"},
    "150": {"geom": 200, "line_no": 3, "perp_m": 0.0, "shots": "auto"},
    "415_B": {"geom": 200, "line_no": 4, "perp_m": 0.0, "shots": "auto"},
}

# Use geometry/config shot positions by default.
# Set to True only if you want SEG2 SOURCE_LOCATION to override.
USE_SEG2_SHOT_POSITION: bool = False

# Filter setup (Ormsby f1-f2-f3-f4)
BP_F1: float = 2.0
BP_F2: float = 4.0
BP_F3: float = 140.0
BP_F4: float = 180.0
BP_FFT_PAD: float = 0.25
BP_REAPPLY: bool = False
BUTTER_ORDER: int = 4
FILTER_DEBOUNCE_MS: int = 40

# Display setup
THEME: str = "light"
T_MAX_MS: float = 150.0
T_DISPLAY_PRE_MS: float = -10.0
CLIP_FACTOR: float = 2.0
FILTER_ON: bool = True
DEFAULT_FILTER_MODE: str = "butter"  # "none" | "butter" | "ormsby" | "cutpass"
DISPLAY_MODE: str = "wiggle"  # "wiggle" | "vd" | "both"
SHOW_TIMELINES: bool = True
WIGGLE_STRETCH: float = 0.85

# Gain and autopicker
GAIN_MODE: str = "norm"  # "none" | "norm" | "agc"
AGC_WINDOW_MS: float = 200.0
AGC_STAT: str = "rms"  # "mean" | "rms"
STA_MS: float = 3.0
LTA_MS: float = 20.0
STALTA_TRIG: float = 3.0
AUTO_PICK_MODE: str = "hilbert_env"  # "stalta" | "maxdiff_zero" | "hilbert_env"
ZERO_X_SEARCH_DIRECTION: str = "backward"  # "backward" | "forward"
HILBERT_ONSET_PCT: float = 0.06 # Fraction of peak envelope amplitude for onset pick (0.0-1.0)
HILBERT_TARGET: str = "onset"  # "onset" | "peak"
PICK_PROCESS_ORDER: str = "F>G>X"  # "F>G>X" | "G>F>X" | "RAW>X"
MANUAL_SNAP_WIN_MS: float = 6.0 # Manual pick snap window for zero-crossing search (ms)
# Auto-pick first-arrival plausibility gate (absolute time, ms)
AUTO_USE_VELOCITY_GATE: bool = True
AUTO_VMIN_M_S: float = 120.0
AUTO_VMAX_M_S: float = 3500.0
AUTO_GATE_PAD_MS: float = 20.0
# Trigger/radio correction (ms).
# Static term is added to SEG2 delay when data are loaded, so display and picks share the same time basis.
# Offset-dependent term remains as residual correction on picks.
TRIGGER_STATIC_MS: float = 6.0
TRIGGER_MS_PER_M: float = 0.0


# ---------------------------------------------------------------------------
# UI/theme helpers
# ---------------------------------------------------------------------------

def _tc() -> dict:
    if THEME == "dark":
        return {
            "fig_bg": "#0f1117",
            "ax_bg": "#1a1d2e",
            "trace": "#4a8fc1",
            "fill_alpha": 0.20,
            "pick_clr": "#e63946",
            "text": "white",
            "label": "#aaa",
            "grid": "#888",
            "tick": "#888",
            "spine": "#333344",
            "leg_face": "#1e2130",
            "leg_edge": "#555",
        }
    return {
        "fig_bg": "white",
        "ax_bg": "#f5f5f5",
        "trace": "#111111",
        "fill_alpha": 0.18,
        "pick_clr": "#cc0000",
        "text": "#111111",
        "label": "#333",
        "grid": "#bbb",
        "tick": "#333",
        "spine": "#aaaaaa",
        "leg_face": "#f0f0f0",
        "leg_edge": "#aaa",
    }


def _ensure_interactive_backend() -> bool:
    """Ensure matplotlib runs an interactive backend so picker windows can open."""

    def _is_non_interactive(name: str) -> bool:
        n = str(name).strip().lower()
        return n in {
            "agg",
            "pdf",
            "ps",
            "svg",
            "template",
            "cairo",
            "module://matplotlib_inline.backend_inline",
            "module://ipykernel.pylab.backend_inline",
        }

    try:
        backend = str(plt.get_backend()).lower()
    except Exception:
        backend = ""

    if backend and not _is_non_interactive(backend):
        return True

    for candidate in ("TkAgg", "QtAgg", "Qt5Agg"):
        try:
            plt.switch_backend(candidate)
            break
        except Exception:
            continue

    try:
        new_backend = str(plt.get_backend()).lower()
    except Exception:
        new_backend = ""

    if (not new_backend) or _is_non_interactive(new_backend):
        print("  [WARN] Matplotlib backend is non-interactive (Agg).")
        return False

    print(f"  [INFO] Using interactive backend: {plt.get_backend()}")
    return True


# ---------------------------------------------------------------------------
# Geometry and SEG2 IO
# ---------------------------------------------------------------------------

def load_geometry(geom_type: int) -> Any:
    """Load receiver x positions from geometry file."""
    path = GEOM_FILES[geom_type]
    positions: list = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.lower().startswith("trace"):
                continue
            parts = s.replace(",", ".").split()
            if len(parts) >= 2:
                try:
                    positions.append(float(parts[1]))
                except ValueError:
                    pass
    return np.array(positions, dtype=float)


def auto_shot_positions(recv_pos: Any) -> dict:
    """Infer standard 3-shot geometry from receiver spread."""
    n = len(recv_pos)
    dx_start = float(recv_pos[1] - recv_pos[0])
    dx_end = float(recv_pos[-1] - recv_pos[-2])
    mid_lo = float(recv_pos[n // 2 - 1])
    mid_hi = float(recv_pos[n // 2])
    return {
        1: round(float(recv_pos[0]) - dx_start, 4),
        2: round((mid_lo + mid_hi) / 2.0, 4),
        3: round(float(recv_pos[-1]) + dx_end, 4),
    }


def _normalize_profile_token(value: Any) -> str:
    """Normalize profile strings for robust matching (LVL150, 150, 214_A, etc.)."""
    s = str(value or "").strip().upper().replace(" ", "")
    s = s.replace("_", "").replace("-", "")
    if s.startswith("LVL"):
        s = s[3:]
    return s


def discover_profile_folders(data_dir: Path) -> list:
    """Return profile folders under data/ that contain SEG2 files."""
    out: list = []
    for p in sorted(data_dir.iterdir() if data_dir.exists() else []):
        if not p.is_dir():
            continue
        if any(p.glob("*.seg2")) or any(p.glob("*.SEG2")):
            out.append(p.name)
    return out


def discover_lvl_geometry_excels(data_dir: Path) -> list:
    """Find likely LVL geometry excel files, newest first."""
    pats = ["LVL*.xls", "LVL*.xlsx", "*lvl*.xls", "*lvl*.xlsx"]
    out: list = []
    seen: set[str] = set()
    for pat in pats:
        for p in data_dir.glob(pat):
            key = str(p.resolve()).lower()
            if key in seen:
                continue
            seen.add(key)
            out.append(p)
    out.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return out


def _to_float_or_none(v: Any) -> float | None:
    try:
        if pd.isna(v):
            return None
    except Exception:
        pass
    s = str(v).strip().replace(",", ".")
    if not s:
        return None
    try:
        return float(s)
    except Exception:
        return None


def _detect_station_table_columns(df: Any) -> tuple | None:
    """Detect header row and key columns for LVL Number / Station / optional X,Y."""
    if df is None or df.empty:
        return None

    def _norm(v: Any) -> str:
        return " ".join(str(v or "").strip().lower().replace("_", " ").replace("-", " ").split())

    max_scan = min(40, int(df.shape[0]))
    for ridx in range(max_scan):
        row = [_norm(v) for v in df.iloc[ridx].tolist()]
        i_prof = None
        i_sta = None
        i_x = None
        i_y = None
        for ci, txt in enumerate(row):
            if not txt:
                continue
            if i_prof is None and (("lvl" in txt and "number" in txt) or txt in ("lvl", "profile", "line")):
                i_prof = ci
            if i_sta is None and "station" in txt:
                i_sta = ci
            if i_x is None and txt in ("x", "east", "easting"):
                i_x = ci
            if i_y is None and txt in ("y", "north", "northing"):
                i_y = ci
        if None not in (i_prof, i_sta):
            return ridx, i_prof, i_sta, i_x, i_y
    return None


def infer_geometry_from_lvl_excel(data_dir: Path, profile_name: str) -> tuple[int | None, Path | None, float | None]:
    """
    Infer geometry from LVL geometry excel by matching profile + stations.

    Preference:
    - distance between station 0.5 and 48.5 derived from X/Y if available,
    - otherwise station numeric difference.
    """
    target = _normalize_profile_token(profile_name)
    files = discover_lvl_geometry_excels(data_dir)
    for path in files:
        try:
            sheets = pd.read_excel(path, sheet_name=None, header=None, dtype=object)
        except Exception:
            continue

        for _sheet_name, df in sheets.items():
            cols = _detect_station_table_columns(df)
            if cols is None:
                continue
            hdr_row, i_prof, i_sta, i_x, i_y = cols

            rows: list[tuple[float, float | None, float | None]] = []
            for ridx in range(hdr_row + 1, int(df.shape[0])):
                p_raw = df.iat[ridx, i_prof] if i_prof < df.shape[1] else None
                if _normalize_profile_token(p_raw) != target:
                    continue
                s_raw = df.iat[ridx, i_sta] if i_sta < df.shape[1] else None
                sta = _to_float_or_none(s_raw)
                if sta is None:
                    continue
                xv = _to_float_or_none(df.iat[ridx, i_x]) if i_x is not None and i_x < df.shape[1] else None
                yv = _to_float_or_none(df.iat[ridx, i_y]) if i_y is not None and i_y < df.shape[1] else None
                rows.append((float(sta), xv, yv))

            if len(rows) < 2:
                continue

            s0 = min(rows, key=lambda r: abs(r[0] - 0.5))
            s1 = min(rows, key=lambda r: abs(r[0] - 48.5))

            dist = None
            if (s0[1] is not None and s0[2] is not None and s1[1] is not None and s1[2] is not None):
                dist = float(np.hypot(float(s1[1]) - float(s0[1]), float(s1[2]) - float(s0[2])))
            else:
                dist = abs(float(s1[0]) - float(s0[0]))

            if abs(dist - 94.0) <= 12.0:
                return 100, path, dist
            if abs(dist - 192.5) <= 16.0:
                return 200, path, dist

    return None, None, None


def read_seg2(path: Path) -> tuple:
    """
    Read one SEG2 shot.

    Returns:
      data, dt_s, n_traces, n_samples, shot_pos, ffid, delay_ms, recv_locs_m
    """
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        st = _read_obspy(str(path), format="SEG2")

    n_traces = len(st)
    n_samp = max(tr.stats.npts for tr in st)
    dt_s = float(st[0].stats.delta)

    data = np.zeros((n_traces, n_samp), dtype=np.float32)
    for i, tr in enumerate(st):
        npts = tr.stats.npts
        data[i, :npts] = tr.data.astype(np.float32)

    shot_pos: Any = None
    ffid: int = 0
    delay_ms: float = 0.0
    recv_locs_m: Any = None

    try:
        hdr0 = dict(st[0].stats.seg2)

        sloc = str(hdr0.get("SOURCE_LOCATION", "")).strip()
        if sloc:
            shot_pos = float(sloc.replace(",", ".").split()[0])

        ssn = str(hdr0.get("SHOT_SEQUENCE_NUMBER", "")).strip()
        if ssn.isdigit():
            ffid = int(ssn)

        delay_s = str(hdr0.get("DELAY", "0")).strip().replace(",", ".")
        delay_ms = float(delay_s) * 1000.0 + float(TRIGGER_STATIC_MS)

        locs = []
        for tr in st:
            h = dict(tr.stats.seg2)
            rl = str(h.get("RECEIVER_LOCATION", "")).strip().replace(",", ".")
            locs.append(float(rl) if rl else None)
        if all(v is not None for v in locs):
            recv_locs_m = np.array(locs, dtype=float)
    except Exception:
        pass

    if ffid == 0:
        digits = "".join(c for c in path.stem if c.isdigit())
        ffid = int(digits) if digits else 0

    return data, dt_s, n_traces, n_samp, shot_pos, ffid, delay_ms, recv_locs_m


# ---------------------------------------------------------------------------
# Signal processing
# ---------------------------------------------------------------------------

def _ormsby_response(freqs: Any, f1: float, f2: float, f3: float, f4: float) -> Any:
    """Trapezoidal Ormsby amplitude response."""
    h = np.zeros(len(freqs), dtype=float)
    m_lo = (freqs > f1) & (freqs <= f2)
    m_pb = (freqs > f2) & (freqs <= f3)
    m_hi = (freqs > f3) & (freqs < f4)
    h[m_lo] = (freqs[m_lo] - f1) / (f2 - f1)
    h[m_pb] = 1.0
    h[m_hi] = (f4 - freqs[m_hi]) / (f4 - f3)
    return h


def ormsby(
    trace: Any,
    dt_s: float,
    f1: float = BP_F1,
    f2: float = BP_F2,
    f3: float = BP_F3,
    f4: float = BP_F4,
    pad_pct: float = BP_FFT_PAD,
    reapply: bool = BP_REAPPLY,
) -> Any:
    """Zero-phase frequency-domain Ormsby filter."""
    n = len(trace)
    n_pad = int(n * (1.0 + pad_pct))
    n_fft = 1 << (n_pad - 1).bit_length()
    spec = np.fft.rfft(trace.astype(np.float64), n=n_fft)
    freqs = np.fft.rfftfreq(n_fft, d=dt_s)
    h = _ormsby_response(freqs, f1, f2, f3, f4)
    if reapply:
        h = h**2
    filtered = np.fft.irfft(spec * h, n=n_fft)
    return filtered[:n].astype(np.float32)


def apply_ormsby_all(data: Any, dt_s: float) -> Any:
    out = np.empty_like(data)
    for i in range(data.shape[0]):
        out[i] = ormsby(data[i], dt_s)
    return out


def apply_ormsby_all_params(data: Any, dt_s: float, f1: float, f2: float, f3: float, f4: float) -> Any:
    out = np.empty_like(data)
    for i in range(data.shape[0]):
        out[i] = ormsby(data[i], dt_s, f1=f1, f2=f2, f3=f3, f4=f4)
    return out


def butterworth_bandpass(
    trace: Any,
    dt_s: float,
    low_hz: float,
    high_hz: float,
    order: int = BUTTER_ORDER,
) -> Any:
    """Zero-phase Butterworth bandpass using scipy SOS filtering."""
    x = np.asarray(trace, dtype=np.float64)
    if x.size < 8:
        return np.asarray(trace, dtype=np.float32)

    fs = 1.0 / max(float(dt_s), 1e-12)
    nyq = 0.5 * fs
    lo = max(0.001, float(low_hz))
    hi = min(float(high_hz), nyq * 0.999)
    if not (0.0 < lo < hi < nyq):
        return np.asarray(trace, dtype=np.float32)

    sos = _scipy_butter(int(max(1, order)), [lo, hi], btype="bandpass", fs=fs, output="sos")
    try:
        y = _scipy_sosfiltfilt(sos, x)
    except Exception:
        y = x
    return np.asarray(y, dtype=np.float32)


def apply_butterworth_all_params(
    data: Any,
    dt_s: float,
    low_hz: float,
    high_hz: float,
    order: int = BUTTER_ORDER,
) -> Any:
    """Apply Butterworth bandpass to all traces in a shot gather."""
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


def _zero_crossing_from_extremum_samples(
    samples: Any,
    dt_s: float,
    win_start_s: float,
    win_end_s: float,
    search_direction: str = "backward",
    use_abs_peak: bool = True,
) -> float | None:
    """Return zero-crossing time (s) near local extremum in a window."""
    x = np.asarray(samples, dtype=float)
    if x.size < 2:
        return None

    i0 = max(0, int(np.floor(win_start_s / dt_s)))
    i1 = min(x.size - 1, int(np.ceil(win_end_s / dt_s)))
    if i1 <= i0:
        return None

    seg = x[i0:i1 + 1]
    if seg.size == 0:
        return None

    if use_abs_peak:
        rel_idx = int(np.argmax(np.abs(seg)))
    else:
        rel_idx = int(np.argmax(seg))
    peak_idx = i0 + rel_idx

    direction = str(search_direction).strip().lower()
    if direction not in {"backward", "forward"}:
        direction = "backward"

    if direction == "backward":
        if x[peak_idx] == 0.0:
            return peak_idx * dt_s
        for right in range(peak_idx, i0, -1):
            left = right - 1
            y0 = float(x[left])
            y1 = float(x[right])
            if y0 == 0.0:
                return left * dt_s
            if y1 == 0.0:
                return right * dt_s
            if y0 * y1 < 0.0:
                frac = -y0 / (y1 - y0)
                return (left + frac) * dt_s
    else:
        if x[peak_idx] == 0.0:
            return peak_idx * dt_s
        for left in range(peak_idx, i1):
            right = left + 1
            y0 = float(x[left])
            y1 = float(x[right])
            if y0 == 0.0:
                return left * dt_s
            if y1 == 0.0:
                return right * dt_s
            if y0 * y1 < 0.0:
                frac = -y0 / (y1 - y0)
                return (left + frac) * dt_s

    return peak_idx * dt_s


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

        # Use threshold between noise floor and peak for better robustness.
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


def apply_gain(
    data: Any,
    dt_s: float,
    mode: str = GAIN_MODE,
    window_ms: float = AGC_WINDOW_MS,
    stat: str = AGC_STAT,
) -> Any:
    """Display-only gain. Does not alter stored pick times."""
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


# ---------------------------------------------------------------------------
# Picker UI
# ---------------------------------------------------------------------------

class FirstBreakPicker:
    """
    Interactive first-break picker for one shot record.

    Controls
    --------
    Left click / drag: place picks
    Right click / drag: delete picks
    Shift + Left click: range-fill picks
    a: auto-pick using selected mode
    f: cycle filter mode (none/butter/ormsby)
    g: cycle gain (none -> norm -> agc)
    o: cycle auto-pick mode (stalta/maxdiff_zero/hilbert_env)
    v: toggle display polarity
    l: toggle timeline guides
    c: top axis channel/offset
    p: save and go previous shot
    n: save and go next shot
    s: save and finalize profile
    q / Esc: quit
    """

    def __init__(
        self,
        data_raw: Any,
        data_filt: Any,
        dt_s: float,
        recv_abs: Any,
        shot_id: int,
        profile_name: str,
        shot_pos_m: float = 0.0,
        delay_ms: float = 0.0,
        existing_picks: dict | None = None,
        qc_dir: Path | None = None,
        save_callback: Any = None,
        header_info: dict | None = None,
    ):
        self.data_raw = data_raw
        self.data_filt = data_filt
        self.dt_s = dt_s
        self.recv_abs = recv_abs
        self.shot_pos_m = shot_pos_m
        self.delay_ms = delay_ms
        self.shot_id = shot_id
        self.profile = profile_name
        self.qc_dir = qc_dir or (OUTPUT_DIR / profile_name)

        self._picks: dict = {}
        self._saved = False
        self._cancelled = False
        self._done = False
        if FILTER_ON:
            mode = str(DEFAULT_FILTER_MODE).strip().lower()
            self._filter_mode = mode if mode in {"none", "butter", "ormsby", "cutpass"} else "butter"
        else:
            self._filter_mode = "none"
        self._inverted = True
        self._top_mode = "channel"
        self.ax_top: Any = None
        self._save_callback = save_callback
        self._header_info = header_info or {}
        self._nav_action = "stay"
        self._drag_pick = False
        self._drag_delete = False
        self._last_drag_idx: Any = None
        self._last_drag_t: float | None = None
        self._range_anchor: Any = None
        self._gain_mode = GAIN_MODE
        self._agc_window_ms = AGC_WINDOW_MS
        self._agc_stat = AGC_STAT
        self._display_mode = DISPLAY_MODE
        self._show_timelines = SHOW_TIMELINES
        self._wiggle_stretch = WIGGLE_STRETCH
        self._f1 = BP_F1
        self._f2 = BP_F2
        self._f3 = BP_F3
        self._f4 = BP_F4
        self._butter_order = BUTTER_ORDER
        self._filter_debounce_ms = max(0, int(FILTER_DEBOUNCE_MS))
        self._filter_timer: Any = None
        self._auto_pick_mode = str(AUTO_PICK_MODE).strip().lower()
        self._pick_order = str(PICK_PROCESS_ORDER).strip().upper()
        self._hilbert_target = str(HILBERT_TARGET).strip().lower()
        self._hover_help_map: dict[Any, str] = {}

        self.data_ormsby = self.data_filt
        self.data_butter = apply_butterworth_all_params(
            self.data_raw,
            self.dt_s,
            low_hz=self._f2,
            high_hz=self._f3,
            order=self._butter_order,
        )
        self.data_cutpass = apply_cutpass_all_params(
            self.data_raw,
            self.dt_s,
            low_cut_hz=self._f1,
            high_cut_hz=self._f4,
            order=self._butter_order,
        )

        if existing_picks:
            self._picks = {int(k): float(v) for k, v in existing_picks.items()}

        self.n_traces = data_raw.shape[0]
        self.n_samp = data_raw.shape[1]
        dt_ms = dt_s * 1000.0
        self.times_ms = delay_ms + np.arange(self.n_samp) * dt_ms

        self.t_disp_start = max(delay_ms, T_DISPLAY_PRE_MS)
        self.t_disp_end = T_MAX_MS
        self._t_view_start = float(self.t_disp_start)
        self._t_view_end = float(self.t_disp_end)
        self._t_full_start = float(delay_ms)
        self._t_full_end = float(delay_ms + (self.n_samp - 1) * dt_ms)
        self._t_view_end = min(self._t_view_end, self._t_full_end)

        off_range = float(np.ptp(recv_abs)) if len(recv_abs) > 1 else 1.0
        self._dx = off_range / max(len(recv_abs) - 1, 1)

        self._build_figure()

    def _build_figure(self):
        c = _tc()
        self.fig, self.ax = plt.subplots(figsize=(16, 8), constrained_layout=False)
        self.fig.subplots_adjust(left=0.06, right=0.76, bottom=0.13, top=0.91)
        self.fig.patch.set_facecolor(c["fig_bg"])
        self.ax.set_facecolor(c["ax_bg"])
        try:
            self.fig.canvas.manager.set_window_title(
                f"LVL Picker -- {self.profile} Shot {self.shot_id}"
            )
        except Exception:
            pass

        self.fig.canvas.mpl_connect("button_press_event", self._on_click)
        self.fig.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.fig.canvas.mpl_connect("motion_notify_event", self._on_hover_help)
        self.fig.canvas.mpl_connect("button_release_event", self._on_release)
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)

        self._build_controls()
        self._redraw()

    def _build_controls(self):
        c = _tc()
        x0 = 0.78
        w = 0.20
        row_h = 0.032
        gap = 0.005
        row3_w = (w - 2 * gap) / 3.0

        self.ax_btn_prev = self.fig.add_axes([x0 + 0 * (row3_w + gap), 0.88, row3_w, row_h])
        self.ax_btn_next = self.fig.add_axes([x0 + 1 * (row3_w + gap), 0.88, row3_w, row_h])
        self.ax_btn_save = self.fig.add_axes([x0 + 2 * (row3_w + gap), 0.88, row3_w, row_h])

        self.ax_btn_auto = self.fig.add_axes([x0 + 0 * (row3_w + gap), 0.84, row3_w, row_h])
        self.ax_btn_inv = self.fig.add_axes([x0 + 1 * (row3_w + gap), 0.84, row3_w, row_h])
        self.ax_btn_tl = self.fig.add_axes([x0 + 2 * (row3_w + gap), 0.84, row3_w, row_h])
        self.ax_btn_apply = self.fig.add_axes([x0 + 0.125, 0.408, 0.03, 0.024])

        btn_face = c["ax_bg"]
        btn_hover = "#e8e8e8" if THEME == "light" else "#2a2d3d"

        self.btn_prev = Button(self.ax_btn_prev, "Prev", color=btn_face, hovercolor=btn_hover)
        self.btn_next = Button(self.ax_btn_next, "Next", color=btn_face, hovercolor=btn_hover)
        self.btn_save = Button(self.ax_btn_save, "Save/Close", color=btn_face, hovercolor=btn_hover)
        self.btn_auto = Button(self.ax_btn_auto, "Auto Picker", color=btn_face, hovercolor=btn_hover)
        self.btn_inv = Button(self.ax_btn_inv, "", color=btn_face, hovercolor=btn_hover)
        self.btn_tl = Button(self.ax_btn_tl, "Timeline", color=btn_face, hovercolor=btn_hover)
        self.btn_apply = Button(self.ax_btn_apply, "Apply", color=btn_face, hovercolor=btn_hover)
        for btn in (
            self.btn_prev,
            self.btn_save,
            self.btn_next,
            self.btn_auto,
            self.btn_inv,
            self.btn_tl,
            self.btn_apply,
        ):
            btn.label.set_color(c["text"])
            btn.label.set_fontsize(8)
            for sp in btn.ax.spines.values():
                sp.set_edgecolor(c["spine"])

        self._update_invert_button_label()
        self._refresh_toggle_buttons()

        self.btn_prev.on_clicked(lambda _e: self._go_prev())
        self.btn_save.on_clicked(lambda _e: self._save_and_finish("finalize"))
        self.btn_next.on_clicked(lambda _e: self._save_and_finish("next"))
        self.btn_auto.on_clicked(lambda _e: self._auto_then_redraw())
        self.btn_inv.on_clicked(lambda _e: self._toggle_invert())
        self.btn_tl.on_clicked(lambda _e: self._toggle_timelines())
        self.btn_apply.on_clicked(lambda _e: self._apply_current_settings_to_picks())

        self._register_hover(self.ax_btn_prev, "Previous shot")
        self._register_hover(self.ax_btn_next, "Next shot")
        self._register_hover(self.ax_btn_save, "Finalize and save picks")
        self._register_hover(self.ax_btn_auto, "Run auto picker on all traces")
        self._register_hover(self.ax_btn_inv, "Toggle display polarity")
        self._register_hover(self.ax_btn_tl, "Toggle timeline guides")
        self._register_hover(self.ax_btn_apply, "Re-apply current mode/settings to existing picks")

        label_fs = 8
        self.fig.text(x0, 0.80, "Gain Control:", fontsize=label_fs, fontweight="bold", color=c["text"], ha="left", va="bottom")
        self.fig.text(x0, 0.76, "AGC stat:", fontsize=label_fs, fontweight="bold", color=c["text"], ha="left", va="bottom")
        self.fig.text(x0, 0.72, "Display:", fontsize=label_fs, fontweight="bold", color=c["text"], ha="left", va="bottom")
        self.fig.text(x0, 0.620, "Filter mode:", fontsize=label_fs, fontweight="bold", color=c["text"], ha="left", va="bottom")
        self.fig.text(x0, 0.580, "Butter order:", fontsize=label_fs, fontweight="bold", color=c["text"], ha="left", va="bottom")
        self.fig.text(x0, 0.495, "Pick mode:", fontsize=label_fs, fontweight="bold", color=c["text"], ha="left", va="bottom")
        self.fig.text(x0, 0.455, "Reihenfolge:", fontsize=label_fs, fontweight="bold", color=c["text"], ha="left", va="bottom")
        self.fig.text(x0, 0.415, "Hilbert:", fontsize=label_fs, fontweight="bold", color=c["text"], ha="left", va="bottom")

        self._gain_labels = ("none", "norm", "agc")
        self._gain_btns = self._make_mode_buttons(
            x0=x0 + 0.055, y=0.792, w=0.145, h=0.026,
            labels=self._gain_labels, callback=self._on_gain_mode,
            tooltips=("No gain", "Per-trace normalize", "Automatic gain control"),
        )

        self._stat_labels = ("rms", "mean")
        self._stat_btns = self._make_mode_buttons(
            x0=x0 + 0.055, y=0.752, w=0.095, h=0.026,
            labels=self._stat_labels, callback=self._on_agc_stat,
            tooltips=("RMS AGC envelope", "Mean-abs AGC envelope"),
        )

        self._disp_labels = ("wiggle", "vd", "both")
        self._disp_btns = self._make_mode_buttons(
            x0=x0 + 0.055, y=0.712, w=0.145, h=0.026,
            labels=self._disp_labels, callback=self._on_display_mode,
            tooltips=("Wiggle only", "Variable-density only", "Wiggle + variable-density"),
        )

        self._filt_labels = ("none", "butter", "ormsby", "cutpass")
        self._filt_btns = self._make_mode_buttons(
            x0=x0 + 0.055, y=0.620, w=0.145, h=0.026,
            labels=self._filt_labels, callback=self._on_filter_mode,
            tooltips=("No filter", "Butterworth bandpass", "Ormsby bandpass", "Low-cut f1 + high-cut f4"),
        )

        self._butter_order_labels = ("2", "4", "6")
        self._butter_order_btns = self._make_mode_buttons(
            x0=x0 + 0.055, y=0.580, w=0.05, h=0.026,
            labels=self._butter_order_labels, callback=self._on_butter_order,
            tooltips=("Butterworth order 2", "Butterworth order 4", "Butterworth order 6"),
        )

        self._pick_mode_labels = ("stalta", "maxdiff_zero", "hilbert_env")
        self._pick_mode_btns = self._make_mode_buttons(
            x0=x0 + 0.055, y=0.488, w=0.145, h=0.026,
            labels=self._pick_mode_labels, callback=self._on_auto_pick_mode,
            tooltips=("Classic STA/LTA", "Peak then zero-cross onset", "Hilbert envelope onset/peak"),
        )

        self._pick_order_labels = ("F>G>X", "G>F>X", "RAW>X")
        self._pick_order_btns = self._make_mode_buttons(
            x0=x0 + 0.055, y=0.448, w=0.145, h=0.026,
            labels=self._pick_order_labels, callback=self._on_pick_order,
            tooltips=("Filter, then gain, then pick", "Gain, then filter, then pick", "Pick on raw trace"),
        )

        self._hilbert_target_labels = ("onset", "peak")
        self._hilbert_target_btns = self._make_mode_buttons(
            x0=x0 + 0.055, y=0.408, w=0.065, h=0.026,
            labels=self._hilbert_target_labels, callback=self._on_hilbert_target,
            tooltips=("Use first threshold crossing", "Use envelope peak"),
        )

        self.ax_agc = self.fig.add_axes([x0 + 0.02, 0.67, 0.25 * w, 0.016], facecolor=c["ax_bg"])
        self.sl_agc = Slider(
            self.ax_agc,
            "AGC",
            5.0,
            500.0,
            valinit=float(self._agc_window_ms),
            valstep=1.0,
        )
        self.sl_agc.label.set_fontweight("bold")
        self.sl_agc.label.set_fontsize(8)
        self.sl_agc.valtext.set_fontsize(8)
        self.sl_agc.on_changed(self._on_agc_window)

        self.ax_wig = self.fig.add_axes([x0 + 0.125, 0.67, 0.25 * w, 0.016], facecolor=c["ax_bg"])
        self.sl_wig = Slider(
            self.ax_wig,
            "Scale",
            0.20,
            5.00,
            valinit=float(self._wiggle_stretch),
            valstep=0.01,
        )
        self.sl_wig.label.set_fontweight("bold")
        self.sl_wig.label.set_fontsize(8)
        self.sl_wig.valtext.set_fontsize(8)
        self.sl_wig.on_changed(self._on_wiggle_stretch)

        self.fig.text(x0, 0.365, "Time View (ms)", fontsize=8, fontweight="bold", color=c["text"], ha="left", va="bottom")
        tmin_hi = max(self._t_full_start + 1.0, self._t_full_end - 1.0)
        tmax_lo = min(self._t_full_end - 1.0, self._t_full_start + 1.0)
        self.ax_tmin = self.fig.add_axes([x0 + 0.02, 0.335, 0.25 * w, 0.014], facecolor=c["ax_bg"])
        self.ax_tmax = self.fig.add_axes([x0 + 0.125, 0.335, 0.25 * w, 0.014], facecolor=c["ax_bg"])
        self.sl_tmin = Slider(self.ax_tmin, "T min", self._t_full_start, tmin_hi,
                      valinit=float(self._t_view_start), valstep=1.0)
        self.sl_tmax = Slider(self.ax_tmax, "T max", tmax_lo, self._t_full_end,
                      valinit=float(self._t_view_end), valstep=1.0)
        self.sl_tmin.label.set_fontweight("bold")
        self.sl_tmax.label.set_fontweight("bold")
        self.sl_tmin.label.set_fontsize(8)
        self.sl_tmax.label.set_fontsize(8)
        self.sl_tmin.valtext.set_fontsize(8)
        self.sl_tmax.valtext.set_fontsize(8)
        self.sl_tmin.on_changed(self._on_tmin)
        self.sl_tmax.on_changed(self._on_tmax)

        self.ax_f1 = self.fig.add_axes([x0 + 0.0125, 0.550, 0.25 * w, 0.014], facecolor=c["ax_bg"])
        self.ax_f2 = self.fig.add_axes([x0 + 0.125, 0.550, 0.25 * w, 0.014], facecolor=c["ax_bg"])
        self.ax_f3 = self.fig.add_axes([x0 + 0.0125, 0.525, 0.25 * w, 0.014], facecolor=c["ax_bg"])
        self.ax_f4 = self.fig.add_axes([x0 + 0.125, 0.525, 0.25 * w, 0.014], facecolor=c["ax_bg"])

        self.sl_f1 = Slider(self.ax_f1, "f1 ", 0.0, 50.0, valinit=float(self._f1), valstep=1)
        self.sl_f2 = Slider(self.ax_f2, "f2 ", 0.5, 80.0, valinit=float(self._f2), valstep=1)
        self.sl_f3 = Slider(self.ax_f3, "f3 ", 20.0, 220.0, valinit=float(self._f3), valstep=1.0)
        self.sl_f4 = Slider(self.ax_f4, "f4 ", 40.0, 260.0, valinit=float(self._f4), valstep=1.0)
        self.sl_f1.on_changed(self._on_filter_sliders)
        self.sl_f2.on_changed(self._on_filter_sliders)
        self.sl_f3.on_changed(self._on_filter_sliders)
        self.sl_f4.on_changed(self._on_filter_sliders)

        self.ax_info = self.fig.add_axes([x0, 0.025, w, 0.30], facecolor=c["ax_bg"])
        self.ax_info.set_xticks([])
        self.ax_info.set_yticks([])
        for sp in self.ax_info.spines.values():
            sp.set_edgecolor(c["spine"])
        self.info_text = self.ax_info.text(
            0.02,
            0.98,
            "",
            va="top",
            ha="left",
            fontsize=8,
            color=c["text"],
            transform=self.ax_info.transAxes,
        )
        self.hover_text = self.fig.text(
            x0,
            0.925,
            "Hover help: move cursor over a button",
            va="center",
            ha="left",
            fontsize=8,
            color=c["label"],
        )

        self._refresh_mode_buttons()

    def _make_mode_buttons(self, x0: float, y: float, w: float, h: float,
                           labels: tuple, callback: Any,
                           tooltips: tuple | None = None) -> list:
        c = _tc()
        gap = 0.004
        n = max(1, len(labels))
        bw = (w - gap * (n - 1)) / n
        buttons: list = []
        for i, lbl in enumerate(labels):
            axb = self.fig.add_axes([x0 + i * (bw + gap), y, bw, h])
            btn = Button(
                axb,
                str(lbl),
                color=c["ax_bg"],
                hovercolor="#e8e8e8" if THEME == "light" else "#2a2d3d",
            )
            btn.label.set_color(c["text"])
            btn.label.set_fontsize(8)
            for sp in btn.ax.spines.values():
                sp.set_edgecolor(c["spine"])
            btn.on_clicked(lambda _e, _lbl=str(lbl): callback(_lbl))
            buttons.append(btn)
            if tooltips is not None and i < len(tooltips):
                self._register_hover(axb, str(tooltips[i]))
        return buttons

    def _register_hover(self, ax: Any, text: str):
        self._hover_help_map[ax] = str(text)

    def _on_hover_help(self, event: Any):
        msg = ""
        if event is not None and event.inaxes is not None:
            msg = self._hover_help_map.get(event.inaxes, "")
        if not msg:
            msg = "Hover help: move cursor over a button"
        if hasattr(self, "hover_text"):
            self.hover_text.set_text(msg)
            self.fig.canvas.draw_idle()

    def _style_mode_buttons(self, buttons: list, labels: tuple, active_label: str):
        c = _tc()
        for btn, lbl in zip(buttons, labels):
            is_active = str(lbl).lower() == str(active_label).lower()
            active_face = "#d9e8ff" if THEME == "light" else "#2a3b5c"
            inactive_face = c["ax_bg"]
            btn.ax.set_facecolor(active_face if is_active else inactive_face)
            btn.color = active_face if is_active else inactive_face
            btn.hovercolor = "#cfe2ff" if THEME == "light" else "#324972"
            btn.label.set_color(c["text"])

    def _refresh_mode_buttons(self):
        if hasattr(self, "_gain_btns"):
            self._style_mode_buttons(self._gain_btns, self._gain_labels, self._gain_mode)
        if hasattr(self, "_stat_btns"):
            self._style_mode_buttons(self._stat_btns, self._stat_labels, self._agc_stat)
        if hasattr(self, "_disp_btns"):
            self._style_mode_buttons(self._disp_btns, self._disp_labels, self._display_mode)
        if hasattr(self, "_filt_btns"):
            self._style_mode_buttons(self._filt_btns, self._filt_labels, self._filter_mode)
        if hasattr(self, "_butter_order_btns"):
            self._style_mode_buttons(self._butter_order_btns, self._butter_order_labels, str(int(self._butter_order)))
        if hasattr(self, "_pick_mode_btns"):
            self._style_mode_buttons(self._pick_mode_btns, self._pick_mode_labels, self._auto_pick_mode)
        if hasattr(self, "_pick_order_btns"):
            self._style_mode_buttons(self._pick_order_btns, self._pick_order_labels, self._pick_order)
        if hasattr(self, "_hilbert_target_btns"):
            self._style_mode_buttons(self._hilbert_target_btns, self._hilbert_target_labels, self._hilbert_target)

    def _refresh_toggle_buttons(self):
        c = _tc()
        if hasattr(self, "btn_tl"):
            active_face = "#d9e8ff" if THEME == "light" else "#2a3b5c"
            inactive_face = c["ax_bg"]
            self.btn_tl.ax.set_facecolor(active_face if self._show_timelines else inactive_face)
            self.btn_tl.color = active_face if self._show_timelines else inactive_face
            self.btn_tl.hovercolor = "#cfe2ff" if THEME == "light" else "#324972"
            self.btn_tl.label.set_color(c["text"])
        if hasattr(self, "btn_inv"):
            active_face = "#d9e8ff" if THEME == "light" else "#2a3b5c"
            inactive_face = c["ax_bg"]
            self.btn_inv.ax.set_facecolor(active_face if self._inverted else inactive_face)
            self.btn_inv.color = active_face if self._inverted else inactive_face
            self.btn_inv.hovercolor = "#cfe2ff" if THEME == "light" else "#324972"
            self.btn_inv.label.set_color(c["text"])

    def _update_invert_button_label(self):
        if hasattr(self, "btn_inv"):
            self.btn_inv.label.set_text("Normal Pol" if self._inverted else "Inverse Pol")
            self._refresh_toggle_buttons()

    def _on_gain_mode(self, label: str):
        self._gain_mode = str(label)
        self._refresh_mode_buttons()
        self._redraw()

    def _on_tmin(self, val: float):
        v = float(val)
        if v >= self._t_view_end - 1.0:
            v = self._t_view_end - 1.0
            self.sl_tmin.set_val(v)
            return
        self._t_view_start = v
        self._redraw()

    def _on_tmax(self, val: float):
        v = float(val)
        if v <= self._t_view_start + 1.0:
            v = self._t_view_start + 1.0
            self.sl_tmax.set_val(v)
            return
        self._t_view_end = v
        self._redraw()

    def _on_agc_stat(self, label: str):
        self._agc_stat = str(label)
        self._refresh_mode_buttons()
        self._redraw()

    def _on_agc_window(self, val: float):
        self._agc_window_ms = float(val)
        if self._gain_mode == "agc":
            self._redraw()

    def _on_wiggle_stretch(self, val: float):
        self._wiggle_stretch = float(val)
        self._redraw()

    def _on_display_mode(self, label: str):
        self._display_mode = str(label)
        self._refresh_mode_buttons()
        self._redraw()

    def _on_auto_pick_mode(self, label: str):
        self._auto_pick_mode = str(label).strip().lower()
        self._refresh_mode_buttons()
        self._redraw()

    def _on_pick_order(self, label: str):
        self._pick_order = str(label).strip().upper()
        self._refresh_mode_buttons()
        self._redraw()

    def _on_hilbert_target(self, label: str):
        self._hilbert_target = str(label).strip().lower()
        self._refresh_mode_buttons()
        self._redraw()

    def _on_filter_mode(self, label: str):
        self._filter_mode = str(label).lower()
        self._refresh_mode_buttons()
        self._redraw()

    def _on_butter_order(self, label: str):
        try:
            order = int(str(label).strip())
        except Exception:
            return
        if order < 1 or order == int(self._butter_order):
            return
        self._butter_order = int(order)
        self._schedule_filter_update()
        self._refresh_mode_buttons()

    def _toggle_invert(self):
        self._inverted = not self._inverted
        self._update_invert_button_label()
        self._redraw()

    def _toggle_timelines(self):
        self._show_timelines = not self._show_timelines
        self._refresh_toggle_buttons()
        self._redraw()

    def _schedule_filter_update(self):
        if self._filter_debounce_ms <= 0:
            self._apply_filter_update_now()
            return

        try:
            if self._filter_timer is not None:
                self._filter_timer.stop()
        except Exception:
            pass

        self._filter_timer = self.fig.canvas.new_timer(interval=int(self._filter_debounce_ms))
        self._filter_timer.single_shot = True
        self._filter_timer.add_callback(self._apply_filter_update_now)
        self._filter_timer.start()

    def _apply_filter_update_now(self):
        self._recompute_filter()
        self._redraw()

    def _apply_current_settings_to_picks(self):
        """Re-snap all current picks using active pick mode/order and settings."""
        if not self._picks:
            print("     Apply: no picks to update.")
            return
        updated = 0
        new_map: dict[int, float] = {}
        for idx, t_old in self._picks.items():
            t_new = round(float(self._snap_pick_time_ms(int(idx), float(t_old))), 2)
            if self._t_view_start <= t_new <= self._t_view_end:
                new_map[int(idx)] = t_new
            else:
                new_map[int(idx)] = float(t_old)
            if abs(float(new_map[int(idx)]) - float(t_old)) > 1e-6:
                updated += 1
        self._picks = new_map
        print(f"     Apply: updated {updated}/{len(self._picks)} picks using mode={self._auto_pick_mode}, order={self._pick_order}.")
        self._redraw()

    def _on_filter_sliders(self, _val: float):
        f1 = float(self.sl_f1.val)
        f2 = float(self.sl_f2.val)
        f3 = float(self.sl_f3.val)
        f4 = float(self.sl_f4.val)
        if not (f1 < f2 < f3 < f4):
            return
        self._f1, self._f2, self._f3, self._f4 = f1, f2, f3, f4
        self._schedule_filter_update()

    def _active_data(self) -> Any:
        if self._filter_mode == "ormsby":
            d = self.data_ormsby
        elif self._filter_mode == "butter":
            d = self.data_butter
        elif self._filter_mode == "cutpass":
            d = self.data_cutpass
        else:
            d = self.data_raw
        d = apply_gain(
            d,
            self.dt_s,
            mode=self._gain_mode,
            window_ms=self._agc_window_ms,
            stat=self._agc_stat,
        )
        return -d if self._inverted else d

    def _recompute_filter(self):
        self.data_ormsby = apply_ormsby_all_params(
            self.data_raw,
            self.dt_s,
            self._f1,
            self._f2,
            self._f3,
            self._f4,
        )
        self.data_butter = apply_butterworth_all_params(
            self.data_raw,
            self.dt_s,
            low_hz=self._f2,
            high_hz=self._f3,
            order=self._butter_order,
        )
        self.data_cutpass = apply_cutpass_all_params(
            self.data_raw,
            self.dt_s,
            low_cut_hz=self._f1,
            high_cut_hz=self._f4,
            order=self._butter_order,
        )

    def _apply_gain_single(self, trace: Any) -> Any:
        arr = np.asarray(trace, dtype=np.float32)[None, :]
        return apply_gain(
            arr,
            self.dt_s,
            mode=self._gain_mode,
            window_ms=self._agc_window_ms,
            stat=self._agc_stat,
        )[0]

    def _apply_filter_single(self, trace: Any) -> Any:
        x = np.asarray(trace, dtype=np.float32)
        if self._filter_mode == "ormsby":
            return ormsby(x, self.dt_s, f1=self._f1, f2=self._f2, f3=self._f3, f4=self._f4)
        if self._filter_mode == "butter":
            return butterworth_bandpass(x, self.dt_s, low_hz=self._f2, high_hz=self._f3, order=self._butter_order)
        if self._filter_mode == "cutpass":
            y = butterworth_highpass(x, self.dt_s, low_hz=self._f1, order=self._butter_order)
            y = butterworth_lowpass(y, self.dt_s, high_hz=self._f4, order=self._butter_order)
            return y
        return x

    def _prepare_trace_for_pick(self, trace: Any) -> Any:
        x = np.asarray(trace, dtype=np.float32)
        if self._pick_order == "G>F>X":
            x = self._apply_gain_single(x)
            x = self._apply_filter_single(x)
        elif self._pick_order == "RAW>X":
            x = x
        else:
            x = self._apply_filter_single(x)
            x = self._apply_gain_single(x)
        if self._inverted:
            x = -x
        return np.asarray(x, dtype=np.float32)

    def _timing_correction_ms(self, idx: int) -> float:
        """Return additional pick-time correction (ms): offset-dependent only."""
        off = abs(float(self.recv_abs[int(idx)]) - float(self.shot_pos_m)) if 0 <= int(idx) < len(self.recv_abs) else 0.0
        return float(TRIGGER_MS_PER_M) * off

    def _auto_abs_gate_ms_for_trace(self, idx: int) -> tuple[float, float]:
        """Return absolute-time gate [ms] for likely first arrivals on one trace."""
        t_abs_min = max(0.0, float(self._t_view_start))
        t_abs_max = min(float(self._t_view_end), float(self._t_full_end))
        if (not AUTO_USE_VELOCITY_GATE) or t_abs_max <= t_abs_min:
            return t_abs_min, t_abs_max

        off = abs(float(self.recv_abs[int(idx)]) - float(self.shot_pos_m)) if 0 <= int(idx) < len(self.recv_abs) else 0.0
        vmin = max(1.0, float(AUTO_VMIN_M_S))
        vmax = max(vmin + 1.0, float(AUTO_VMAX_M_S))
        pad = max(0.0, float(AUTO_GATE_PAD_MS))

        gate_lo = max(0.0, (off / vmax) * 1000.0 - pad)
        gate_hi = (off / vmin) * 1000.0 + pad

        lo = max(t_abs_min, gate_lo)
        hi = min(t_abs_max, gate_hi)
        if hi <= lo:
            return t_abs_min, t_abs_max
        return lo, hi

    def _snap_pick_time_ms(self, idx: int, t_hint_ms: float) -> float:
        mode = str(self._auto_pick_mode).strip().lower()
        if mode == "stalta":
            return float(t_hint_ms) + self._timing_correction_ms(idx)

        tr = self._prepare_trace_for_pick(self.data_raw[idx])
        t_rel_hint_s = (float(t_hint_ms) - float(self.delay_ms)) / 1000.0
        # For manual correction, prefer onset BEFORE the click to avoid snapping to wavelet tail.
        back_s = max(self.dt_s, 1.25 * float(MANUAL_SNAP_WIN_MS) / 1000.0)
        fwd_s = max(self.dt_s, 0.25 * float(MANUAL_SNAP_WIN_MS) / 1000.0)
        t0 = max(0.0, t_rel_hint_s - back_s)
        t1 = min((self.n_samp - 1) * self.dt_s, t_rel_hint_s + fwd_s)
        if t1 <= t0:
            return float(t_hint_ms) + self._timing_correction_ms(idx)

        if mode == "maxdiff_zero":
            ts = _zero_crossing_from_extremum_samples(
                samples=tr,
                dt_s=self.dt_s,
                win_start_s=t0,
                win_end_s=t1,
                search_direction=ZERO_X_SEARCH_DIRECTION,
                use_abs_peak=True,
            )
            if ts is None:
                return float(t_hint_ms) + self._timing_correction_ms(idx)
            return float(self.delay_ms + 1000.0 * ts + self._timing_correction_ms(idx))

        if mode == "hilbert_env":
            tr_obj = Trace(data=np.asarray(tr, dtype=np.float32))
            tr_obj.stats.delta = float(self.dt_s)
            _env, t_peak, t_onset = hilbert_envelope_pick(
                trace=tr_obj,
                win_start_s=t0,
                win_end_s=t1,
                onset_pct=HILBERT_ONSET_PCT,
            )
            if self._hilbert_target == "peak":
                ts = t_peak if t_peak is not None else t_onset
            else:
                ts = t_onset if t_onset is not None else t_peak
            if ts is None:
                return float(t_hint_ms) + self._timing_correction_ms(idx)
            return float(self.delay_ms + 1000.0 * float(ts) + self._timing_correction_ms(idx))

        return float(t_hint_ms) + self._timing_correction_ms(idx)

    def _draw_traces(self):
        c = _tc()
        data = self._active_data()

        dt_ms = self.dt_s * 1000.0
        i0 = max(0, int((self._t_view_start - self.delay_ms) / dt_ms))
        i1 = min(self.n_samp, int((self._t_view_end - self.delay_ms) / dt_ms) + 2)
        if i1 <= i0:
            i1 = min(self.n_samp, i0 + 2)
        ts = slice(i0, i1)
        t_ms = self.times_ms[ts]

        if self._display_mode in ("vd", "both"):
            d_img = data[:, ts]
            vmax = float(np.percentile(np.abs(d_img), 98)) if d_img.size else 1.0
            if vmax < 1e-12:
                vmax = 1.0
            self.ax.imshow(
                d_img.T,
                aspect="auto",
                cmap="seismic",
                vmin=-vmax,
                vmax=vmax,
                extent=[self.recv_abs.min(), self.recv_abs.max(), t_ms[-1], t_ms[0]],
                alpha=0.45,
                interpolation="nearest",
                zorder=1,
            )

        if self._display_mode == "vd":
            return

        stds = np.array([data[i, ts].astype(float).std() for i in range(self.n_traces)])
        valid = stds[stds > 1e-20]
        med = float(np.median(valid)) if len(valid) else 1.0
        norms = np.clip(stds, med * 0.3, med * 3.0)
        norms = np.where(norms > 1e-20, norms, med)
        scale = norms * CLIP_FACTOR

        for i, off in enumerate(self.recv_abs):
            tr = data[i, ts].astype(float)
            tr_n = np.clip(tr / scale[i], -1.0, 1.0)
            x_w = off + tr_n * self._dx * self._wiggle_stretch
            self.ax.plot(x_w, t_ms, color=c["trace"], lw=0.5, alpha=0.85)
            pos = np.where(tr_n > 0.0, tr_n, 0.0)
            self.ax.fill_betweenx(
                t_ms,
                off,
                off + pos * self._dx * self._wiggle_stretch,
                color=c["trace"],
                alpha=c["fill_alpha"],
            )

    def _draw_picks(self):
        c = _tc()
        for idx, t in self._picks.items():
            if idx < len(self.recv_abs):
                self.ax.plot(
                    self.recv_abs[idx],
                    t,
                    "v",
                    color=c["pick_clr"],
                    ms=7,
                    zorder=10,
                    markeredgecolor=c["pick_clr"],
                    markeredgewidth=0.6,
                )

    def _redraw(self):
        if getattr(self, "ax_top", None) is not None:
            try:
                self.ax_top.remove()
            except Exception:
                pass
            self.ax_top = None

        c = _tc()
        self.ax.clear()
        self.ax.set_facecolor(c["ax_bg"])
        self._draw_traces()
        self._draw_picks()

        self.ax.axvline(self.shot_pos_m, color="#ffcc00", lw=1.2, ls="--", alpha=0.70, zorder=3)

        if self._show_timelines:
            self.ax.axhline(0.0, color="#00aa66", lw=1.0, ls=":", alpha=0.75, zorder=3)
            for ms in range(10, int(self._t_view_end) + 1, 10):
                if self._t_view_start <= float(ms) <= self._t_view_end:
                    major = ms % 20 == 0
                    self.ax.axhline(
                        float(ms),
                        color="#666666" if THEME == "light" else "#aaaaaa",
                        lw=0.50 if major else 0.30,
                        ls=":" if major else "--",
                        alpha=0.35 if major else 0.20,
                        zorder=2,
                    )

        margin = self._dx * 0.5
        self.ax.set_xlim(self.recv_abs.min() - margin, self.recv_abs.max() + margin)
        self.ax.set_ylim(self._t_view_end, self._t_view_start)

        dt_ms = self.dt_s * 1000.0
        n_samp = self.n_samp
        t_end = self.delay_ms + (n_samp - 1) * dt_ms
        if self._filter_mode == "ormsby":
            filt_str = f"Ormsby {self._f1:.0f}-{self._f2:.0f}-{self._f3:.0f}-{self._f4:.0f} Hz"
        elif self._filter_mode == "butter":
            filt_str = f"Butterworth order {self._butter_order} ({self._f2:.0f}-{self._f3:.0f} Hz)"
        elif self._filter_mode == "cutpass":
            filt_str = f"Cut/Pass order {self._butter_order} (low f1={self._f1:.0f} Hz, high f4={self._f4:.0f} Hz)"
        else:
            filt_str = "Filter none"
        if self._inverted:
            filt_str += " [INV]"
        title = (
            f"Profile {self.profile} | Shot {self.shot_id} | {n_samp} smp "
            f"dt={dt_ms:.4f} ms delay={self.delay_ms:.1f} ms end={t_end:.1f} ms "
            f"| {len(self._picks)}/{self.n_traces} picks | {filt_str} | gain={self._gain_mode} | mode={self._auto_pick_mode} | order={self._pick_order}"
        )
        self.ax.set_title(title, color=c["text"], fontsize=9)
        self.ax.set_xlabel("Receiver position (m)", color=c["label"], fontsize=8)
        self.ax.set_ylabel("Time (ms)", color=c["label"])

        self.ax.text(
            0.0,
            -0.14,
            "L:pick R:delete Shift+L:range a:auto x:apply f:filter g:gain o:mode r:order h:hilbert v:polarity l:timeline c:top axis s/n:save+next p:prev q:quit",
            transform=self.ax.transAxes,
            ha="left",
            va="top",
            color=c["label"],
            fontsize=7,
            clip_on=False,
        )

        self.ax.grid(True, lw=0.3, alpha=0.30, color=c["grid"])
        self.ax.tick_params(colors=c["tick"])
        for sp in self.ax.spines.values():
            sp.set_edgecolor(c["spine"])

        handle = Line2D([0], [0], marker="v", linestyle="none", color=c["pick_clr"], markersize=7, label="First-break pick")
        self.ax.legend(
            handles=[handle],
            fontsize=8,
            facecolor=c["leg_face"],
            edgecolor=c["leg_edge"],
            labelcolor=c["text"],
            loc="lower right",
        )

        ax_top = self.ax.twiny()
        ax_top.set_xlim(self.ax.get_xlim())
        n_recv = len(self.recv_abs)
        step = max(1, n_recv // 8)
        idxs = list(range(0, n_recv, step))
        if (n_recv - 1) not in idxs:
            idxs.append(n_recv - 1)

        if self._top_mode == "channel":
            ax_top.set_xticks(self.recv_abs[idxs])
            ax_top.set_xticklabels([str(i + 1) for i in idxs], fontsize=7)
            ax_top.set_xlabel("Channel [c -> offset]", color=c["label"], fontsize=8)
        else:
            spm = self.shot_pos_m
            ax_top.set_xticks(self.recv_abs[idxs])
            ax_top.set_xticklabels([f"{self.recv_abs[i] - spm:+.0f}" for i in idxs], fontsize=7)
            ax_top.set_xlabel("Signed offset from shot (m) [c -> channel]", color=c["label"], fontsize=8)

        ax_top.tick_params(colors=c["tick"], labelsize=7)
        for sp in ax_top.spines.values():
            sp.set_edgecolor(c["spine"])

        self.ax_top = ax_top

        ffid = self._header_info.get("ffid", "?")
        ntr = self._header_info.get("n_tr", self.n_traces)
        nsamp = self._header_info.get("n_samp", self.n_samp)
        info_block = (
            "Header Info\n"
            f"FFID           : {ffid}\n"
            f"Traces         : {ntr}\n"
            f"Samples        : {nsamp}\n"
            f"dt (ms)        : {dt_ms:.4f}\n"
            f"Delay (ms)     : {self.delay_ms:.1f}\n"
            f"Shot (m)       : {self.shot_pos_m:.2f}\n"
            f"Gain           : {self._gain_mode} ({self._agc_stat}, {self._agc_window_ms:.0f} ms)\n"
            f"Display / Scale: {self._display_mode} / {self._wiggle_stretch:.2f}\n"
            f"Auto mode      : {self._auto_pick_mode}\n"
            f"Hilbert target  : {self._hilbert_target}\n"
            f"Pick order      : {self._pick_order}\n"
            f"Filter mode    : {self._filter_mode}\n"
            f"Ormsby (Hz)    : {self._f1:.1f}-{self._f2:.1f}-{self._f3:.1f}-{self._f4:.1f}\n"
            f"Butter (Hz)    : {self._f2:.1f}-{self._f3:.1f} (order {self._butter_order})\n"
            f"Cut/Pass (Hz)  : low=f1={self._f1:.1f}, high=f4={self._f4:.1f} (order {self._butter_order})\n"
            f"Time corr (ms) : delay includes +{TRIGGER_STATIC_MS:.2f}; residual +{TRIGGER_MS_PER_M:.5f}*offset(m)\n"
            f"Time view (ms) : {self._t_view_start:.1f} .. {self._t_view_end:.1f}\n"
            f"Polarity       : {'inverse' if self._inverted else 'normal'}\n"
            f"Timelines      : {'on' if self._show_timelines else 'off'}"
        )
        if hasattr(self, "info_text"):
            self.info_text.set_text(info_block)

        self._refresh_mode_buttons()
        self._refresh_toggle_buttons()
        self.fig.patch.set_facecolor(c["fig_bg"])
        self.fig.canvas.draw_idle()

    def _nearest_idx(self, x: float) -> int:
        return int(np.abs(self.recv_abs - x).argmin())

    def _on_click(self, event: Any):
        ax_top_ref = self.ax_top
        valid_axes = (self.ax,) if ax_top_ref is None else (self.ax, ax_top_ref)
        if event.inaxes not in valid_axes or event.x is None or event.y is None:
            return

        try:
            xdata, ydata = self.ax.transData.inverted().transform((event.x, event.y))
        except Exception:
            if event.xdata is None or event.ydata is None:
                return
            xdata, ydata = event.xdata, event.ydata

        idx = self._nearest_idx(xdata)
        if event.key == "shift" and event.button == 1:
            self._range_fill(idx, float(ydata))
            self._redraw()
            return

        if event.button == 1:
            t = round(float(self._snap_pick_time_ms(idx, float(ydata))), 2)
            if self._t_view_start <= t <= self._t_view_end:
                self._picks[idx] = t
                self._drag_pick = True
                self._last_drag_idx = idx
                self._last_drag_t = t
        elif event.button == 3:
            self._picks.pop(idx, None)
            self._drag_delete = True
            self._last_drag_idx = idx
            self._last_drag_t = None

        self._redraw()

    def _on_motion(self, event: Any):
        if not (self._drag_pick or self._drag_delete):
            return

        ax_top_ref = self.ax_top
        valid_axes = (self.ax,) if ax_top_ref is None else (self.ax, ax_top_ref)
        if event.inaxes not in valid_axes or event.x is None or event.y is None:
            return

        try:
            xdata, ydata = self.ax.transData.inverted().transform((event.x, event.y))
        except Exception:
            return

        idx = self._nearest_idx(xdata)

        if self._drag_pick:
            t = round(float(self._snap_pick_time_ms(idx, float(ydata))), 2)
            if self._last_drag_idx is None:
                if self._t_view_start <= t <= self._t_view_end:
                    self._picks[idx] = t
                self._last_drag_idx = idx
                self._last_drag_t = t
                self._redraw()
                return

            prev_idx = int(self._last_drag_idx)
            prev_t = float(self._last_drag_t if self._last_drag_t is not None else t)
            if idx == prev_idx:
                return

            step = 1 if idx > prev_idx else -1
            span = abs(idx - prev_idx)
            for n, ii in enumerate(range(prev_idx + step, idx + step, step), start=1):
                frac = n / float(span)
                ti = round(float(prev_t + frac * (t - prev_t)), 2)
                if self._t_view_start <= ti <= self._t_view_end:
                    self._picks[ii] = ti

            self._last_drag_idx = idx
            self._last_drag_t = t

        elif self._drag_delete:
            if self._last_drag_idx is None:
                self._picks.pop(idx, None)
                self._last_drag_idx = idx
                self._redraw()
                return

            prev_idx = int(self._last_drag_idx)
            if idx == prev_idx:
                return

            lo, hi = (prev_idx, idx) if prev_idx < idx else (idx, prev_idx)
            for ii in range(lo, hi + 1):
                self._picks.pop(ii, None)

            self._last_drag_idx = idx
            self._last_drag_t = None

        self._redraw()

    def _on_release(self, _event: Any):
        self._drag_pick = False
        self._drag_delete = False
        self._last_drag_idx = None
        self._last_drag_t = None

    def _range_fill(self, idx: int, ydata: float):
        t = round(float(ydata), 2)
        if not (self._t_view_start <= t <= self._t_view_end):
            return
        if self._range_anchor is None:
            self._range_anchor = (idx, t)
            print(f"     Range anchor set at trace {idx + 1}, t={t:.2f} ms")
            return

        i0, t0 = self._range_anchor
        i1, t1 = idx, t
        if i0 == i1:
            self._picks[i0] = t1
            self._range_anchor = None
            return

        lo, hi = (i0, i1) if i0 < i1 else (i1, i0)
        for ii in range(lo, hi + 1):
            frac = (ii - i0) / float(i1 - i0)
            tt = round(float(t0 + frac * (t1 - t0)), 2)
            if self._t_view_start <= tt <= self._t_view_end:
                self._picks[ii] = tt

        self._range_anchor = None
        print(f"     Range fill: traces {lo + 1}-{hi + 1}")

    def _save_and_finish(self, nav: str):
        if self._save_callback is not None:
            self._save_callback(self._picks)
        try:
            self._save_qc_image()
        except Exception:
            pass
        self._saved = True
        self._nav_action = nav
        self._done = True
        try:
            self.fig.canvas.stop_event_loop()
        except Exception:
            pass

    def _go_prev(self):
        self._save_and_finish("prev")

    def _auto_then_redraw(self):
        self._auto_pick()
        self._redraw()

    def _on_key(self, event: Any):
        global THEME
        key = event.key
        if key in ("q", "escape"):
            self._cancelled = True
            self._nav_action = "quit"
            self._done = True
            try:
                self.fig.canvas.stop_event_loop()
            except Exception:
                pass
        elif key == "s":
            self._save_and_finish("finalize")
        elif key == "n":
            self._save_and_finish("next")
        elif key == "p":
            self._go_prev()
        elif key == "k":
            self._save_qc_image()
        elif key == "f":
            self._filter_mode = {"none": "butter", "butter": "ormsby", "ormsby": "cutpass", "cutpass": "none"}.get(self._filter_mode, "none")
            self._refresh_mode_buttons()
            self._redraw()
        elif key == "u":
            self._f2 = max(self._f1 + 0.5, self._f2 + 1.0)
            self._schedule_filter_update()
        elif key == "j":
            self._f2 = max(self._f1 + 0.5, self._f2 - 1.0)
            self._schedule_filter_update()
        elif key == "i":
            self._f3 = min(self._f4 - 1.0, self._f3 + 5.0)
            self._schedule_filter_update()
        elif key == "m":
            self._f3 = max(self._f2 + 1.0, self._f3 - 5.0)
            self._schedule_filter_update()
        elif key == "g":
            self._gain_mode = {"none": "norm", "norm": "agc", "agc": "none"}.get(self._gain_mode, "none")
            self._redraw()
        elif key == "o":
            self._auto_pick_mode = {
                "stalta": "maxdiff_zero",
                "maxdiff_zero": "hilbert_env",
                "hilbert_env": "stalta",
            }.get(self._auto_pick_mode, "stalta")
            self._refresh_mode_buttons()
            self._redraw()
        elif key == "r":
            self._pick_order = {
                "F>G>X": "G>F>X",
                "G>F>X": "RAW>X",
                "RAW>X": "F>G>X",
            }.get(self._pick_order, "F>G>X")
            self._refresh_mode_buttons()
            self._redraw()
        elif key == "h":
            self._hilbert_target = {"onset": "peak", "peak": "onset"}.get(self._hilbert_target, "onset")
            self._refresh_mode_buttons()
            self._redraw()
        elif key == "x":
            self._apply_current_settings_to_picks()
        elif key == "v":
            self._toggle_invert()
        elif key == "l":
            self._toggle_timelines()
        elif key == "t":
            THEME = "light" if THEME == "dark" else "dark"
            self._redraw()
        elif key == "a":
            self._auto_pick()
            self._redraw()
        elif key == "c":
            self._top_mode = "offset" if self._top_mode == "channel" else "channel"
            self._redraw()

    def _auto_pick(self):
        mode = str(self._auto_pick_mode).strip().lower()
        if mode == "maxdiff_zero":
            self._auto_pick_maxdiff_zero()
            return
        if mode == "hilbert_env":
            self._auto_pick_hilbert_env()
            return
        self._auto_pick_stalta()

    def _auto_pick_stalta(self):
        data = self._active_data()
        count = 0
        n_sta = max(1, int(STA_MS / 1000.0 / self.dt_s))
        n_lta = max(3, int(LTA_MS / 1000.0 / self.dt_s))
        start_idx = max(0, int((0.0 - self.delay_ms) / (self.dt_s * 1000.0)))
        end_idx = min(self.n_samp - n_sta - 1, int((self._t_view_end - self.delay_ms) / (self.dt_s * 1000.0)))

        for i in range(self.n_traces):
            tr = data[i].astype(float)
            tr = tr - np.mean(tr[: max(1, n_lta)])
            char = np.maximum(tr, 0.0)
            search_lo = max(n_lta, start_idx)
            search_hi = max(search_lo + 2, end_idx)
            local_char = char[search_lo:search_hi]
            noise_floor = float(np.median(local_char)) if local_char.size else 0.0
            amp_floor = max(1e-12, 4.0 * noise_floor)
            picked_idx = None

            for k in range(search_lo, search_hi - 1):
                lta = char[k - n_lta : k].mean()
                if lta < 1e-30:
                    continue
                sta = char[k : k + n_sta].mean()
                ratio = sta / lta
                next_lta = char[k + 1 - n_lta : k + 1].mean()
                if next_lta < 1e-30:
                    continue
                next_sta = char[k + 1 : k + 1 + n_sta].mean()
                next_ratio = next_sta / next_lta
                if ratio >= STALTA_TRIG and next_ratio >= STALTA_TRIG and sta >= amp_floor:
                    picked_idx = k
                    break

            if picked_idx is None:
                continue

            t_abs = round(self.delay_ms + picked_idx * self.dt_s * 1000.0, 2)
            t_abs = float(t_abs) + self._timing_correction_ms(i)
            if self._t_view_start <= t_abs <= self._t_view_end:
                self._picks[i] = t_abs
                count += 1

        print(f"     Auto-pick (stalta): {count}/{self.n_traces} placed")

    def _auto_pick_maxdiff_zero(self):
        count = 0
        for i in range(self.n_traces):
            tr = self._prepare_trace_for_pick(self.data_raw[i])
            abs_lo_ms, abs_hi_ms = self._auto_abs_gate_ms_for_trace(i)
            win_start_s = max(0.0, (abs_lo_ms - float(self.delay_ms)) / 1000.0)
            win_end_s = max(win_start_s + self.dt_s, (abs_hi_ms - float(self.delay_ms)) / 1000.0)
            t_s = _zero_crossing_from_extremum_samples(
                samples=tr,
                dt_s=self.dt_s,
                win_start_s=win_start_s,
                win_end_s=win_end_s,
                search_direction=ZERO_X_SEARCH_DIRECTION,
                use_abs_peak=True,
            )
            if t_s is None:
                continue
            t_abs = round(self.delay_ms + (t_s * 1000.0) + self._timing_correction_ms(i), 2)
            if self._t_view_start <= t_abs <= self._t_view_end:
                self._picks[i] = t_abs
                count += 1
        print(f"     Auto-pick (maxdiff_zero): {count}/{self.n_traces} placed")

    def _auto_pick_hilbert_env(self):
        count = 0
        for i in range(self.n_traces):
            abs_lo_ms, abs_hi_ms = self._auto_abs_gate_ms_for_trace(i)
            win_start_s = max(0.0, (abs_lo_ms - float(self.delay_ms)) / 1000.0)
            win_end_s = max(win_start_s + self.dt_s, (abs_hi_ms - float(self.delay_ms)) / 1000.0)
            tr = Trace(data=np.asarray(self._prepare_trace_for_pick(self.data_raw[i]), dtype=np.float32))
            tr.stats.delta = float(self.dt_s)
            _, peak_time_s, onset_time_s = hilbert_envelope_pick(
                trace=tr,
                win_start_s=win_start_s,
                win_end_s=win_end_s,
                onset_pct=HILBERT_ONSET_PCT,
            )
            if self._hilbert_target == "peak":
                t_s = peak_time_s if peak_time_s is not None else onset_time_s
            else:
                t_s = onset_time_s if onset_time_s is not None else peak_time_s
            if t_s is None:
                continue
            t_abs = round(self.delay_ms + (float(t_s) * 1000.0) + self._timing_correction_ms(i), 2)
            if self._t_view_start <= t_abs <= self._t_view_end:
                self._picks[i] = t_abs
                count += 1
        print(f"     Auto-pick (hilbert_env): {count}/{self.n_traces} placed")

    def _save_qc_image(self):
        self.qc_dir.mkdir(parents=True, exist_ok=True)
        out = self.qc_dir / f"{self.profile}_shot{self.shot_id:02d}_picks.png"
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        bbox = self.ax.get_tightbbox(renderer).expanded(1.02, 1.04)
        bbox_inches = bbox.transformed(self.fig.dpi_scale_trans.inverted())
        self.fig.savefig(str(out), dpi=180, bbox_inches=bbox_inches, facecolor=self.ax.get_facecolor())
        print(f"\n     QC image -> {out.name}")

    def run(self) -> Any:
        self._done = False

        def _on_close(_evt: Any):
            if not (self._saved or self._cancelled):
                self._cancelled = True
            self._done = True
            try:
                self.fig.canvas.stop_event_loop()
            except Exception:
                pass

        self.fig.canvas.mpl_connect("close_event", _on_close)
        plt.show(block=False)

        while not self._done:
            try:
                self.fig.canvas.start_event_loop(0.05)
            except Exception:
                self._cancelled = True
                break

        try:
            if plt.fignum_exists(self.fig.number):
                plt.close(self.fig)
        except Exception:
            pass

        if self._cancelled:
            return {"status": "quit", "picks": self._picks}
        return {"status": self._nav_action or "next", "picks": self._picks}


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def _picks_json_path(profile_name: str) -> Path:
    return OUTPUT_DIR / profile_name / "picks.json"


def _session_picks_json_path(profile_name: str) -> Path:
    return OUTPUT_DIR / profile_name / "picks.session.json"


def load_picks_json(profile_name: str) -> dict:
    p = _picks_json_path(profile_name)
    if not p.exists():
        return {}
    try:
        with open(p, encoding="utf-8") as fh:
            raw = json.load(fh)
        return {int(k): {int(ti): float(tv) for ti, tv in v.items()} for k, v in raw.items()}
    except Exception as exc:
        print(f"  [WARN] Could not load picks.json: {exc}")
        return {}


def load_session_picks_json(profile_name: str) -> dict:
    p = _session_picks_json_path(profile_name)
    if not p.exists():
        return {}
    try:
        with open(p, encoding="utf-8") as fh:
            raw = json.load(fh)
        return {int(k): {int(ti): float(tv) for ti, tv in v.items()} for k, v in raw.items()}
    except Exception as exc:
        print(f"  [WARN] Could not load picks.session.json: {exc}")
        return {}


def save_picks_json(profile_name: str, all_picks: dict):
    p = _picks_json_path(profile_name)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump({str(k): {str(ti): tv for ti, tv in v.items()} for k, v in all_picks.items()}, fh, indent=2)
    print(f"     Picks saved -> {p.relative_to(CWD.parent)}")


def save_session_picks_json(profile_name: str, all_picks: dict):
    p = _session_picks_json_path(profile_name)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump({str(k): {str(ti): tv for ti, tv in v.items()} for k, v in all_picks.items()}, fh, indent=2)
    print(f"     Session picks saved -> {p.relative_to(CWD.parent)}")


def clear_session_picks_json(profile_name: str):
    p = _session_picks_json_path(profile_name)
    try:
        if p.exists():
            p.unlink()
            print(f"     Session file cleared -> {p.relative_to(CWD.parent)}")
    except Exception as exc:
        print(f"  [WARN] Could not clear session file: {exc}")


# ---------------------------------------------------------------------------
# Minimal picking workflow (no analysis/export)
# ---------------------------------------------------------------------------

def process_profile(profile_name: str, geom_override: int | None = None):
    cfg = PROFILES.get(profile_name)
    has_profile_cfg = cfg is not None
    if cfg is None:
        cfg = {"geom": 100, "line_no": profile_name, "perp_m": 0.0, "shots": "auto"}
        print(f"  [INFO] Profile '{profile_name}' not in PROFILES; using dynamic defaults.")

    data_dir = DATA_DIR / profile_name
    if not data_dir.exists():
        print(f"[ERROR] Data folder not found: {data_dir}")
        return

    if not _ensure_interactive_backend():
        return

    geom_type: int
    geom_src: str
    if geom_override is not None:
        geom_type = int(geom_override)
        geom_src = "CLI"
    elif has_profile_cfg and cfg.get("geom") in (100, 200):
        geom_type = int(cfg.get("geom"))
        geom_src = "profile config"
    else:
        geom_xls, xls_path, xls_dist = infer_geometry_from_lvl_excel(DATA_DIR, profile_name)
        if geom_xls in (100, 200):
            geom_type = int(geom_xls)
            if xls_path is not None and xls_dist is not None:
                geom_src = f"excel ({xls_path.name}, dist~{xls_dist:.2f} m)"
            else:
                geom_src = "excel"
        else:
            geom_type = 100
            geom_src = "default"

    recv_positions = load_geometry(geom_type)

    shots_cfg = cfg.get("shots", "auto")
    if shots_cfg == "auto" or not isinstance(shots_cfg, dict):
        shots_cfg = auto_shot_positions(recv_positions)

    print(
        f"  Geometry {geom_type} m : {len(recv_positions)} receivers, "
        f"{recv_positions[0]:.2f} - {recv_positions[-1]:.2f} m ({geom_src})"
    )
    print(
        "  Shot positions: "
        + "  ".join(f"Shot{k}={v:.3f}m" for k, v in shots_cfg.items())
    )

    seg2_candidates = list(data_dir.glob("*.seg2")) + list(data_dir.glob("*.SEG2"))
    seg2_unique = {str(p.resolve()).lower(): p for p in seg2_candidates}

    def _ffid(path: Path) -> int:
        digits = "".join(c for c in path.stem if c.isdigit())
        return int(digits) if digits else 0

    seg2_files = sorted(seg2_unique.values(), key=_ffid)
    if not seg2_files:
        print(f"[ERROR] No .seg2 files found in {data_dir}")
        return

    print(f"  Found {len(seg2_files)} SEG2 file(s)")

    final_picks: dict = load_picks_json(profile_name)
    session_picks: dict = load_session_picks_json(profile_name)
    all_picks: dict = session_picks if session_picks else {sid: dict(vals) for sid, vals in final_picks.items()}
    if session_picks:
        print("  Session picks found: resuming from picks.session.json")

    shot_cache: list = []
    for file_idx, seg2_path in enumerate(seg2_files):
        shot_id = file_idx + 1
        data_raw, dt_s, n_tr, n_samp, shot_pos_hdr, ffid_hdr, delay_ms, _ = read_seg2(seg2_path)

        n_show = min(n_tr, len(recv_positions))
        shot_pos_cfg = shots_cfg.get(shot_id)
        shot_pos_source = "default"
        shot_pos_m = 0.0

        if shot_pos_cfg is not None:
            shot_pos_m = float(shot_pos_cfg)
            shot_pos_source = "config"
        elif shot_pos_hdr is not None:
            shot_pos_m = float(shot_pos_hdr)
            shot_pos_source = "header"

        if USE_SEG2_SHOT_POSITION and shot_pos_hdr is not None:
            shot_pos_m = float(shot_pos_hdr)
            shot_pos_source = "header"

        shot_cache.append(
            {
                "shot_id": shot_id,
                "seg2_path": seg2_path,
                "data_raw": data_raw,
                "dt_s": dt_s,
                "n_tr": n_tr,
                "n_samp": n_samp,
                "shot_pos_hdr": shot_pos_hdr,
                "ffid_hdr": ffid_hdr,
                "delay_ms": delay_ms,
                "n_show": n_show,
                "shot_pos_m": shot_pos_m,
                "shot_pos_source": shot_pos_source,
            }
        )

    idx = 0
    finalized = False
    qc_dir = OUTPUT_DIR / profile_name

    while 0 <= idx < len(shot_cache):
        shot = shot_cache[idx]
        shot_id = int(shot["shot_id"])
        seg2_path = shot["seg2_path"]

        data_raw = shot["data_raw"]
        dt_s = float(shot["dt_s"])
        n_tr = int(shot["n_tr"])
        n_samp = int(shot["n_samp"])
        ffid_hdr = shot["ffid_hdr"]
        delay_ms = float(shot["delay_ms"])
        n_show = int(shot["n_show"])
        shot_pos_m = float(shot["shot_pos_m"])
        shot_pos_source = shot.get("shot_pos_source", "default")

        dt_ms = dt_s * 1000.0
        t_end = delay_ms + (n_samp - 1) * dt_ms

        print(f"\n  -- Shot {shot_id} ({seg2_path.name}) --")
        print(
            f"     FFID={ffid_hdr} | {n_tr} traces | dt={dt_ms:.4f} ms | "
            f"delay={delay_ms:.1f} ms | {n_samp} smp ({delay_ms:.1f} to {t_end:.1f} ms)"
        )
        print(f"     Shot pos : {shot_pos_m:.2f} m ({shot_pos_source.upper()})")

        geom_slice = recv_positions[:n_show]
        data_slice = data_raw[:n_show]
        data_filt = apply_ormsby_all(data_slice, dt_s)

        def _save_cb(picks_for_shot: dict, _sid: int = shot_id) -> None:
            all_picks[_sid] = picks_for_shot
            save_session_picks_json(profile_name, all_picks)
            print(f"     Saved {len(picks_for_shot)} pick(s) to session.")

        picker = FirstBreakPicker(
            data_slice,
            data_filt,
            dt_s,
            geom_slice,
            shot_id,
            profile_name,
            shot_pos_m=shot_pos_m,
            delay_ms=delay_ms,
            existing_picks=all_picks.get(shot_id, {}),
            qc_dir=qc_dir,
            save_callback=_save_cb,
            header_info={"ffid": ffid_hdr, "n_tr": n_tr, "n_samp": n_samp, "shot_pos_hdr": shot.get("shot_pos_hdr")},
        )

        result = picker.run() or {"status": "quit", "picks": all_picks.get(shot_id, {})}
        status = result.get("status", "next")
        all_picks[shot_id] = dict(result.get("picks", {}))
        save_session_picks_json(profile_name, all_picks)

        if status == "prev":
            idx = max(0, idx - 1)
            continue
        if status == "next":
            idx = min(len(shot_cache) - 1, idx + 1)
            if idx == len(shot_cache) - 1 and shot_id == len(shot_cache):
                print("     Already at last shot; press s to finalize or p to review.")
            continue
        if status == "finalize":
            finalized = True
            break
        if status == "quit":
            break

    total_picks = sum(len(v) for v in all_picks.values())
    print(f"\n  Total picks in memory: {total_picks}")

    if finalized:
        save_picks_json(profile_name, all_picks)
        clear_session_picks_json(profile_name)
        print(f"  Finalized picks -> {(OUTPUT_DIR / profile_name).relative_to(CWD.parent)}")
    else:
        print("  Session left open. Continue later to finalize picks.json.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="fb_picker.py",
        description="Clean first-break picker only (no refraction analysis).",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python fb_picker.py 150\n"
            "  python fb_picker.py 150 200\n"
            "  python fb_picker.py --all\n"
        ),
    )
    parser.add_argument("profile", nargs="?", default=None, help="Profile folder name, e.g. 120")
    parser.add_argument("geometry", nargs="?", default=None, help="Optional geometry override: 100 or 200")
    parser.add_argument("--all", action="store_true", help="Process all profile folders detected under data/")
    args = parser.parse_args()

    geom_override: int | None = None
    if args.geometry is not None:
        g = str(args.geometry).strip().lower().replace(" ", "")
        if g in ("100", "geometry100"):
            geom_override = 100
        elif g in ("200", "geometry200"):
            geom_override = 200
        else:
            print("[ERROR] Invalid geometry override. Use 100 or 200.")
            return

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.all:
        targets = discover_profile_folders(DATA_DIR)
        if not targets:
            print("[ERROR] No profile folders with SEG2 files found under data/.")
            return
    elif args.profile:
        targets = [args.profile]
    else:
        chosen_profile = None
        try:
            root = tk.Tk()
            root.withdraw()
            sel = filedialog.askdirectory(title="Select profile data folder", initialdir=str(DATA_DIR), mustexist=True)
            root.destroy()
            if sel:
                p = Path(sel)
                if p.parent.resolve() == DATA_DIR.resolve():
                    chosen_profile = p.name
        except Exception:
            chosen_profile = None

        if chosen_profile:
            targets = [chosen_profile]
        else:
            print("\nAvailable profiles:")
            print(f"  {'Name':<12}  {'Geom':>6}  Data folder")
            print(f"  {'-' * 12}  {'-' * 6}  {'-' * 30}")
            dynamic_profiles = discover_profile_folders(DATA_DIR)
            if not dynamic_profiles:
                print("  (none found)")
            for pname in dynamic_profiles:
                pcfg = PROFILES.get(pname, {"geom": 100})
                folder = DATA_DIR / pname
                n_seg2 = len(list(folder.glob("*.seg2"))) + len(list(folder.glob("*.SEG2")))
                print(f"  {pname:<12}  {int(pcfg.get('geom', 100)):>5}m  found ({n_seg2} SEG2 files)")
            print("\nUsage: python fb_picker.py <profile> [geometry]")
            return

    for pname in targets:
        print(f"\n{'=' * 60}\nProfile : {pname}\n{'=' * 60}")
        process_profile(pname, geom_override=geom_override)

    print("\nAll done.")


if __name__ == "__main__":
    main()
