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

from obspy import read as _read_obspy  # type: ignore[import-untyped]
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.widgets import Button, Slider


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

# Display setup
THEME: str = "light"
T_MAX_MS: float = 150.0
T_DISPLAY_PRE_MS: float = -10.0
CLIP_FACTOR: float = 2.0
FILTER_ON: bool = True
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
        delay_ms = float(delay_s) * 1000.0

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
    a: STA/LTA auto-pick
    f: filter on/off
    g: cycle gain (none -> norm -> agc)
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
        self._filter_on = FILTER_ON
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

        self.ax_btn_apply = self.fig.add_axes([x0 + 0.145, 0.505, 0.04, 0.022])

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
        self.btn_apply.on_clicked(lambda _e: self._apply_filter_controls())

        self.fig.text(x0, 0.80, "Gain:", fontsize=8, fontweight="bold", color=c["text"])
        self.fig.text(x0, 0.76, "AGC stat:", fontsize=8, fontweight="bold", color=c["text"])
        self.fig.text(x0, 0.72, "Display:", fontsize=8, fontweight="bold", color=c["text"])
        self.fig.text(x0, 0.602, "Filter (Hz)", fontsize=8, fontweight="bold", color=c["text"])

        self._gain_labels = ("none", "norm", "agc")
        self._gain_btns = self._make_mode_buttons(
            x0=x0 + 0.055, y=0.792, w=0.145, h=0.026,
            labels=self._gain_labels, callback=self._on_gain_mode,
        )

        self._stat_labels = ("rms", "mean")
        self._stat_btns = self._make_mode_buttons(
            x0=x0 + 0.055, y=0.752, w=0.095, h=0.026,
            labels=self._stat_labels, callback=self._on_agc_stat,
        )

        self._disp_labels = ("wiggle", "vd", "both")
        self._disp_btns = self._make_mode_buttons(
            x0=x0 + 0.055, y=0.712, w=0.145, h=0.026,
            labels=self._disp_labels, callback=self._on_display_mode,
        )

        self.ax_agc = self.fig.add_axes([x0 + 0.025, 0.67, 0.5 * w, 0.016], facecolor=c["ax_bg"])
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

        self.ax_wig = self.fig.add_axes([x0 + 0.025, 0.64, 0.5 * w, 0.016], facecolor=c["ax_bg"])
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

        self.fig.text(x0, 0.472, "Time View (ms)", fontsize=8, fontweight="bold", color=c["text"], ha="left", va="bottom")
        tmin_hi = max(self._t_full_start + 1.0, self._t_full_end - 1.0)
        tmax_lo = min(self._t_full_end - 1.0, self._t_full_start + 1.0)
        self.ax_tmin = self.fig.add_axes([x0 + 0.025, 0.440, 0.1, 0.014], facecolor=c["ax_bg"])
        self.ax_tmax = self.fig.add_axes([x0 + 0.025, 0.415, 0.1, 0.014], facecolor=c["ax_bg"])
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

        self.ax_f1 = self.fig.add_axes([x0 + 0.0125, 0.580, 0.1, 0.014], facecolor=c["ax_bg"])
        self.ax_f2 = self.fig.add_axes([x0 + 0.0125, 0.555, 0.1, 0.014], facecolor=c["ax_bg"])
        self.ax_f3 = self.fig.add_axes([x0 + 0.0125, 0.530, 0.1, 0.014], facecolor=c["ax_bg"])
        self.ax_f4 = self.fig.add_axes([x0 + 0.0125, 0.505, 0.1, 0.014], facecolor=c["ax_bg"])

        self.sl_f1 = Slider(self.ax_f1, "f1", 0.0, 50.0, valinit=float(self._f1), valstep=0.5)
        self.sl_f2 = Slider(self.ax_f2, "f2", 0.5, 80.0, valinit=float(self._f2), valstep=0.5)
        self.sl_f3 = Slider(self.ax_f3, "f3", 20.0, 220.0, valinit=float(self._f3), valstep=1.0)
        self.sl_f4 = Slider(self.ax_f4, "f4", 40.0, 260.0, valinit=float(self._f4), valstep=1.0)

        self.ax_info = self.fig.add_axes([x0, 0.125, w, 0.265], facecolor=c["ax_bg"])
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

        self._refresh_mode_buttons()

    def _make_mode_buttons(self, x0: float, y: float, w: float, h: float, labels: tuple, callback: Any) -> list:
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
        return buttons

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

    def _toggle_invert(self):
        self._inverted = not self._inverted
        self._update_invert_button_label()
        self._redraw()

    def _toggle_timelines(self):
        self._show_timelines = not self._show_timelines
        self._refresh_toggle_buttons()
        self._redraw()

    def _apply_filter_controls(self):
        f1 = float(self.sl_f1.val)
        f2 = float(self.sl_f2.val)
        f3 = float(self.sl_f3.val)
        f4 = float(self.sl_f4.val)
        if not (f1 < f2 < f3 < f4):
            print("     [WARN] Need f1 < f2 < f3 < f4")
            return
        self._f1, self._f2, self._f3, self._f4 = f1, f2, f3, f4
        self._recompute_filter()
        self._redraw()

    def _active_data(self) -> Any:
        d = self.data_filt if self._filter_on else self.data_raw
        d = apply_gain(
            d,
            self.dt_s,
            mode=self._gain_mode,
            window_ms=self._agc_window_ms,
            stat=self._agc_stat,
        )
        return -d if self._inverted else d

    def _recompute_filter(self):
        self.data_filt = apply_ormsby_all_params(self.data_raw, self.dt_s, self._f1, self._f2, self._f3, self._f4)

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
        filt_str = (
            f"Ormsby {self._f1:.0f}-{self._f2:.0f}-{self._f3:.0f}-{self._f4:.0f} Hz "
            + ("ON" if self._filter_on else "OFF")
            + (" [INV]" if self._inverted else "")
        )
        title = (
            f"Profile {self.profile} | Shot {self.shot_id} | {n_samp} smp "
            f"dt={dt_ms:.4f} ms delay={self.delay_ms:.1f} ms end={t_end:.1f} ms "
            f"| {len(self._picks)}/{self.n_traces} picks | {filt_str} | gain={self._gain_mode}"
        )
        self.ax.set_title(title, color=c["text"], fontsize=9)
        self.ax.set_xlabel("Receiver position (m)", color=c["label"], fontsize=8)
        self.ax.set_ylabel("Time (ms)", color=c["label"])

        self.ax.text(
            0.0,
            -0.14,
            "L:pick R:delete Shift+L:range a:auto f:filter g:gain v:polarity l:timeline c:top axis s/n:save+next p:prev q:quit",
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
            f"Filter (Hz)    : {self._f1:.1f}-{self._f2:.1f}-{self._f3:.1f}-{self._f4:.1f}\n"
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
            t = round(float(ydata), 2)
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
            t = round(float(ydata), 2)
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
            self._filter_on = not self._filter_on
            self._redraw()
        elif key == "u":
            self._f2 = max(self._f1 + 0.5, self._f2 + 1.0)
            self._recompute_filter()
            self._redraw()
        elif key == "j":
            self._f2 = max(self._f1 + 0.5, self._f2 - 1.0)
            self._recompute_filter()
            self._redraw()
        elif key == "i":
            self._f3 = min(self._f4 - 1.0, self._f3 + 5.0)
            self._recompute_filter()
            self._redraw()
        elif key == "m":
            self._f3 = max(self._f2 + 1.0, self._f3 - 5.0)
            self._recompute_filter()
            self._redraw()
        elif key == "g":
            self._gain_mode = {"none": "norm", "norm": "agc", "agc": "none"}.get(self._gain_mode, "none")
            self._redraw()
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
            if self._t_view_start <= t_abs <= self._t_view_end:
                self._picks[i] = t_abs
                count += 1

        print(f"     Auto-pick: {count}/{self.n_traces} placed")

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
    if cfg is None:
        print(f"[ERROR] Profile '{profile_name}' not in PROFILES. Known: {list(PROFILES)}")
        return

    data_dir = DATA_DIR / profile_name
    if not data_dir.exists():
        print(f"[ERROR] Data folder not found: {data_dir}")
        return

    if not _ensure_interactive_backend():
        return

    geom_type = int(geom_override) if geom_override is not None else int(cfg["geom"])
    recv_positions = load_geometry(geom_type)

    shots_cfg = cfg.get("shots", "auto")
    if shots_cfg == "auto" or not isinstance(shots_cfg, dict):
        shots_cfg = auto_shot_positions(recv_positions)

    print(
        f"  Geometry {geom_type} m : {len(recv_positions)} receivers, "
        f"{recv_positions[0]:.2f} - {recv_positions[-1]:.2f} m"
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
    parser.add_argument("--all", action="store_true", help="Process all profiles in PROFILES")
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
        targets = list(PROFILES)
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
            for pname, pcfg in PROFILES.items():
                folder = DATA_DIR / pname
                status = "found" if folder.exists() else "MISSING"
                n_seg2 = len(list(folder.glob("*.seg2"))) if folder.exists() else 0
                print(f"  {pname:<12}  {pcfg['geom']:>5}m  {status} ({n_seg2} SEG2 files)")
            print("\nUsage: python fb_picker.py <profile> [geometry]")
            return

    for pname in targets:
        print(f"\n{'=' * 60}\nProfile : {pname}\n{'=' * 60}")
        process_profile(pname, geom_override=geom_override)

    print("\nAll done.")


if __name__ == "__main__":
    main()
