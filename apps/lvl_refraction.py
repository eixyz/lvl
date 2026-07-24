"""
lvl_refraction.py  --  First-break picking + LVL refraction analysis
====================================================================
Project : Low Velocity Layer (LVL) refraction seismic surveys

Workflow (mirrors ProMAX flow)
------------------------------
1.  Read SEG2 files -- one file per shot, sorted by FFID (file number)
2.  Assign receiver geometry from standard geometry files
    (geometry100.txt = 100 m spread,  geometry200.txt = 200 m spread)
3.  Apply bulk time static  (like ProMAX: bulk shift static -> add)
4.  Apply zero-phase Ormsby bandpass filter
    (frequency domain, 25 % zero-padding for FFT, 4 corner frequencies)
5.  Interactive first-break picking  (matplotlib GUI)
6.  Automatic refraction analysis:
    -- AIC-based breakpoint detection (or manual)
    -- T-X segment fitting -> velocities
    -- Perpendicular offset correction  (true offset = sqrt(inline^2 + perp^2))
    -- Intercept-time depth inversion (2- or 3-layer)
7.  Export:
    -- <profile>_picks.xlsx     Config + per-shot formulas + Analysis sheet
    -- <profile>_picks_clean.txt  compatible with lvl.ipynb
    -- <profile>_tx_picks.png     T-X summary with fitted lines, time downward
    -- <profile>_shot_XX_picks.png  per-shot QC image (saved on shot navigation/finalize)

Interactive picker controls
---------------------------
    Left-click / drag       place picks quickly (interpolates across skipped traces)
    Right-click / drag      delete picks quickly (range delete while dragging)
    Shift + Left-click      range-fill picks between two traces
    n                       save to session and go to next shot
    p                       save to session and go to previous shot
    s                       Save/Close (finalize picks + exports)
    a                       STA/LTA auto-pick (gain-aware)
    f                       toggle Ormsby filter on / off
    v                       invert trace polarity (display only)
    g                       cycle gain mode: norm/agc/none
    l                       toggle timeline guides
    c                       toggle top axis: channel <-> signed offset
    k                       save per-shot QC screenshot (PNG)
    q / Esc                 quit picking

Directory layout
----------------
  data/
    geometry100.txt
    geometry200.txt
    120/       Rec_00001.seg2  Rec_00002.seg2  Rec_00003.seg2
    150/  ...
  scripts/
    lvl_refraction.py   <- this file
  output/
    120/       <auto-created>

How to run
----------
  conda activate seiseng
  cd d:\\Daten\\seismic\\lvl\\scripts

  python lvl_refraction.py 120              # interactive picking for profile 120
  python lvl_refraction.py 120 --export-only  # re-export without picking
  python lvl_refraction.py --all            # pick all profiles in sequence
  python lvl_refraction.py                  # list available profiles
"""

from __future__ import annotations

import sys
from pathlib import Path

# CWD = Path(__file__).resolve().parent
# PROJECT_DIR = CWD.parent
# sys.path.append(str(PROJECT_DIR))

PROJECT_DIR = Path(__file__).resolve().parents[1]

if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))
# Importing directories
from src.common.paths import *
from src.common import paths as core_paths
from src.io import project_io as pio
from src.common.settings import APP_NAME, APP_VERSION
from src.picker.picker import FirstBreakPicker as FirstBreakPickerV2
from src.picker.settings import PickerSettings as PickerSettingsV2

import argparse
import json
import math
import datetime
import importlib
import subprocess 
import tkinter as tk
from tkinter import filedialog

from typing import Any
from zoneinfo import ZoneInfo
from src.communication.control_bridge import read_latest_picker_command

# ---------------------------------------------------------------------------
# Auto-install missing packages into the active environment
# ---------------------------------------------------------------------------
def _ensure(pkg: str, mod: str = "") -> Any:
    m = mod or pkg
    try:
        return importlib.import_module(m)
    except ImportError:
        print(f"  Installing '{pkg}' ...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "-q"])
        return importlib.import_module(m)

np = _ensure("numpy")
pd = _ensure("pandas")
_ensure("obspy")
_ensure("matplotlib")
_ensure("openpyxl")
_ensure("scipy")
_ensure("xlrd")

from obspy import Trace, read as _read_obspy          # type: ignore[import-untyped]
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.widgets import Button, Slider, TextBox
from scipy.signal import butter as _scipy_butter, hilbert as _scipy_hilbert, sosfiltfilt as _scipy_sosfiltfilt
import openpyxl
from openpyxl.styles import PatternFill, Font, Alignment
from openpyxl.utils import get_column_letter
import pandas as pd


# ===========================================================================
# CONFIG  --  edit this section for each project
# ===========================================================================

# ===========================================================================
# CONFIG  --  single source of truth in src/common/settings.py + paths.py
# ===========================================================================

from src.common.settings import (
    GEOM_FILES, PROFILES, USE_SEG2_SHOT_POSITION,
    BP_F1, BP_F2, BP_F3, BP_F4, BP_FFT_PAD, BP_REAPPLY, BUTTER_ORDER,
    FILTER_DEBOUNCE_MS, BULK_SHIFT_MS,
    THEME, T_MAX_MS, T_DISPLAY_PRE_MS, CLIP_FACTOR, FILTER_ON,
    DEFAULT_FILTER_MODE, DISPLAY_MODE, SHOW_TIMELINES, WIGGLE_STRETCH,
    GAIN_MODE, AGC_WINDOW_MS, AGC_STAT,
    STA_MS, LTA_MS, STALTA_TRIG, AUTO_PICK_MODE, ZERO_X_SEARCH_DIRECTION,
    HILBERT_ONSET_PCT, HILBERT_TARGET, PICK_PROCESS_ORDER, MANUAL_SNAP_WIN_MS,
    AUTO_USE_VELOCITY_GATE, AUTO_VMIN_M_S, AUTO_VMAX_M_S, AUTO_GATE_PAD_MS,
    TRIGGER_STATIC_MS, TRIGGER_MS_PER_M,
    N_LAYERS, AUTO_BREAKPOINTS, MANUAL_BREAKS,
    theme_colors, _ensure_interactive_backend,
)

# Processing engine (moved out of this file; see src/picker, src/refraction, src/io).
from src.picker.preprocessing import (
    ormsby, apply_ormsby_all, apply_ormsby_all_params,
    butterworth_bandpass, apply_butterworth_all_params,
    butterworth_highpass, butterworth_lowpass, apply_cutpass_all_params,
    apply_bulk_static, apply_gain,
)
from src.picker.features import hilbert_envelope_pick
from src.picker.refinement import _zero_crossing_from_extremum_samples
from src.io.seg2_reader import read_seg2, read_seg2_acquisition_time_de, read_seg2_mid_xyz
from src.io.pick_reader import (
    load_layer_json, load_layer_session_json, save_layer_json,
    save_layer_session_json, clear_layer_session_json,
    load_picks_json, load_session_picks_json, save_picks_json,
    save_session_picks_json, clear_session_picks_json,
)
from src.utils.geometry import (
    load_geometry, auto_shot_positions, discover_field_report_excels,
    discover_profile_folders, infer_geometry_from_spread_length,
    infer_geometry_from_seg2_file, load_midpoint_xyz_from_geometry_excels,
    load_profile_geometry_from_excels, infer_geometry_from_field_report,
    load_profile_offsets_from_excel, read_perpendicular_offsets,
    _normalize_profile_token,
)
from src.utils.math_utils import true_offset
from src.refraction.travel_times import fit_line, depth_2layer, depth_3layer
from src.refraction.velocity_model import (
    compute_layer_averages, _predict_time_from_fit, _choose_shot_fit,
)
from src.refraction.layer_analysis import (
    build_corrected_pick_data, _fit_layers_from_windows, _compute_fit_rms,
    _recompute_fit_derived, _drop_fit_layer, build_analysis_from_layers,
)
from src.refraction.processing import parse_shot_value_map, resolve_perp_by_shot
from src.io.exporters import (
    export_velocity_summary_excel, export_excel, export_picks_txt,
    export_tx_plot, export_fit_plot, export_corrected_qc_plot,
    export_arrivals_observed_computed_plot, export_layer_fit_rms_plot,
)

# ===========================================================================
# END OF CONFIG
# ===========================================================================


# ---------------------------------------------------------------------------
# Theme colour palette
# ---------------------------------------------------------------------------
def _tc() -> dict:
    """Return a colour dict for the current THEME (toggled live via the 't' key)."""
    return theme_colors(THEME)


def prompt_offset_model_by_shot(perp_by_shot: dict,
                                inline_shift_by_shot: dict | None = None,
                                shots_info: list | None = None,
                                all_picks: dict | None = None,
                                recv_positions: Any | None = None) -> tuple:
    """
    Interactive UI to edit per-shot geometric corrections.

    - PO (m): perpendicular source-to-line distance.
    - X-shift (m): optional inline distance shift for far-offset acquisition geometry.

    A live preview panel redraws all shots together when values change.
    """
    out_po = {int(k): float(v) for k, v in (perp_by_shot or {}).items()}
    out_shift = {int(k): float(v) for k, v in (inline_shift_by_shot or {}).items()}
    for sid in out_po:
        out_shift.setdefault(int(sid), 0.0)

    if not out_po:
        return out_po, out_shift

    try:
        backend = str(plt.get_backend()).lower()
    except Exception:
        backend = ""
    if "agg" in backend and not _ensure_interactive_backend():
        print("  [INFO] Using default PO values (interactive PO window unavailable).")
        return out_po, out_shift

    c = _tc()
    shot_ids = sorted(out_po)
    n = len(shot_ids)
    fig_h = max(5.6, min(11.0, 3.3 + 0.55 * n))
    fig = plt.figure(figsize=(11.2, fig_h))
    fig.patch.set_facecolor(c["fig_bg"])
    fig.text(0.04, 0.95, "Shot Geometry Setup (PO + X-shift)",
             color=c["text"], fontsize=12, fontweight="bold", va="top")
    fig.text(0.04, 0.90,
             "Adjust PO and optional inline shift. Press Enter in a field for live update.",
             color=c["label"], fontsize=9, va="top")

    top = 0.82
    row_h = 0.055
    textboxes_po: dict = {}
    textboxes_shift: dict = {}
    for i, sid in enumerate(shot_ids):
        y = top - i * row_h
        fig.text(0.04, y + 0.017, f"Shot {sid}", color=c["text"], fontsize=9,
                 ha="left", va="center")
        ax_tb_po = fig.add_axes([0.11, y, 0.12, 0.038])
        tb_po = TextBox(ax_tb_po, "PO", initial=f"{out_po[sid]:.3f}")
        tb_po.label.set_color(c["text"])
        tb_po.label.set_fontsize(8)
        tb_po.text_disp.set_color(c["text"])
        tb_po.text_disp.set_fontsize(9)
        ax_tb_po.set_facecolor(c["ax_bg"])
        for sp in ax_tb_po.spines.values():
            sp.set_edgecolor(c["spine"])

        ax_tb_shift = fig.add_axes([0.25, y, 0.12, 0.038])
        tb_shift = TextBox(ax_tb_shift, "X-shift", initial=f"{out_shift.get(sid, 0.0):.3f}")
        tb_shift.label.set_color(c["text"])
        tb_shift.label.set_fontsize(8)
        tb_shift.text_disp.set_color(c["text"])
        tb_shift.text_disp.set_fontsize(9)
        ax_tb_shift.set_facecolor(c["ax_bg"])
        for sp in ax_tb_shift.spines.values():
            sp.set_edgecolor(c["spine"])

        textboxes_po[sid] = tb_po
        textboxes_shift[sid] = tb_shift

    ax_prev = fig.add_axes([0.41, 0.16, 0.55, 0.70])
    ax_prev.set_facecolor(c["ax_bg"])
    for sp in ax_prev.spines.values():
        sp.set_edgecolor(c["spine"])
    ax_prev.grid(True, lw=0.3, alpha=0.30, color=c["grid"])
    ax_prev.tick_params(colors=c["tick"])
    ax_prev.set_xlabel("Signed inline x (m, +right / -left)", color=c["label"], fontsize=8)
    ax_prev.set_ylabel("FB time (ms)", color=c["label"], fontsize=8)
    ax_prev.set_title("Live preview: all shots", color=c["text"], fontsize=9)

    def _read_float(tb: Any, fallback: float) -> float:
        try:
            txt = str(tb.text).strip()
        except Exception:
            txt = ""
        if not txt:
            return float(fallback)
        try:
            return float(txt.replace(",", "."))
        except Exception:
            return float(fallback)

    def _read_ui_values() -> tuple:
        po_new = dict(out_po)
        sh_new = dict(out_shift)
        for sid in shot_ids:
            po_new[sid] = _read_float(textboxes_po[sid], po_new.get(sid, 0.0))
            sh_new[sid] = _read_float(textboxes_shift[sid], sh_new.get(sid, 0.0))
        return po_new, sh_new

    def _redraw_preview(po_map: dict, sh_map: dict):
        ax_prev.clear()
        ax_prev.set_facecolor(c["ax_bg"])
        ax_prev.grid(True, lw=0.3, alpha=0.30, color=c["grid"])
        ax_prev.tick_params(colors=c["tick"])
        ax_prev.set_xlabel("Signed inline x (m, +right / -left)", color=c["label"], fontsize=8)
        ax_prev.set_ylabel("FB time (ms)", color=c["label"], fontsize=8)
        ax_prev.set_title("Live preview: all shots", color=c["text"], fontsize=9)

        if not (shots_info and all_picks and recv_positions is not None):
            ax_prev.text(0.5, 0.5, "Preview unavailable (missing picks/geometry)",
                         transform=ax_prev.transAxes, ha="center", va="center",
                         color=c["label"], fontsize=8)
            fig.canvas.draw_idle()
            return

        pal = ["#e63946", "#2a9d8f", "#e9c46a", "#457b9d", "#f4a261", "#8d99ae"]
        plotted = 0
        for i, (sid, shot_pos) in enumerate(shots_info):
            sid_i = int(sid)
            picks = apply_bulk_static((all_picks or {}).get(sid_i, {}))
            if not picks:
                continue
            po = float(po_map.get(sid_i, 0.0))
            x_shift = float(sh_map.get(sid_i, 0.0))
            xvals: list = []
            tvals: list = []
            for tr in sorted(picks):
                if tr >= len(recv_positions):
                    continue
                inline_abs = abs(float(recv_positions[tr]) - float(shot_pos))
                inline_signed = float(recv_positions[tr]) - float(shot_pos)
                inline_corr = max(0.0, inline_abs + x_shift)
                x_signed = float(np.sign(inline_signed) * inline_corr)
                xvals.append(x_signed)
                tvals.append(float(picks[tr]))
            if len(xvals) < 2:
                continue
            col = pal[i % len(pal)]
            order = np.argsort(np.asarray(xvals, dtype=float))
            xs = np.asarray(xvals, dtype=float)[order]
            ts = np.asarray(tvals, dtype=float)[order]
            ax_prev.plot(xs, ts, "o-", ms=3.5, lw=1.0, alpha=0.9, color=col,
                         label=f"Shot {sid_i}  PO={po:.2f}m  dX={x_shift:.2f}m")
            plotted += 1

        if plotted > 0:
            ax_prev.legend(fontsize=7, facecolor=c["leg_face"], edgecolor=c["leg_edge"],
                           labelcolor=c["text"], loc="best")
            ax_prev.set_ylim(T_MAX_MS, 0.0)
        else:
            ax_prev.text(0.5, 0.5, "No picks available for preview", transform=ax_prev.transAxes,
                         ha="center", va="center", color=c["label"], fontsize=8)
        fig.canvas.draw_idle()

    def _apply_preview(_evt: Any = None):
        po_map, sh_map = _read_ui_values()
        _redraw_preview(po_map, sh_map)

    _apply_preview()

    status_ax = fig.add_axes([0.04, 0.14, 0.33, 0.045])
    status_ax.set_facecolor(c["ax_bg"])
    status_ax.set_xticks([])
    status_ax.set_yticks([])
    for sp in status_ax.spines.values():
        sp.set_edgecolor(c["spine"])
    status_txt = status_ax.text(0.02, 0.5, "Ready", color=c["text"],
                                fontsize=8, va="center", ha="left",
                                transform=status_ax.transAxes)

    # Keep controls below input fields so preview plot stays visually clean.
    ax_ok = fig.add_axes([0.04, 0.08, 0.105, 0.045])
    ax_ref = fig.add_axes([0.155, 0.08, 0.105, 0.045])
    ax_def = fig.add_axes([0.27, 0.08, 0.105, 0.045])
    btn_ok = Button(ax_ok, "Use Values", color=c["ax_bg"],
                    hovercolor="#e8e8e8" if THEME == "light" else "#2a2d3d")
    btn_ref = Button(ax_ref, "Refresh", color=c["ax_bg"],
                     hovercolor="#e8e8e8" if THEME == "light" else "#2a2d3d")
    btn_def = Button(ax_def, "Use Defaults", color=c["ax_bg"],
                     hovercolor="#e8e8e8" if THEME == "light" else "#2a2d3d")
    for btn in (btn_ok, btn_ref, btn_def):
        btn.label.set_color(c["text"])
        btn.label.set_fontsize(8)
        for sp in btn.ax.spines.values():
            sp.set_edgecolor(c["spine"])

    state = {"done": False, "use_defaults": False}

    def _finish(use_defaults: bool):
        state["use_defaults"] = use_defaults
        state["done"] = True
        try:
            fig.canvas.stop_event_loop()
        except Exception:
            pass

    btn_ok.on_clicked(lambda _e: _finish(False))
    btn_ref.on_clicked(_apply_preview)
    btn_def.on_clicked(lambda _e: _finish(True))

    for sid in shot_ids:
        textboxes_po[sid].on_submit(lambda _txt, _sid=sid: _apply_preview())
        textboxes_shift[sid].on_submit(lambda _txt, _sid=sid: _apply_preview())
        try:
            textboxes_po[sid].on_text_change(lambda _txt, _sid=sid: _apply_preview())
            textboxes_shift[sid].on_text_change(lambda _txt, _sid=sid: _apply_preview())
        except Exception:
            pass

    def _on_key(evt: Any):
        if evt.key == "enter":
            _finish(False)
        elif evt.key in ("escape", "q"):
            _finish(True)

    def _on_close(_evt: Any):
        _finish(True)

    fig.canvas.mpl_connect("key_press_event", _on_key)
    fig.canvas.mpl_connect("close_event", _on_close)
    plt.show(block=False)
    while not state["done"] and plt.fignum_exists(fig.number):
        try:
            fig.canvas.start_event_loop(0.05)
        except Exception:
            break

    if state["use_defaults"]:
        try:
            if plt.fignum_exists(fig.number):
                plt.close(fig)
        except Exception:
            pass
        return out_po, out_shift

    for sid in shot_ids:
        cur_po = out_po[sid]
        cur_shift = out_shift.get(sid, 0.0)
        out_po[sid] = _read_float(textboxes_po[sid], cur_po)
        out_shift[sid] = _read_float(textboxes_shift[sid], cur_shift)
        if not np.isfinite(out_po[sid]):
            status_txt.set_text(f"Invalid PO for Shot {sid}; keeping {cur_po:.3f}")
            out_po[sid] = cur_po
        if not np.isfinite(out_shift[sid]):
            status_txt.set_text(f"Invalid X-shift for Shot {sid}; keeping {cur_shift:.3f}")
            out_shift[sid] = cur_shift

    try:
        if plt.fignum_exists(fig.number):
            plt.close(fig)
    except Exception:
        pass
    return out_po, out_shift




# AnalysisWorkflow moved to src/refraction/pipeline.py (shared with lvl_studio.py).
from src.refraction.pipeline import AnalysisWorkflow

def _pick_layer_windows_from_plot(x_vals: Any, t_vals: Any,
                                  title: str) -> dict:
    """
        Let user pick up to 6 x-boundaries on a T-X plot:
            x1-x2 -> layer 1
            x3-x4 -> layer 2
            x5-x6 -> layer 3
        Dedicated UI controls:
            LMB = add point, MMB = delete nearest point, RMB-hold = preview trendline/apparent velocity,
            Undo = remove last, Save = continue, Skip = ignore side
    """
    try:
        backend = str(plt.get_backend()).lower()
    except Exception:
        backend = ""
    if "agg" in backend and not _ensure_interactive_backend():
        return {"windows": [None, None, None], "picked_points": []}

    c = _tc()
    fig, ax = plt.subplots(figsize=(8.6, 5.4), constrained_layout=False)
    fig.subplots_adjust(left=0.08, right=0.74, bottom=0.12, top=0.90)
    fig.patch.set_facecolor(c["fig_bg"])
    ax.set_facecolor(c["ax_bg"])
    ax.plot(x_vals, t_vals, "o-", ms=4, lw=1.2, color="#1f77b4")
    ax.set_xlabel("Corrected true offset XO (m)", color=c["label"])
    ax.set_ylabel("Corrected travel-time (ms)", color=c["label"])
    ax.set_title(title + "\nPick x1..x6 (LMB), delete near point (MMB), hold RMB for trend preview.",
                 color=c["text"], fontsize=10)
    ax.grid(True, lw=0.3, alpha=0.3, color=c["grid"])
    ax.tick_params(colors=c["tick"])
    for sp in ax.spines.values():
        sp.set_edgecolor(c["spine"])

    pick_names = ["x1", "x2", "x3", "x4", "x5", "x6"]
    state = {
        "xs": [],
        "ts": [],
        "done": False,
        "skip": False,
        "preview_active": False,
        "preview_anchor_idx": None,
    }
    click_lines: list = []
    click_labels: list = []

    preview_line, = ax.plot([], [], "--", lw=1.2, color="#ff7f0e", alpha=0.9, zorder=6, visible=False)
    preview_txt = ax.text(
        0.98,
        0.02,
        "",
        transform=ax.transAxes,
        va="bottom",
        ha="right",
        fontsize=8,
        color=c["text"],
        bbox=dict(boxstyle="round,pad=0.30", facecolor=c["ax_bg"], edgecolor=c["spine"], alpha=0.90),
        visible=False,
    )

    info_ax = fig.add_axes([0.76, 0.62, 0.22, 0.28])
    info_ax.set_facecolor(c["ax_bg"])
    info_ax.set_xticks([])
    info_ax.set_yticks([])
    for sp in info_ax.spines.values():
        sp.set_edgecolor(c["spine"])
    info_txt = info_ax.text(
        0.05, 0.95,
        "Layer windows\n"
        "x1-x2: layer1\n"
        "x3-x4: layer2\n"
        "x5-x6: layer3\n\n"
        "Selected points:\nnone",
        transform=info_ax.transAxes,
        color=c["text"], fontsize=8, va="top", ha="left",
    )

    legend_ax = fig.add_axes([0.76, 0.12, 0.22, 0.16])
    legend_ax.set_facecolor(c["ax_bg"])
    legend_ax.set_xticks([])
    legend_ax.set_yticks([])
    for sp in legend_ax.spines.values():
        sp.set_edgecolor(c["spine"])
    legend_ax.text(
        0.05,
        0.95,
        "Mouse controls\n"
        "LMB: add point\n"
        "MMB: delete nearest\n"
        "RMB hold: trend preview\n"
        "RMB release: clear preview",
        transform=legend_ax.transAxes,
        color=c["text"],
        fontsize=8,
        va="top",
        ha="left",
    )

    ax_undo = fig.add_axes([0.76, 0.50, 0.22, 0.07])
    ax_save = fig.add_axes([0.76, 0.41, 0.22, 0.07])
    ax_skip = fig.add_axes([0.76, 0.32, 0.22, 0.07])

    hover = "#e8e8e8" if THEME == "light" else "#2a2d3d"
    btn_undo = Button(ax_undo, "Undo last", color=c["ax_bg"], hovercolor=hover)
    btn_save = Button(ax_save, "Save side", color=c["ax_bg"], hovercolor=hover)
    btn_skip = Button(ax_skip, "Skip side", color=c["ax_bg"], hovercolor=hover)
    for btn in (btn_undo, btn_save, btn_skip):
        btn.label.set_color(c["text"])
        btn.label.set_fontsize(8)
        for sp in btn.ax.spines.values():
            sp.set_edgecolor(c["spine"])

    def _clear_preview() -> None:
        state["preview_active"] = False
        state["preview_anchor_idx"] = None
        preview_line.set_visible(False)
        preview_txt.set_visible(False)

    def _nearest_point_idx(xv: float, tv: float) -> int | None:
        if not state["xs"]:
            return None
        xr = max(1e-9, float(np.max(x_vals) - np.min(x_vals)))
        tr = max(1e-9, float(np.max(t_vals) - np.min(t_vals)))
        best_idx = None
        best_d2 = float("inf")
        for i, (px, pt) in enumerate(zip(state["xs"], state["ts"])):
            dx = (float(px) - float(xv)) / xr
            dt = (float(pt) - float(tv)) / tr
            d2 = dx * dx + dt * dt
            if d2 < best_d2:
                best_d2 = d2
                best_idx = i
        if best_idx is None:
            return None
        return int(best_idx) if best_d2 <= (0.06 ** 2) else None

    def _redraw_selected_markers() -> None:
        while click_lines:
            ln = click_lines.pop()
            try:
                ln.remove()
            except Exception:
                pass
        while click_labels:
            tx = click_labels.pop()
            try:
                tx.remove()
            except Exception:
                pass

        for idx, (xv, tv) in enumerate(zip(state["xs"], state["ts"])):
            nm = pick_names[idx] if idx < len(pick_names) else f"x{idx + 1}"
            ln = ax.axvline(x=float(xv), color="#e63946", lw=1.0, ls="--", alpha=0.85)
            tx = ax.text(float(xv), float(tv), f"{nm}\n{float(tv):.1f}ms", color="#e63946", fontsize=8,
                         rotation=90, va="bottom", ha="center")
            click_lines.append(ln)
            click_labels.append(tx)

    def _update_preview(xc: float, tc: float) -> None:
        if not state["preview_active"]:
            return
        idx = state.get("preview_anchor_idx", None)
        if idx is None or idx < 0 or idx >= len(state["xs"]):
            _clear_preview()
            return
        xa = float(state["xs"][idx])
        ta = float(state["ts"][idx])
        xb = float(xc)
        tb = float(tc)
        preview_line.set_data([xa, xb], [ta, tb])
        preview_line.set_visible(True)
        dt_ms = float(tb - ta)
        dx_m = float(xb - xa)
        if abs(dt_ms) > 1e-9:
            v_app = 1000.0 * abs(dx_m) / abs(dt_ms)
            v_txt = f"{v_app:.1f} m/s"
        else:
            v_txt = "inf"
        preview_txt.set_text(
            f"Preview\n"
            f"from x={xa:.2f}, t={ta:.2f}\n"
            f"to   x={xb:.2f}, t={tb:.2f}\n"
            f"Vapp ~ {v_txt}"
        )
        preview_txt.set_visible(True)

    def _update_info():
        if not state["xs"]:
            sel = "none"
        else:
            parts = []
            for i, xv in enumerate(state["xs"]):
                nm = pick_names[i] if i < len(pick_names) else f"x{i}"
                tv = state["ts"][i] if i < len(state["ts"]) else float("nan")
                parts.append(f"{nm}={xv:.2f}, t{i+1}={tv:.2f}")
            sel = "\n".join(parts)
        info_txt.set_text(
            "Layer windows\n"
            "x1-x2: layer1\n"
            "x3-x4: layer2\n"
            "x5-x6: layer3\n\n"
            f"Selected points:\n{sel}"
        )

    def _finish(skip: bool):
        if state["done"]:
            return
        state["skip"] = skip
        state["done"] = True
        try:
            fig.canvas.stop_event_loop()
        except Exception:
            pass

    def _undo(_evt: Any):
        if not state["xs"]:
            return
        state["xs"].pop()
        if state["ts"]:
            state["ts"].pop()
        _redraw_selected_markers()
        _clear_preview()
        _update_info()
        fig.canvas.draw_idle()

    def _save(_evt: Any):
        _finish(False)

    def _skip(_evt: Any):
        _finish(True)

    btn_undo.on_clicked(_undo)
    btn_save.on_clicked(_save)
    btn_skip.on_clicked(_skip)

    def _on_click(evt: Any):
        if evt.inaxes is not ax or evt.xdata is None or evt.ydata is None:
            return
        if evt.button == 2:
            idx_near = _nearest_point_idx(float(evt.xdata), float(evt.ydata))
            if idx_near is None:
                return
            state["xs"].pop(idx_near)
            state["ts"].pop(idx_near)
            _redraw_selected_markers()
            _clear_preview()
            _update_info()
            fig.canvas.draw_idle()
            return
        if evt.button == 3:
            if not state["xs"]:
                return
            idx_near = _nearest_point_idx(float(evt.xdata), float(evt.ydata))
            state["preview_anchor_idx"] = idx_near if idx_near is not None else (len(state["xs"]) - 1)
            state["preview_active"] = True
            _update_preview(float(evt.xdata), float(evt.ydata))
            fig.canvas.draw_idle()
            return
        if evt.button != 1:
            return
        if len(state["xs"]) >= 6:
            return
        xv = float(evt.xdata)
        tv = float(evt.ydata)
        state["xs"].append(xv)
        state["ts"].append(tv)
        _redraw_selected_markers()
        _clear_preview()
        _update_info()
        fig.canvas.draw_idle()

    def _on_motion(evt: Any):
        if not state["preview_active"]:
            return
        if evt.inaxes is not ax or evt.xdata is None or evt.ydata is None:
            return
        _update_preview(float(evt.xdata), float(evt.ydata))
        fig.canvas.draw_idle()

    def _on_release(evt: Any):
        if evt.button == 3 and state["preview_active"]:
            _clear_preview()
            fig.canvas.draw_idle()

    def _on_leave(_evt: Any):
        if state["preview_active"]:
            _clear_preview()
            fig.canvas.draw_idle()

    def _on_key(evt: Any):
        if evt.key == "enter":
            _finish(False)
        elif evt.key in ("escape", "q"):
            _finish(True)
        elif evt.key in ("backspace", "delete"):
            _undo(evt)

    def _on_close(_evt: Any):
        _finish(True)

    fig.canvas.mpl_connect("button_press_event", _on_click)
    fig.canvas.mpl_connect("motion_notify_event", _on_motion)
    fig.canvas.mpl_connect("button_release_event", _on_release)
    fig.canvas.mpl_connect("axes_leave_event", _on_leave)
    fig.canvas.mpl_connect("key_press_event", _on_key)
    fig.canvas.mpl_connect("close_event", _on_close)
    _update_info()
    plt.show(block=False)
    while not state["done"] and plt.fignum_exists(fig.number):
        try:
            fig.canvas.start_event_loop(0.05)
        except Exception:
            break

    xs = state["xs"]
    try:
        if plt.fignum_exists(fig.number):
            plt.close(fig)
    except Exception:
        pass

    if state["skip"] or not xs:
        return {"windows": [None, None, None], "picked_points": []}

    windows: list = []
    for i in range(0, 6, 2):
        if i + 1 < len(xs):
            a, b = xs[i], xs[i + 1]
            windows.append((min(a, b), max(a, b)))
        else:
            windows.append(None)
    picked_points = [(float(xv), float(tv)) for xv, tv in zip(state.get("xs", []), state.get("ts", []))]
    return {"windows": windows[:3], "picked_points": picked_points}










def _review_layer_fit_interactive(x_vals: Any, t_vals: Any,
                                  fit_res: dict, title: str) -> str:
    """
    Review layer fit and choose action.

    Returns one of: "accept", "repick", "skip", "drop1", "drop2", "drop3".
    """
    try:
        backend = str(plt.get_backend()).lower()
    except Exception:
        backend = ""
    if "agg" in backend and not _ensure_interactive_backend():
        return "accept"

    c = _tc()
    x = np.asarray(x_vals, dtype=float)
    t = np.asarray(t_vals, dtype=float)
    if x.size < 2 or t.size < 2:
        return "accept"

    rms = _compute_fit_rms(x, t, fit_res)
    nfit = int(sum(_predict_time_from_fit(float(xx), fit_res) is not None for xx in x))

    fig, (ax_tx, ax_cmp) = plt.subplots(1, 2, figsize=(12.2, 5.0), constrained_layout=False)
    fig.subplots_adjust(left=0.06, right=0.82, bottom=0.12, top=0.88, wspace=0.22)
    fig.patch.set_facecolor(c["fig_bg"])
    ax_tx.set_facecolor(c["ax_bg"])
    ax_cmp.set_facecolor(c["ax_bg"])

    ax_tx.plot(x, t, "o", ms=4, color="#1f77b4", alpha=0.9, label="Observed")
    wins = list((fit_res or {}).get("windows", []) or [])
    segs = list((fit_res or {}).get("segments", []) or [])
    vel_lines: list = []
    for i, seg in enumerate(segs):
        sl = float(seg.get("slope_ms_m", 0.0))
        ic = float(seg.get("intercept_ms", 0.0))
        vel_val = float(seg.get("velocity_m_s", 0.0) or 0.0)
        vel_lines.append(f"L{i+1}: {vel_val:.0f} m/s")
        if sl <= 0.0:
            continue
        w = wins[i] if i < len(wins) else None
        if w is not None:
            x0, x1 = float(w[0]), float(w[1])
        else:
            x0, x1 = float(seg.get("x0", 0.0)), float(seg.get("x1", 0.0))
        if x1 <= x0:
            continue
        xl = np.linspace(x0, x1, 40)
        tl = sl * xl + ic
        ax_tx.plot(xl, tl, "-", lw=1.6, alpha=0.9, color="#e63946")

    if vel_lines:
        ax_tx.text(
            0.98, 0.02,
            "Apparent velocity\n" + "\n".join(vel_lines),
            transform=ax_tx.transAxes,
            va="bottom", ha="right",
            fontsize=8,
            color=c["text"],
            bbox=dict(boxstyle="round,pad=0.35",
                      facecolor=c["ax_bg"], edgecolor=c["spine"], alpha=0.90),
        )

    pred_all: list = []
    obs_all: list = []
    for xv, tv in zip(x, t):
        tp = _predict_time_from_fit(float(xv), fit_res)
        if tp is None:
            continue
        obs_all.append(float(tv))
        pred_all.append(float(tp))

    if obs_all:
        oa = np.asarray(obs_all, dtype=float)
        pa = np.asarray(pred_all, dtype=float)
        ax_cmp.scatter(oa, pa, s=24, color="#2a9d8f", alpha=0.9)
        lo = float(min(oa.min(), pa.min()))
        hi = float(max(oa.max(), pa.max()))
        pad = max(1.0, 0.03 * (hi - lo))
        lo -= pad
        hi += pad
        ax_cmp.plot([lo, hi], [lo, hi], "--", lw=1.1, color="#555555")
        ax_cmp.set_xlim(lo, hi)
        ax_cmp.set_ylim(lo, hi)
    else:
        ax_cmp.text(0.5, 0.5, "No predicted points", transform=ax_cmp.transAxes,
                    ha="center", va="center", color=c["label"], fontsize=8)

    ax_tx.set_xlabel("Corrected true offset XO (m)", color=c["label"], fontsize=8)
    ax_tx.set_ylabel("Corrected travel-time (ms)", color=c["label"], fontsize=8)
    ax_tx.set_title("Observed with fitted segments", color=c["text"], fontsize=9)
    ax_tx.grid(True, lw=0.3, alpha=0.3, color=c["grid"])
    ax_tx.tick_params(colors=c["tick"], labelsize=8)

    ax_cmp.set_xlabel("Observed (ms)", color=c["label"], fontsize=8)
    ax_cmp.set_ylabel("Computed (ms)", color=c["label"], fontsize=8)
    ax_cmp.set_title("Observed vs Computed", color=c["text"], fontsize=9)
    ax_cmp.grid(True, lw=0.3, alpha=0.3, color=c["grid"])
    ax_cmp.tick_params(colors=c["tick"], labelsize=8)

    for ax in (ax_tx, ax_cmp):
        for sp in ax.spines.values():
            sp.set_edgecolor(c["spine"])

    fig.suptitle(
        f"{title} | RMS={rms:.3f} ms (Nfit={nfit})",
        color=c["text"], fontsize=10,
    )

    panel = fig.add_axes([0.84, 0.28, 0.14, 0.40])
    panel.set_facecolor(c["ax_bg"])
    panel.set_xticks([])
    panel.set_yticks([])
    for sp in panel.spines.values():
        sp.set_edgecolor(c["spine"])
    panel.text(0.06, 0.95,
               "Fit review\n\n"
               "Accept: keep this fit\n"
               "Repick: choose windows again\n"
               "Skip: use previous/empty\n"
               "Drop L1/L2/L3: ignore one layer\n\n"
               "Keys: Enter=Accept, r=Repick, 1/2/3=Drop, Esc=Skip",
               transform=panel.transAxes, va="top", ha="left",
               fontsize=8, color=c["text"])

    ax_acc = fig.add_axes([0.84, 0.24, 0.14, 0.06])
    ax_rep = fig.add_axes([0.84, 0.17, 0.14, 0.06])
    ax_d1 = fig.add_axes([0.84, 0.10, 0.045, 0.06])
    ax_d2 = fig.add_axes([0.8875, 0.10, 0.045, 0.06])
    ax_d3 = fig.add_axes([0.935, 0.10, 0.045, 0.06])
    ax_skp = fig.add_axes([0.84, 0.03, 0.14, 0.06])
    hover = "#e8e8e8" if THEME == "light" else "#2a2d3d"
    b_acc = Button(ax_acc, "Accept", color=c["ax_bg"], hovercolor=hover)
    b_rep = Button(ax_rep, "Repick", color=c["ax_bg"], hovercolor=hover)
    b_d1 = Button(ax_d1, "L1", color=c["ax_bg"], hovercolor=hover)
    b_d2 = Button(ax_d2, "L2", color=c["ax_bg"], hovercolor=hover)
    b_d3 = Button(ax_d3, "L3", color=c["ax_bg"], hovercolor=hover)
    b_skp = Button(ax_skp, "Skip", color=c["ax_bg"], hovercolor=hover)
    for b in (b_acc, b_rep, b_d1, b_d2, b_d3, b_skp):
        b.label.set_color(c["text"])
        b.label.set_fontsize(8)
        for sp in b.ax.spines.values():
            sp.set_edgecolor(c["spine"])

    state = {"done": False, "choice": "accept"}

    def _finish(choice: str):
        if state["done"]:
            return
        state["done"] = True
        state["choice"] = str(choice)
        try:
            fig.canvas.stop_event_loop()
        except Exception:
            pass

    b_acc.on_clicked(lambda _e: _finish("accept"))
    b_rep.on_clicked(lambda _e: _finish("repick"))
    b_d1.on_clicked(lambda _e: _finish("drop1"))
    b_d2.on_clicked(lambda _e: _finish("drop2"))
    b_d3.on_clicked(lambda _e: _finish("drop3"))
    b_skp.on_clicked(lambda _e: _finish("skip"))

    def _on_key(evt: Any):
        if evt.key == "enter":
            _finish("accept")
        elif evt.key in ("r", "R"):
            _finish("repick")
        elif evt.key == "1":
            _finish("drop1")
        elif evt.key == "2":
            _finish("drop2")
        elif evt.key == "3":
            _finish("drop3")
        elif evt.key in ("escape", "q"):
            _finish("skip")

    def _on_close(_evt: Any):
        _finish("skip")

    fig.canvas.mpl_connect("key_press_event", _on_key)
    fig.canvas.mpl_connect("close_event", _on_close)
    plt.show(block=False)
    while not state["done"] and plt.fignum_exists(fig.number):
        try:
            fig.canvas.start_event_loop(0.05)
        except Exception:
            break
    try:
        if plt.fignum_exists(fig.number):
            plt.close(fig)
    except Exception:
        pass
    return str(state.get("choice", "accept"))


def pick_layer_windows_interactive(profile_name: str,
                                   corrected_by_shot: dict,
                                   existing_results: dict | None = None) -> dict:
    """
    Interactive layer picking for each shot side (L/R).

    Parameters
    ----------
    profile_name : str
        Profile identifier used in plot titles.
    corrected_by_shot : dict
        Output of build_corrected_pick_data(), keyed by shot id.
    existing_results : dict | None, optional
        Previously saved layer picks. If a side is skipped, previous values
        can be kept from this mapping.

    Returns
    -------
    dict
        Layer fit payload by shot and side:
        {shot_id: {"L"|"R"|"ALL": fit_result}}.
    """
    existing = existing_results or {}
    results: dict = {}
    print("\n  -- Layer windows picking (x1..x6) --")
    print("     For each shot side: pick x1,x2,x3,x4,x5,x6; Enter to skip side.")

    for shot_id in sorted(corrected_by_shot):
        rows = corrected_by_shot.get(shot_id, [])
        if not rows:
            continue
        per_side: dict = {}
        side_candidates: list = []
        left_rows = [r for r in rows if r.get("side") == "L"]
        right_rows = [r for r in rows if r.get("side") == "R"]
        min_side_points = 2
        if len(left_rows) >= min_side_points:
            side_candidates.append(("L", left_rows))
        if len(right_rows) >= min_side_points:
            side_candidates.append(("R", right_rows))
        if side_candidates and len(side_candidates) < 2 and (len(left_rows) > 0 and len(right_rows) > 0):
            print(f"     Shot {shot_id}: one side has too few picks for standalone fit "
                  f"(L={len(left_rows)}, R={len(right_rows)}).")

        if not side_candidates and len(rows) >= 3:
            side_candidates.append(("ALL", rows))
            print(f"     Shot {shot_id}: using combined side (ALL) for layer-window picking.")

        if not side_candidates:
            print(f"     Shot {shot_id}: skipped layer-window plot (need >=3 corrected picks).")
            continue

        for side, side_rows in side_candidates:
            side_rows = sorted(side_rows, key=lambda r: r.get("true_off_m", 0.0))
            x = np.array([r.get("true_off_m", 0.0) for r in side_rows], dtype=float)
            t = np.array([r["fb_interp_inline_ms"] for r in side_rows], dtype=float)
            title = f"Profile {profile_name} | Shot {shot_id} | Side {side}"
            print(f"     Opening layer-window picker: Shot {shot_id} Side {side}")

            while True:
                pick_payload = _pick_layer_windows_from_plot(x, t, title)
                windows = pick_payload.get("windows", [None, None, None])
                picked_points = pick_payload.get("picked_points", [])

                if all(w is None for w in windows):
                    prev = ((existing.get(int(shot_id), {}) or {}).get(side))
                    if prev:
                        per_side[side] = prev
                        print(f"     Shot {shot_id} Side {side}: skipped by user, kept previous layer picks.")
                    else:
                        fit_res = _fit_layers_from_windows(x, t, windows)
                        fit_res["windows"] = windows
                        fit_res["picked_points"] = picked_points
                        fit_res["n_points"] = int(len(x))
                        fit_res["skipped"] = True
                        fit_res["rms_ms"] = 0.0
                        per_side[side] = fit_res
                        print(f"     Shot {shot_id} Side {side}: skipped by user, saved empty layer selection.")
                    break

                fit_res = _fit_layers_from_windows(x, t, windows)
                fit_res["windows"] = windows
                fit_res["picked_points"] = picked_points
                fit_res["n_points"] = int(len(x))
                fit_res["skipped"] = False
                fit_res["rms_ms"] = _compute_fit_rms(x, t, fit_res)

                choice = _review_layer_fit_interactive(x, t, fit_res, title)
                if choice == "repick":
                    print(f"     Shot {shot_id} Side {side}: re-pick requested.")
                    continue
                if choice in ("drop1", "drop2", "drop3"):
                    drop_idx = {"drop1": 0, "drop2": 1, "drop3": 2}[choice]
                    fit_res = _drop_fit_layer(fit_res, drop_idx)
                    fit_res["rms_ms"] = _compute_fit_rms(x, t, fit_res)
                    print(f"     Shot {shot_id} Side {side}: dropped layer {drop_idx + 1} and reviewing again.")
                    continue
                if choice == "skip":
                    prev = ((existing.get(int(shot_id), {}) or {}).get(side))
                    if prev:
                        per_side[side] = prev
                        print(f"     Shot {shot_id} Side {side}: skipped in review, kept previous layer picks.")
                    else:
                        fit_res_empty = _fit_layers_from_windows(x, t, [None, None, None])
                        fit_res_empty["windows"] = [None, None, None]
                        fit_res_empty["picked_points"] = []
                        fit_res_empty["n_points"] = int(len(x))
                        fit_res_empty["skipped"] = True
                        fit_res_empty["rms_ms"] = 0.0
                        per_side[side] = fit_res_empty
                        print(f"     Shot {shot_id} Side {side}: skipped in review.")
                    break

                per_side[side] = fit_res
                print(f"     Shot {shot_id} Side {side}: accepted fit, RMS={fit_res.get('rms_ms', 0.0):.3f} ms")
                break

        if per_side:
            results[int(shot_id)] = per_side

    return results














# ---------------------------------------------------------------------------
# Interactive first-break picker
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
        control_file: Path | None = None,
        show_plot_controls: bool = True,
    ):
        self.data_raw = data_raw
        self.data_filt = data_filt
        self.dt_s = dt_s
        self.recv_abs = recv_abs
        self.shot_pos_m = shot_pos_m
        self.delay_ms = delay_ms
        self.shot_id = shot_id
        self.profile = profile_name
        self.qc_dir = qc_dir or core_paths.require_active_project().plots_dir_for(profile_name)

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
        self._control_file = Path(control_file) if control_file else None
        self._last_control_id = 0
        self._show_plot_controls = bool(show_plot_controls)
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
        if self._show_plot_controls:
            self.fig.subplots_adjust(left=0.06, right=0.76, bottom=0.13, top=0.91)
        else:
            self.fig.subplots_adjust(left=0.06, right=0.98, bottom=0.13, top=0.91)
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

        if self._show_plot_controls:
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

    def _poll_external_control(self):
        self._last_control_id, payload = read_latest_picker_command(self._control_file, self._last_control_id)
        if payload is None:
            return

        action = str(payload.get("action", "")).strip().lower()
        try:
            if action == "prev":
                self._go_prev()
                return
            if action == "next":
                self._save_and_finish("next")
                return
            if action == "finalize":
                self._save_and_finish("finalize")
                return
            if action == "quit":
                self._cancelled = True
                self._nav_action = "quit"
                self._done = True
                try:
                    self.fig.canvas.stop_event_loop()
                except Exception:
                    pass
                return
            if action == "auto":
                self._auto_then_redraw()
                return
            if action == "invert":
                self._toggle_invert()
                return
            if action == "timeline":
                self._toggle_timelines()
                return
            if action == "filter_toggle":
                self._filter_mode = {"none": "butter", "butter": "ormsby", "ormsby": "cutpass", "cutpass": "none"}.get(self._filter_mode, "none")
                self._refresh_mode_buttons()
                self._redraw()
                return
            if action == "save_image":
                self._save_qc_image()
                return
            if action == "gain":
                self._on_gain_mode(str(payload.get("value", "norm")))
                return
            if action == "agc_stat":
                self._on_agc_stat(str(payload.get("value", "rms")))
                return
            if action == "display":
                self._on_display_mode(str(payload.get("value", "both")))
                return
            if action == "agc_window":
                self._on_agc_window(float(payload.get("value", self._agc_window_ms)))
                self._redraw()
                return
            if action == "wiggle_scale":
                self._on_wiggle_stretch(float(payload.get("value", self._wiggle_stretch)))
                self._redraw()
                return
            if action == "set_filter":
                f1 = float(payload.get("f1", self._f1))
                f2 = float(payload.get("f2", self._f2))
                f3 = float(payload.get("f3", self._f3))
                f4 = float(payload.get("f4", self._f4))
                if f1 < f2 < f3 < f4:
                    self._f1, self._f2, self._f3, self._f4 = f1, f2, f3, f4
                    self._schedule_filter_update()
                return
        except Exception as exc:
            print(f"     [WARN] control command failed: {action} ({exc})")

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
                self._poll_external_control()
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
# Picks persistence  (JSON -- always RAW, no static applied)
# ---------------------------------------------------------------------------































# ---------------------------------------------------------------------------
# Export: picks_clean.txt  (compatible with lvl.ipynb)
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# Export: Excel workbook  (Config + per-shot formulas + Analysis)
# ---------------------------------------------------------------------------

_FILL_HDR  = PatternFill("solid", fgColor="2E4057")
_FILL_IN   = PatternFill("solid", fgColor="FFFACD")   # lemon  = editable input
_FILL_FORM = PatternFill("solid", fgColor="E8F4F8")   # blue   = Excel formula
_FILL_OK   = PatternFill("solid", fgColor="DDFFDD")   # green  = depth result
_FONT_HDR  = Font(bold=True, color="FFFFFF")
_FONT_BOLD = Font(bold=True)
_FONT_ITA  = Font(italic=True, color="888888")








# ---------------------------------------------------------------------------
# Export: T-X picks plot  (time increases downward + analysis overlay)
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# Main processing pipeline
# ---------------------------------------------------------------------------

# process_profile moved to src/refraction/pipeline.py (shared with lvl_studio.py).
from src.refraction.pipeline import process_profile






# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    """CLI entry point: parse arguments, resolve targets, and process profiles."""
    from src.common.logging_setup import setup_file_logging
    log_path = setup_file_logging("lvl_refraction")
    print(f"{APP_NAME} (refraction engine) v{APP_VERSION}")
    print(f"Logging to: {log_path}")
    parser = argparse.ArgumentParser(
        prog="lvl_refraction.py",
        description="LVL refraction seismic: interactive picking + analysis.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python lvl_refraction.py 120\n"
            "  python lvl_refraction.py 120 --export-only\n"
            "  python lvl_refraction.py --all\n"
        ),
    )
    parser.add_argument("profile", nargs="?", default=None,
                        help="Profile folder name (e.g. 120 or 214_A)")
    parser.add_argument("geometry", nargs="?", default=None,
                        help="Legacy geometry override (positional): 100, 200, geometry100, geometry200")
    parser.add_argument("--geom", default=None,
                        help="Geometry override: 100, 200, geometry100, geometry200")
    parser.add_argument("--all", action="store_true",
                        help="Process all profiles in PROFILES")
    parser.add_argument("--export-only", action="store_true",
                        help="Skip picking; re-analyse and re-export existing picks")
    parser.add_argument("--perp-by-shot", default=None,
                        help="Override perpendicular offset per shot, e.g. 1:0,2:3.5,3:0")
    parser.add_argument("--inline-shift-by-shot", default=None,
                        help="Optional inline shift per shot (m), e.g. 1:50,2:0,3:50")
    parser.add_argument("--perp-excel", default=None,
                        help="Excel file with FFID/profile/perpendicular offsets defaults")
    parser.add_argument("--perp-sheet", default=None,
                        help="Excel sheet name/index for --perp-excel (default: first sheet)")
    parser.add_argument("--perp-col-ffid", default="A",
                        help="Excel column for FFID (default: A)")
    parser.add_argument("--perp-col-perp", default="D",
                        help="Excel column for perpendicular offset (default: D)")
    parser.add_argument("--perp-col-profile", default="F",
                        help="Excel column for profile token, e.g. LVL150 (default: F)")
    parser.add_argument("--perp-col-inline-shift", default=None,
                        help="Optional Excel column for inline X-shift in meters")
    parser.add_argument("--coord-excel", action="append", default=None,
                        help="Optional LVL coordinate workbook path (can be passed multiple times)")
    parser.add_argument("--device-type", default="sw_maps", choices=["sw_maps", "geomax"],
                        help="Device/source type hint for loaders (default: sw_maps)")
    parser.add_argument("--control-file", default=None,
                        help="Path to JSON command file for external GUI live control")
    parser.add_argument("--minimal-plot-controls", action="store_true",
                        help="Hide in-plot buttons/widgets and use external GUI controls")
    parser.add_argument("--no-layer-pick", action="store_true",
                        help="Skip interactive x0..x5 layer-window picking")
    parser.add_argument("--auto-pick", action="store_true",
                        help="Pick every shot automatically with Picker Engine V2 "
                             "instead of opening the interactive picker "
                             "(preprocessing -> features -> likelihood -> coherence "
                             "-> path optimization -> confidence/quality). Runs "
                             "straight through to layer analysis and export.")
    parser.add_argument("--project", default=None,
                        help="Path to a project folder (created with lvl_studio or "
                             "src.io.project_io.create_project). If omitted, you'll be "
                             "prompted to pick or create one.")
    args = parser.parse_args()

    # -- resolve the active project -----------------------------------------
    project = None
    if args.project:
        try:
            project = pio.open_project(args.project)
        except pio.ProjectNotFoundError:
            print(f"[ERROR] No project.json found at: {args.project}")
            return
    else:
        try:
            root = tk.Tk()
            root.withdraw()
            sel = filedialog.askdirectory(
                title="Select or create a project folder",
                initialdir=str(core_paths.PROJECTS_ROOT) if core_paths.PROJECTS_ROOT.exists() else str(Path.cwd()),
                mustexist=False,
            )
            root.destroy()
        except Exception:
            sel = None
        if not sel:
            print("[ERROR] No project selected. Re-run with --project <path>.")
            return
        if pio.is_project_dir(sel):
            project = pio.open_project(sel)
        else:
            name = Path(sel).name
            project = pio.create_project(name, parent_dir=Path(sel).parent)
            print(f"  Created new project '{project.name}' at {project.root}")

    project.activate()
    print(f"  Project : {project.name}  ({project.root})")
    if project.raw_folder is None:
        print(
            "  [WARN] Project has no raw_folder set. Set one via lvl_studio's "
            "'Project Settings...', or edit project.json directly."
        )

    try:
        perp_override = parse_shot_value_map(args.perp_by_shot)
    except Exception as exc:
        print(f"[ERROR] Invalid --perp-by-shot value: {exc}")
        return

    try:
        inline_shift_override = parse_shot_value_map(args.inline_shift_by_shot)
    except Exception as exc:
        print(f"[ERROR] Invalid --inline-shift-by-shot value: {exc}")
        return

    geom_override: int | None = None
    geom_raw = args.geom if args.geom is not None else args.geometry
    if geom_raw is not None:
        g = str(geom_raw).strip().lower().replace(" ", "")
        if g in ("100", "geometry100"):
            geom_override = 100
        elif g in ("200", "geometry200"):
            geom_override = 200
        else:
            print("[ERROR] Invalid geometry override. Use one of:")
            print("        100, 200, geometry100, geometry200")
            return

    project.paths.results_dir.mkdir(parents=True, exist_ok=True)

    sheet_val: str | int | None = args.perp_sheet
    if sheet_val is not None:
        s_try = str(sheet_val).strip()
        if s_try.isdigit():
            sheet_val = int(s_try)

    report_paths: list = []
    if args.perp_excel:
        report_paths = [str(Path(args.perp_excel))]
    elif project.metadata_folder is not None:
        report_paths = [str(p) for p in discover_field_report_excels(project.metadata_folder)]
        if report_paths:
            print("  Auto field reports detected:")
            for p in report_paths:
                try:
                    print(f"    - {Path(p).relative_to(project.metadata_folder.parent)}")
                except Exception:
                    print(f"    - {p}")

    geometry_paths: list = []
    if args.coord_excel:
        for cp in args.coord_excel:
            try:
                geometry_paths.append(str(Path(cp)))
            except Exception:
                pass
    for p in getattr(project, "manual_geometry_files", []):
        sp = str(p)
        if Path(p).exists() and sp not in geometry_paths:
            geometry_paths.append(sp)
    if project.geometry_folder is not None:
        for pattern in ("LVL*.xls*", "LVL*.csv", "LVL*.txt", "LVL*.dat"):
            for p in project.geometry_folder.glob(pattern):
                name = p.name.lower()
                if "field_report" in name or "fieldreport" in name:
                    continue
                sp = str(p)
                if sp not in geometry_paths:
                    geometry_paths.append(sp)

    print(f"  Device type: {args.device_type}")

    perp_excel_cfg = {
        "paths": report_paths,
        "sheet": sheet_val,
        "ffid_col": args.perp_col_ffid,
        "perp_col": args.perp_col_perp,
        "profile_col": args.perp_col_profile,
        "inline_shift_col": args.perp_col_inline_shift,
        "geometry_paths": geometry_paths,
        "device_type": args.device_type,
    }

    raw_dir = project.raw_folder
    if raw_dir is None:
        print("[ERROR] Project has no raw_folder set (see 'Project Settings...' in lvl_studio).")
        return

    if args.all:
        targets = discover_profile_folders(raw_dir)
        if not targets:
            print(f"[ERROR] No profile folders with SEG2 files found under {raw_dir}.")
            return
    elif args.profile:
        targets = [args.profile]
    else:
        chosen_profile = None
        try:
            root = tk.Tk()
            root.withdraw()
            sel = filedialog.askdirectory(
                title="Select profile data folder",
                initialdir=str(raw_dir),
                mustexist=True,
            )
            root.destroy()
            if sel:
                p = Path(sel)
                has_seg2 = any(p.glob("*.seg2")) or any(p.glob("*.SEG2"))
                if has_seg2:
                    chosen_profile = p.name
        except Exception:
            chosen_profile = None

        if chosen_profile:
            print(f"  Selected folder -> profile '{chosen_profile}'")
            targets = [chosen_profile]
        else:
            print("\nAvailable profiles:")
            print(f"  {'Name':<16}  {'Geom':>6}  {'Data folder':<30}")
            print(f"  {'-'*16}  {'-'*6}  {'-'*30}")
            dynamic_profiles = discover_profile_folders(raw_dir)
            if not dynamic_profiles:
                print("  (none found)")
            for pname in dynamic_profiles:
                pcfg = PROFILES.get(pname, {"geom": 200})
                folder = raw_dir / pname
                n_seg2 = len(list(folder.glob("*.seg2"))) + len(list(folder.glob("*.SEG2")))
                print(f"  {pname:<16}  {int(pcfg.get('geom', 200)):>5}m  "
                      f"found ({n_seg2} SEG2 files)")
            print(f"\nUsage:  python lvl_refraction.py <profile> --project <path>  [--export-only]")
            return

    for pname in targets:
        print(f"\n{'='*60}\nProfile : {pname}\n{'='*60}")
        process_profile(pname, pick_mode=not args.export_only,
                        geom_override=geom_override,
                        perp_by_shot_override=perp_override,
                        inline_shift_by_shot_override=inline_shift_override,
                        perp_excel_cfg=perp_excel_cfg,
                        enable_layer_pick=not args.no_layer_pick,
                        control_file=(Path(args.control_file) if args.control_file else None),
                        show_plot_controls=not args.minimal_plot_controls,
                        raw_dir=raw_dir,
                        auto_pick=args.auto_pick)
        project.register_profile(pname)

    pio.save_project(project)
    print("\nAll done.")


if __name__ == "__main__":
    main()