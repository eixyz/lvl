from __future__ import annotations

import json
import os
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from obspy import read as obspy_read

try:
    from PyQt6 import QtCore, QtGui, QtWidgets
    QT_API = "PyQt6"
except Exception:
    from PyQt5 import QtCore, QtGui, QtWidgets
    QT_API = "PyQt5"

from matplotlib.backend_bases import MouseButton
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


def _dock_area_left() -> Any:
    area = getattr(QtCore.Qt, "DockWidgetArea", None)
    return area.LeftDockWidgetArea if area is not None else QtCore.Qt.LeftDockWidgetArea


def _dock_area_right() -> Any:
    area = getattr(QtCore.Qt, "DockWidgetArea", None)
    return area.RightDockWidgetArea if area is not None else QtCore.Qt.RightDockWidgetArea


def _dock_area_bottom() -> Any:
    area = getattr(QtCore.Qt, "DockWidgetArea", None)
    return area.BottomDockWidgetArea if area is not None else QtCore.Qt.BottomDockWidgetArea


def _orient_horizontal() -> Any:
    orient = getattr(QtCore.Qt, "Orientation", None)
    return orient.Horizontal if orient is not None else QtCore.Qt.Horizontal


def _orient_vertical() -> Any:
    orient = getattr(QtCore.Qt, "Orientation", None)
    return orient.Vertical if orient is not None else QtCore.Qt.Vertical


def _dock_features(*feature_names: str) -> Any:
    feat_enum = getattr(QtWidgets.QDockWidget, "DockWidgetFeature", None)
    if feat_enum is not None:
        value = feat_enum.NoDockWidgetFeatures
        for name in feature_names:
            value = value | getattr(feat_enum, name)
        return value
    value = 0
    for name in feature_names:
        value |= getattr(QtWidgets.QDockWidget, name)
    return value


DATA_DIR = Path(__file__).resolve().parent.parent / "data"
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "results"


@dataclass
class ShotGather:
    source_file: Path
    shot_label: str
    ffid: int
    dt_ms: float
    delay_ms: float
    data: np.ndarray
    header: dict[str, Any]

    @property
    def n_traces(self) -> int:
        return int(self.data.shape[0])

    @property
    def n_samples(self) -> int:
        return int(self.data.shape[1])


@dataclass
class StudioSettings:
    clip_pct: float = 99.0
    max_traces_render: int = 384
    cmap: str = "gray"
    display_mode: str = "wiggle"
    wiggle_scale: float = 1.0
    clip_factor: float = 2.0
    use_seg2_delay: bool = True
    trigger_static_ms: float = 6.0
    extra_delay_ms: float = 0.0
    polarity: str = "invert"
    gain_mode: str = "norm"
    agc_window_ms: float = 200.0
    agc_stat: str = "rms"
    filter_mode: str = "butter"
    f1: float = 2.0
    f2: float = 4.0
    f3: float = 140.0
    f4: float = 180.0
    butter_order: int = 4
    auto_pick_mode: str = "hilbert_env"
    hilbert_target: str = "onset"
    hilbert_onset_pct: float = 0.08
    pick_order: str = "F>G>X"
    manual_snap_win_ms: float = 6.0
    geometry_override: str = "auto"
    x_axis_mode: str = "geom_x"
    invert_y_axis: bool = True
    t_min_ms: float = -10.0
    t_max_ms: float = 150.0
    show_grid: bool = False
    show_timelines: bool = True
    theme: str = "light"
    use_seg2_shot_position: bool = False


class SeismicDisplay(QtWidgets.QWidget):
    def __init__(self, parent: QtWidgets.QWidget | None = None):
        super().__init__(parent)
        self.fig = Figure(figsize=(11, 7), dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.fig.add_subplot(111)
        self.last_stride = 1
        self.last_n_view = 0
        self.last_view_to_trace: list[int] = []
        self.last_view_x: list[float] = []
        self._current_shot: ShotGather | None = None

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.canvas)

    def draw_shot(
        self,
        shot: ShotGather,
        settings: StudioSettings,
        data_view_full: np.ndarray,
        t0_ms: float,
        x_positions: np.ndarray,
        x_label: str,
        shot_x: float | None = None,
        layers: list[tuple[str, str, dict[int, float], bool]] | None = None,
    ):
        self._current_shot = shot
        self.ax.clear()

        dark = str(settings.theme).lower() == "dark"
        fig_bg = "#0f1117" if dark else "white"
        ax_bg = "#1a1d2e" if dark else "#f5f5f5"
        trace_col = "#e6e6e6" if dark else "#111111"
        tick_col = "#cccccc" if dark else "#333333"
        self.fig.patch.set_facecolor(fig_bg)
        self.ax.set_facecolor(ax_bg)

        data = np.asarray(data_view_full, dtype=np.float32)
        n_tr = data.shape[0]
        stride = max(1, int(np.ceil(n_tr / float(max(1, settings.max_traces_render)))))
        self.last_stride = stride
        self.last_view_to_trace = list(range(0, n_tr, stride))
        self.last_view_x = [float(x_positions[i]) for i in self.last_view_to_trace]
        data_view = data[::stride, :]
        self.last_n_view = int(data_view.shape[0])

        clip = float(np.percentile(np.abs(data_view), settings.clip_pct))
        if not np.isfinite(clip) or clip <= 0:
            clip = 1.0

        dt = float(shot.dt_ms)
        t0 = float(t0_ms)
        t = t0 + np.arange(data_view.shape[1], dtype=np.float32) * dt
        t_min = float(settings.t_min_ms)
        t_max = float(settings.t_max_ms)
        if t_max <= t_min:
            t_max = t_min + 1.0
        tmask = (t >= t_min) & (t <= t_max)
        if not np.any(tmask):
            tmask[:] = True
        data_view = data_view[:, tmask]
        t = t[tmask]
        t0v = float(t[0])
        t1v = float(t[-1])
        xv = np.asarray(self.last_view_x, dtype=float)
        if xv.size == 0:
            xv = np.arange(1, data_view.shape[0] + 1, dtype=float)

        if xv.size >= 2:
            x_step = float(np.median(np.diff(xv)))
            if abs(x_step) < 1e-9:
                x_step = 1.0
        else:
            x_step = 1.0
        dx = abs(x_step)
        x0 = float(xv[0] - 0.5 * x_step)
        x1 = float(xv[-1] + 0.5 * x_step)

        mode = str(settings.display_mode).strip().lower()
        if mode in ("vd", "both"):
            vmax = float(np.percentile(np.abs(data_view), 98)) if data_view.size else 1.0
            if vmax < 1e-12:
                vmax = 1.0
            self.ax.imshow(
                data_view.T,
                cmap=settings.cmap,
                aspect="auto",
                interpolation="nearest",
                vmin=-vmax,
                vmax=vmax,
                extent=[x0, x1, t1v, t0v],
                alpha=0.45 if mode == "both" else 1.0,
                zorder=1,
            )

        if mode in ("wiggle", "both"):
            # Script-compatible per-trace normalization (median-std clip).
            stds = data_view.std(axis=1).astype(float)
            valid = stds[stds > 1e-20]
            med = float(np.median(valid)) if valid.size else 1.0
            norms = np.clip(stds, med * 0.3, med * 3.0)
            norms = np.where(norms > 1e-20, norms, med)
            scale = norms * float(settings.clip_factor)
            defl = dx * float(settings.wiggle_scale)
            for i in range(data_view.shape[0]):
                base_x = float(xv[i])
                tr = data_view[i].astype(float)
                tr_n = np.clip(tr / scale[i], -1.0, 1.0)
                x = base_x + tr_n * defl
                self.ax.plot(x, t, color=trace_col, linewidth=0.5, alpha=0.85, zorder=4)
                pos = np.where(tr_n > 0.0, tr_n, 0.0)
                self.ax.fill_betweenx(t, base_x, base_x + pos * defl, color=trace_col, alpha=0.18, zorder=3)

        # Shot-position (orange dashed) and zero-time (green dashed) guides.
        if shot_x is not None and np.isfinite(float(shot_x)):
            self.ax.axvline(float(shot_x), color="#ffcc00", lw=1.2, ls="--", alpha=0.75, zorder=5)
        if t0v <= 0.0 <= t1v:
            self.ax.axhline(0.0, color="#00aa66", lw=1.0, ls=":", alpha=0.8, zorder=5)

        if settings.show_timelines:
            start_ms = int(np.ceil(max(t0v, 0.0) / 10.0) * 10)
            for ms in range(start_ms, int(t1v) + 1, 10):
                if not (t0v <= ms <= t1v) or ms == 0:
                    continue
                major = (ms % 20 == 0)
                self.ax.axhline(
                    float(ms),
                    color="#aaaaaa" if dark else "#666666",
                    lw=0.5 if major else 0.3,
                    ls=":" if major else "--",
                    alpha=0.35 if major else 0.2,
                    zorder=2,
                )

        legend_handles = []
        for name, color, layer_picks, is_active in (layers or []):
            if not layer_picks:
                continue
            xs2: list[float] = []
            ys2: list[float] = []
            for view_idx, tr_idx in enumerate(self.last_view_to_trace):
                if tr_idx in layer_picks:
                    xs2.append(float(xv[view_idx]))
                    ys2.append(float(layer_picks[tr_idx]))
            if not xs2:
                continue
            lw = 1.4 if is_active else 1.0
            msz = 3.2 if is_active else 2.4
            alpha = 1.0 if is_active else 0.85
            zorder = 7 if is_active else 6
            label = f"{name} (active)" if is_active else name
            (line,) = self.ax.plot(
                xs2, ys2, color=color, linewidth=lw, marker="o", markersize=msz,
                alpha=alpha, zorder=zorder, label=label,
            )
            legend_handles.append(line)
        if len(legend_handles) > 1:
            leg = self.ax.legend(handles=legend_handles, loc="upper right", fontsize=7, framealpha=0.85)
            for txt in leg.get_texts():
                txt.set_color(tick_col)

        self.ax.set_title(
            f"{shot.shot_label} | FFID {shot.ffid} | traces {shot.n_traces} (view {data_view.shape[0]}) | samples {shot.n_samples}",
            fontsize=10,
        )
        self.ax.set_xlabel(x_label)
        self.ax.set_ylabel("Time (ms)")
        self.ax.set_xlim(x0, x1)
        self.ax.set_ylim(t0v, t1v)
        self.ax.tick_params(colors=tick_col)
        self.ax.xaxis.label.set_color(tick_col)
        self.ax.yaxis.label.set_color(tick_col)
        self.ax.title.set_color(tick_col)
        for spine in self.ax.spines.values():
            spine.set_edgecolor(tick_col)
        if settings.invert_y_axis:
            self.ax.invert_yaxis()
        self.ax.margins(x=0.0, y=0.0)
        self.ax.grid(bool(settings.show_grid), alpha=0.15)
        self.fig.tight_layout()
        self.canvas.draw_idle()


class LvlStudioWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LVL Studio - Independent Seismic Workspace")
        self.resize(1500, 900)

        self.settings = StudioSettings()
        self.shots: list[ShotGather] = []
        # Multi-layer picks: "Mine" is the default, always-active editable layer.
        # Imported layers are overlaid read-only until made active.
        self.pick_layers: dict[str, dict[int, dict[int, float]]] = {"Mine": {}}
        self.layer_colors: dict[str, str] = {"Mine": "#d7191c"}
        self.layer_visible: dict[str, bool] = {"Mine": True}
        self.active_layer: str = "Mine"
        self._layer_color_cycle = ["#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a", "#b15928", "#a6cee3", "#e31a1c", "#008080"]
        self.current_idx = -1
        self.pick_enabled = True
        self._drag_pick = False
        self._drag_delete = False
        self._range_anchor: tuple[int, float] | None = None
        self._last_drag_trace: int | None = None
        self._last_drag_t: float | None = None
        self._display_window: QtWidgets.QMainWindow | None = None
        self._current_folder: Path | None = None
        self._last_analysis: dict[str, Any] | None = None
        self._custom_geom_positions: np.ndarray | None = None
        self._custom_geom_path: str | None = None
        self._view_xlim: tuple[float, float] | None = None
        self._view_ylim: tuple[float, float] | None = None

        self._build_ui()
        self._build_actions()
        self._build_shortcuts()

    @property
    def picks_by_shot(self) -> dict[int, dict[int, float]]:
        """Picks of the currently active layer (kept as a property so all
        existing pick-editing code transparently operates on the active layer)."""
        return self.pick_layers.setdefault(self.active_layer, {})

    @picks_by_shot.setter
    def picks_by_shot(self, value: dict[int, dict[int, float]]):
        self.pick_layers[self.active_layer] = value

    def _reset_pick_layers(self):
        self.pick_layers = {"Mine": {}}
        self.layer_colors = {"Mine": "#d7191c"}
        self.layer_visible = {"Mine": True}
        self.active_layer = "Mine"
        self._refresh_layers_list()

    def _build_ui(self):
        self.display = SeismicDisplay(self)
        self.setCentralWidget(self.display)

        self.btn_open = QtWidgets.QPushButton("Open SEG2 Folder")
        self.btn_prev = QtWidgets.QPushButton("Prev Shot")
        self.btn_next = QtWidgets.QPushButton("Next Shot")

        self.spin_clip = QtWidgets.QDoubleSpinBox()
        self.spin_clip.setRange(90.0, 100.0)
        self.spin_clip.setDecimals(2)
        self.spin_clip.setSingleStep(0.25)
        self.spin_clip.setValue(self.settings.clip_pct)

        self.spin_max_tr = QtWidgets.QSpinBox()
        self.spin_max_tr.setRange(32, 5000)
        self.spin_max_tr.setSingleStep(16)
        self.spin_max_tr.setValue(self.settings.max_traces_render)

        self.spin_wiggle = QtWidgets.QDoubleSpinBox()
        self.spin_wiggle.setRange(0.1, 4.0)
        self.spin_wiggle.setDecimals(2)
        self.spin_wiggle.setSingleStep(0.1)
        self.spin_wiggle.setValue(self.settings.wiggle_scale)

        self.cmb_mode = QtWidgets.QComboBox()
        self.cmb_mode.addItems(["wiggle", "vd", "both"])
        self.cmb_mode.setCurrentText(self.settings.display_mode)

        self.cmb_x_axis = QtWidgets.QComboBox()
        self.cmb_x_axis.addItems(["geom_x", "offset", "channel"])
        self.cmb_x_axis.setCurrentText(self.settings.x_axis_mode)

        self.chk_invert_y = QtWidgets.QCheckBox("Invert Y axis")
        self.chk_invert_y.setChecked(self.settings.invert_y_axis)
        self.chk_grid = QtWidgets.QCheckBox("Show grid")
        self.chk_grid.setChecked(self.settings.show_grid)
        self.chk_timelines = QtWidgets.QCheckBox("Timelines (10/20 ms)")
        self.chk_timelines.setChecked(self.settings.show_timelines)
        self.cmb_theme = QtWidgets.QComboBox()
        self.cmb_theme.addItems(["light", "dark"])
        self.cmb_theme.setCurrentText(self.settings.theme)

        self.spin_tmin = QtWidgets.QDoubleSpinBox()
        self.spin_tmin.setRange(-1000.0, 5000.0)
        self.spin_tmin.setDecimals(2)
        self.spin_tmin.setValue(self.settings.t_min_ms)

        self.spin_tmax = QtWidgets.QDoubleSpinBox()
        self.spin_tmax.setRange(-1000.0, 5000.0)
        self.spin_tmax.setDecimals(2)
        self.spin_tmax.setValue(self.settings.t_max_ms)

        self.cmb_cmap = QtWidgets.QComboBox()
        self.cmb_cmap.addItems(["gray", "seismic", "viridis", "magma", "cividis", "turbo", "plasma"])
        self.cmb_cmap.setCurrentText(self.settings.cmap)

        self.btn_pick_toggle = QtWidgets.QPushButton("Picker: ON")
        self.btn_pick_toggle.setCheckable(True)
        self.btn_pick_toggle.setChecked(True)
        self.btn_auto_pick = QtWidgets.QPushButton("Auto Pick (Current Shot)")
        self.btn_save_picks = QtWidgets.QPushButton("Save Picks")
        self.btn_clear_picks = QtWidgets.QPushButton("Clear Shot Picks")
        self.chk_layer_pick = QtWidgets.QCheckBox("Layer pick UI")
        self.chk_layer_pick.setChecked(False)
        self.btn_run_compute = QtWidgets.QPushButton("Run Full Computation")
        self.btn_run_analysis = QtWidgets.QPushButton("Run Interactive Analysis")

        self.chk_use_seg2_delay = QtWidgets.QCheckBox("Use SEG2 delay")
        self.chk_use_seg2_delay.setChecked(self.settings.use_seg2_delay)
        self.chk_use_seg2_shotpos = QtWidgets.QCheckBox("Use SEG2 shot position")
        self.chk_use_seg2_shotpos.setChecked(self.settings.use_seg2_shot_position)
        self.spin_trigger_static = QtWidgets.QDoubleSpinBox()
        self.spin_trigger_static.setRange(-200.0, 200.0)
        self.spin_trigger_static.setDecimals(2)
        self.spin_trigger_static.setSingleStep(0.5)
        self.spin_trigger_static.setValue(self.settings.trigger_static_ms)
        self.spin_extra_delay = QtWidgets.QDoubleSpinBox()
        self.spin_extra_delay.setRange(-200.0, 200.0)
        self.spin_extra_delay.setDecimals(2)
        self.spin_extra_delay.setSingleStep(0.5)
        self.spin_extra_delay.setValue(self.settings.extra_delay_ms)

        self.cmb_polarity = QtWidgets.QComboBox()
        self.cmb_polarity.addItems(["normal", "invert"])
        self.cmb_polarity.setCurrentText(self.settings.polarity)

        self.cmb_gain_mode = QtWidgets.QComboBox()
        self.cmb_gain_mode.addItems(["none", "norm", "agc"])
        self.cmb_gain_mode.setCurrentText(self.settings.gain_mode)
        self.spin_agc_window = QtWidgets.QDoubleSpinBox()
        self.spin_agc_window.setRange(1.0, 2000.0)
        self.spin_agc_window.setDecimals(1)
        self.spin_agc_window.setSingleStep(10.0)
        self.spin_agc_window.setValue(self.settings.agc_window_ms)
        self.cmb_agc_stat = QtWidgets.QComboBox()
        self.cmb_agc_stat.addItems(["rms", "mean"])
        self.cmb_agc_stat.setCurrentText(self.settings.agc_stat)

        self.cmb_filter_mode = QtWidgets.QComboBox()
        self.cmb_filter_mode.addItems(["none", "butter", "ormsby", "cutpass"])
        self.cmb_filter_mode.setCurrentText(self.settings.filter_mode)
        self.spin_f1 = QtWidgets.QDoubleSpinBox()
        self.spin_f1.setRange(0.1, 500.0)
        self.spin_f1.setDecimals(2)
        self.spin_f1.setValue(self.settings.f1)
        self.spin_f2 = QtWidgets.QDoubleSpinBox()
        self.spin_f2.setRange(0.1, 500.0)
        self.spin_f2.setDecimals(2)
        self.spin_f2.setValue(self.settings.f2)
        self.spin_f3 = QtWidgets.QDoubleSpinBox()
        self.spin_f3.setRange(0.1, 500.0)
        self.spin_f3.setDecimals(2)
        self.spin_f3.setValue(self.settings.f3)
        self.spin_f4 = QtWidgets.QDoubleSpinBox()
        self.spin_f4.setRange(0.1, 500.0)
        self.spin_f4.setDecimals(2)
        self.spin_f4.setValue(self.settings.f4)
        self.spin_butter_order = QtWidgets.QSpinBox()
        self.spin_butter_order.setRange(1, 12)
        self.spin_butter_order.setValue(self.settings.butter_order)

        self.cmb_auto_mode = QtWidgets.QComboBox()
        self.cmb_auto_mode.addItems(["stalta", "maxdiff_zero", "hilbert_env"])
        self.cmb_auto_mode.setCurrentText(self.settings.auto_pick_mode)
        self.cmb_hilbert_target = QtWidgets.QComboBox()
        self.cmb_hilbert_target.addItems(["onset", "peak"])
        self.cmb_hilbert_target.setCurrentText(self.settings.hilbert_target)
        self.spin_hilbert_onset_pct = QtWidgets.QDoubleSpinBox()
        self.spin_hilbert_onset_pct.setRange(0.0, 1.0)
        self.spin_hilbert_onset_pct.setDecimals(3)
        self.spin_hilbert_onset_pct.setSingleStep(0.01)
        self.spin_hilbert_onset_pct.setValue(self.settings.hilbert_onset_pct)
        self.cmb_pick_order = QtWidgets.QComboBox()
        self.cmb_pick_order.addItems(["F>G>X", "G>F>X", "RAW>X"])
        self.cmb_pick_order.setCurrentText(self.settings.pick_order)
        self.spin_manual_snap = QtWidgets.QDoubleSpinBox()
        self.spin_manual_snap.setRange(0.1, 50.0)
        self.spin_manual_snap.setDecimals(2)
        self.spin_manual_snap.setSingleStep(0.5)
        self.spin_manual_snap.setValue(self.settings.manual_snap_win_ms)

        self.cmb_geom_view = QtWidgets.QComboBox()
        self.cmb_geom_view.addItems(["auto", "100", "200", "custom"])
        self.cmb_geom_view.setCurrentText(self.settings.geometry_override)

        self.layers_list = QtWidgets.QListWidget(self)
        self.layers_list.setToolTip("Checked = visible on plot. Select a layer and click 'Set Active' to edit it.")
        self.btn_import_layer = QtWidgets.QPushButton("Import Picks File...")
        self.btn_set_active_layer = QtWidgets.QPushButton("Set Active Layer")
        self.btn_remove_layer = QtWidgets.QPushButton("Remove Layer")
        self.lbl_active_layer = QtWidgets.QLabel("Active layer: Mine")

        controls = QtWidgets.QToolBox(self)

        page_nav = QtWidgets.QWidget(self)
        nav_form = QtWidgets.QFormLayout(page_nav)
        nav_form.setContentsMargins(6, 6, 6, 6)
        nav_form.addRow(self.btn_open)
        nav_form.addRow(self.btn_prev, self.btn_next)

        page_disp = QtWidgets.QWidget(self)
        disp_form = QtWidgets.QFormLayout(page_disp)
        disp_form.setContentsMargins(6, 6, 6, 6)
        disp_form.addRow("Display mode", self.cmb_mode)
        disp_form.addRow("X axis", self.cmb_x_axis)
        disp_form.addRow(self.chk_invert_y)
        disp_form.addRow(self.chk_grid)
        disp_form.addRow(self.chk_timelines)
        disp_form.addRow("Theme", self.cmb_theme)
        disp_form.addRow("Time min ms", self.spin_tmin)
        disp_form.addRow("Time max ms", self.spin_tmax)
        disp_form.addRow("Clip percentile", self.spin_clip)
        disp_form.addRow("Max traces in view", self.spin_max_tr)
        disp_form.addRow("Wiggle scale", self.spin_wiggle)
        disp_form.addRow("Colormap", self.cmb_cmap)

        page_proc = QtWidgets.QWidget(self)
        proc_form = QtWidgets.QFormLayout(page_proc)
        proc_form.setContentsMargins(6, 6, 6, 6)
        proc_form.addRow("Geometry", self.cmb_geom_view)
        proc_form.addRow(self.chk_use_seg2_delay)
        proc_form.addRow(self.chk_use_seg2_shotpos)
        proc_form.addRow("TRIGGER_STATIC_MS", self.spin_trigger_static)
        proc_form.addRow("Extra delay ms", self.spin_extra_delay)
        proc_form.addRow("Polarity", self.cmb_polarity)
        proc_form.addRow("Gain", self.cmb_gain_mode)
        proc_form.addRow("AGC window ms", self.spin_agc_window)
        proc_form.addRow("AGC stat", self.cmb_agc_stat)
        proc_form.addRow("Filter", self.cmb_filter_mode)
        proc_form.addRow("Butter order", self.spin_butter_order)
        proc_form.addRow("f1 low", self.spin_f1)
        proc_form.addRow("f2 low pass", self.spin_f2)
        proc_form.addRow("f3 high pass", self.spin_f3)
        proc_form.addRow("f4 high", self.spin_f4)

        page_pick = QtWidgets.QWidget(self)
        pick_form = QtWidgets.QFormLayout(page_pick)
        pick_form.setContentsMargins(6, 6, 6, 6)
        pick_form.addRow("Auto pick mode", self.cmb_auto_mode)
        pick_form.addRow("Hilbert target", self.cmb_hilbert_target)
        pick_form.addRow("HILBERT_ONSET_PCT", self.spin_hilbert_onset_pct)
        pick_form.addRow("Pick order", self.cmb_pick_order)
        pick_form.addRow("Manual snap win ms", self.spin_manual_snap)
        pick_form.addRow(self.btn_pick_toggle)
        pick_form.addRow(self.btn_auto_pick)
        pick_form.addRow(self.btn_save_picks)
        pick_form.addRow(self.btn_clear_picks)

        page_compute = QtWidgets.QWidget(self)
        cmp_form = QtWidgets.QFormLayout(page_compute)
        cmp_form.setContentsMargins(6, 6, 6, 6)
        cmp_form.addRow(self.chk_layer_pick)
        cmp_form.addRow(self.btn_run_analysis)
        cmp_form.addRow(self.btn_run_compute)

        page_layers = QtWidgets.QWidget(self)
        layers_v = QtWidgets.QVBoxLayout(page_layers)
        layers_v.setContentsMargins(6, 6, 6, 6)
        layers_v.addWidget(self.lbl_active_layer)
        layers_v.addWidget(self.layers_list)
        layers_v.addWidget(self.btn_import_layer)
        layers_v.addWidget(self.btn_set_active_layer)
        layers_v.addWidget(self.btn_remove_layer)

        controls.addItem(page_nav, "Navigation")
        controls.addItem(page_disp, "Display")
        controls.addItem(page_proc, "Processing")
        controls.addItem(page_layers, "Pick Layers")
        controls.addItem(page_pick, "Picking")
        controls.addItem(page_compute, "Compute")

        self.control_dock = QtWidgets.QDockWidget("Controls", self)
        self.control_dock.setWidget(controls)
        self.control_dock.setFeatures(_dock_features("DockWidgetMovable", "DockWidgetFloatable"))
        self.control_dock.setMinimumWidth(180)
        self.control_dock.setMaximumWidth(340)
        self.addDockWidget(_dock_area_left(), self.control_dock)

        self.shot_list = QtWidgets.QListWidget(self)
        self.shot_dock = QtWidgets.QDockWidget("Shot Browser", self)
        self.shot_dock.setWidget(self.shot_list)
        self.shot_dock.setFeatures(_dock_features("DockWidgetMovable", "DockWidgetFloatable", "DockWidgetClosable"))
        self.addDockWidget(_dock_area_left(), self.shot_dock)

        self.summary_table = QtWidgets.QTableWidget(self)
        self.summary_table.setColumnCount(2)
        self.summary_table.setHorizontalHeaderLabels(["Header / Summary", "Value"])
        self.summary_dock = QtWidgets.QDockWidget("Shot Summary", self)
        self.summary_dock.setWidget(self.summary_table)
        self.summary_dock.setFeatures(_dock_features("DockWidgetMovable", "DockWidgetFloatable", "DockWidgetClosable"))
        self.addDockWidget(_dock_area_bottom(), self.summary_dock)

        results_container = QtWidgets.QWidget(self)
        results_layout = QtWidgets.QVBoxLayout(results_container)
        results_layout.setContentsMargins(4, 4, 4, 4)
        self.results_table = QtWidgets.QTableWidget(self)
        self.results_table.setColumnCount(2)
        self.results_table.setHorizontalHeaderLabels(["Result", "Value"])
        self.btn_save_excel = QtWidgets.QPushButton("Save Excel Summary")
        self.btn_save_excel.setEnabled(False)
        results_layout.addWidget(self.results_table)
        results_layout.addWidget(self.btn_save_excel)
        self.results_dock = QtWidgets.QDockWidget("Analysis Results", self)
        self.results_dock.setWidget(results_container)
        self.results_dock.setFeatures(_dock_features("DockWidgetMovable", "DockWidgetFloatable", "DockWidgetClosable"))
        self.addDockWidget(_dock_area_right(), self.results_dock)

        self.splitDockWidget(self.control_dock, self.shot_dock, _orient_vertical())
        self.resizeDocks([self.control_dock, self.shot_dock], [230, 330], _orient_vertical())
        self.resizeDocks([self.control_dock], [250], _orient_horizontal())
        self.resizeDocks([self.control_dock, self.shot_dock], [230, 260], _orient_horizontal())

        self.statusBar().showMessage("Open a SEG2 folder to start.")

        self.btn_open.clicked.connect(self.open_folder)
        self.btn_prev.clicked.connect(self.prev_shot)
        self.btn_next.clicked.connect(self.next_shot)
        self.btn_pick_toggle.clicked.connect(self._toggle_pick_mode)
        self.btn_auto_pick.clicked.connect(self._auto_pick_current_shot)
        self.btn_save_picks.clicked.connect(self._save_picks_current_profile)
        self.btn_clear_picks.clicked.connect(self._clear_current_shot_picks)
        self.btn_run_compute.clicked.connect(self._run_full_computation)
        self.btn_run_analysis.clicked.connect(self._run_interactive_analysis)
        self.btn_save_excel.clicked.connect(self._save_excel_summary)
        self.btn_import_layer.clicked.connect(self._import_pick_layer)
        self.btn_set_active_layer.clicked.connect(self._set_active_layer_from_selection)
        self.btn_remove_layer.clicked.connect(self._remove_selected_layer)
        self.layers_list.itemChanged.connect(self._on_layer_item_changed)
        self.shot_list.currentRowChanged.connect(self.set_shot)
        self.spin_clip.valueChanged.connect(self._on_view_settings_changed)
        self.spin_max_tr.valueChanged.connect(self._on_view_settings_changed)
        self.spin_wiggle.valueChanged.connect(self._on_view_settings_changed)
        self.cmb_mode.currentTextChanged.connect(self._on_view_settings_changed)
        self.cmb_x_axis.currentTextChanged.connect(self._on_view_settings_changed)
        self.chk_invert_y.toggled.connect(self._on_view_settings_changed)
        self.chk_grid.toggled.connect(self._on_view_settings_changed)
        self.chk_timelines.toggled.connect(self._on_view_settings_changed)
        self.cmb_theme.currentTextChanged.connect(self._on_view_settings_changed)
        self.spin_tmin.valueChanged.connect(self._on_view_settings_changed)
        self.spin_tmax.valueChanged.connect(self._on_view_settings_changed)
        self.cmb_cmap.currentTextChanged.connect(self._on_view_settings_changed)
        self.cmb_geom_view.currentTextChanged.connect(self._on_geometry_changed)
        self.chk_use_seg2_delay.toggled.connect(self._on_view_settings_changed)
        self.chk_use_seg2_shotpos.toggled.connect(self._on_view_settings_changed)
        self.spin_trigger_static.valueChanged.connect(self._on_view_settings_changed)
        self.spin_extra_delay.valueChanged.connect(self._on_view_settings_changed)
        self.cmb_polarity.currentTextChanged.connect(self._on_view_settings_changed)
        self.cmb_gain_mode.currentTextChanged.connect(self._on_view_settings_changed)
        self.spin_agc_window.valueChanged.connect(self._on_view_settings_changed)
        self.cmb_agc_stat.currentTextChanged.connect(self._on_view_settings_changed)
        self.cmb_filter_mode.currentTextChanged.connect(self._on_view_settings_changed)
        self.spin_butter_order.valueChanged.connect(self._on_view_settings_changed)
        self.spin_f1.valueChanged.connect(self._on_view_settings_changed)
        self.spin_f2.valueChanged.connect(self._on_view_settings_changed)
        self.spin_f3.valueChanged.connect(self._on_view_settings_changed)
        self.spin_f4.valueChanged.connect(self._on_view_settings_changed)
        self.cmb_auto_mode.currentTextChanged.connect(self._on_view_settings_changed)
        self.cmb_hilbert_target.currentTextChanged.connect(self._on_view_settings_changed)
        self.spin_hilbert_onset_pct.valueChanged.connect(self._on_view_settings_changed)
        self.cmb_pick_order.currentTextChanged.connect(self._on_view_settings_changed)
        self.spin_manual_snap.valueChanged.connect(self._on_view_settings_changed)
        self.display.canvas.mpl_connect("button_press_event", self._on_canvas_press)
        self.display.canvas.mpl_connect("motion_notify_event", self._on_canvas_motion)
        self.display.canvas.mpl_connect("button_release_event", self._on_canvas_release)
        self.display.canvas.mpl_connect("scroll_event", self._on_canvas_scroll)

    def _build_actions(self):
        QAction = QtGui.QAction
        menu_file = self.menuBar().addMenu("File")
        act_open = QAction("Open SEG2 Folder", self)
        act_open.triggered.connect(self.open_folder)
        menu_file.addAction(act_open)

        menu_file.addSeparator()
        act_exit = QAction("Exit", self)
        act_exit.setShortcut("Ctrl+Q")
        act_exit.triggered.connect(self.close)
        menu_file.addAction(act_exit)

        menu_view = self.menuBar().addMenu("View")

        act_undock = QAction("Undock Display", self)
        act_undock.triggered.connect(self.undock_display)
        menu_view.addAction(act_undock)

        act_dock = QAction("Dock Display", self)
        act_dock.triggered.connect(self.dock_display)
        menu_view.addAction(act_dock)

        act_zoom_reset = QAction("Reset Zoom", self)
        act_zoom_reset.triggered.connect(self.reset_zoom)
        menu_view.addAction(act_zoom_reset)

        menu_pick = self.menuBar().addMenu("Picker")
        act_toggle = QAction("Enable Picker", self)
        act_toggle.setCheckable(True)
        act_toggle.setChecked(True)
        act_toggle.triggered.connect(lambda checked: self._set_pick_mode(bool(checked)))
        menu_pick.addAction(act_toggle)

        act_clear = QAction("Clear Current Shot Picks", self)
        act_clear.triggered.connect(self._clear_current_shot_picks)
        menu_pick.addAction(act_clear)

        act_save = QAction("Save Picks JSON", self)
        act_save.triggered.connect(self._save_picks_json)
        menu_pick.addAction(act_save)

        menu_compute = self.menuBar().addMenu("Compute")
        act_run_compute = QAction("Run Full Computation", self)
        act_run_compute.triggered.connect(self._run_full_computation)
        menu_compute.addAction(act_run_compute)

        nav_menu = self.menuBar().addMenu("Navigate")
        act_prev = QAction("Previous Shot", self)
        act_prev.setShortcut("Left")
        act_prev.triggered.connect(self.prev_shot)
        nav_menu.addAction(act_prev)

        act_next = QAction("Next Shot", self)
        act_next.setShortcut("Right")
        act_next.triggered.connect(self.next_shot)
        nav_menu.addAction(act_next)

    def _on_view_settings_changed(self):
        self.settings.clip_pct = float(self.spin_clip.value())
        self.settings.max_traces_render = int(self.spin_max_tr.value())
        self.settings.wiggle_scale = float(self.spin_wiggle.value())
        self.settings.display_mode = self.cmb_mode.currentText()
        self.settings.x_axis_mode = self.cmb_x_axis.currentText().strip().lower()
        self.settings.invert_y_axis = bool(self.chk_invert_y.isChecked())
        self.settings.show_grid = bool(self.chk_grid.isChecked())
        self.settings.show_timelines = bool(self.chk_timelines.isChecked())
        self.settings.theme = self.cmb_theme.currentText().strip().lower()
        self.settings.use_seg2_shot_position = bool(self.chk_use_seg2_shotpos.isChecked())
        self.settings.t_min_ms = float(self.spin_tmin.value())
        self.settings.t_max_ms = float(self.spin_tmax.value())
        self.settings.cmap = self.cmb_cmap.currentText()
        self.settings.geometry_override = self.cmb_geom_view.currentText().strip().lower()
        self.settings.use_seg2_delay = bool(self.chk_use_seg2_delay.isChecked())
        self.settings.trigger_static_ms = float(self.spin_trigger_static.value())
        self.settings.extra_delay_ms = float(self.spin_extra_delay.value())
        self.settings.polarity = self.cmb_polarity.currentText().strip().lower()
        self.settings.gain_mode = self.cmb_gain_mode.currentText().strip().lower()
        self.settings.agc_window_ms = float(self.spin_agc_window.value())
        self.settings.agc_stat = self.cmb_agc_stat.currentText().strip().lower()
        self.settings.filter_mode = self.cmb_filter_mode.currentText().strip().lower()
        self.settings.butter_order = int(self.spin_butter_order.value())
        self.settings.f1 = float(self.spin_f1.value())
        self.settings.f2 = float(self.spin_f2.value())
        self.settings.f3 = float(self.spin_f3.value())
        self.settings.f4 = float(self.spin_f4.value())
        self.settings.auto_pick_mode = self.cmb_auto_mode.currentText().strip().lower()
        self.settings.hilbert_target = self.cmb_hilbert_target.currentText().strip().lower()
        self.settings.hilbert_onset_pct = float(self.spin_hilbert_onset_pct.value())
        self.settings.pick_order = self.cmb_pick_order.currentText().strip().upper()
        self.settings.manual_snap_win_ms = float(self.spin_manual_snap.value())
        self._render_current()

    def _backend(self):
        import lvl_refraction as lr
        return lr

    def _on_geometry_changed(self, text: str):
        t = str(text).strip().lower()
        if t == "custom":
            start_dir = str(DATA_DIR if DATA_DIR.exists() else Path.cwd())
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self, "Select geometry file (same format as geometry100.txt)",
                start_dir, "Geometry (*.txt);;All files (*.*)"
            )
            if not path:
                # Revert to the previous selection.
                prev = self.settings.geometry_override if self.settings.geometry_override != "custom" else "auto"
                self.cmb_geom_view.blockSignals(True)
                self.cmb_geom_view.setCurrentText(prev)
                self.cmb_geom_view.blockSignals(False)
                return
            try:
                self._custom_geom_positions = self._parse_geometry_file(Path(path))
                self._custom_geom_path = path
                self.statusBar().showMessage(
                    f"Custom geometry loaded: {Path(path).name} "
                    f"({self._custom_geom_positions.size} pts, span "
                    f"{abs(self._custom_geom_positions[-1] - self._custom_geom_positions[0]):.1f} m)"
                )
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self, "Geometry", f"Could not read geometry file:\n{exc}")
                prev = self.settings.geometry_override if self.settings.geometry_override != "custom" else "auto"
                self.cmb_geom_view.blockSignals(True)
                self.cmb_geom_view.setCurrentText(prev)
                self.cmb_geom_view.blockSignals(False)
                return
        self._on_view_settings_changed()

    def _parse_geometry_file(self, path: Path) -> np.ndarray:
        """Parse a geometry file (same format as geometry100.txt): 'trace  position'."""
        positions: list[float] = []
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
                elif len(parts) == 1:
                    try:
                        positions.append(float(parts[0]))
                    except ValueError:
                        pass
        arr = np.array(positions, dtype=float)
        if arr.size < 2:
            raise ValueError("Geometry file needs at least 2 positions.")
        return arr

    def _nearest_geom_by_span(self, span: float, margin: float | None = None) -> int | None:
        """Return known geometry type whose spread span is closest to `span` (within margin)."""
        lr = self._backend()
        best = None
        best_diff = 1e18
        for gt in sorted(getattr(lr, "GEOM_FILES", {100: None, 200: None}).keys()):
            try:
                r = np.asarray(lr.load_geometry(int(gt)), dtype=float)
            except Exception:
                continue
            gs = abs(float(r[-1] - r[0]))
            d = abs(gs - float(span))
            if d < best_diff:
                best_diff = d
                best = int(gt)
        if margin is not None and best_diff > float(margin):
            return None
        return best

    def _resolve_geom_type(self, shot: ShotGather) -> int | None:
        g = str(self.settings.geometry_override).lower()
        if g in ("100", "200"):
            return int(g)
        if g == "custom":
            if self._custom_geom_positions is not None and self._custom_geom_positions.size >= 2:
                span = abs(float(self._custom_geom_positions[-1] - self._custom_geom_positions[0]))
                return self._nearest_geom_by_span(span)
            return None
        lr = self._backend()
        # 1) LVL coordinate workbook distance -> nearest geometry within ~10 m.
        try:
            profile = self._profile_name_from_open_folder()
            if profile:
                geometry_paths = [str(p) for p in self._discover_geometry_excels()]
                if geometry_paths:
                    res = lr.load_profile_geometry_from_excels(profile_name=profile, excel_paths=geometry_paths)
                    length = res[0] if res else None
                    if length and float(length) > 0.0:
                        gt = self._nearest_geom_by_span(float(length), margin=10.0)
                        if gt in (100, 200):
                            return int(gt)
        except Exception:
            pass
        # 2) SEG2 receiver span inference.
        try:
            geom, _length = lr.infer_geometry_from_seg2_file(shot.source_file)
            if geom in (100, 200):
                return int(geom)
        except Exception:
            pass
        # 3) Trace-count fallback.
        try:
            best = None
            best_diff = 1e18
            for gt in (100, 200):
                sz = int(np.asarray(lr.load_geometry(gt)).size)
                d = abs(sz - shot.n_traces)
                if d < best_diff:
                    best_diff = d
                    best = gt
            return best
        except Exception:
            return None

    def _receiver_positions_for_shot(self, shot: ShotGather) -> np.ndarray:
        g = str(self.settings.geometry_override).lower()
        recv = None
        try:
            if g == "custom" and self._custom_geom_positions is not None:
                recv = np.asarray(self._custom_geom_positions, dtype=float)
            else:
                gt = self._resolve_geom_type(shot)
                if gt in (100, 200):
                    recv = np.asarray(self._backend().load_geometry(int(gt)), dtype=float)
        except Exception:
            recv = None
        if recv is None or recv.size < 2:
            return np.arange(1, shot.n_traces + 1, dtype=float)
        if recv.size >= shot.n_traces:
            return recv[: shot.n_traces]
        dx = float(np.median(np.diff(recv)))
        extra = np.arange(shot.n_traces - recv.size, dtype=float) * dx + recv[-1] + dx
        return np.concatenate([recv, extra])

    def _shot_source_position(self, shot: ShotGather, recv_positions: np.ndarray) -> float:
        src_loc = str((shot.header or {}).get("SOURCE_LOCATION", "")).strip().replace(",", ".")
        if src_loc:
            try:
                return float(src_loc.split()[0])
            except Exception:
                pass
        if recv_positions.size >= 2:
            return float(0.5 * (recv_positions[0] + recv_positions[-1]))
        return 0.0

    def _x_axis_values(self, shot: ShotGather, recv_positions: np.ndarray) -> tuple[np.ndarray, str]:
        mode = str(self.settings.x_axis_mode).lower()
        if mode == "channel":
            return np.arange(1, shot.n_traces + 1, dtype=float), "Channel"
        if mode == "offset":
            sp = self._shot_source_position(shot, recv_positions)
            return np.asarray(recv_positions, dtype=float) - float(sp), "Offset (m)"
        return np.asarray(recv_positions, dtype=float), "Geometry X (m)"

    def _auto_shot_positions_for_all(self, recv_positions: np.ndarray) -> dict[int, float]:
        """Derive per-shot inline positions from geometry, matching the script.

        For 3 shots this uses the backend `auto_shot_positions` exactly
        (shot1 before spread, shot2 spread midpoint, shot3 after spread).
        """
        rp = np.asarray(recv_positions, dtype=float)
        n = len(self.shots)
        out: dict[int, float] = {}
        if n == 0:
            return out
        if rp.size < 2:
            for i in range(n):
                out[i] = float(rp[0]) if rp.size else float(i)
            return out
        if n == 3:
            try:
                auto = self._backend().auto_shot_positions(rp)
                return {0: float(auto[1]), 1: float(auto[2]), 2: float(auto[3])}
            except Exception:
                pass
        dx0 = float(rp[1] - rp[0])
        dxe = float(rp[-1] - rp[-2])
        if n == 1:
            out[0] = float(0.5 * (rp[0] + rp[-1]))
            return out
        for i in range(n):
            if i == 0:
                out[i] = float(rp[0] - dx0)
            elif i == n - 1:
                out[i] = float(rp[-1] + dxe)
            else:
                frac = i / float(n - 1)
                out[i] = float(rp[0] + frac * (rp[-1] - rp[0]))
        return out

    def _shot_geom_position(self, shot: ShotGather, recv_positions: np.ndarray) -> float:
        # Match the script: geometry-derived positions by default; SEG2
        # SOURCE_LOCATION only when explicitly enabled.
        if self.settings.use_seg2_shot_position:
            src = str((shot.header or {}).get("SOURCE_LOCATION", "")).strip().replace(",", ".")
            if src:
                try:
                    return float(src.split()[0])
                except Exception:
                    pass
        positions = self._auto_shot_positions_for_all(recv_positions)
        rp = np.asarray(recv_positions, dtype=float)
        return positions.get(self.current_idx, float(np.mean(rp)) if rp.size else 0.0)

    def _shot_x_in_axis(self, shot: ShotGather, recv_positions: np.ndarray) -> float | None:
        mode = str(self.settings.x_axis_mode).lower()
        sp = self._shot_geom_position(shot, recv_positions)
        if mode == "offset":
            return 0.0
        if mode == "channel":
            rp = np.asarray(recv_positions, dtype=float)
            if rp.size < 2:
                return None
            ch = np.arange(1, rp.size + 1, dtype=float)
            try:
                order = np.argsort(rp)
                return float(np.interp(sp, rp[order], ch[order]))
            except Exception:
                return None
        return float(sp)

    def _processed_trace_matrix(self, shot: ShotGather) -> np.ndarray:
        lr = self._backend()
        data = np.asarray(shot.data, dtype=np.float32)
        filt = str(self.settings.filter_mode).lower()
        if filt == "butter":
            lo = min(self.settings.f2, self.settings.f3 - 0.1)
            hi = max(self.settings.f3, lo + 0.2)
            data = lr.apply_butterworth_all_params(
                data,
                shot.dt_ms / 1000.0,
                low_hz=lo,
                high_hz=hi,
                order=int(self.settings.butter_order),
            )
        elif filt == "ormsby":
            f1, f2, f3, f4 = self.settings.f1, self.settings.f2, self.settings.f3, self.settings.f4
            if f1 < f2 < f3 < f4:
                data = lr.apply_ormsby_all_params(data, shot.dt_ms / 1000.0, f1=f1, f2=f2, f3=f3, f4=f4)
        elif filt == "cutpass":
            low_cut = min(self.settings.f2, self.settings.f3 - 0.1)
            high_cut = max(self.settings.f3, low_cut + 0.2)
            data = lr.apply_cutpass_all_params(
                data,
                shot.dt_ms / 1000.0,
                low_cut_hz=low_cut,
                high_cut_hz=high_cut,
                order=int(self.settings.butter_order),
            )

        data = lr.apply_gain(
            data,
            shot.dt_ms / 1000.0,
            mode=self.settings.gain_mode,
            window_ms=self.settings.agc_window_ms,
            stat=self.settings.agc_stat,
        )
        if self.settings.polarity == "invert":
            data = -data
        return np.asarray(data, dtype=np.float32)

    def _filtered_trace_matrix(self, shot: ShotGather) -> np.ndarray:
        lr = self._backend()
        data = np.asarray(shot.data, dtype=np.float32)
        filt = str(self.settings.filter_mode).lower()
        if filt == "butter":
            lo = min(self.settings.f2, self.settings.f3 - 0.1)
            hi = max(self.settings.f3, lo + 0.2)
            data = lr.apply_butterworth_all_params(
                data,
                shot.dt_ms / 1000.0,
                low_hz=lo,
                high_hz=hi,
                order=int(self.settings.butter_order),
            )
        elif filt == "ormsby":
            f1, f2, f3, f4 = self.settings.f1, self.settings.f2, self.settings.f3, self.settings.f4
            if f1 < f2 < f3 < f4:
                data = lr.apply_ormsby_all_params(data, shot.dt_ms / 1000.0, f1=f1, f2=f2, f3=f3, f4=f4)
        elif filt == "cutpass":
            low_cut = min(self.settings.f2, self.settings.f3 - 0.1)
            high_cut = max(self.settings.f3, low_cut + 0.2)
            data = lr.apply_cutpass_all_params(
                data,
                shot.dt_ms / 1000.0,
                low_cut_hz=low_cut,
                high_cut_hz=high_cut,
                order=int(self.settings.butter_order),
            )
        if self.settings.polarity == "invert":
            data = -data
        return np.asarray(data, dtype=np.float32)

    def _gain_trace_matrix(self, shot: ShotGather) -> np.ndarray:
        lr = self._backend()
        data = self._filtered_trace_matrix(shot)
        data = lr.apply_gain(
            data,
            shot.dt_ms / 1000.0,
            mode=self.settings.gain_mode,
            window_ms=self.settings.agc_window_ms,
            stat=self.settings.agc_stat,
        )
        return np.asarray(data, dtype=np.float32)

    def _trace_for_manual_snap(self, shot: ShotGather, trace_idx: int) -> np.ndarray:
        order = str(self.settings.pick_order).upper()
        raw = np.asarray(shot.data[trace_idx], dtype=np.float32)
        fil = self._filtered_trace_matrix(shot)[trace_idx]
        gain = self._gain_trace_matrix(shot)[trace_idx]
        if order == "RAW>X":
            return raw
        if order == "G>F>X":
            return gain
        return fil

    def _snap_pick_time_ms(self, shot: ShotGather, trace_idx: int, click_t_ms: float) -> float:
        lr = self._backend()
        mode = str(self.settings.auto_pick_mode).lower()
        dt_s = shot.dt_ms / 1000.0
        t0 = self._time_zero_ms(shot)
        # STA/LTA: honor cursor position exactly (no snap).
        if mode == "stalta":
            return float(click_t_ms)

        tr = self._trace_for_manual_snap(shot, trace_idx)
        t_rel = (float(click_t_ms) - t0) / 1000.0
        win = float(self.settings.manual_snap_win_ms)
        back_s = max(dt_s, 1.25 * win / 1000.0)
        fwd_s = max(dt_s, 0.25 * win / 1000.0)
        a = max(0.0, t_rel - back_s)
        b = min((len(tr) - 1) * dt_s, t_rel + fwd_s)
        if b <= a:
            return float(click_t_ms)

        if mode == "maxdiff_zero":
            ts = lr._zero_crossing_from_extremum_samples(
                samples=tr,
                dt_s=dt_s,
                win_start_s=a,
                win_end_s=b,
                search_direction="backward",
                use_abs_peak=True,
            )
            if ts is None:
                return float(click_t_ms)
            return float(t0 + 1000.0 * ts)

        tro = lr.Trace(data=np.asarray(tr, dtype=np.float32))
        tro.stats.delta = float(dt_s)
        _env, t_peak, t_onset = lr.hilbert_envelope_pick(
            trace=tro,
            win_start_s=a,
            win_end_s=b,
            onset_pct=float(self.settings.hilbert_onset_pct),
        )
        if self.settings.hilbert_target == "peak":
            ts = t_peak if t_peak is not None else t_onset
        else:
            ts = t_onset if t_onset is not None else t_peak
        if ts is None:
            return float(click_t_ms)
        return float(t0 + 1000.0 * float(ts))

    def _time_zero_ms(self, shot: ShotGather) -> float:
        t0 = 0.0
        if self.settings.use_seg2_delay:
            t0 += float(shot.delay_ms)
        t0 += float(self.settings.trigger_static_ms)
        t0 += float(self.settings.extra_delay_ms)
        return t0

    def _auto_gate_local_s(self, i: int, rp: np.ndarray, shot_pos: float,
                           t0: float, dt_s: float, n_samp: int) -> tuple[float, float]:
        """Velocity-gated search window (seconds from sample 0), like the script."""
        vmin, vmax, pad = 120.0, 3500.0, 20.0
        t_lo_view = float(self.settings.t_min_ms)
        t_hi_view = float(self.settings.t_max_ms)
        off = abs(float(rp[i]) - float(shot_pos)) if i < len(rp) else 0.0
        lo = max(t_lo_view, (off / vmax) * 1000.0 - pad)
        hi = min(t_hi_view, (off / vmin) * 1000.0 + pad)
        if hi <= lo:
            lo, hi = t_lo_view, t_hi_view
        a = max(0.0, (lo - t0) / 1000.0)
        b = min((n_samp - 1) * dt_s, (hi - t0) / 1000.0)
        if b <= a:
            a, b = 0.0, (n_samp - 1) * dt_s
        return a, b

    def _auto_pick_current_shot(self):
        if self.current_idx < 0 or self.current_idx >= len(self.shots):
            return
        shot = self.shots[self.current_idx]
        lr = self._backend()
        prepared = self._gain_trace_matrix(shot)
        dt_s = shot.dt_ms / 1000.0
        n_samp = prepared.shape[1]
        t0 = self._time_zero_ms(shot)
        mode = str(self.settings.auto_pick_mode).lower()
        picks = self._shot_pick_map(self.current_idx)
        rp = self._receiver_positions_for_shot(shot)
        shot_pos = self._shot_geom_position(shot, rp)
        count = 0

        if mode == "stalta":
            n_sta = max(1, int(0.003 / dt_s))
            n_lta = max(3, int(0.020 / dt_s))
            trig = 3.0
            for i in range(prepared.shape[0]):
                a, b = self._auto_gate_local_s(i, rp, shot_pos, t0, dt_s, n_samp)
                lo_k = max(n_lta, int(a / dt_s))
                hi_k = min(n_samp - n_sta - 1, int(b / dt_s))
                tr = prepared[i].astype(float)
                tr = tr - np.mean(tr[: max(1, n_lta)])
                char = np.maximum(tr, 0.0)
                found = None
                for k in range(lo_k, max(lo_k + 1, hi_k)):
                    lta = char[k - n_lta:k].mean()
                    if lta <= 1e-30:
                        continue
                    sta = char[k:k + n_sta].mean()
                    if sta / lta >= trig:
                        found = k
                        break
                if found is not None:
                    picks[i] = round(t0 + found * shot.dt_ms, 2)
                    count += 1
        elif mode == "maxdiff_zero":
            for i in range(prepared.shape[0]):
                a, b = self._auto_gate_local_s(i, rp, shot_pos, t0, dt_s, n_samp)
                zt = lr._zero_crossing_from_extremum_samples(
                    samples=prepared[i].astype(float),
                    dt_s=dt_s,
                    win_start_s=a,
                    win_end_s=b,
                    search_direction="backward",
                    use_abs_peak=True,
                )
                if zt is not None:
                    picks[i] = round(t0 + zt * 1000.0, 2)
                    count += 1
        else:
            for i in range(prepared.shape[0]):
                a, b = self._auto_gate_local_s(i, rp, shot_pos, t0, dt_s, n_samp)
                tr = lr.Trace(data=np.asarray(prepared[i], dtype=np.float32))
                tr.stats.delta = float(dt_s)
                _env, peak_s, onset_s = lr.hilbert_envelope_pick(
                    trace=tr,
                    win_start_s=a,
                    win_end_s=b,
                    onset_pct=float(self.settings.hilbert_onset_pct),
                )
                if self.settings.hilbert_target == "peak":
                    ts = peak_s if peak_s is not None else onset_s
                else:
                    ts = onset_s if onset_s is not None else peak_s
                if ts is not None:
                    picks[i] = round(t0 + float(ts) * 1000.0, 2)
                    count += 1

        self.statusBar().showMessage(f"Auto-pick ({mode}): {count}/{shot.n_traces} placed")
        self._render_current()

    def _set_pick_mode(self, enabled: bool):
        self.pick_enabled = bool(enabled)
        self.btn_pick_toggle.setChecked(bool(enabled))
        self.btn_pick_toggle.setText("Picker: ON" if enabled else "Picker: OFF")

    def _toggle_pick_mode(self):
        self._set_pick_mode(bool(self.btn_pick_toggle.isChecked()))

    def _shot_pick_map(self, idx: int) -> dict[int, float]:
        return self.picks_by_shot.setdefault(int(idx), {})

    # ------------------------------------------------------------------
    # Multi-layer pick management (own picks + imported reference picks)
    # ------------------------------------------------------------------

    def _next_layer_color(self) -> str:
        used = set(self.layer_colors.values())
        for c in self._layer_color_cycle:
            if c not in used:
                return c
        return self._layer_color_cycle[len(self.layer_colors) % len(self._layer_color_cycle)]

    def _refresh_layers_list(self):
        self.layers_list.blockSignals(True)
        self.layers_list.clear()
        for name in self.pick_layers:
            item = QtWidgets.QListWidgetItem(name)
            color = self.layer_colors.get(name, "#888888")
            pix = QtGui.QPixmap(12, 12)
            pix.fill(QtGui.QColor(color))
            item.setIcon(QtGui.QIcon(pix))
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable
                          if hasattr(QtCore.Qt, "ItemFlag") else item.flags() | QtCore.Qt.ItemIsUserCheckable)
            checked = self.layer_visible.get(name, True)
            check_state = getattr(QtCore.Qt, "CheckState", None)
            item.setCheckState(
                (check_state.Checked if checked else check_state.Unchecked)
                if check_state is not None
                else (QtCore.Qt.Checked if checked else QtCore.Qt.Unchecked)
            )
            if name == self.active_layer:
                font = item.font()
                font.setBold(True)
                item.setFont(font)
                item.setText(f"{name}  [active]")
            self.layers_list.addItem(item)
        self.layers_list.blockSignals(False)
        self.lbl_active_layer.setText(f"Active layer: {self.active_layer}")

    def _selected_layer_name(self) -> str | None:
        item = self.layers_list.currentItem()
        if item is None:
            return None
        return str(item.text()).replace("  [active]", "").strip()

    def _on_layer_item_changed(self, item: Any):
        name = str(item.text()).replace("  [active]", "").strip()
        check_state = getattr(QtCore.Qt, "CheckState", None)
        checked_val = check_state.Checked if check_state is not None else QtCore.Qt.Checked
        self.layer_visible[name] = (item.checkState() == checked_val)
        self._render_current()

    def _set_active_layer_from_selection(self):
        name = self._selected_layer_name()
        if not name or name not in self.pick_layers:
            QtWidgets.QMessageBox.information(self, "Pick Layers", "Select a layer in the list first.")
            return
        self.active_layer = name
        self.layer_visible[name] = True
        self._refresh_layers_list()
        self._render_current()
        self.statusBar().showMessage(f"Active layer set to '{name}'. Editing/auto-pick/analysis now use this layer.")

    def _remove_selected_layer(self):
        name = self._selected_layer_name()
        if not name or name not in self.pick_layers:
            return
        if name == "Mine":
            QtWidgets.QMessageBox.information(self, "Pick Layers", "The default 'Mine' layer cannot be removed.")
            return
        del self.pick_layers[name]
        self.layer_colors.pop(name, None)
        self.layer_visible.pop(name, None)
        if self.active_layer == name:
            self.active_layer = "Mine"
        self._refresh_layers_list()
        self._render_current()

    def _import_pick_layer(self):
        if not self.shots:
            QtWidgets.QMessageBox.information(self, "Import Picks", "Open a SEG2 folder first.")
            return
        start_dir = str(self._current_folder or DATA_DIR if DATA_DIR.exists() else Path.cwd())
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Import picks file", start_dir,
            "Picks (*.txt *.json);;Text picks (*.txt);;JSON picks (*.json);;All files (*.*)",
        )
        if not path:
            return
        try:
            p = Path(path)
            if p.suffix.lower() == ".json":
                layer_data = self._parse_json_picks_file(p)
            else:
                layer_data = self._parse_txt_picks_file(p)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Import Picks Error", str(exc))
            return

        n_shots_with_picks = sum(1 for v in layer_data.values() if v)
        if n_shots_with_picks == 0:
            QtWidgets.QMessageBox.warning(self, "Import Picks", "No usable picks found in this file.")
            return

        base_name = Path(path).stem
        name = base_name
        suffix = 2
        while name in self.pick_layers:
            name = f"{base_name} ({suffix})"
            suffix += 1

        self.pick_layers[name] = layer_data
        self.layer_colors[name] = self._next_layer_color()
        self.layer_visible[name] = True
        self._refresh_layers_list()
        self._render_current()
        self.statusBar().showMessage(
            f"Imported layer '{name}' from {Path(path).name} ({n_shots_with_picks} shot(s) with picks). "
            "Use 'Set Active Layer' to edit/use it for computation."
        )

    def _parse_txt_picks_file(self, path: Path) -> dict[int, dict[int, float]]:
        """Parse Ensemble/#/SOURCE/CHAN/OFFSET/FB_PICK text (as exported by lvl_refraction.py)."""
        out: dict[int, dict[int, float]] = {}
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                s = line.strip()
                if not s:
                    continue
                parts = s.replace(",", ".").split()
                if not parts or not parts[0].lstrip("-").isdigit():
                    continue  # header or non-data line
                if len(parts) < 6:
                    continue
                try:
                    source = int(float(parts[2]))
                    chan = int(float(parts[3]))
                    fb_pick = float(parts[5])
                except Exception:
                    continue
                shot_idx = source - 1
                trace_idx = chan - 1
                out.setdefault(shot_idx, {})[trace_idx] = fb_pick
        return out

    def _parse_json_picks_file(self, path: Path) -> dict[int, dict[int, float]]:
        """Parse picks.json/picks.session.json-style ({shot_id: {trace_idx: ms}})
        or the studio's own Save Picks JSON export format."""
        with open(path, encoding="utf-8") as fh:
            raw = json.load(fh)
        out: dict[int, dict[int, float]] = {}
        if isinstance(raw, dict) and "shots" in raw and isinstance(raw["shots"], list):
            for entry in raw["shots"]:
                idx = int(entry.get("index", 0))
                tr_map = entry.get("picks_ms_by_trace_index", {}) or {}
                d: dict[int, float] = {}
                for k, v in tr_map.items():
                    try:
                        d[int(k) - 1] = float(v)
                    except Exception:
                        continue
                if d:
                    out[idx] = d
            return out
        for sid_str, tr_map in (raw or {}).items():
            try:
                idx = int(sid_str) - 1
            except Exception:
                continue
            d = {}
            for ti_str, tv in (tr_map or {}).items():
                try:
                    d[int(ti_str)] = float(tv)
                except Exception:
                    continue
            if d:
                out[idx] = d
        return out

    def _pick_trace_from_event(self, event: Any) -> tuple[int, float] | None:
        if not self.pick_enabled:
            return None
        if self.current_idx < 0 or self.current_idx >= len(self.shots):
            return None
        if event.inaxes != self.display.ax:
            return None
        if event.xdata is None or event.ydata is None:
            return None

        if not self.display.last_view_x:
            return None
        x_evt = float(event.xdata)
        view_idx = int(np.argmin(np.abs(np.asarray(self.display.last_view_x, dtype=float) - x_evt)))
        if view_idx < 0 or view_idx >= len(self.display.last_view_to_trace):
            return None
        trace_idx = int(self.display.last_view_to_trace[view_idx])
        t_ms = float(event.ydata)
        return trace_idx, t_ms

    def _on_canvas_press(self, event: Any):
        pick = self._pick_trace_from_event(event)
        if pick is None:
            return
        trace_idx, t_ms = pick

        shot_picks = self._shot_pick_map(self.current_idx)
        if event.button == MouseButton.LEFT:
            # Shift+Left = anchor-based range fill between two traces.
            if str(getattr(event, "key", "") or "").lower() == "shift":
                self._range_fill(trace_idx, t_ms)
                self._render_current()
                return
            shot = self.shots[self.current_idx]
            t_snap = self._snap_pick_time_ms(shot, trace_idx, t_ms)
            shot_picks[trace_idx] = round(t_snap, 2)
            self._drag_pick = True
            self._last_drag_trace = trace_idx
            self._last_drag_t = float(t_snap)
        elif event.button == MouseButton.RIGHT:
            # Right-click deletes; right-drag deletes a range while moving.
            shot_picks.pop(trace_idx, None)
            self._drag_delete = True
            self._last_drag_trace = trace_idx
            self._last_drag_t = None
        else:
            return
        self._render_current()

    def _range_fill(self, trace_idx: int, t_ms: float):
        shot_picks = self._shot_pick_map(self.current_idx)
        t = round(float(t_ms), 2)
        if self._range_anchor is None:
            self._range_anchor = (int(trace_idx), t)
            self.statusBar().showMessage(f"Range anchor set at trace {trace_idx + 1}, t={t:.2f} ms")
            return
        i0, t0 = self._range_anchor
        i1, t1 = int(trace_idx), t
        if i0 == i1:
            shot_picks[i0] = t1
            self._range_anchor = None
            return
        lo, hi = (i0, i1) if i0 < i1 else (i1, i0)
        for ii in range(lo, hi + 1):
            frac = (ii - i0) / float(i1 - i0)
            shot_picks[ii] = round(float(t0 + frac * (t1 - t0)), 2)
        self._range_anchor = None
        self.statusBar().showMessage(f"Range fill: traces {lo + 1}-{hi + 1}")

    def _on_canvas_motion(self, event: Any):
        # Live hover readout (trace / time under cursor).
        hov = self._pick_trace_from_event(event)
        if hov is not None and not (self._drag_pick or self._drag_delete):
            hv_idx, hv_t = hov
            self.statusBar().showMessage(
                f"Trace {hv_idx + 1} | t={hv_t:.2f} ms | picks {len(self._shot_pick_map(self.current_idx))}"
            )

        if self._drag_delete:
            pick = self._pick_trace_from_event(event)
            if pick is None:
                return
            trace_idx, _t = pick
            shot_picks = self._shot_pick_map(self.current_idx)
            if self._last_drag_trace is None:
                shot_picks.pop(trace_idx, None)
                self._last_drag_trace = trace_idx
                self._render_current()
                return
            prev = int(self._last_drag_trace)
            if trace_idx == prev:
                return
            lo, hi = (prev, trace_idx) if prev < trace_idx else (trace_idx, prev)
            for ii in range(lo, hi + 1):
                shot_picks.pop(ii, None)
            self._last_drag_trace = trace_idx
            self._render_current()
            return

        if not self._drag_pick:
            return
        pick = self._pick_trace_from_event(event)
        if pick is None:
            return
        trace_idx, t_ms = pick
        shot = self.shots[self.current_idx]
        t_ms = self._snap_pick_time_ms(shot, trace_idx, t_ms)
        shot_picks = self._shot_pick_map(self.current_idx)

        if self._last_drag_trace is None or self._last_drag_t is None:
            shot_picks[trace_idx] = round(t_ms, 2)
            self._last_drag_trace = trace_idx
            self._last_drag_t = float(t_ms)
            self._render_current()
            return

        i0 = int(self._last_drag_trace)
        t0 = float(self._last_drag_t)
        i1 = int(trace_idx)
        t1 = float(t_ms)

        if i0 == i1:
            shot_picks[i1] = round(t1, 2)
        else:
            lo, hi = (i0, i1) if i0 < i1 else (i1, i0)
            for ii in range(lo, hi + 1):
                frac = (ii - i0) / float(i1 - i0)
                shot_picks[ii] = round(float(t0 + frac * (t1 - t0)), 2)

        self._last_drag_trace = i1
        self._last_drag_t = t1
        self._render_current()

    def _on_canvas_release(self, _event: Any):
        self._drag_pick = False
        self._drag_delete = False
        self._last_drag_trace = None
        self._last_drag_t = None

    def keyPressEvent(self, event: Any):
        key = str(getattr(event, "text", lambda: "")()).lower()
        if self._handle_shortcut_key(key):
            return
        super().keyPressEvent(event)

    def _build_shortcuts(self):
        QShortcut = getattr(QtGui, "QShortcut", None) or getattr(QtWidgets, "QShortcut", None)
        if QShortcut is None:
            return
        ctx = None
        sc_ctx = getattr(QtCore.Qt, "ShortcutContext", None)
        if sc_ctx is not None:
            ctx = sc_ctx.ApplicationShortcut
        else:
            ctx = getattr(QtCore.Qt, "ApplicationShortcut", None)

        keys = ["n", "p", "a", "s", "k", "v", "g", "f", "l", "c", "t"]
        self._shortcuts = []
        for key in keys:
            sc = QShortcut(QtGui.QKeySequence(key), self)
            if ctx is not None:
                try:
                    sc.setContext(ctx)
                except Exception:
                    pass
            sc.activated.connect(lambda k=key: self._handle_shortcut_key(k))
            self._shortcuts.append(sc)

    def _handle_shortcut_key(self, key: str) -> bool:
        key = str(key).lower()
        if key == "n":
            self.next_shot()
        elif key == "p":
            self.prev_shot()
        elif key == "a":
            self._auto_pick_current_shot()
        elif key in ("s", "k"):
            self._save_picks_current_profile()
        elif key == "v":
            self.cmb_polarity.setCurrentText("normal" if self.settings.polarity == "invert" else "invert")
        elif key == "g":
            cur = self.cmb_gain_mode.currentText().strip().lower()
            self.cmb_gain_mode.setCurrentText({"none": "norm", "norm": "agc", "agc": "none"}.get(cur, "norm"))
        elif key == "f":
            cur = self.cmb_filter_mode.currentText().strip().lower()
            self.cmb_filter_mode.setCurrentText(
                {"none": "butter", "butter": "ormsby", "ormsby": "cutpass", "cutpass": "none"}.get(cur, "butter")
            )
        elif key == "l":
            self.chk_grid.setChecked(not self.chk_grid.isChecked())
        elif key == "c":
            cur = self.cmb_x_axis.currentText().strip().lower()
            self.cmb_x_axis.setCurrentText({"channel": "geom_x", "geom_x": "offset", "offset": "channel"}.get(cur, "geom_x"))
        elif key == "t":
            self.cmb_theme.setCurrentText("dark" if self.settings.theme == "light" else "light")
        else:
            return False
        return True

    def _save_picks_json(self):
        if not self.shots:
            return
        dst, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save picks JSON", "picks_studio.json", "JSON (*.json)")
        if not dst:
            return
        payload: dict[str, Any] = {
            "shots": [],
        }
        for idx, shot in enumerate(self.shots):
            picks = self.picks_by_shot.get(idx, {})
            payload["shots"].append(
                {
                    "index": idx,
                    "ffid": int(shot.ffid),
                    "label": shot.shot_label,
                    "source_file": str(shot.source_file),
                    "picks_ms_by_trace_index": {str(k + 1): float(v) for k, v in sorted(picks.items())},
                }
            )
        with open(dst, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)
        self.statusBar().showMessage(f"Picks saved: {dst}")

    def _clear_current_shot_picks(self):
        if self.current_idx < 0:
            return
        self.picks_by_shot[self.current_idx] = {}
        self._render_current()

    def undock_display(self):
        if self._display_window is not None:
            self._display_window.raise_()
            return
        w = QtWidgets.QMainWindow(self)
        w.setWindowTitle("LVL Studio Display")
        w.setCentralWidget(self.display)
        w.resize(1400, 900)
        w.showMaximized()
        self._display_window = w

    def dock_display(self):
        if self._display_window is None:
            return
        disp = self._display_window.takeCentralWidget()
        if disp is not None:
            self.setCentralWidget(disp)
        self._display_window.close()
        self._display_window = None

    def reset_zoom(self):
        self._view_xlim = None
        self._view_ylim = None
        self._render_current()

    def _on_canvas_scroll(self, event: Any):
        # Cursor-centered scroll zoom that does NOT block left-click picking
        # (promax/geotomo-style). Zoom persists across re-renders.
        ax = self.display.ax
        if event.inaxes != ax or event.xdata is None or event.ydata is None:
            return
        base = 1.2
        scale = (1.0 / base) if str(event.button) == "up" else base
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        xd, yd = float(event.xdata), float(event.ydata)
        nx0 = xd - (xd - x0) * scale
        nx1 = xd + (x1 - xd) * scale
        ny0 = yd - (yd - y0) * scale
        ny1 = yd + (y1 - yd) * scale
        ax.set_xlim(nx0, nx1)
        ax.set_ylim(ny0, ny1)
        self._view_xlim = (nx0, nx1)
        self._view_ylim = (ny0, ny1)
        self.display.canvas.draw_idle()

    def open_folder(self):
        start_dir = str(DATA_DIR if DATA_DIR.exists() else Path.cwd())
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Select SEG2 folder", start_dir)
        if not folder:
            return
        self._current_folder = Path(folder)
        try:
            self.shots = self._load_seg2_folder(Path(folder))
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Load Error", str(exc))
            return

        self._reset_pick_layers()
        self.shot_list.clear()
        for s in self.shots:
            self.shot_list.addItem(
                f"{s.shot_label} | FFID {s.ffid} | traces {s.n_traces} | samples {s.n_samples}"
            )

        loaded_msg = ""
        profile_name = self._profile_name_from_open_folder()
        if profile_name:
            n_loaded = self._autoload_picks(profile_name)
            if n_loaded:
                loaded_msg = f" | resumed {n_loaded} shot picks from output/{profile_name}"

        if self.shots:
            self.shot_list.setCurrentRow(0)
            self.statusBar().showMessage(f"Loaded {len(self.shots)} shots from {folder}{loaded_msg}")
        else:
            self.current_idx = -1
            self.summary_table.setRowCount(0)
            self.statusBar().showMessage("No SEG2 files found.")

    def _autoload_picks(self, profile_name: str) -> int:
        """Load existing picks (session first, then finalized) into the studio."""
        base = OUTPUT_DIR / profile_name
        for name in ("picks.session.json", "picks.json"):
            p = base / name
            if not p.exists():
                continue
            try:
                with open(p, encoding="utf-8") as fh:
                    raw = json.load(fh)
            except Exception:
                continue
            count = 0
            for sid_str, tr_map in (raw or {}).items():
                try:
                    sid = int(sid_str)
                except Exception:
                    continue
                idx = sid - 1
                if idx < 0 or idx >= len(self.shots):
                    continue
                shot_picks = self._shot_pick_map(idx)
                for ti_str, tv in (tr_map or {}).items():
                    try:
                        shot_picks[int(ti_str)] = float(tv)
                    except Exception:
                        continue
                if shot_picks:
                    count += 1
            if count:
                return count
        return 0

    def _save_picks_current_profile(self):
        profile_name = self._profile_name_from_open_folder()
        if not profile_name or not self.shots:
            # Fall back to explicit JSON export when not in a data/<profile> folder.
            self._save_picks_json()
            return
        base = OUTPUT_DIR / profile_name
        base.mkdir(parents=True, exist_ok=True)
        payload: dict[str, dict[str, float]] = {}
        for idx in range(len(self.shots)):
            picks = self.picks_by_shot.get(idx, {})
            if not picks:
                continue
            payload[str(idx + 1)] = {str(int(k)): float(v) for k, v in sorted(picks.items())}
        for name in ("picks.json", "picks.session.json"):
            with open(base / name, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2)
        self.statusBar().showMessage(f"Picks saved to output/{profile_name}/picks.json")

    def _profile_name_from_open_folder(self) -> str | None:
        if self._current_folder is None:
            return None
        folder = self._current_folder.resolve()
        try:
            if folder.parent.resolve() == DATA_DIR.resolve():
                return folder.name
        except Exception:
            return None
        return None

    def _export_picks_to_lvl_refraction_session(self, profile_name: str):
        from lvl_refraction import save_session_picks_json

        all_picks: dict[int, dict[int, float]] = {}
        for idx in range(len(self.shots)):
            shot_id = idx + 1
            picks = self.picks_by_shot.get(idx, {})
            all_picks[shot_id] = {int(k): float(v) for k, v in picks.items()}
        save_session_picks_json(profile_name, all_picks)

    def _run_full_computation(self):
        profile_name = self._profile_name_from_open_folder()
        if not profile_name:
            QtWidgets.QMessageBox.warning(
                self,
                "Computation Setup",
                "Open a profile folder under lvl/data/<profile> to run full computation with lvl_refraction backend.",
            )
            return
        if not self.shots:
            QtWidgets.QMessageBox.warning(self, "Computation Setup", "No shots loaded.")
            return

        geom_txt = self.cmb_geom_view.currentText().strip().lower()
        geom_override = 100 if geom_txt == "100" else 200 if geom_txt == "200" else None
        enable_layer_pick = bool(self.chk_layer_pick.isChecked())

        try:
            self._export_picks_to_lvl_refraction_session(profile_name)

            from lvl_refraction import process_profile, discover_field_report_excels

            report_paths = [str(p) for p in discover_field_report_excels(DATA_DIR)]
            geometry_paths: list[str] = []
            for p in DATA_DIR.glob("LVL*.xls*"):
                name = p.name.lower()
                if "field_report" in name or "fieldreport" in name:
                    continue
                geometry_paths.append(str(p))
            perp_excel_cfg = {
                "paths": report_paths,
                "sheet": None,
                "ffid_col": "A",
                "perp_col": "D",
                "profile_col": "F",
                "inline_shift_col": None,
                "geometry_paths": geometry_paths,
                "device_type": "sw_maps",
            }

            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor if hasattr(QtCore.Qt, "CursorShape") else QtCore.Qt.WaitCursor)
            self.statusBar().showMessage(f"Running full computation for profile {profile_name}...")

            process_profile(
                profile_name=profile_name,
                pick_mode=False,
                geom_override=geom_override,
                perp_excel_cfg=perp_excel_cfg,
                enable_layer_pick=enable_layer_pick,
                control_file=None,
                show_plot_controls=True,
            )

            self.statusBar().showMessage(f"Computation finished for {profile_name}. See lvl/output/{profile_name}")
            QtWidgets.QMessageBox.information(
                self,
                "Computation Finished",
                f"Full computation finished for profile '{profile_name}'.\nOutputs are in lvl/output/{profile_name}.",
            )
        except Exception as exc:
            tb = traceback.format_exc(limit=8)
            self.statusBar().showMessage(f"Computation failed: {exc}")
            QtWidgets.QMessageBox.critical(self, "Computation Error", f"{exc}\n\n{tb}")
        finally:
            try:
                QtWidgets.QApplication.restoreOverrideCursor()
            except Exception:
                pass

    def _derive_shot_positions(self, recv_positions: np.ndarray) -> list[tuple[int, float]]:
        # Geometry-derived positions (script parity). This keeps the middle
        # shot at the spread midpoint so its picks split cleanly into L/R.
        auto = self._auto_shot_positions_for_all(recv_positions)
        rp = np.asarray(recv_positions, dtype=float)
        out: list[tuple[int, float]] = []
        for idx in range(len(self.shots)):
            sid = idx + 1
            shot = self.shots[idx]
            shot_pos = None
            if self.settings.use_seg2_shot_position:
                src_loc = str(shot.header.get("SOURCE_LOCATION", "")).strip().replace(",", ".")
                if src_loc:
                    try:
                        shot_pos = float(src_loc.split()[0])
                    except Exception:
                        shot_pos = None
            if shot_pos is None:
                shot_pos = auto.get(idx)
            if shot_pos is None:
                shot_pos = float(np.mean(rp)) if rp.size else float(idx + 1)
            out.append((sid, float(shot_pos)))
        return out

    def _discover_report_paths(self) -> list[Path]:
        """Field-report Excel files, searching data/ and data/<profile> like the script."""
        lr = self._backend()
        out: list[Path] = []
        seen: set[str] = set()
        search_dirs = [DATA_DIR]
        if self._current_folder is not None:
            search_dirs.append(self._current_folder)
        for d in search_dirs:
            try:
                for p in lr.discover_field_report_excels(d):
                    key = os.path.normcase(str(Path(p).resolve()))
                    if key not in seen:
                        seen.add(key)
                        out.append(Path(p))
            except Exception:
                continue
        return out

    def _load_perp_overrides(self, profile_name: str) -> tuple[dict, dict, list[str]]:
        """Load per-shot PO (and inline shift) from field reports; returns (perp, shift, sources)."""
        lr = self._backend()
        ffid_by_shot = {idx + 1: int(self.shots[idx].ffid) for idx in range(len(self.shots))}
        perp: dict[int, float] = {}
        shift: dict[int, float] = {}
        sources: list[str] = []
        for rp in self._discover_report_paths():
            try:
                po_map, sh_map = lr.load_profile_offsets_from_excel(
                    path=rp,
                    profile_name=profile_name,
                    ffid_by_shot=ffid_by_shot,
                )
            except Exception as exc:
                self.statusBar().showMessage(f"[WARN] PO parse failed in {rp.name}: {exc}")
                continue
            if po_map:
                perp.update(po_map)
                sources.append(rp.name)
            if sh_map:
                shift.update(sh_map)
        return perp, shift, sources

    def _run_interactive_analysis(self):
        if not self.shots:
            QtWidgets.QMessageBox.warning(self, "Analysis", "No shots loaded.")
            return
        if not any(self.picks_by_shot.get(i) for i in range(len(self.shots))):
            QtWidgets.QMessageBox.warning(self, "Analysis", "No picks available. Pick arrivals first.")
            return

        lr = self._backend()
        recv_positions = self._receiver_positions_for_shot(self.shots[0])
        shots_info = self._derive_shot_positions(recv_positions)
        shot_label_pos = {sid: pos for sid, pos in shots_info}
        all_picks = {idx + 1: {int(k): float(v) for k, v in (self.picks_by_shot.get(idx, {}) or {}).items()} for idx in range(len(self.shots))}
        profile = self._profile_name_from_open_folder() or "studio"
        geom_type = self._resolve_geom_type(self.shots[0]) or 200
        cfg = {"perp_m": 0.0, "line_no": profile, "geom": int(geom_type)}

        perp_override, shift_override, po_sources = self._load_perp_overrides(profile)
        if po_sources:
            self.statusBar().showMessage(
                "PO from field report(s): " + ", ".join(po_sources)
                + " | " + ", ".join(f"S{k}={v:.2f}m" for k, v in sorted(perp_override.items()))
            )
        else:
            self.statusBar().showMessage("No field-report PO found; using PO=0.0 (editable in the PO window).")

        try:
            wf = lr.AnalysisWorkflow(
                profile_name=profile,
                cfg=cfg,
                shots_info=shots_info,
                shot_label_pos=shot_label_pos,
                all_picks=all_picks,
                recv_positions=recv_positions,
                perp_override=perp_override or None,
                inline_shift_override=shift_override or None,
                enable_layer_pick=True,
                existing_layer_results={},
            )
            out = wf.run()
            layer_results = out.get("layer_results", {})
            corrected_by_shot = out.get("corrected_by_shot", {})
            avg = lr.compute_layer_averages(layer_results, shot_label_pos)

            self._last_analysis = {
                "profile": profile,
                "cfg": cfg,
                "recv_positions": recv_positions,
                "shots_info": shots_info,
                "layer_results": layer_results,
                "corrected_by_shot": corrected_by_shot,
                "avg": avg,
                "perp_override": perp_override,
            }
            self._populate_results(avg, perp_override, po_sources)
            self.btn_save_excel.setEnabled(True)
            v0 = float((avg.get("V0", {}) or {}).get("avg", 0.0) or 0.0)
            v1 = float((avg.get("V1", {}) or {}).get("avg", 0.0) or 0.0)
            v2 = float((avg.get("V2", {}) or {}).get("avg", 0.0) or 0.0)
            self.statusBar().showMessage(f"Analysis complete. V0={v0:.1f} V1={v1:.1f} V2={v2:.1f} m/s. Results shown at right.")
        except Exception as exc:
            tb = traceback.format_exc(limit=8)
            self.statusBar().showMessage(f"Interactive analysis failed: {exc}")
            QtWidgets.QMessageBox.critical(self, "Interactive Analysis Error", f"{exc}\n\n{tb}")

    def _populate_results(self, avg: dict, perp_override: dict, po_sources: list[str]):
        def g(section: str, key: str = "avg") -> float:
            return float((avg.get(section, {}) or {}).get(key, 0.0) or 0.0)

        rows: list[tuple[str, str]] = [
            ("V0 avg (m/s)", f"{g('V0'):.2f}"),
            ("V1 avg (m/s)", f"{g('V1'):.2f}"),
            ("V2 avg (m/s)", f"{g('V2'):.2f}"),
            ("TI1 avg (ms)", f"{g('ti1_ms'):.3f}"),
            ("TI2 avg (ms)", f"{g('ti2_ms'):.3f}"),
            ("H0 depth (m)", f"{float(avg.get('h1_m', 0.0) or 0.0):.3f}"),
            ("H1 depth (m)", f"{float(avg.get('h2_m', 0.0) or 0.0):.3f}"),
            ("DR total depth (m)", f"{float(avg.get('h1_m', 0.0) or 0.0) + float(avg.get('h2_m', 0.0) or 0.0):.3f}"),
        ]
        if perp_override:
            rows.append(("PO by shot (m)", ", ".join(f"S{k}={v:.2f}" for k, v in sorted(perp_override.items()))))
        rows.append(("PO source", ", ".join(po_sources) if po_sources else "none (PO=0)"))

        self.results_table.setRowCount(len(rows))
        for i, (k, v) in enumerate(rows):
            self.results_table.setItem(i, 0, QtWidgets.QTableWidgetItem(k))
            self.results_table.setItem(i, 1, QtWidgets.QTableWidgetItem(v))
        self.results_table.resizeColumnsToContents()

    def _save_excel_summary(self):
        if not self._last_analysis:
            QtWidgets.QMessageBox.information(self, "Save Excel", "Run Interactive Analysis first.")
            return
        lr = self._backend()
        a = self._last_analysis
        try:
            analysis = lr.build_analysis_from_layers(a["corrected_by_shot"], a["layer_results"])
            geometry_paths = [str(p) for p in self._discover_geometry_excels()]
            out_dir = OUTPUT_DIR / a["profile"]
            out_dir.mkdir(parents=True, exist_ok=True)
            path = lr.export_velocity_summary_excel(
                profile_name=a["profile"],
                cfg=a["cfg"],
                recv_positions=a["recv_positions"],
                shots_info=a["shots_info"],
                layer_results=a["layer_results"],
                analysis=analysis,
                output_dir=OUTPUT_DIR,
                geometry_excel_paths=geometry_paths,
            )
            self.statusBar().showMessage(f"Excel summary saved: {path}")
            QtWidgets.QMessageBox.information(self, "Save Excel", f"Velocity summary saved to:\n{path}")
        except Exception as exc:
            tb = traceback.format_exc(limit=8)
            QtWidgets.QMessageBox.critical(self, "Save Excel Error", f"{exc}\n\n{tb}")

    def _discover_geometry_excels(self) -> list[Path]:
        out: list[Path] = []
        for p in DATA_DIR.glob("LVL*.xls*"):
            name = p.name.lower()
            if "field_report" in name or "fieldreport" in name:
                continue
            out.append(p)
        return out

    def _load_seg2_folder(self, folder: Path) -> list[ShotGather]:
        files_raw = list(folder.glob("*.seg2")) + list(folder.glob("*.SEG2"))
        seen: set[str] = set()
        files: list[Path] = []
        for fp in files_raw:
            try:
                key = os.path.normcase(str(fp.resolve()))
            except Exception:
                key = os.path.normcase(str(fp))
            if key in seen:
                continue
            seen.add(key)
            files.append(fp)
        files = sorted(files, key=lambda p: p.name.lower())
        if not files:
            return []

        shots: list[ShotGather] = []
        for fp in files:
            st = obspy_read(str(fp), format="SEG2")
            if len(st) == 0:
                continue

            n_tr = len(st)
            n_samp = min(int(tr.stats.npts) for tr in st)
            data = np.zeros((n_tr, n_samp), dtype=np.float32)

            for i, tr in enumerate(st):
                d = np.asarray(tr.data, dtype=np.float32)
                data[i, :] = d[:n_samp]

            dt_ms = float(st[0].stats.delta) * 1000.0

            seg2 = getattr(st[0].stats, "seg2", {}) or {}
            delay_raw = seg2.get("DELAY", 0.0)

            # FFID: SEG2 SHOT_SEQUENCE_NUMBER first, then filename digits, else shot number.
            ssn = str(seg2.get("SHOT_SEQUENCE_NUMBER", "")).strip()
            digits = "".join(ch for ch in fp.stem if ch.isdigit())
            if ssn.isdigit():
                ffid = int(ssn)
            elif digits:
                ffid = int(digits)
            else:
                ffid = len(shots) + 1
            try:
                delay_ms = float(str(delay_raw).strip().replace(",", ".")) * 1000.0
            except Exception:
                delay_ms = 0.0

            shots.append(
                ShotGather(
                    source_file=fp,
                    shot_label=fp.stem,
                    ffid=ffid,
                    dt_ms=dt_ms,
                    delay_ms=delay_ms,
                    data=data,
                    header=dict(seg2),
                )
            )

        shots.sort(key=lambda s: s.ffid)
        return shots

    def set_shot(self, idx: int):
        if idx < 0 or idx >= len(self.shots):
            self.current_idx = -1
            return
        self.current_idx = idx
        self._render_current()

    def _render_current(self):
        if self.current_idx < 0 or self.current_idx >= len(self.shots):
            return

        shot = self.shots[self.current_idx]
        x_positions = self._receiver_positions_for_shot(shot)
        x_axis, x_label = self._x_axis_values(shot, x_positions)
        shot_x = self._shot_x_in_axis(shot, x_positions)
        data_proc = self._processed_trace_matrix(shot)
        t0 = self._time_zero_ms(shot)
        picks = self.picks_by_shot.get(self.current_idx, {})
        layers = [
            (name, self.layer_colors.get(name, "#888888"), layer_picks.get(self.current_idx, {}), name == self.active_layer)
            for name, layer_picks in self.pick_layers.items()
            if self.layer_visible.get(name, True)
        ]
        self.display.draw_shot(
            shot,
            self.settings,
            data_view_full=data_proc,
            t0_ms=t0,
            x_positions=x_axis,
            x_label=x_label,
            shot_x=shot_x,
            layers=layers,
        )
        # Preserve an active scroll-zoom across re-renders (e.g. after each pick).
        if self._view_xlim is not None and self._view_ylim is not None:
            try:
                self.display.ax.set_xlim(self._view_xlim)
                self.display.ax.set_ylim(self._view_ylim)
                self.display.canvas.draw_idle()
            except Exception:
                pass
        self._fill_shot_summary(shot, len(picks))
        profile = self._profile_name_from_open_folder() or "(custom)"
        self.statusBar().showMessage(
            f"Profile {profile} | Shot {self.current_idx + 1}/{len(self.shots)} (FFID {shot.ffid}) | dt={shot.dt_ms:.4f} ms | delay={shot.delay_ms:.2f} ms | picks {len(picks)}/{shot.n_traces}"
        )

    def _summary_rows_from_header(self, header: dict[str, Any]) -> Iterable[tuple[str, str]]:
        keys_preferred = [
            "FILE_NUMBER",
            "DELAY",
            "SAMPLE_INTERVAL",
            "STACK",
            "SHOT_SEQUENCE_NUMBER",
            "SOURCE_LOCATION",
            "RECEIVER_LOCATION",
            "ACQUISITION_DATE",
            "ACQUISITION_TIME",
            "DATE",
            "TIME",
        ]
        seen = set()
        for k in keys_preferred:
            if k in header:
                seen.add(k)
                yield k, str(header.get(k))
        for k in sorted(header.keys()):
            if k in seen:
                continue
            val = header.get(k)
            if isinstance(val, (str, int, float)):
                yield str(k), str(val)

    def _fill_shot_summary(self, shot: ShotGather, pick_count: int):
        data = np.asarray(shot.data, dtype=np.float32)
        rms = np.sqrt(np.mean(data * data, axis=1))
        peak = np.max(np.abs(data), axis=1)

        rows: list[tuple[str, str]] = [
            ("Shot label", shot.shot_label),
            ("Source file", str(shot.source_file.name)),
            ("FFID", str(shot.ffid)),
            ("Trace count", str(shot.n_traces)),
            ("Sample count", str(shot.n_samples)),
            ("Sample interval (ms)", f"{shot.dt_ms:.6g}"),
            ("SEG2 delay (ms)", f"{shot.delay_ms:.6g}"),
            ("TRIGGER_STATIC_MS", f"{self.settings.trigger_static_ms:.6g}"),
            ("Extra delay (ms)", f"{self.settings.extra_delay_ms:.6g}"),
            ("Time zero used (ms)", f"{self._time_zero_ms(shot):.6g}"),
            ("Polarity", self.settings.polarity),
            ("Gain", self.settings.gain_mode),
            ("Filter", self.settings.filter_mode),
            ("Display mode", self.settings.display_mode),
            ("RMS amplitude mean", f"{float(np.mean(rms)):.6g}"),
            ("Peak abs amplitude max", f"{float(np.max(peak)):.6g}"),
            ("Picks on this shot", str(int(pick_count))),
        ]

        for k, v in self._summary_rows_from_header(shot.header):
            rows.append((f"SEG2 {k}", v))
            if len(rows) >= 26:
                break

        self.summary_table.setRowCount(len(rows))
        for i, (k, v) in enumerate(rows):
            self.summary_table.setItem(i, 0, QtWidgets.QTableWidgetItem(k))
            self.summary_table.setItem(i, 1, QtWidgets.QTableWidgetItem(v))
        self.summary_table.resizeColumnsToContents()

    def prev_shot(self):
        if not self.shots:
            return
        idx = max(0, self.current_idx - 1)
        self.shot_list.setCurrentRow(idx)

    def next_shot(self):
        if not self.shots:
            return
        idx = min(len(self.shots) - 1, self.current_idx + 1)
        self.shot_list.setCurrentRow(idx)


def main() -> int:
    app = QtWidgets.QApplication(sys.argv)
    win = LvlStudioWindow()
    win.statusBar().showMessage(f"{QT_API} active. Open a SEG2 folder to start.")
    win.show()
    if hasattr(app, "exec"):
        return app.exec()
    return app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())
