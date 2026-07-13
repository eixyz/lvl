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
    display_mode: str = "both"
    wiggle_scale: float = 0.45


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
        self._current_shot: ShotGather | None = None

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.canvas)

    def draw_shot(self, shot: ShotGather, settings: StudioSettings, picks: dict[int, float] | None = None):
        self._current_shot = shot
        self.ax.clear()

        data = np.asarray(shot.data, dtype=np.float32)
        n_tr = data.shape[0]
        stride = max(1, int(np.ceil(n_tr / float(max(1, settings.max_traces_render)))))
        self.last_stride = stride
        self.last_view_to_trace = list(range(0, n_tr, stride))
        data_view = data[::stride, :]
        self.last_n_view = int(data_view.shape[0])

        clip = float(np.percentile(np.abs(data_view), settings.clip_pct))
        if not np.isfinite(clip) or clip <= 0:
            clip = 1.0

        dt = float(shot.dt_ms)
        t0 = float(shot.delay_ms)
        t1 = t0 + (data_view.shape[1] - 1) * dt

        mode = str(settings.display_mode).strip().lower()
        if mode in ("vd", "both"):
            self.ax.imshow(
                data_view.T,
                cmap=settings.cmap,
                aspect="auto",
                interpolation="nearest",
                vmin=-clip,
                vmax=clip,
                extent=[0.5, data_view.shape[0] + 0.5, t1, t0],
            )

        if mode in ("wiggle", "both"):
            ns = data_view.shape[1]
            t = t0 + np.arange(ns, dtype=np.float32) * dt
            for i in range(data_view.shape[0]):
                tr = data_view[i].astype(np.float32)
                x = (i + 1) + (tr / clip) * float(settings.wiggle_scale)
                self.ax.plot(x, t, color="black", linewidth=0.5, alpha=0.95)

        if picks:
            xs: list[float] = []
            ys: list[float] = []
            for view_idx, tr_idx in enumerate(self.last_view_to_trace, start=1):
                if tr_idx in picks:
                    xs.append(float(view_idx))
                    ys.append(float(picks[tr_idx]))
            if xs:
                self.ax.plot(xs, ys, color="#d7191c", linewidth=1.2, marker="o", markersize=3)

        self.ax.set_title(
            f"{shot.shot_label} | FFID {shot.ffid} | traces {shot.n_traces} (view {data_view.shape[0]}) | samples {shot.n_samples}",
            fontsize=10,
        )
        self.ax.set_xlabel("Station index (possibly decimated for view)")
        self.ax.set_ylabel("Time (ms)")
        self.ax.grid(alpha=0.15)
        self.fig.tight_layout()
        self.canvas.draw_idle()


class LvlStudioWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LVL Studio - Independent Seismic Workspace")
        self.resize(1500, 900)

        self.settings = StudioSettings()
        self.shots: list[ShotGather] = []
        self.picks_by_shot: dict[int, dict[int, float]] = {}
        self.current_idx = -1
        self.pick_enabled = True
        self._display_window: QtWidgets.QMainWindow | None = None
        self._current_folder: Path | None = None

        self._build_ui()
        self._build_actions()

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
        self.spin_wiggle.setRange(0.1, 1.5)
        self.spin_wiggle.setDecimals(2)
        self.spin_wiggle.setSingleStep(0.05)
        self.spin_wiggle.setValue(self.settings.wiggle_scale)

        self.cmb_mode = QtWidgets.QComboBox()
        self.cmb_mode.addItems(["wiggle", "vd", "both"])
        self.cmb_mode.setCurrentText(self.settings.display_mode)

        self.cmb_cmap = QtWidgets.QComboBox()
        self.cmb_cmap.addItems(["gray", "seismic", "viridis", "magma", "cividis", "turbo", "plasma"])
        self.cmb_cmap.setCurrentText(self.settings.cmap)

        self.btn_pick_toggle = QtWidgets.QPushButton("Picker: ON")
        self.btn_pick_toggle.setCheckable(True)
        self.btn_pick_toggle.setChecked(True)
        self.btn_save_picks = QtWidgets.QPushButton("Save Picks")
        self.btn_clear_picks = QtWidgets.QPushButton("Clear Shot Picks")
        self.cmb_compute_geom = QtWidgets.QComboBox()
        self.cmb_compute_geom.addItems(["auto", "100", "200"])
        self.chk_layer_pick = QtWidgets.QCheckBox("Layer pick UI")
        self.chk_layer_pick.setChecked(False)
        self.btn_run_compute = QtWidgets.QPushButton("Run Full Computation")

        controls = QtWidgets.QWidget(self)
        form = QtWidgets.QFormLayout(controls)
        form.setContentsMargins(6, 6, 6, 6)
        form.addRow(self.btn_open)
        form.addRow(self.btn_prev, self.btn_next)
        form.addRow("Display mode", self.cmb_mode)
        form.addRow("Clip percentile", self.spin_clip)
        form.addRow("Max traces in view", self.spin_max_tr)
        form.addRow("Wiggle scale", self.spin_wiggle)
        form.addRow("Colormap", self.cmb_cmap)
        form.addRow(self.btn_pick_toggle)
        form.addRow(self.btn_save_picks)
        form.addRow(self.btn_clear_picks)
        form.addRow("Computation geom", self.cmb_compute_geom)
        form.addRow(self.chk_layer_pick)
        form.addRow(self.btn_run_compute)

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

        self.splitDockWidget(self.control_dock, self.shot_dock, _orient_vertical())
        self.resizeDocks([self.control_dock, self.shot_dock], [230, 330], _orient_vertical())
        self.resizeDocks([self.control_dock], [250], _orient_horizontal())
        self.resizeDocks([self.control_dock, self.shot_dock], [230, 260], _orient_horizontal())

        self.statusBar().showMessage("Open a SEG2 folder to start.")

        self.btn_open.clicked.connect(self.open_folder)
        self.btn_prev.clicked.connect(self.prev_shot)
        self.btn_next.clicked.connect(self.next_shot)
        self.btn_pick_toggle.clicked.connect(self._toggle_pick_mode)
        self.btn_save_picks.clicked.connect(self._save_picks_json)
        self.btn_clear_picks.clicked.connect(self._clear_current_shot_picks)
        self.btn_run_compute.clicked.connect(self._run_full_computation)
        self.shot_list.currentRowChanged.connect(self.set_shot)
        self.spin_clip.valueChanged.connect(self._on_view_settings_changed)
        self.spin_max_tr.valueChanged.connect(self._on_view_settings_changed)
        self.spin_wiggle.valueChanged.connect(self._on_view_settings_changed)
        self.cmb_mode.currentTextChanged.connect(self._on_view_settings_changed)
        self.cmb_cmap.currentTextChanged.connect(self._on_view_settings_changed)
        self.display.canvas.mpl_connect("button_press_event", self._on_canvas_click)

    def _build_actions(self):
        QAction = QtGui.QAction
        menu_file = self.menuBar().addMenu("File")
        act_open = QAction("Open SEG2 Folder", self)
        act_open.triggered.connect(self.open_folder)
        menu_file.addAction(act_open)

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
        self.settings.cmap = self.cmb_cmap.currentText()
        self._render_current()

    def _set_pick_mode(self, enabled: bool):
        self.pick_enabled = bool(enabled)
        self.btn_pick_toggle.setChecked(bool(enabled))
        self.btn_pick_toggle.setText("Picker: ON" if enabled else "Picker: OFF")

    def _toggle_pick_mode(self):
        self._set_pick_mode(bool(self.btn_pick_toggle.isChecked()))

    def _shot_pick_map(self, idx: int) -> dict[int, float]:
        return self.picks_by_shot.setdefault(int(idx), {})

    def _on_canvas_click(self, event: Any):
        if not self.pick_enabled:
            return
        if self.current_idx < 0 or self.current_idx >= len(self.shots):
            return
        if event.inaxes != self.display.ax:
            return
        if event.xdata is None or event.ydata is None:
            return

        view_idx = int(round(float(event.xdata))) - 1
        if view_idx < 0 or view_idx >= len(self.display.last_view_to_trace):
            return
        trace_idx = int(self.display.last_view_to_trace[view_idx])
        t_ms = float(event.ydata)

        shot_picks = self._shot_pick_map(self.current_idx)
        if event.button == MouseButton.LEFT:
            shot_picks[trace_idx] = round(t_ms, 2)
        elif event.button == MouseButton.RIGHT:
            if trace_idx in shot_picks:
                del shot_picks[trace_idx]
            else:
                nearest = None
                nearest_dist = 1e9
                for tr_idx in shot_picks:
                    dist = abs(tr_idx - trace_idx)
                    if dist < nearest_dist:
                        nearest = tr_idx
                        nearest_dist = dist
                if nearest is not None and nearest_dist <= max(1, self.display.last_stride):
                    del shot_picks[nearest]
        else:
            return
        self._render_current()

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
        self._render_current()

    def open_folder(self):
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Select SEG2 folder")
        if not folder:
            return
        self._current_folder = Path(folder)
        try:
            self.shots = self._load_seg2_folder(Path(folder))
        except Exception as exc:
            QtWidgets.QMessageBox.critical(self, "Load Error", str(exc))
            return

        self.shot_list.clear()
        for s in self.shots:
            self.shot_list.addItem(
                f"{s.shot_label} | FFID {s.ffid} | traces {s.n_traces} | samples {s.n_samples}"
            )

        if self.shots:
            self.shot_list.setCurrentRow(0)
            self.statusBar().showMessage(f"Loaded {len(self.shots)} shots from {folder}")
        else:
            self.current_idx = -1
            self.summary_table.setRowCount(0)
            self.statusBar().showMessage("No SEG2 files found.")

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

        geom_txt = self.cmb_compute_geom.currentText().strip().lower()
        geom_override = 100 if geom_txt == "100" else 200 if geom_txt == "200" else None
        enable_layer_pick = bool(self.chk_layer_pick.isChecked())

        try:
            self._export_picks_to_lvl_refraction_session(profile_name)

            from lvl_refraction import process_profile

            QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor if hasattr(QtCore.Qt, "CursorShape") else QtCore.Qt.WaitCursor)
            self.statusBar().showMessage(f"Running full computation for profile {profile_name}...")

            process_profile(
                profile_name=profile_name,
                pick_mode=False,
                geom_override=geom_override,
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
            ffid_raw = seg2.get("FILE_NUMBER", 0)
            delay_raw = seg2.get("DELAY", 0.0)

            try:
                ffid = int(ffid_raw)
            except Exception:
                ffid = len(shots) + 1
            try:
                delay_ms = float(delay_raw)
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
        picks = self.picks_by_shot.get(self.current_idx, {})
        self.display.draw_shot(shot, self.settings, picks=picks)
        self._fill_shot_summary(shot, len(picks))

        self.statusBar().showMessage(
            f"Shot {self.current_idx + 1}/{len(self.shots)} | FFID {shot.ffid} | traces {shot.n_traces} | samples {shot.n_samples}"
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
            ("Delay (ms)", f"{shot.delay_ms:.6g}"),
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
