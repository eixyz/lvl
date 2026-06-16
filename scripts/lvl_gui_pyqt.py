from __future__ import annotations

import ast
import json
import os
import subprocess
import sys
import importlib
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
LVL_SCRIPT = SCRIPT_DIR / "lvl_refraction_gui_support.py"


def _ensure(pkg: str, mod: str = "") -> Any:
    module_name = mod or pkg
    try:
        return importlib.import_module(module_name)
    except ImportError:
        print(f"Installing '{pkg}' ...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg, "-q"])
        return importlib.import_module(module_name)


def _load_qt_modules() -> tuple[Any, Any, Any, bool]:
    try:
        return (
            importlib.import_module("PyQt6.QtCore"),
            importlib.import_module("PyQt6.QtGui"),
            importlib.import_module("PyQt6.QtWidgets"),
            True,
        )
    except ImportError:
        return (
            _ensure("PyQt5", "PyQt5.QtCore"),
            _ensure("PyQt5", "PyQt5.QtGui"),
            _ensure("PyQt5", "PyQt5.QtWidgets"),
            False,
        )


QtCore, QtGui, QtWidgets, IS_QT6 = _load_qt_modules()

QProcess = QtCore.QProcess
Qt = QtCore.Qt
QApplication = QtWidgets.QApplication
QCheckBox = QtWidgets.QCheckBox
QComboBox = QtWidgets.QComboBox
QDoubleSpinBox = QtWidgets.QDoubleSpinBox
QGridLayout = QtWidgets.QGridLayout
QGroupBox = QtWidgets.QGroupBox
QHBoxLayout = QtWidgets.QHBoxLayout
QLabel = QtWidgets.QLabel
QLineEdit = QtWidgets.QLineEdit
QMainWindow = QtWidgets.QMainWindow
QMessageBox = QtWidgets.QMessageBox
QPlainTextEdit = QtWidgets.QPlainTextEdit
QPushButton = QtWidgets.QPushButton
QVBoxLayout = QtWidgets.QVBoxLayout
QWidget = QtWidgets.QWidget


def load_profiles(script_path: Path) -> list[str]:
    if not script_path.exists():
        return []
    try:
        src = script_path.read_text(encoding="utf-8")
        tree = ast.parse(src)
        for node in tree.body:
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name) and target.id == "PROFILES":
                        value = ast.literal_eval(node.value)
                        if isinstance(value, dict):
                            return sorted(str(k) for k in value.keys())
    except Exception:
        return []
    return []


class LvlControlWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LVL Refraction Control (PyQt)")
        self.resize(1080, 760)

        self.control_file = SCRIPT_DIR.parent / "output" / "_picker_control.json"
        self._cmd_seq = 0

        self.process = QProcess(self)
        self.process.setProgram(sys.executable)
        self.process.setWorkingDirectory(str(SCRIPT_DIR))
        self.process.readyReadStandardOutput.connect(self._on_stdout)
        self.process.readyReadStandardError.connect(self._on_stderr)
        self.process.started.connect(self._on_started)
        self.process.finished.connect(self._on_finished)

        self._build_ui()
        self._load_profiles()

    def _build_ui(self):
        central = QWidget(self)
        root = QVBoxLayout(central)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(8)

        cfg_group = QGroupBox("Run Configuration")
        cfg_layout = QGridLayout(cfg_group)

        self.profile_combo = QComboBox()
        self.profile_combo.setEditable(True)
        self.profile_combo.setMinimumWidth(220)

        self.geometry_combo = QComboBox()
        self.geometry_combo.addItems(["auto", "100", "200"])

        self.cb_all = QCheckBox("Process all profiles (--all)")
        self.cb_export_only = QCheckBox("Export only (--export-only)")

        self.cb_all.toggled.connect(self._on_all_toggled)

        cfg_layout.addWidget(QLabel("Profile:"), 0, 0)
        cfg_layout.addWidget(self.profile_combo, 0, 1)
        cfg_layout.addWidget(QLabel("Geometry override:"), 0, 2)
        cfg_layout.addWidget(self.geometry_combo, 0, 3)
        cfg_layout.addWidget(self.cb_all, 1, 0, 1, 2)
        cfg_layout.addWidget(self.cb_export_only, 1, 2, 1, 2)

        btn_row = QHBoxLayout()
        self.btn_run = QPushButton("Run")
        self.btn_stop = QPushButton("Stop")
        self.btn_clear = QPushButton("Clear Info")
        self.btn_open_output = QPushButton("Open Output Folder")

        self.btn_run.clicked.connect(self.run_process)
        self.btn_stop.clicked.connect(self.stop_process)
        self.btn_clear.clicked.connect(self._clear_log)
        self.btn_open_output.clicked.connect(self.open_output_folder)

        self.btn_stop.setEnabled(False)

        for btn in (self.btn_run, self.btn_stop, self.btn_open_output, self.btn_clear):
            btn.setMinimumHeight(26)

        btn_row.addWidget(self.btn_run)
        btn_row.addWidget(self.btn_stop)
        btn_row.addWidget(self.btn_open_output)
        btn_row.addStretch(1)
        btn_row.addWidget(self.btn_clear)

        live_group = QGroupBox("Live Picker Controls (GUI -> plot)")
        live_layout = QGridLayout(live_group)

        self.btn_prev_shot = QPushButton("Prev Shot")
        self.btn_save_next = QPushButton("Save/Next")
        self.btn_auto = QPushButton("Auto Pick")
        self.btn_invert = QPushButton("Invert")
        self.btn_timeline = QPushButton("Timelines")
        self.btn_filter_toggle = QPushButton("Filter ON/OFF")
        self.btn_save_png = QPushButton("Save PNG")
        self.btn_quit_pick = QPushButton("Quit Picking")

        for btn in (
            self.btn_prev_shot, self.btn_save_next, self.btn_auto, self.btn_invert,
            self.btn_timeline, self.btn_filter_toggle, self.btn_save_png, self.btn_quit_pick,
        ):
            btn.setMinimumHeight(24)
            btn.setMaximumHeight(28)

        self.gain_combo = QComboBox()
        self.gain_combo.addItems(["norm", "agc", "none"])
        self.stat_combo = QComboBox()
        self.stat_combo.addItems(["rms", "mean"])
        self.display_combo = QComboBox()
        self.display_combo.addItems(["both", "wiggle", "vd"])

        for combo in (self.gain_combo, self.stat_combo, self.display_combo):
            combo.setMaximumWidth(120)

        self.agc_spin = QDoubleSpinBox()
        self.agc_spin.setRange(5.0, 500.0)
        self.agc_spin.setSingleStep(1.0)
        self.agc_spin.setValue(200.0)
        self.agc_spin.setSuffix(" ms")
        self.agc_spin.setMaximumWidth(130)

        self.scale_spin = QDoubleSpinBox()
        self.scale_spin.setRange(0.20, 5.00)
        self.scale_spin.setSingleStep(0.01)
        self.scale_spin.setValue(0.85)
        self.scale_spin.setMaximumWidth(130)

        self.f1_edit = QLineEdit("2.0")
        self.f2_edit = QLineEdit("4.0")
        self.f3_edit = QLineEdit("140.0")
        self.f4_edit = QLineEdit("180.0")
        self.btn_apply_filter = QPushButton("Apply f1-f4")
        self.btn_apply_filter.setMinimumHeight(24)
        self.btn_apply_filter.setMaximumHeight(28)

        for le in (self.f1_edit, self.f2_edit, self.f3_edit, self.f4_edit):
            le.setMaximumWidth(90)

        self.btn_prev_shot.clicked.connect(lambda: self._send_cmd({"action": "prev"}))
        self.btn_save_next.clicked.connect(lambda: self._send_cmd({"action": "next"}))
        self.btn_auto.clicked.connect(lambda: self._send_cmd({"action": "auto"}))
        self.btn_invert.clicked.connect(lambda: self._send_cmd({"action": "invert"}))
        self.btn_timeline.clicked.connect(lambda: self._send_cmd({"action": "timeline"}))
        self.btn_filter_toggle.clicked.connect(lambda: self._send_cmd({"action": "filter_toggle"}))
        self.btn_save_png.clicked.connect(lambda: self._send_cmd({"action": "save_image"}))
        self.btn_quit_pick.clicked.connect(lambda: self._send_cmd({"action": "quit"}))

        self.gain_combo.currentTextChanged.connect(
            lambda v: self._send_cmd({"action": "gain", "value": str(v).strip().lower()}))
        self.stat_combo.currentTextChanged.connect(
            lambda v: self._send_cmd({"action": "agc_stat", "value": str(v).strip().lower()}))
        self.display_combo.currentTextChanged.connect(
            lambda v: self._send_cmd({"action": "display", "value": str(v).strip().lower()}))
        self.agc_spin.valueChanged.connect(
            lambda v: self._send_cmd({"action": "agc_window", "value": float(v)}))
        self.scale_spin.valueChanged.connect(
            lambda v: self._send_cmd({"action": "wiggle_scale", "value": float(v)}))
        self.btn_apply_filter.clicked.connect(self._apply_filter_from_gui)

        live_layout.addWidget(self.btn_prev_shot, 0, 0)
        live_layout.addWidget(self.btn_save_next, 0, 1)
        live_layout.addWidget(self.btn_auto, 0, 2)
        live_layout.addWidget(self.btn_quit_pick, 0, 3)

        live_layout.addWidget(self.btn_invert, 1, 0)
        live_layout.addWidget(self.btn_timeline, 1, 1)
        live_layout.addWidget(self.btn_filter_toggle, 1, 2)
        live_layout.addWidget(self.btn_save_png, 1, 3)

        live_layout.addWidget(QLabel("Gain"), 2, 0)
        live_layout.addWidget(self.gain_combo, 2, 1)
        live_layout.addWidget(QLabel("AGC stat"), 2, 2)
        live_layout.addWidget(self.stat_combo, 2, 3)

        live_layout.addWidget(QLabel("Display"), 3, 0)
        live_layout.addWidget(self.display_combo, 3, 1)
        live_layout.addWidget(QLabel("AGC window"), 3, 2)
        live_layout.addWidget(self.agc_spin, 3, 3)

        live_layout.addWidget(QLabel("Scale"), 4, 0)
        live_layout.addWidget(self.scale_spin, 4, 1)
        live_layout.addWidget(QLabel("f1"), 4, 2)
        live_layout.addWidget(self.f1_edit, 4, 3)

        live_layout.addWidget(QLabel("f2"), 5, 0)
        live_layout.addWidget(self.f2_edit, 5, 1)
        live_layout.addWidget(QLabel("f3"), 5, 2)
        live_layout.addWidget(self.f3_edit, 5, 3)

        live_layout.addWidget(QLabel("f4"), 6, 0)
        live_layout.addWidget(self.f4_edit, 6, 1)
        live_layout.addWidget(self.btn_apply_filter, 6, 2, 1, 2)

        log_group = QGroupBox("Status / Info")
        log_layout = QVBoxLayout(log_group)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(180)
        if IS_QT6:
            self.log.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        else:
            self.log.setLineWrapMode(QPlainTextEdit.NoWrap)
        log_layout.addWidget(self.log)

        self.lbl_state = QLabel("State: idle")
        self.lbl_last_cmd = QLabel("Last command: -")
        log_layout.addWidget(self.lbl_state)
        log_layout.addWidget(self.lbl_last_cmd)

        root.addWidget(cfg_group)
        root.addLayout(btn_row)
        root.addWidget(live_group)
        root.addWidget(log_group, stretch=1)

        self.setCentralWidget(central)
        self._set_live_controls_enabled(False)

    def _load_profiles(self):
        profiles = load_profiles(LVL_SCRIPT)
        if profiles:
            self.profile_combo.addItems(profiles)
            self.profile_combo.setCurrentIndex(0)
        else:
            self.profile_combo.addItem("120")

    def _on_all_toggled(self, checked: bool):
        self.profile_combo.setEnabled(not checked)

    def _append_log(self, text: str):
        if not text:
            return
        end_cursor = (QtGui.QTextCursor.MoveOperation.End
                      if IS_QT6 else QtGui.QTextCursor.End)
        self.log.moveCursor(end_cursor)
        try:
            self.log.document().setMaximumBlockCount(400)
        except Exception:
            pass
        self.log.insertPlainText(text)
        self.log.moveCursor(end_cursor)

    def _clear_log(self):
        self.log.clear()

    def _build_args(self) -> list[str]:
        args: list[str] = [str(LVL_SCRIPT)]

        if self.cb_all.isChecked():
            args.append("--all")
        else:
            profile = self.profile_combo.currentText().strip()
            if not profile:
                raise ValueError("Profile is required unless --all is enabled.")
            args.append(profile)
            geom = self.geometry_combo.currentText().strip().lower()
            if geom in ("100", "200"):
                args.append(geom)

        if self.cb_export_only.isChecked():
            args.append("--export-only")

        args.extend([
            "--control-file", str(self.control_file),
            "--minimal-plot-controls",
        ])

        return args

    def _set_live_controls_enabled(self, enabled: bool):
        widgets = [
            self.btn_prev_shot, self.btn_save_next, self.btn_auto, self.btn_invert,
            self.btn_timeline, self.btn_filter_toggle, self.btn_save_png, self.btn_quit_pick,
            self.gain_combo, self.stat_combo, self.display_combo,
            self.agc_spin, self.scale_spin,
            self.f1_edit, self.f2_edit, self.f3_edit, self.f4_edit,
            self.btn_apply_filter,
        ]
        for w in widgets:
            w.setEnabled(enabled)

    def _send_cmd(self, cmd: dict):
        not_running = (QProcess.ProcessState.NotRunning
                       if IS_QT6 else QProcess.NotRunning)
        if self.process.state() == not_running:
            return
        self.control_file.parent.mkdir(parents=True, exist_ok=True)
        self._cmd_seq += 1
        payload = dict(cmd)
        payload["id"] = self._cmd_seq
        with open(self.control_file, "w", encoding="utf-8") as fh:
            json.dump(payload, fh)
        self.lbl_last_cmd.setText(f"Last command: {payload}")

    def _apply_filter_from_gui(self):
        try:
            f1 = float(self.f1_edit.text().strip())
            f2 = float(self.f2_edit.text().strip())
            f3 = float(self.f3_edit.text().strip())
            f4 = float(self.f4_edit.text().strip())
        except Exception:
            QMessageBox.warning(self, "Invalid filter", "f1/f2/f3/f4 must be numeric.")
            return
        if not (f1 < f2 < f3 < f4):
            QMessageBox.warning(self, "Invalid filter", "Need f1 < f2 < f3 < f4.")
            return
        self._send_cmd({"action": "set_filter", "f1": f1, "f2": f2, "f3": f3, "f4": f4})

    def run_process(self):
        not_running = (QProcess.ProcessState.NotRunning
                       if IS_QT6 else QProcess.NotRunning)
        if self.process.state() != not_running:
            QMessageBox.information(self, "Busy", "A process is already running.")
            return

        if not LVL_SCRIPT.exists():
            QMessageBox.critical(self, "Missing Script", f"Cannot find: {LVL_SCRIPT}")
            return

        try:
            args = self._build_args()
        except ValueError as exc:
            QMessageBox.warning(self, "Invalid Input", str(exc))
            return

        self._append_log("\n=== Starting ===\n")
        self._append_log(f"Command: {sys.executable} {' '.join(args)}\n\n")
        self.lbl_state.setText("State: starting")

        try:
            if self.control_file.exists():
                self.control_file.unlink()
        except Exception:
            pass

        self.process.setArguments(args)
        self.process.start()

    def stop_process(self):
        not_running = (QProcess.ProcessState.NotRunning
                       if IS_QT6 else QProcess.NotRunning)
        if self.process.state() == not_running:
            return
        self._append_log("\n[INFO] Stop requested...\n")
        self.process.terminate()
        if not self.process.waitForFinished(2000):
            self.process.kill()

    def open_output_folder(self):
        out_dir = SCRIPT_DIR.parent / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        if sys.platform.startswith("win"):
            os.startfile(str(out_dir))
        else:
            QMessageBox.information(self, "Output", str(out_dir))

    def _on_started(self):
        self.btn_run.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self._set_live_controls_enabled(True)
        self.lbl_state.setText("State: running")

    def _on_finished(self, exit_code: int, _exit_status: QProcess.ExitStatus):
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        self._set_live_controls_enabled(False)
        self._append_log(f"\n=== Finished (exit={exit_code}) ===\n")
        self.lbl_state.setText(f"State: finished (exit={exit_code})")

    def _on_stdout(self):
        data = bytes(self.process.readAllStandardOutput()).decode("utf-8", errors="replace")
        self._append_log(data)

    def _on_stderr(self):
        data = bytes(self.process.readAllStandardError()).decode("utf-8", errors="replace")
        self._append_log(data)


def main():
    app = QApplication(sys.argv)
    win = LvlControlWindow()
    win.show()
    exec_fn = getattr(app, "exec", None) or getattr(app, "exec_", None)
    if exec_fn is None:
        raise RuntimeError("Could not find Qt application exec method.")
    sys.exit(exec_fn())


if __name__ == "__main__":
    main()
