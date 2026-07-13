from __future__ import annotations

import subprocess
import sys
import threading
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, ttk

from lvl_modules.app_paths import DATA_DIR, OUTPUT_DIR
from lvl_modules.control_bridge import write_picker_command
from lvl_modules.run_config import RefractionRunRequest, build_refraction_command


SCRIPT_DIR = DATA_DIR.parent / "scripts"


class CommandCenter:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("LVL Command Center")
        self.root.geometry("980x680")

        self.proc: subprocess.Popen[str] | None = None
        self.control_file = OUTPUT_DIR / "_picker_control.json"
        self._cmd_seq = 0

        self.profile_var = tk.StringVar(value="")
        self.geom_var = tk.StringVar(value="auto")
        self.device_var = tk.StringVar(value="sw_maps")
        self.seg2_dir_var = tk.StringVar(value="")
        self.field_report_var = tk.StringVar(value="")
        self.coord_file_var = tk.StringVar(value="")
        self.export_only_var = tk.BooleanVar(value=False)

        self._build_ui()
        self._load_profiles()

    def _build_ui(self) -> None:
        top = ttk.Frame(self.root, padding=10)
        top.pack(fill=tk.X)

        ttk.Label(top, text="Profile").grid(row=0, column=0, sticky="w", padx=4, pady=4)
        self.profile_combo = ttk.Combobox(top, textvariable=self.profile_var, state="normal", width=16)
        self.profile_combo.grid(row=0, column=1, sticky="w", padx=4, pady=4)

        ttk.Label(top, text="Geometry").grid(row=0, column=2, sticky="w", padx=4, pady=4)
        geom_combo = ttk.Combobox(top, textvariable=self.geom_var, values=["auto", "100", "200"], state="readonly", width=10)
        geom_combo.grid(row=0, column=3, sticky="w", padx=4, pady=4)

        ttk.Label(top, text="Device Type").grid(row=0, column=4, sticky="w", padx=4, pady=4)
        dev_combo = ttk.Combobox(top, textvariable=self.device_var, values=["sw_maps", "geomax"], state="readonly", width=12)
        dev_combo.grid(row=0, column=5, sticky="w", padx=4, pady=4)

        ttk.Checkbutton(top, text="Export only", variable=self.export_only_var).grid(row=0, column=6, sticky="w", padx=8, pady=4)

        row1 = ttk.Frame(self.root, padding=(10, 0, 10, 0))
        row1.pack(fill=tk.X)
        ttk.Label(row1, text="SEG2 Folder").pack(side=tk.LEFT, padx=4)
        ttk.Entry(row1, textvariable=self.seg2_dir_var).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=4)
        ttk.Button(row1, text="Browse", command=self._pick_seg2_folder).pack(side=tk.LEFT, padx=4)

        row2 = ttk.Frame(self.root, padding=(10, 0, 10, 0))
        row2.pack(fill=tk.X)
        ttk.Label(row2, text="Field Report").pack(side=tk.LEFT, padx=4)
        ttk.Entry(row2, textvariable=self.field_report_var).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=4)
        ttk.Button(row2, text="Browse", command=self._pick_field_report).pack(side=tk.LEFT, padx=4)

        row3 = ttk.Frame(self.root, padding=(10, 0, 10, 0))
        row3.pack(fill=tk.X)
        ttk.Label(row3, text="Coordinate File").pack(side=tk.LEFT, padx=4)
        ttk.Entry(row3, textvariable=self.coord_file_var).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=4)
        ttk.Button(row3, text="Browse", command=self._pick_coord_file).pack(side=tk.LEFT, padx=4)

        btns = ttk.Frame(self.root, padding=10)
        btns.pack(fill=tk.X)
        ttk.Button(btns, text="Run", command=self._run).pack(side=tk.LEFT, padx=4)
        ttk.Button(btns, text="Stop", command=self._stop).pack(side=tk.LEFT, padx=4)
        ttk.Button(btns, text="Open Output", command=self._open_output).pack(side=tk.LEFT, padx=4)
        ttk.Button(btns, text="Clear Log", command=self._clear_log).pack(side=tk.LEFT, padx=4)

        live = ttk.LabelFrame(self.root, text="Live Picker Controls (JSON bridge)", padding=8)
        live.pack(fill=tk.X, padx=10, pady=(0, 8))

        ttk.Button(live, text="Prev Shot", command=lambda: self._send_cmd({"action": "prev"})).grid(row=0, column=0, padx=3, pady=3)
        ttk.Button(live, text="Save/Next", command=lambda: self._send_cmd({"action": "next"})).grid(row=0, column=1, padx=3, pady=3)
        ttk.Button(live, text="Finalize", command=lambda: self._send_cmd({"action": "finalize"})).grid(row=0, column=2, padx=3, pady=3)
        ttk.Button(live, text="Auto Pick", command=lambda: self._send_cmd({"action": "auto"})).grid(row=0, column=3, padx=3, pady=3)
        ttk.Button(live, text="Quit Picking", command=lambda: self._send_cmd({"action": "quit"})).grid(row=0, column=4, padx=3, pady=3)

        ttk.Button(live, text="Invert", command=lambda: self._send_cmd({"action": "invert"})).grid(row=1, column=0, padx=3, pady=3)
        ttk.Button(live, text="Timelines", command=lambda: self._send_cmd({"action": "timeline"})).grid(row=1, column=1, padx=3, pady=3)
        ttk.Button(live, text="Filter On/Off", command=lambda: self._send_cmd({"action": "filter_toggle"})).grid(row=1, column=2, padx=3, pady=3)
        ttk.Button(live, text="Save PNG", command=lambda: self._send_cmd({"action": "save_image"})).grid(row=1, column=3, padx=3, pady=3)

        ttk.Label(live, text="Gain").grid(row=2, column=0, sticky="e", padx=3, pady=3)
        self.gain_combo = ttk.Combobox(live, values=["none", "norm", "agc"], state="readonly", width=10)
        self.gain_combo.set("norm")
        self.gain_combo.grid(row=2, column=1, sticky="w", padx=3, pady=3)
        self.gain_combo.bind("<<ComboboxSelected>>", lambda _e: self._send_cmd({"action": "gain", "value": self.gain_combo.get()}))

        ttk.Label(live, text="AGC stat").grid(row=2, column=2, sticky="e", padx=3, pady=3)
        self.stat_combo = ttk.Combobox(live, values=["rms", "mean"], state="readonly", width=10)
        self.stat_combo.set("rms")
        self.stat_combo.grid(row=2, column=3, sticky="w", padx=3, pady=3)
        self.stat_combo.bind("<<ComboboxSelected>>", lambda _e: self._send_cmd({"action": "agc_stat", "value": self.stat_combo.get()}))

        ttk.Label(live, text="Display").grid(row=2, column=4, sticky="e", padx=3, pady=3)
        self.display_combo = ttk.Combobox(live, values=["wiggle", "vd", "both"], state="readonly", width=10)
        self.display_combo.set("both")
        self.display_combo.grid(row=2, column=5, sticky="w", padx=3, pady=3)
        self.display_combo.bind("<<ComboboxSelected>>", lambda _e: self._send_cmd({"action": "display", "value": self.display_combo.get()}))

        self.log = tk.Text(self.root, wrap="none", height=28)
        self.log.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    def _load_profiles(self) -> None:
        profiles = sorted(p.name for p in DATA_DIR.iterdir() if p.is_dir()) if DATA_DIR.exists() else []
        self.profile_combo["values"] = profiles
        if profiles and not self.profile_var.get():
            self.profile_var.set(profiles[0])

    def _pick_seg2_folder(self) -> None:
        sel = filedialog.askdirectory(title="Select SEG2 folder", initialdir=str(DATA_DIR))
        if not sel:
            return
        self.seg2_dir_var.set(sel)
        p = Path(sel)
        if p.parent.resolve() == DATA_DIR.resolve():
            self.profile_var.set(p.name)

    def _pick_field_report(self) -> None:
        sel = filedialog.askopenfilename(
            title="Select Field Report Excel",
            filetypes=[("Excel", "*.xls *.xlsx"), ("All files", "*.*")],
            initialdir=str(DATA_DIR),
        )
        if sel:
            self.field_report_var.set(sel)

    def _pick_coord_file(self) -> None:
        sel = filedialog.askopenfilename(
            title="Select Coordinate Excel",
            filetypes=[("Excel", "*.xls *.xlsx"), ("All files", "*.*")],
            initialdir=str(DATA_DIR),
        )
        if sel:
            self.coord_file_var.set(sel)

    def _append(self, text: str) -> None:
        self.log.insert(tk.END, text)
        self.log.see(tk.END)

    def _clear_log(self) -> None:
        self.log.delete("1.0", tk.END)

    def _build_command(self) -> list[str]:
        profile = self.profile_var.get().strip()
        if not profile:
            raise ValueError("Profile is required")

        request = RefractionRunRequest(
            profile=profile,
            geometry=self.geom_var.get().strip().lower() or "auto",
            device_type=self.device_var.get().strip().lower() or "sw_maps",
            export_only=bool(self.export_only_var.get()),
            perp_excel=self.field_report_var.get().strip(),
            coord_excel=self.coord_file_var.get().strip(),
        )
        return build_refraction_command(
            request=request,
            control_file=self.control_file,
            minimal_plot_controls=True,
        )

    def _send_cmd(self, payload: dict) -> None:
        if not (self.proc and self.proc.poll() is None):
            self._append("[INFO] No active run for live command.\n")
            return
        self._cmd_seq += 1
        msg = write_picker_command(self.control_file, payload, self._cmd_seq)
        self._append(f"[CMD] {msg}\n")

    def _run(self) -> None:
        if self.proc and self.proc.poll() is None:
            self._append("[INFO] Process is already running.\n")
            return
        try:
            cmd = self._build_command()
        except Exception as exc:
            self._append(f"[ERROR] {exc}\n")
            return

        self._append("\n=== Starting ===\n")
        self._append("Command: " + " ".join(cmd) + "\n\n")
        try:
            if self.control_file.exists():
                self.control_file.unlink()
        except Exception:
            pass
        self.proc = subprocess.Popen(
            cmd,
            cwd=str(SCRIPT_DIR),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        threading.Thread(target=self._pump_output, daemon=True).start()

    def _pump_output(self) -> None:
        assert self.proc is not None
        for line in self.proc.stdout or []:
            self.root.after(0, self._append, line)
        code = self.proc.wait()
        self.root.after(0, self._append, f"\n=== Finished (exit={code}) ===\n")

    def _stop(self) -> None:
        if not self.proc or self.proc.poll() is not None:
            return
        self.proc.terminate()
        self._append("[INFO] Stop requested.\n")

    def _open_output(self) -> None:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        if sys.platform.startswith("win"):
            import os
            os.startfile(str(OUTPUT_DIR))
        else:
            self._append(str(OUTPUT_DIR) + "\n")


def main() -> None:
    root = tk.Tk()
    CommandCenter(root)
    root.mainloop()


if __name__ == "__main__":
    main()
