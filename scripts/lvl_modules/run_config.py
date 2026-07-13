from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

from lvl_modules.app_paths import REFRACTION_SCRIPT


@dataclass
class RefractionRunRequest:
    profile: str
    geometry: str = "auto"
    device_type: str = "sw_maps"
    export_only: bool = False
    perp_excel: str = ""
    coord_excel: str = ""


def build_refraction_command(
    request: RefractionRunRequest,
    control_file: Path | None = None,
    minimal_plot_controls: bool = True,
    python_executable: str | None = None,
    script_path: Path = REFRACTION_SCRIPT,
) -> list[str]:
    profile = str(request.profile).strip()
    if not profile:
        raise ValueError("Profile is required")

    pyexe = python_executable or sys.executable
    cmd = [pyexe, str(script_path), profile]

    geom = str(request.geometry).strip().lower()
    if geom in ("100", "200"):
        cmd.extend(["--geom", geom])

    if request.export_only:
        cmd.append("--export-only")

    perp_excel = str(request.perp_excel).strip()
    if perp_excel:
        cmd.extend(["--perp-excel", perp_excel])

    coord_excel = str(request.coord_excel).strip()
    if coord_excel:
        cmd.extend(["--coord-excel", coord_excel])

    device = str(request.device_type).strip().lower() or "sw_maps"
    cmd.extend(["--device-type", device])

    if control_file is not None:
        cmd.extend(["--control-file", str(control_file)])

    if minimal_plot_controls:
        cmd.append("--minimal-plot-controls")

    return cmd
