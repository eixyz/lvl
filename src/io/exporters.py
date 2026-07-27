"""Excel/text/plot exporters for LVL refraction profile results.

All functions here are pure I/O + rendering: they take already-computed
picks/layer-fit/analysis data structures and write files (xlsx/txt/png).
No PyQt or interactive matplotlib event loops are used.
"""
from __future__ import annotations

import datetime
import math
from pathlib import Path
from typing import Any

import numpy as np
import openpyxl
from openpyxl.styles import PatternFill, Font, Alignment
from openpyxl.utils import get_column_letter
import matplotlib.pyplot as plt

from src.common.paths import require_active_project
from src.common.settings import (
    T_MAX_MS, theme_colors, _ensure_interactive_backend,
    BP_F1, BP_F2, BP_F3, BP_F4, BP_FFT_PAD, BULK_SHIFT_MS,
    STA_MS, LTA_MS, STALTA_TRIG, N_LAYERS,
)
from src.utils.geometry import load_profile_geometry_from_excels, load_midpoint_xyz_from_geometry_excels
from src.utils.math_utils import true_offset
from src.picker.preprocessing import apply_bulk_static
from src.refraction.travel_times import fit_line, depth_2layer, depth_3layer
from src.refraction.velocity_model import compute_layer_averages, _predict_time_from_fit, _choose_shot_fit

# Workbook style constants shared by the Excel export helpers below.
_FILL_HDR  = PatternFill("solid", fgColor="2E4057")
_FILL_IN   = PatternFill("solid", fgColor="FFFACD")   # lemon  = editable input
_FILL_FORM = PatternFill("solid", fgColor="E8F4F8")   # blue   = Excel formula
_FILL_OK   = PatternFill("solid", fgColor="DDFFDD")   # green  = depth result
_FONT_HDR  = Font(bold=True, color="FFFFFF")
_FONT_BOLD = Font(bold=True)
_FONT_ITA  = Font(italic=True, color="888888")

def _project():
    return require_active_project()

def report_directory(profile: str) -> Path:
    return _project().ensure_dir(_project().reports_dir_for(profile))


def plots_directory(profile: str) -> Path:
    p = _project().plots_dir_for(profile)
    return _project().ensure_dir(p)


def velocity_directory() -> Path:
    p = _project().velocity_dir
    return _project().ensure_dir(p)


def export_processing_report(
    profile_name: str,
    cfg: dict,
    analysis: dict | None = None,
    layer_results: dict | None = None,
    perp_by_shot: dict | None = None,
    inline_shift_by_shot: dict | None = None,
    po_sources: list | None = None,
    geometry_paths: list | None = None,
    report_paths: list | None = None,
    manual_geometry_paths: list | None = None,
    raw_folder: Path | str | None = None,
    picker_method: str | None = None,
    output_files: list | None = None,
    parameters: dict | None = None,
    filename_suffix: str = "",
) -> Path:
    """Write a human-readable processing_report.md for one profile.

    This exists purely for traceability/handover: which geometry file,
    field-report file, and coordinate file were actually MATCHED and used
    to process this profile (not every candidate file that happened to be
    searched), with what parameters, producing what results - so someone
    looking at `results/reports/<profile>/` six months from now (or a
    colleague who didn't run the processing themselves) can see exactly
    what went into it without re-deriving it from picks.json.
    """
    out_dir = report_directory(profile_name)
    out_path = out_dir / f"{profile_name}_processing_report{filename_suffix}.md"

    def _fmt_paths(paths) -> str:
        paths = [p for p in (paths or []) if p]
        if not paths:
            return "(none matched)"
        return "\n".join(f"- `{p}`" for p in paths)

    def _fmt_by_shot(d: dict, unit: str) -> str:
        if not d:
            return "(none set - default 0.0)"
        return ", ".join(f"S{k}={float(v):.2f}{unit}" for k, v in sorted(d.items()))

    avg = compute_layer_averages(layer_results or {}, {}) if layer_results else {}

    def g(section: str, key: str = "avg") -> float:
        return float((avg.get(section, {}) or {}).get(key, 0.0) or 0.0)

    lines = [
        f"# Processing report - profile {profile_name}",
        "",
        f"Generated: {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"Picker method: {picker_method or '(not recorded)'}",
        "",
        "## Parameters used",
        f"- Geometry type: {cfg.get('geom', '?')}",
        f"- Line number: {cfg.get('line_no', profile_name)}",
    ]
    for key, val in (parameters or {}).items():
        if key in ("geom_type", "line_no"):
            continue  # already shown above
        lines.append(f"- {key}: {val}")
    lines += [
        "",
        "## Input files actually used",
        "(files that were searched but did not match this profile are not listed)",
        "",
    ]
    if raw_folder:
        lines += [f"Raw SEG2 folder: `{raw_folder}`", ""]
    lines += [
        "Geometry/coordinate file(s) matched:",
        _fmt_paths(geometry_paths),
        "",
        "Manually uploaded geometry files consulted (subset of the above, if any matched):",
        _fmt_paths(manual_geometry_paths),
        "",
        "Field-report file(s) that provided a perpendicular offset:",
        _fmt_paths(report_paths),
        "",
        "## Perpendicular offset (PO) by shot",
        _fmt_by_shot(perp_by_shot, "m"),
        f"Source: {', '.join(po_sources) if po_sources else 'none (PO=0.0)'}",
        "",
        "## Inline shift by shot",
        _fmt_by_shot(inline_shift_by_shot, "m"),
        "",
    ]

    if avg:
        lines += [
            "## Results summary",
            f"- V0 avg: {g('V0'):.1f} m/s",
            f"- V1 avg: {g('V1'):.1f} m/s",
            f"- V2 avg: {g('V2'):.1f} m/s",
            f"- H0 depth: {float(avg.get('h1_m', 0.0) or 0.0):.2f} m",
            f"- H1 depth: {float(avg.get('h2_m', 0.0) or 0.0):.2f} m",
            "",
        ]

    if output_files:
        lines += ["## Output files written this run"] + [f"- `{p}`" for p in output_files] + [""]

    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  processing_report -> {_project().relative(out_path)}")
    return out_path

def export_velocity_summary_excel(profile_name: str,
                                  cfg: dict,
                                  recv_positions: Any,
                                  shots_info: list,
                                  layer_results: dict,
                                  analysis: dict,
                                  output_dir: Path,
                                  geometry_excel_paths: list | None = None,
                                  filename: str = "lvl_velocity_summary.xlsx",
                                  acquisition_time_de: str | None = None,
                                  seg2_mid_xyz: tuple | None = None) -> Path:
    """Write/append one profile row to a consolidated velocity summary workbook."""
    headers = [
        "NUM", "NAME", "SP", "SPREAD", "X_UTM", "Y_UTM", "ELEV",
        "V0E", "V0C", "V0", "H0", "V1", "TI1", "H1", "V2", "TI2", "DR",
        "LINE_LENGTH_SP_M", "LINE_LENGTH_M", "MIDDLE_STATION", "ACQ_DATETIME_DE",
        "V0_LEFT", "V0_RIGHT", "V1_LEFT", "V1_CENTER", "V1_RIGHT", "V2_LEFT", "V2_CENTER", "V2_RIGHT",
        "PROCESSED_ON",
    ]

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / filename
    if out_path.exists():
        wb = openpyxl.load_workbook(str(out_path))
        ws = wb.active
        for ci, h in enumerate(headers, start=1):
            old_h = ws.cell(row=1, column=ci).value
            if str(old_h).strip() != str(h):
                _chdr(ws, 1, ci, h)
    else:
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = "Summary"
        for ci, h in enumerate(headers, start=1):
            _chdr(ws, 1, ci, h)

    def _highlight_summary_headers() -> None:
        # Emphasize average/depth/intercept-time fields in the summary header row.
        target_headers = {
            "V0", "H0", "V1", "TI1", "H1", "V2", "TI2", "DR",
        }
        fill = PatternFill("solid", fgColor="FFD966")
        font = Font(bold=True, color="000000")
        for ci in range(1, ws.max_column + 1):
            hv = ws.cell(row=1, column=ci).value
            if str(hv).strip() in target_headers:
                cell = ws.cell(row=1, column=ci)
                cell.fill = fill
                cell.font = font
                cell.alignment = Alignment(horizontal="center")

    _highlight_summary_headers()

    shot_pos_by_id = {int(sid): float(pos) for sid, pos in shots_info}
    shot_ids = sorted(shot_pos_by_id, key=lambda sid: shot_pos_by_id[sid])
    sid_l = shot_ids[0] if shot_ids else None
    sid_m = shot_ids[len(shot_ids) // 2] if shot_ids else None
    sid_r = shot_ids[-1] if shot_ids else None

    def _get_v(shot_id: int | None, key: str) -> float:
        if shot_id is None:
            return 0.0
        per_side = (layer_results or {}).get(int(shot_id), {}) or {}
        _, fit = _choose_shot_fit(per_side)
        if not fit:
            return 0.0
        return float((fit or {}).get(key, 0.0) or 0.0)

    avg = compute_layer_averages(layer_results or {}, shot_pos_by_id)
    v0a = float((avg.get("V0", {}) or {}).get("avg", 0.0) or 0.0)
    v1a = float((avg.get("V1", {}) or {}).get("avg", 0.0) or 0.0)
    v2a = float((avg.get("V2", {}) or {}).get("avg", 0.0) or 0.0)
    t1 = float((avg.get("ti1_ms", {}) or {}).get("avg", 0.0) or 0.0)
    t2 = float((avg.get("ti2_ms", {}) or {}).get("avg", 0.0) or 0.0)
    d1 = float(avg.get("h1_m", 0.0) or 0.0)
    d2 = float(avg.get("h2_m", 0.0) or 0.0)

    name_base = f"LVL{str(profile_name).strip()}"
    name_new = f"{name_base}_new"
    sp_count = int(len(shots_info or []))
    profile_len_shot_points_m = float(abs(float(recv_positions[-1]) - float(recv_positions[0]))) if len(recv_positions) >= 2 else 0.0
    line_length_m = float(profile_len_shot_points_m)
    x_mid = y_mid = z_mid = None
    src = None
    acq_dt_final = acquisition_time_de

    g_len, gx, gy, gz, gsrc, s_min, s_mid, s_max, acq_dt_excel = load_profile_geometry_from_excels(
        profile_name=profile_name,
        excel_paths=(geometry_excel_paths or []),
    )
    excel_has_profile = any(v is not None for v in (g_len, gx, gy, gz, s_mid))
    if g_len is not None and g_len > 0.0:
        line_length_m = float(g_len)
    if gx is not None and gy is not None and gz is not None:
        # Coordinate workbook is primary source for UTM + ortho height.
        x_mid, y_mid, z_mid = float(gx), float(gy), float(gz)
        src = gsrc
    elif not excel_has_profile and seg2_mid_xyz is not None and len(seg2_mid_xyz) >= 3:
        # Only if profile is missing in coordinate workbook, fallback to SEG2-derived coordinates.
        sx, sy, sz = seg2_mid_xyz[0], seg2_mid_xyz[1], seg2_mid_xyz[2]
        if sx is not None and sy is not None and sz is not None:
            x_mid, y_mid, z_mid = float(sx), float(sy), float(sz)
            src = "SEG2(mid-shot)"
    if acq_dt_excel:
        acq_dt_final = acq_dt_excel
    elif len(recv_positions) >= 1:
        station_mid = (0.5 + (float(len(recv_positions)) + 0.5)) / 2.0
        x_mid, y_mid, z_mid, src = load_midpoint_xyz_from_geometry_excels(
            profile_name,
            station_mid=station_mid,
            excel_paths=(geometry_excel_paths or []),
        )

    if src is not None:
        src_rel = _project().relative(Path(src))
        if s_min is not None and s_mid is not None and s_max is not None:
            print(
                f"  Geometry loaded from: {src_rel} | "
                f"Stations {s_min:.1f}-{s_mid:.1f}-{s_max:.1f} | "
                f"LineLength={line_length_m:.3f} m"
            )
        else:
            print(f"  Midpoint XYZ loaded from: {src_rel}")

    v0e_vals = [v for v in (float(_get_v(sid_l, "V0_m_s")), float(_get_v(sid_r, "V0_m_s"))) if v > 0.0]
    v0e = float(np.mean(v0e_vals)) if v0e_vals else 0.0
    v0c = float(_get_v(sid_m, "V0_m_s")) if sid_m is not None else 0.0

    row_template = [
        "",  # NUM assigned after target row selection
        name_base,
        sp_count,
        f"{profile_len_shot_points_m:.1f} m",
        x_mid if x_mid is not None else "",
        y_mid if y_mid is not None else "",
        z_mid if z_mid is not None else "",
        round(v0e, 3),
        round(v0c, 3),
        round(v0a, 3),
        round(d1, 3),
        round(v1a, 3),
        round(t1, 3),
        round(d2, 3),
        round(v2a, 3),
        round(t2, 3),
        round(d1 + d2, 3),
        round(profile_len_shot_points_m, 3),
        round(line_length_m, 3),
        round(float(s_mid), 3) if s_mid is not None else "",
        acq_dt_final or "",
        round(_get_v(sid_l, "V0_m_s"), 3),
        round(_get_v(sid_r, "V0_m_s"), 3),
        round(_get_v(sid_l, "V1_m_s"), 3),
        round(_get_v(sid_m, "V1_m_s"), 3),
        round(_get_v(sid_r, "V1_m_s"), 3),
        round(_get_v(sid_l, "V2_m_s"), 3),
        round(_get_v(sid_m, "V2_m_s"), 3),
        round(_get_v(sid_r, "V2_m_s"), 3),
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    ]
    name_col = 2
    legacy_line_val = cfg.get("line_no", None)
    row_base = None
    row_new = None
    row_legacy = None
    for rr in range(2, ws.max_row + 1):
        val = ws.cell(row=rr, column=name_col).value
        val_s = str(val).strip()
        if val_s == name_base:
            row_base = rr
        elif val_s == name_new:
            row_new = rr
        elif legacy_line_val is not None and val_s == str(legacy_line_val).strip():
            row_legacy = rr

    def _next_num() -> int:
        nums = []
        for rr in range(2, ws.max_row + 1):
            v = ws.cell(row=rr, column=1).value
            try:
                nums.append(int(float(str(v).strip())))
            except Exception:
                continue
        return (max(nums) + 1) if nums else 1

    base_exists = (row_base is not None) or (row_legacy is not None)
    if base_exists:
        # Keep the original/base row intact; update/create the "_new" row.
        row = list(row_template)
        row[1] = name_new
        target_row = row_new if row_new is not None else (ws.max_row + 1)
    else:
        # No base row yet: write/update canonical row by profile token.
        row = list(row_template)
        row[1] = name_base
        target_row = row_base if row_base is not None else (ws.max_row + 1)

    existing_num = ws.cell(row=target_row, column=1).value if target_row <= ws.max_row else None
    try:
        row[0] = int(float(str(existing_num).strip()))
    except Exception:
        row[0] = _next_num()

    for ci, vv in enumerate(row, start=1):
        ws.cell(row=target_row, column=ci, value=vv)

    _autofit_xl(ws)
    try:
        wb.save(str(out_path))
        print(f"  Velocity summary -> {_project().relative(out_path)}")
        return out_path
    except PermissionError:
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        alt = out_path.with_name(f"{out_path.stem}_{stamp}{out_path.suffix}")
        wb.save(str(alt))
        print(f"  [WARN] Summary file locked; wrote fallback -> {_project().relative(alt)}")
        return alt

def _chdr(ws: Any, row: int, col: int, val: str):
    c = ws.cell(row=row, column=col, value=val)
    c.fill = _FILL_HDR
    c.font = _FONT_HDR
    c.alignment = Alignment(horizontal="center")

def _autofit_xl(ws: Any, min_w: int = 6, max_w: int = 30):
    for col in ws.columns:
        best = max((len(str(c.value)) for c in col if c.value is not None),
                   default=min_w)
        ws.column_dimensions[col[0].column_letter].width = min(best + 3, max_w)

def export_excel(profile_name: str, shots_info: list,
                 all_picks: dict, recv_positions: Any,
                 analysis: dict, cfg: dict,
                 corrected_by_shot: dict | None = None,
                 layer_results: dict | None = None,
                 perp_by_shot: dict | None = None,
                 filename_suffix: str = "") -> Path:
    """
    Write <profile>_picks.xlsx with four sheets.

    Config    Yellow cells = user input.
              B5 = perp_m, B6 = bulk_static_ms.
              All per-shot formula cells reference these two cells so
              changing perp_m or the static updates everything instantly.

    Shot_N    Per-shot picks with Excel formulas:
                True_off  = SQRT(Inline_off^2 + Config!$B$5^2)
                FB_corr   = FB_raw + Config!$B$6
              Summary rows: SLOPE / INTERCEPT / RSQ on corrected picks.

    Combined  All shots concatenated.

    Analysis  Python-computed velocities, intercept times, depths (green),
              R^2, RMS.  Depth formula reference below the table.
    """
    out_dir = _project().ensure_dir(_project().picks_dir_for(profile_name))
    xl_path = out_dir / f"{profile_name}_picks{filename_suffix}.xlsx"
    wb      = openpyxl.Workbook()

    # -- Config sheet -------------------------------------------------------
    ws_cfg = wb.active
    ws_cfg.title = "Config"
    ws_cfg["A1"] = f"LVL Refraction Analysis  --  Profile {profile_name}"
    ws_cfg["A1"].font = Font(bold=True, size=13)
    ws_cfg.merge_cells("A1:D1")

    def _inp(row: int, name: str, value: Any, unit: str = "", note: str = ""):
        ws_cfg.cell(row=row, column=1, value=name).font = _FONT_BOLD
        c = ws_cfg.cell(row=row, column=2, value=value)
        c.fill = _FILL_IN
        c.alignment = Alignment(horizontal="right")
        ws_cfg.cell(row=row, column=3, value=unit)
        if note:
            ws_cfg.cell(row=row, column=4, value=note).font = _FONT_ITA

    _inp(3,  "Profile name",      profile_name,           "",   "Folder name under data/")
    _inp(4,  "Geometry type",     cfg.get("geom", "?"),   "m",  "100 = 100m spread,  200 = 200m spread")
    _inp(5,  "Perp. distance",    cfg.get("perp_m", 0.0), "m",
             "Perpendicular shot-to-line offset.  "
             "True_off = SQRT(inline^2 + perp^2)  <- referenced in Shot sheets col G")
    _inp(6,  "Bulk static",       BULK_SHIFT_MS,          "ms",
             "Added to raw FB times.  ProMAX: bulk shift static -> add.  "
             "Negative = earlier.  <- referenced in Shot sheets col I")
    _inp(7,  "N layers",          N_LAYERS,                "",   "2 or 3")
    _inp(8,  "Ormsby f1",         BP_F1,                  "Hz", "Low-cut ramp start")
    _inp(9,  "Ormsby f2",         BP_F2,                  "Hz", "Low-cut ramp end / pass start")
    _inp(10, "Ormsby f3",         BP_F3,                  "Hz", "Pass end / high-cut ramp start")
    _inp(11, "Ormsby f4",         BP_F4,                  "Hz", "High-cut ramp end")
    _inp(12, "FFT zero-padding",  int(BP_FFT_PAD * 100),  "%",  "25 = ProMAX default")
    _inp(13, "STA window",        STA_MS,                 "ms", "Short-term average (auto-pick)")
    _inp(14, "LTA window",        LTA_MS,                 "ms", "Long-term average (auto-pick)")
    _inp(15, "STA/LTA threshold", STALTA_TRIG,             "",  "Trigger ratio")
    if perp_by_shot:
        shot_po_str = ", ".join(f"{sid}:{perp_by_shot[sid]:.3f}" for sid in sorted(perp_by_shot))
        _inp(16, "Perp by shot", shot_po_str, "m",
             "Used for corrected offsets Ck=sqrt((xk-SP)^2+PO^2)")

    ws_cfg.column_dimensions["A"].width = 24
    ws_cfg.column_dimensions["B"].width = 14
    ws_cfg.column_dimensions["C"].width = 6
    ws_cfg.column_dimensions["D"].width = 60

    # -- Per-shot sheets ----------------------------------------------------
    combined_rows: list = []

    for shot_id, shot_pos_m in shots_info:
        corr_rows = (corrected_by_shot or {}).get(shot_id, [])
        if not corr_rows:
            raw_picks = all_picks.get(shot_id, {})
            if not raw_picks:
                continue
            local_po = float((perp_by_shot or {}).get(shot_id, cfg.get("perp_m", 0.0)))
            corr_rows = []
            for trace_idx in sorted(raw_picks):
                if trace_idx >= len(recv_positions):
                    continue
                recv_pos = float(recv_positions[trace_idx])
                inline   = recv_pos - shot_pos_m
                raw_fb = float(raw_picks[trace_idx])
                bulk_fb = float(apply_bulk_static({trace_idx: raw_fb})[trace_idx])
                corr_rows.append({
                    "trace_idx": trace_idx,
                    "trace_no": trace_idx + 1,
                    "recv_pos_m": recv_pos,
                    "shot_pos_m": shot_pos_m,
                    "inline_signed_m": inline,
                    "inline_abs_m": abs(inline),
                    "perp_m": local_po,
                    "true_off_m": float(true_offset(abs(inline), local_po)),
                    "fb_raw_ms": raw_fb,
                    "fb_bulk_ms": bulk_fb,
                    "fb_interp_inline_ms": bulk_fb,
                    "fb_interp_geom_ms": bulk_fb,
                    "side": "L" if inline < 0 else "R",
                })

        if not corr_rows:
            continue

        ws = wb.create_sheet(title=f"Shot_{shot_id}")
        col_headers = ["Trace", "FFID", "Recv_pos_m", "Shot_pos_m",
               "Inline_signed_m", "Inline_abs_m", "Inline_shift_m", "Inline_corr_m",
               "Perp_m", "True_off_m",
                   "FB_raw_ms", "FB_bulk_ms", "FB_interp_geom_ms", "Side"]
        for ci, h in enumerate(col_headers, 1):
            _chdr(ws, 1, ci, h)

        for row_i, row_d in enumerate(corr_rows, start=2):
            ws.cell(row=row_i, column=1, value=int(row_d["trace_no"]))
            ws.cell(row=row_i, column=2, value=shot_id)
            ws.cell(row=row_i, column=3, value=round(float(row_d["recv_pos_m"]), 3))
            ws.cell(row=row_i, column=4, value=round(float(row_d["shot_pos_m"]), 3))
            ws.cell(row=row_i, column=5, value=round(float(row_d["inline_signed_m"]), 3))
            ws.cell(row=row_i, column=6, value=round(float(row_d["inline_abs_m"]), 3))
            ws.cell(row=row_i, column=7, value=round(float(row_d.get("inline_shift_m", 0.0)), 3))
            ws.cell(row=row_i, column=8, value=round(float(row_d.get("inline_corr_m", row_d["inline_abs_m"])), 3))
            ws.cell(row=row_i, column=9, value=round(float(row_d["perp_m"]), 3))
            ws.cell(row=row_i, column=10, value=round(float(row_d["true_off_m"]), 3))
            ws.cell(row=row_i, column=11, value=round(float(row_d["fb_raw_ms"]), 3))
            ws.cell(row=row_i, column=12, value=round(float(row_d["fb_bulk_ms"]), 3))
            ws.cell(row=row_i, column=13, value=round(float(row_d.get("fb_interp_geom_ms", row_d["fb_interp_inline_ms"])), 3))
            ws.cell(row=row_i, column=14, value=str(row_d.get("side", "R")))

            combined_rows.append({
                "shot_id": shot_id,
                "trace": int(row_d["trace_no"]),
                "recv_pos": float(row_d["recv_pos_m"]),
                "shot_pos": float(row_d["shot_pos_m"]),
                "inline_abs": float(row_d["inline_abs_m"]),
                "inline_shift": float(row_d.get("inline_shift_m", 0.0)),
                "inline_corr": float(row_d.get("inline_corr_m", row_d["inline_abs_m"])),
                "true_off": float(row_d["true_off_m"]),
                "fb_raw": float(row_d["fb_raw_ms"]),
                "fb_bulk": float(row_d["fb_bulk_ms"]),
                "fb_interp": float(row_d.get("fb_interp_geom_ms", row_d["fb_interp_inline_ms"])),
                "side": str(row_d.get("side", "R")),
            })

        # Summary: SLOPE / INTERCEPT / RSQ on corrected picks vs true offset
        n_data  = len(corr_rows)
        sum_row = n_data + 3
        _chdr(ws, sum_row, 1, "Metric")
        _chdr(ws, sum_row, 2, "Value")
        _chdr(ws, sum_row, 3, "Range used")

        x_rng = f"J2:J{n_data + 1}"
        t_rng = f"M2:M{n_data + 1}"
        for ri, (label, formula, basis) in enumerate([
            ("Velocity (m/s)",
             f"=IFERROR(1000/SLOPE({t_rng},{x_rng}),\"n/a\")",
             "SLOPE(FB_interp_geom, True_off)"),
            ("Intercept (ms)",
             f"=IFERROR(INTERCEPT({t_rng},{x_rng}),\"n/a\")",
             "INTERCEPT(FB_interp_geom, True_off)"),
            ("R^2",
             f"=IFERROR(RSQ({t_rng},{x_rng}),\"n/a\")",
             "RSQ(FB_interp_geom, True_off)"),
        ], start=sum_row + 1):
            ws.cell(row=ri, column=1, value=label).font = _FONT_BOLD
            ws.cell(row=ri, column=2, value=formula).fill = _FILL_FORM
            ws.cell(row=ri, column=3, value=basis).font  = _FONT_ITA

        _autofit_xl(ws)

    # -- Combined sheet -----------------------------------------------------
    if combined_rows:
        ws_c = wb.create_sheet(title="Combined")
        for ci, h in enumerate(["Shot_ID", "Trace", "Recv_pos_m", "Shot_pos_m",
                     "Inline_abs_m", "Inline_shift_m", "Inline_corr_m", "True_off_m", "FB_raw_ms",
                     "FB_bulk_ms", "FB_interp_geom_ms", "Side"], 1):
            _chdr(ws_c, 1, ci, h)
        for ri, row in enumerate(combined_rows, start=2):
            ws_c.cell(row=ri, column=1, value=row["shot_id"])
            ws_c.cell(row=ri, column=2, value=row["trace"])
            ws_c.cell(row=ri, column=3, value=round(row["recv_pos"], 3))
            ws_c.cell(row=ri, column=4, value=round(row["shot_pos"], 3))
            ws_c.cell(row=ri, column=5, value=round(row["inline_abs"], 3))
            ws_c.cell(row=ri, column=6, value=round(row.get("inline_shift", 0.0), 3))
            ws_c.cell(row=ri, column=7, value=round(row.get("inline_corr", row["inline_abs"]), 3))
            ws_c.cell(row=ri, column=8, value=round(row["true_off"], 3))
            ws_c.cell(row=ri, column=9, value=round(row["fb_raw"], 3))
            ws_c.cell(row=ri, column=10, value=round(row["fb_bulk"], 3))
            ws_c.cell(row=ri, column=11, value=round(row["fb_interp"], 3))
            ws_c.cell(row=ri, column=12, value=row["side"])
        _autofit_xl(ws_c)

    # -- Layer picks sheet -------------------------------------------------
    if layer_results:
        ws_l = wb.create_sheet(title="Layer_Picks")
        l_hdrs = ["Shot_ID", "Side", "x1", "x2", "x3", "x4", "x5", "x6",
              "t1 (ms)", "t2 (ms)", "t3 (ms)", "t4 (ms)", "t5 (ms)", "t6 (ms)",
              "V0 (m/s)", "V1 (m/s)", "V2 (m/s)",
              "ti1 (ms)", "ti2 (ms)", "h1 (m)", "h2 (m)", "n_points"]
        for ci, h in enumerate(l_hdrs, 1):
            _chdr(ws_l, 1, ci, h)

        rr = 2
        for shot_id in sorted(layer_results):
            shot_sides = layer_results.get(shot_id, {})
            ordered_sides = [s for s in ("L", "R", "ALL") if s in shot_sides]
            ordered_sides += [s for s in shot_sides if s not in ordered_sides]
            for side in ordered_sides:
                res = shot_sides.get(side)
                if not res:
                    continue
                wins = res.get("windows", [None, None, None])
                pts = list(res.get("picked_points", []))
                pts += [(None, None)] * (6 - len(pts))

                def _w(i: int, j: int) -> Any:
                    w = wins[i] if i < len(wins) else None
                    return round(float(w[j]), 3) if w else ""

                def _px(i: int) -> Any:
                    xv = pts[i][0] if i < len(pts) else None
                    return round(float(xv), 3) if xv is not None else ""

                def _pt(i: int) -> Any:
                    tv = pts[i][1] if i < len(pts) else None
                    return round(float(tv), 3) if tv is not None else ""

                ws_l.cell(row=rr, column=1, value=shot_id)
                ws_l.cell(row=rr, column=2, value=side)
                ws_l.cell(row=rr, column=3, value=_px(0) if _px(0) != "" else _w(0, 0))
                ws_l.cell(row=rr, column=4, value=_px(1) if _px(1) != "" else _w(0, 1))
                ws_l.cell(row=rr, column=5, value=_px(2) if _px(2) != "" else _w(1, 0))
                ws_l.cell(row=rr, column=6, value=_px(3) if _px(3) != "" else _w(1, 1))
                ws_l.cell(row=rr, column=7, value=_px(4) if _px(4) != "" else _w(2, 0))
                ws_l.cell(row=rr, column=8, value=_px(5) if _px(5) != "" else _w(2, 1))
                ws_l.cell(row=rr, column=9, value=_pt(0))
                ws_l.cell(row=rr, column=10, value=_pt(1))
                ws_l.cell(row=rr, column=11, value=_pt(2))
                ws_l.cell(row=rr, column=12, value=_pt(3))
                ws_l.cell(row=rr, column=13, value=_pt(4))
                ws_l.cell(row=rr, column=14, value=_pt(5))
                ws_l.cell(row=rr, column=15, value=round(float(res.get("V0_m_s", 0.0)), 2))
                ws_l.cell(row=rr, column=16, value=round(float(res.get("V1_m_s", 0.0)), 2))
                ws_l.cell(row=rr, column=17, value=round(float(res.get("V2_m_s", 0.0)), 2))
                ws_l.cell(row=rr, column=18, value=round(float(res.get("ti1_ms", 0.0)), 3))
                ws_l.cell(row=rr, column=19, value=round(float(res.get("ti2_ms", 0.0)), 3))
                ws_l.cell(row=rr, column=20, value=round(float(res.get("h1_m", 0.0)), 3))
                ws_l.cell(row=rr, column=21, value=round(float(res.get("h2_m", 0.0)), 3))
                ws_l.cell(row=rr, column=22, value=int(res.get("n_points", 0)))

                for col in (20, 21):
                    cell = ws_l.cell(row=rr, column=col)
                    if cell.value not in ("", 0, 0.0):
                        cell.fill = _FILL_OK
                rr += 1

        _autofit_xl(ws_l)

    # -- Layer averages sheet ---------------------------------------------
    shot_pos_by_id = {int(sid): float(sp) for sid, sp in shots_info}
    layer_avg = compute_layer_averages(layer_results or {}, shot_pos_by_id)
    if layer_avg:
        ws_la = wb.create_sheet(title="Layer_Averages")
        hdrs = ["Layer", "V_off (m/s)", "V_center (m/s)", "V_avg (m/s)",
                "ti_off (ms)", "ti_center (ms)", "ti_avg (ms)"]
        for ci, h in enumerate(hdrs, 1):
            _chdr(ws_la, 1, ci, h)

        layer_map = [("V0", "ti1_ms", "Layer 1"), ("V1", "ti1_ms", "Layer 2"), ("V2", "ti2_ms", "Layer 3")]
        rr = 2
        for vkey, tkey, lname in layer_map:
            vv = layer_avg.get(vkey, {})
            tt = layer_avg.get(tkey, {})
            ws_la.cell(row=rr, column=1, value=lname)
            ws_la.cell(row=rr, column=2, value=round(float(vv.get("off", 0.0)), 3))
            ws_la.cell(row=rr, column=3, value=round(float(vv.get("center", 0.0)), 3))
            ws_la.cell(row=rr, column=4, value=round(float(vv.get("avg", 0.0)), 3))
            ws_la.cell(row=rr, column=5, value=round(float(tt.get("off", 0.0)), 3))
            ws_la.cell(row=rr, column=6, value=round(float(tt.get("center", 0.0)), 3))
            ws_la.cell(row=rr, column=7, value=round(float(tt.get("avg", 0.0)), 3))
            rr += 1

        ws_la.cell(row=rr + 1, column=1, value="Depth h1 (m)").font = _FONT_BOLD
        ws_la.cell(row=rr + 1, column=2, value=round(float(layer_avg.get("h1_m", 0.0)), 3)).fill = _FILL_OK
        ws_la.cell(row=rr + 2, column=1, value="Depth h2 (m)").font = _FONT_BOLD
        ws_la.cell(row=rr + 2, column=2, value=round(float(layer_avg.get("h2_m", 0.0)), 3)).fill = _FILL_OK
        ws_la.cell(row=rr + 4, column=1, value="Off-end shots").font = _FONT_BOLD
        ws_la.cell(row=rr + 4, column=2, value=", ".join(str(s) for s in layer_avg.get("off_shots", [])))
        ws_la.cell(row=rr + 5, column=1, value="Center shots").font = _FONT_BOLD
        ws_la.cell(row=rr + 5, column=2, value=", ".join(str(s) for s in layer_avg.get("center_shots", [])))
        _autofit_xl(ws_la)

    # -- Analysis sheet -----------------------------------------------------
    ws_a = wb.create_sheet(title="Analysis")
    a_hdrs = ["Shot_ID", "FFID", "V1 (m/s)", "V2 (m/s)", "V3 (m/s)",
              "ti_1 (ms)", "ti_2 (ms)", "h1 (m)", "h2 (m)", "RMS (ms)", "n_picks"]
    for ci, h in enumerate(a_hdrs, 1):
        _chdr(ws_a, 1, ci, h)

    row_a = 2
    for shot_id, _ in shots_info:
        res = (analysis or {}).get(shot_id)
        if res is None:
            continue
        deps = res.get("depths", {})
        segs = res.get("segments", [])

        def _vel(i: int) -> Any:
            return round(segs[i]["velocity_m_s"], 1) if i < len(segs) else ""

        def _dep(key: str) -> Any:
            val = deps.get(key)
            return round(val, 2) if val else ""

        ws_a.cell(row=row_a, column=1,  value=shot_id)
        ws_a.cell(row=row_a, column=2,  value=shot_id)
        ws_a.cell(row=row_a, column=3,  value=_vel(0))
        ws_a.cell(row=row_a, column=4,  value=_vel(1))
        ws_a.cell(row=row_a, column=5,  value=_vel(2))
        ws_a.cell(row=row_a, column=6,
                  value=round(deps["ti1_ms"], 3) if deps.get("ti1_ms") else "")
        ws_a.cell(row=row_a, column=7,
                  value=round(deps["ti2_ms"], 3) if deps.get("ti2_ms") else "")
        ws_a.cell(row=row_a, column=8,  value=_dep("h1_m"))
        ws_a.cell(row=row_a, column=9,  value=_dep("h2_m"))
        ws_a.cell(row=row_a, column=10, value=round(res.get("rms_ms", 0.0), 3))
        ws_a.cell(row=row_a, column=11, value=len(all_picks.get(shot_id, {})))

        for col in (8, 9):
            cell = ws_a.cell(row=row_a, column=col)
            if cell.value != "":
                cell.fill = _FILL_OK
        row_a += 1

    note_row = row_a + 2
    notes = [
        "Depth formulas (intercept-time method):",
        "  2-layer:  h1 = (ti1_ms/1000 x V1 x V2) / (2 x SQRT(V2^2 - V1^2))",
        "  3-layer:  h2 = (ti2_eff_ms/1000 x V1 x V3) / (2 x SQRT(V3^2 - V1^2))"
        "    where  ti2_eff = ti2 - 2 x h1 x cos(ic12) / V1 x 1000",
    ]
    for oi, note in enumerate(notes):
        r = note_row + oi
        ws_a.cell(row=r, column=1, value=note).font = (
            _FONT_BOLD if oi == 0 else _FONT_ITA)
        ws_a.merge_cells(f"A{r}:K{r}")

    _autofit_xl(ws_a)
    wb.save(str(xl_path))
    print(f"  Excel     -> {_project().relative(xl_path)}")
    return xl_path

def export_picks_txt(profile_name: str, shots_info: list,
                     all_picks: dict, recv_positions: Any) -> Path:
    """Write Ensemble / SOURCE / CHAN / OFFSET / FB_PICK text file."""
    out_dir  = _project().ensure_dir(_project().picks_dir_for(profile_name))
    txt_path = out_dir / f"{profile_name}_picks_clean.txt"
    header   = (f"{'Ensemble':>10} {'#':>4} {'SOURCE':>8} "
                f"{'CHAN':>6} {'OFFSET':>10} {'FB_PICK':>10}")

    with open(txt_path, "w", encoding="utf-8") as fh:
        fh.write(header + "\n")
        ens = 1
        for shot_idx, (shot_id, shot_pos_m) in enumerate(shots_info):
            corr = apply_bulk_static(all_picks.get(shot_id, {}))
            if not corr:
                continue
            if shot_idx > 0:
                fh.write("\n")
            for trace_idx in sorted(corr):
                if trace_idx >= len(recv_positions):
                    continue
                offset_m = recv_positions[trace_idx] - shot_pos_m
                fh.write(f"{ens:>10} {ens:>4} {shot_id:>8} "
                         f"{trace_idx+1:>6} {offset_m:>10.2f} "
                         f"{corr[trace_idx]:>10.3f}\n")
                ens += 1

    print(f"  picks_txt -> {_project().relative(txt_path)}")
    return txt_path

def export_tx_plot(profile_name: str, shots_info: list,
                   all_picks: dict, recv_positions: Any,
                   analysis: dict, perp_m: float = 0.0,
                   shot_label_pos: dict | None = None,
                   corrected_by_shot: dict | None = None,
                   layer_results: dict | None = None,
                   theme: str = "light") -> Path:
    """
    Final T-X summary plot.
    - Scatter of corrected picks (absolute true offset vs corrected pick time)
    - Dashed fitted-line overlay with velocity labels
    - Y-axis inverted: t = 0 at top  (seismic convention: time grows downward)
    - Theme follows the current THEME setting at export time
    """
    c = theme_colors(theme)
    fig, (ax_ch, ax_tx) = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=False)
    fig.patch.set_facecolor(c["fig_bg"])
    ax_ch.set_facecolor(c["ax_bg"])
    ax_tx.set_facecolor(c["ax_bg"])

    pal = ["#e63946", "#2a9d8f", "#e9c46a", "#a8dadc", "#f4a261", "#457b9d"]
    mks = ["o", "s", "^", "D", "P", "X"]

    for shot_id, shot_pos_m in shots_info:
        corr = apply_bulk_static(all_picks.get(shot_id, {}))
        if not corr:
            continue
        col = pal[(shot_id - 1) % len(pal)]
        mk  = mks[(shot_id - 1) % len(mks)]
        label_pos_m = (shot_label_pos or {}).get(shot_id, shot_pos_m)

        x_geom, t_vals = [], []
        for idx in sorted(corr):
            if idx < len(recv_positions):
                x_geom.append(float(recv_positions[idx]))
                t_vals.append(corr[idx])

        # Left panel: channel-domain view (QC style)
        tr = [idx + 1 for idx in sorted(corr) if idx < len(recv_positions)]
        ax_ch.plot(tr, t_vals, "o-", ms=3, lw=1.0, color=col,
               label=f"Shot {shot_id} (@ {label_pos_m:.1f} m)")

        # Right panel: geometry T-X view
        ax_tx.scatter(x_geom, t_vals, color=col, marker=mk, s=28, zorder=5,
                  label=f"Shot {shot_id} (@ {label_pos_m:.1f} m)")

        # Dashed fitted T-X lines for each available side (L/R/ALL).
        side_map = (layer_results or {}).get(shot_id, {}) or {}
        rows_corr = list((corrected_by_shot or {}).get(shot_id, []))
        ordered_sides = [s for s in ("L", "R", "ALL") if s in side_map]
        ordered_sides += [s for s in side_map if s not in ordered_sides]
        for side in ordered_sides:
            fit_res = side_map.get(side)
            if not fit_res or not rows_corr:
                continue
            if side in ("L", "R"):
                rows_plot = [r for r in rows_corr if r.get("side") == side]
                if not rows_plot:
                    continue
            else:
                rows_plot = rows_corr

            x_pred: list = []
            t_pred: list = []
            for rr in rows_plot:
                xg = float(rr.get("recv_pos_m", shot_pos_m))
                x_abs = float(rr.get("true_off_m", abs(xg - shot_pos_m)))
                tp = _predict_time_from_fit(x_abs, fit_res)
                if tp is None:
                    continue
                x_pred.append(xg)
                t_pred.append(float(tp))

            if not x_pred:
                continue

            order = np.argsort(np.asarray(x_pred, dtype=float))
            xp = np.asarray(x_pred, dtype=float)[order]
            tp = np.asarray(t_pred, dtype=float)[order]
            ls = "--" if side != "R" else ":"
            ax_tx.plot(xp, tp, ls, color=col, lw=1.6, alpha=c["fit_alpha"])

            segs = list((fit_res or {}).get("segments", []) or [])
            for seg in segs:
                vel = float(seg.get("velocity_m_s", 0.0) or 0.0)
                if vel <= 0.0:
                    continue
                sl = float(seg.get("slope_ms_m", 0.0) or 0.0)
                ic = float(seg.get("intercept_ms", 0.0) or 0.0)
                if sl <= 0.0:
                    continue
                xm = float(np.median(xp))
                tm = sl * abs(xm - shot_pos_m) + ic
                ax_tx.text(xm, tm - 3.5, f"{side}:{vel:.0f}", color=col,
                       fontsize=7, ha="center", fontweight="bold")

    # Seismic convention: t = 0 at top, time grows downward
    ax_ch.set_ylim(T_MAX_MS, 0.0)
    ax_tx.set_ylim(T_MAX_MS, 0.0)
    ax_ch.set_xlabel("Channel", color=c["label"])
    ax_ch.set_ylabel("First-break time  (ms)", color=c["label"])
    ax_ch.set_title("Shot picks vs channel", color=c["text"], fontsize=10)

    ax_tx.set_xlabel("Receiver position (m)", color=c["label"])
    ax_tx.set_ylabel("First-break time  (ms)", color=c["label"])
    ax_tx.set_title("T-X picks + fitted arrivals", color=c["text"], fontsize=10)

    # Channel ticks on top of T-X panel
    if len(recv_positions) > 0:
        ax_tx_top = ax_tx.twiny()
        ax_tx_top.set_xlim(ax_tx.get_xlim())
        xs = np.asarray(recv_positions, dtype=float)
        step = max(1, len(xs) // 10)
        idxs = list(range(0, len(xs), step))
        if (len(xs) - 1) not in idxs:
            idxs.append(len(xs) - 1)
        ax_tx_top.set_xticks(xs[idxs])
        ax_tx_top.set_xticklabels([str(i + 1) for i in idxs], fontsize=7)
        ax_tx_top.set_xlabel("Channel", color=c["label"], fontsize=8)
        ax_tx_top.tick_params(colors=c["tick"], labelsize=7)
        for sp in ax_tx_top.spines.values():
            sp.set_edgecolor(c["spine"])

    # Velocity legend text block
    if layer_results:
        txt_lines: list = []
        for sid in sorted(layer_results):
            side_map = layer_results.get(sid, {}) or {}
            for side in [s for s in ("L", "R", "ALL") if s in side_map]:
                lr = side_map.get(side) or {}
                txt_lines.append(
                    f"S{sid}-{side}: V0={lr.get('V0_m_s',0):.0f}, "
                    f"V1={lr.get('V1_m_s',0):.0f}, V2={lr.get('V2_m_s',0):.0f}"
                )
        if txt_lines:
            ax_tx.text(0.01, 0.01, "\n".join(txt_lines[:10]), transform=ax_tx.transAxes,
                       fontsize=7, color=c["text"], va="bottom", ha="left",
                       bbox=dict(facecolor=c["ax_bg"], edgecolor=c["spine"], alpha=0.75))

    for ax in (ax_ch, ax_tx):
        ax.grid(True, lw=0.3, alpha=0.3, color=c["grid"])
        ax.tick_params(colors=c["tick"])
        for sp in ax.spines.values():
            sp.set_edgecolor(c["spine"])

    handles, labels = ax_ch.get_legend_handles_labels()
    if handles:
        ax_ch.legend(fontsize=8, facecolor=c["leg_face"],
                     edgecolor=c["leg_edge"], labelcolor=c["text"], loc="best")
    handles2, labels2 = ax_tx.get_legend_handles_labels()
    if handles2:
        ax_tx.legend(fontsize=8, facecolor=c["leg_face"],
                     edgecolor=c["leg_edge"], labelcolor=c["text"], loc="best")

    fig.suptitle(f"Profile {profile_name} - Picks, T-X, and fit overlays",
                 color=c["text"], fontsize=11)
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.11, top=0.90, wspace=0.12)

    out_dir = _project().ensure_dir(_project().plots_dir_for(profile_name))
    out = out_dir / f"{profile_name}_tx_picks.png"
    fig.savefig(str(out), dpi=180, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  T-X plot  -> {_project().relative(out)}")
    return out

def export_fit_plot(profile_name: str,
                    corrected_by_shot: dict,
                    layer_results: dict,
                    filename_suffix: str = "",
                    theme: str = "light") -> Path:
    """Save observed-vs-computed fit scatter with global RMS annotation."""
    c = theme_colors(theme)
    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    fig.patch.set_facecolor(c["fig_bg"])
    ax.set_facecolor(c["ax_bg"])

    pal = ["#e63946", "#2a9d8f", "#e9c46a", "#457b9d", "#f4a261", "#8d99ae"]
    all_obs: list = []
    all_pred: list = []

    for i, shot_id in enumerate(sorted(corrected_by_shot)):
        rows = corrected_by_shot.get(shot_id, [])
        side_map = (layer_results or {}).get(shot_id, {}) or {}
        fit_side, fit_res = _choose_shot_fit(side_map)
        if not rows or not fit_res:
            continue

        if fit_side == "ALL":
            rows_fit = rows
        else:
            rows_fit = [r for r in rows if r.get("side") == fit_side]
            if not rows_fit:
                rows_fit = rows

        obs: list = []
        pred: list = []
        for r in rows_fit:
            t_obs = float(r.get("fb_interp_inline_ms", 0.0))
            t_pred = _predict_time_from_fit(float(r.get("true_off_m", 0.0)), fit_res)
            if t_pred is None:
                continue
            obs.append(t_obs)
            pred.append(float(t_pred))

        if not obs:
            continue

        oa = np.asarray(obs, dtype=float)
        pa = np.asarray(pred, dtype=float)
        shot_rms = float(np.sqrt(np.mean((oa - pa) ** 2)))
        col = pal[i % len(pal)]
        ax.scatter(oa, pa, s=26, color=col, alpha=0.86,
                   label=f"Shot {shot_id} ({fit_side}) RMS={shot_rms:.2f} ms")
        all_obs.extend(obs)
        all_pred.extend(pred)

    if all_obs:
        obs_a = np.asarray(all_obs, dtype=float)
        pred_a = np.asarray(all_pred, dtype=float)
        rms_all = float(np.sqrt(np.mean((obs_a - pred_a) ** 2)))
        lo = float(min(obs_a.min(), pred_a.min()))
        hi = float(max(obs_a.max(), pred_a.max()))
        pad = max(1.0, 0.03 * (hi - lo))
        lo -= pad
        hi += pad
        ax.plot([lo, hi], [lo, hi], "--", color="#444444", lw=1.2, label="Ideal: y=x")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.text(0.02, 0.98, f"Global RMS = {rms_all:.3f} ms\nN = {len(obs_a)}",
                transform=ax.transAxes, va="top", ha="left",
                fontsize=9, color=c["text"])
    else:
        ax.text(0.5, 0.5, "No fitted points available", transform=ax.transAxes,
                ha="center", va="center", color=c["text"])

    ax.set_xlabel("Observed FB time (ms)", color=c["label"])
    ax.set_ylabel("Computed FB time (ms)", color=c["label"])
    ax.set_title(f"Observed vs Computed Fit - Profile {profile_name}",
                 color=c["text"], fontsize=11)
    ax.grid(True, lw=0.3, alpha=0.3, color=c["grid"])
    ax.tick_params(colors=c["tick"])
    for sp in ax.spines.values():
        sp.set_edgecolor(c["spine"])
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=8, facecolor=c["leg_face"], edgecolor=c["leg_edge"],
                  labelcolor=c["text"], loc="best")
    # tight_layout emits warnings with dense title/legend combinations.
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.12, top=0.92)

    out_dir = _project().ensure_dir(_project().plots_dir_for(profile_name))
    out = out_dir / f"{profile_name}_fit_rms{filename_suffix}.png"
    fig.savefig(str(out), dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Fit plot  -> {_project().relative(out)}")
    return out

def export_corrected_qc_plot(profile_name: str,
                             corrected_by_shot: dict,
                             shot_label_pos: dict | None = None,
                             layer_results: dict | None = None,
                             filename_suffix: str = "",
                             show_plot: bool = False,
                             theme: str = "light") -> Path:
    """QC visualization for corrected picks and true-offset mapping."""
    c = theme_colors(theme)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    fig.patch.set_facecolor(c["fig_bg"])
    ax1.set_facecolor(c["ax_bg"])
    ax2.set_facecolor(c["ax_bg"])

    pal = ["#e63946", "#2a9d8f", "#e9c46a", "#457b9d", "#f4a261", "#8d99ae"]

    for i, shot_id in enumerate(sorted(corrected_by_shot)):
        rows = corrected_by_shot.get(shot_id, [])
        if not rows:
            continue
        col = pal[i % len(pal)]
        label_pos = (shot_label_pos or {}).get(shot_id, 0.0)

        tr = [r["trace_no"] for r in rows]
        x_geom = [r["recv_pos_m"] for r in rows]
        t_interp = [r["fb_interp_inline_ms"] for r in rows]
        ax1.plot(tr, t_interp, "o-", ms=3, lw=1.0, color=col,
                 label=f"Shot {shot_id} (@{label_pos:.1f} m)")

        ax2.plot(x_geom, t_interp, "o-", ms=3, lw=1.0, color=col, alpha=0.9)

    # Secondary top axis on geometry-x panel: channels
    sample_rows: list = []
    for sid in sorted(corrected_by_shot):
        rr = corrected_by_shot.get(sid, [])
        if rr:
            sample_rows = rr
            break
    if sample_rows:
        ch_from_x = sorted([(float(r["recv_pos_m"]), int(r["trace_no"]))
                            for r in sample_rows], key=lambda p: p[0])
        if ch_from_x:
            x_vals = np.array([p[0] for p in ch_from_x], dtype=float)
            ch_vals = [p[1] for p in ch_from_x]
            ax2_top = ax2.twiny()
            ax2_top.set_xlim(ax2.get_xlim())
            step = max(1, len(x_vals) // 10)
            idxs = list(range(0, len(x_vals), step))
            if (len(x_vals) - 1) not in idxs:
                idxs.append(len(x_vals) - 1)
            ax2_top.set_xticks(x_vals[idxs])
            ax2_top.set_xticklabels([str(ch_vals[i]) for i in idxs], fontsize=7)
            ax2_top.set_xlabel("Channel", color=c["label"], fontsize=8)
            ax2_top.tick_params(colors=c["tick"], labelsize=7)
            for sp in ax2_top.spines.values():
                sp.set_edgecolor(c["spine"])

    ax1.set_xlabel("Channel", color=c["label"])
    ax1.set_ylabel("First-break time (ms)", color=c["label"])
    ax1.set_title("Interpolated FB vs Channel", color=c["text"], fontsize=10)

    ax2.set_xlabel("Geometry x (m)", color=c["label"])
    ax2.set_title("Interpolated FB vs geometry x", color=c["text"], fontsize=10)

    for ax in (ax1, ax2):
        ax.grid(True, lw=0.3, alpha=0.3, color=c["grid"])
        ax.tick_params(colors=c["tick"])
        for sp in ax.spines.values():
            sp.set_edgecolor(c["spine"])

    ax1.legend(fontsize=8, facecolor=c["leg_face"], edgecolor=c["leg_edge"],
               labelcolor=c["text"])
    ax1.set_ylim(T_MAX_MS, 0.0)
    ax2.set_ylim(T_MAX_MS, 0.0)

    if layer_results:
        txt_lines: list = []
        for sid in sorted(layer_results):
            shot_sides = layer_results.get(sid, {}) or {}
            ordered_sides = [s for s in ("L", "R", "ALL") if s in shot_sides]
            ordered_sides += [s for s in shot_sides if s not in ordered_sides]
            for side in ordered_sides:
                lr = shot_sides.get(side)
                if not lr:
                    continue
                txt_lines.append(
                    f"S{sid}-{side}: V0={lr.get('V0_m_s',0):.0f}, "
                    f"V1={lr.get('V1_m_s',0):.0f}, V2={lr.get('V2_m_s',0):.0f} m/s"
                )
        if txt_lines:
            ax2.text(0.01, 0.01, "\n".join(txt_lines[:8]), transform=ax2.transAxes,
                     fontsize=7, color=c["text"], va="bottom", ha="left",
                     bbox=dict(facecolor=c["ax_bg"], edgecolor=c["spine"], alpha=0.75))

    fig.suptitle(f"Corrected first-break QC  --  Profile {profile_name}",
                 color=c["text"], fontsize=11)
    fig.tight_layout()

    out_dir = _project().ensure_dir(_project().plots_dir_for(profile_name))
    out = out_dir / f"{profile_name}_corrected_qc{filename_suffix}.png"
    if show_plot:
        try:
            backend = str(plt.get_backend()).lower()
        except Exception:
            backend = ""
        if "agg" in backend:
            show_plot = _ensure_interactive_backend()
    if show_plot:
        try:
            plt.show(block=True)
        except Exception:
            pass
    fig.savefig(str(out), dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  corrected_qc -> {_project().relative(out)}")
    return out

def export_arrivals_observed_computed_plot(profile_name: str,
                                           corrected_by_shot: dict,
                                           layer_results: dict,
                                           filename_suffix: str = "",
                                           theme: str = "light") -> Path:
    """
    Plot observed and computed first arrivals together vs geometry x.

    This complements the observed-vs-computed scatter by keeping x (offset)
    explicit and overlaying both curves for each shot/selected side.
    """
    c = theme_colors(theme)
    fig, axs = plt.subplots(2, 2, figsize=(13.5, 9.0), constrained_layout=False)
    fig.patch.set_facecolor(c["fig_bg"])
    ax_geom = axs[0, 0]
    ax_chan = axs[0, 1]
    ax_scatter = axs[1, 0]
    ax_res = axs[1, 1]
    for ax in (ax_geom, ax_chan, ax_scatter, ax_res):
        ax.set_facecolor(c["ax_bg"])

    pal = ["#e63946", "#2a9d8f", "#e9c46a", "#457b9d", "#f4a261", "#8d99ae"]
    all_obs: list = []
    all_pred: list = []
    all_x_true: list = []
    all_chan: list = []

    for i, shot_id in enumerate(sorted(corrected_by_shot)):
        rows = list(corrected_by_shot.get(shot_id, []))
        side_map = (layer_results or {}).get(shot_id, {}) or {}
        if not rows or not side_map:
            continue

        ordered_sides = [s for s in ("L", "R", "ALL") if s in side_map]
        ordered_sides += [s for s in side_map if s not in ordered_sides]

        for side in ordered_sides:
            fit_res = side_map.get(side)
            if not fit_res:
                continue

            if side in ("L", "R"):
                rows_plot = [r for r in rows if r.get("side") == side]
                if not rows_plot:
                    continue
            else:
                rows_plot = rows

            x_obs: list = []
            ch_obs: list = []
            x_true_obs: list = []
            t_obs: list = []
            t_cmp: list = []
            for rr in rows_plot:
                x = float(rr.get("recv_pos_m", 0.0))
                x_abs = float(rr.get("true_off_m", abs(x)))
                ch = float(rr.get("trace_no", 0.0))
                tobs = float(rr.get("fb_interp_inline_ms", 0.0))
                tpred = _predict_time_from_fit(x_abs, fit_res)
                if tpred is None:
                    continue
                x_obs.append(x)
                x_true_obs.append(x_abs)
                ch_obs.append(ch)
                t_obs.append(tobs)
                t_cmp.append(float(tpred))
                all_obs.append(tobs)
                all_pred.append(float(tpred))
                all_x_true.append(x_abs)
                all_chan.append(ch)

            if not x_obs:
                continue

            order = np.argsort(np.asarray(x_obs, dtype=float))
            xs = np.asarray(x_obs, dtype=float)[order]
            to = np.asarray(t_obs, dtype=float)[order]
            tc = np.asarray(t_cmp, dtype=float)[order]
            col = pal[i % len(pal)]
            xc = np.asarray(ch_obs, dtype=float)[order]
            xt = np.asarray(x_true_obs, dtype=float)[order]
            ls = "--" if side != "R" else ":"
            ax_geom.plot(xs, to, "o", ms=3.8, color=col, alpha=0.95,
                     label=f"S{shot_id}-{side} obs")
            ax_geom.plot(xs, tc, ls, lw=1.3, color=col, alpha=0.9,
                     label=f"S{shot_id}-{side} comp")
            ax_chan.plot(xc, to, "o", ms=3.6, color=col, alpha=0.95)
            ax_chan.plot(xc, tc, ls, lw=1.2, color=col, alpha=0.9)

            # Residuals over true offset for this shot-side
            ax_res.plot(xt, (to - tc), ".", ms=5, color=col, alpha=0.85)

    if all_obs:
        oa = np.asarray(all_obs, dtype=float)
        pa = np.asarray(all_pred, dtype=float)
        rms = float(np.sqrt(np.mean((oa - pa) ** 2)))
        ax_geom.text(0.02, 0.98, f"Global RMS = {rms:.3f} ms\nN = {len(oa)}",
            transform=ax_geom.transAxes, va="top", ha="left",
                fontsize=9, color=c["text"])

        # Observed vs computed scatter panel
        ax_scatter.scatter(oa, pa, s=22, color="#2a9d8f", alpha=0.8)
        lo = float(min(oa.min(), pa.min()))
        hi = float(max(oa.max(), pa.max()))
        pad = max(1.0, 0.03 * (hi - lo))
        lo -= pad
        hi += pad
        ax_scatter.plot([lo, hi], [lo, hi], "--", color="#555555", lw=1.1)
        ax_scatter.set_xlim(lo, hi)
        ax_scatter.set_ylim(lo, hi)
    else:
        ax_geom.text(0.5, 0.5, "No observed/computed pairs", transform=ax_geom.transAxes,
                ha="center", va="center", color=c["text"])

    ax_geom.set_ylim(T_MAX_MS, 0.0)
    ax_chan.set_ylim(T_MAX_MS, 0.0)
    ax_geom.set_xlabel("Geometry x (m)", color=c["label"])
    ax_geom.set_ylabel("First-break time (ms)", color=c["label"])
    ax_geom.set_title("Observed + computed vs geometry", color=c["text"], fontsize=10)
    ax_chan.set_xlabel("Channel", color=c["label"])
    ax_chan.set_ylabel("First-break time (ms)", color=c["label"])
    ax_chan.set_title("Observed + computed vs channel", color=c["text"], fontsize=10)

    ax_scatter.set_xlabel("Observed (ms)", color=c["label"])
    ax_scatter.set_ylabel("Computed (ms)", color=c["label"])
    ax_scatter.set_title("Observed vs computed", color=c["text"], fontsize=10)

    ax_res.axhline(0.0, color="#666666", lw=1.0, ls="--", alpha=0.8)
    ax_res.set_xlabel("True offset XO (m)", color=c["label"])
    ax_res.set_ylabel("Residual (obs-comp) ms", color=c["label"])
    ax_res.set_title("Residuals over offset", color=c["text"], fontsize=10)

    for ax in (ax_geom, ax_chan, ax_scatter, ax_res):
        ax.grid(True, lw=0.3, alpha=0.3, color=c["grid"])
        ax.tick_params(colors=c["tick"])
        for sp in ax.spines.values():
            sp.set_edgecolor(c["spine"])

    handles, labels = ax_geom.get_legend_handles_labels()
    if handles:
        ax_geom.legend(fontsize=8, facecolor=c["leg_face"], edgecolor=c["leg_edge"],
                       labelcolor=c["text"], loc="best", ncol=2)

    # Optional top axis with channel ticks for easier interpretation.
    sample_rows: list = []
    for sid in sorted(corrected_by_shot):
        rr = corrected_by_shot.get(sid, [])
        if rr:
            sample_rows = rr
            break
    if sample_rows:
        pairs = sorted([(float(r.get("recv_pos_m", 0.0)), int(r.get("trace_no", 0)))
                        for r in sample_rows], key=lambda p: p[0])
        xs = np.array([p[0] for p in pairs], dtype=float)
        ch = [p[1] for p in pairs]
        if xs.size:
            ax_top = ax_geom.twiny()
            ax_top.set_xlim(ax_geom.get_xlim())
            step = max(1, len(xs) // 10)
            idxs = list(range(0, len(xs), step))
            if (len(xs) - 1) not in idxs:
                idxs.append(len(xs) - 1)
            ax_top.set_xticks(xs[idxs])
            ax_top.set_xticklabels([str(ch[i]) for i in idxs], fontsize=7)
            ax_top.set_xlabel("Channel", color=c["label"], fontsize=8)
            ax_top.tick_params(colors=c["tick"], labelsize=7)
            for sp in ax_top.spines.values():
                sp.set_edgecolor(c["spine"])

    fig.suptitle(f"Profile {profile_name} - Observed/computed diagnostics", color=c["text"], fontsize=11)
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.08, top=0.93, wspace=0.18, hspace=0.22)

    out_dir = _project().ensure_dir(_project().plots_dir_for(profile_name))
    out = out_dir / f"{profile_name}_arrivals_obs_comp{filename_suffix}.png"
    fig.savefig(str(out), dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  arrivals   -> {_project().relative(out)}")
    return out

def export_layer_fit_rms_plot(profile_name: str,
                              corrected_by_shot: dict,
                              layer_results: dict,
                              filename_suffix: str = "",
                              theme: str = "light") -> Path:
    """
    Plot layer-average computed arrivals and residual diagnostics over true offset.

    Observed picks are plotted as markers; layer lines and the minimum-time
    envelope are plotted as solid curves. RMS is reported globally and by
    controlling layer (which layer produced the minimum-time prediction).
    """
    c = theme_colors(theme)
    fig, (ax_tx, ax_res) = plt.subplots(1, 2, figsize=(13.0, 5.8), constrained_layout=False)
    fig.patch.set_facecolor(c["fig_bg"])
    ax_tx.set_facecolor(c["ax_bg"])
    ax_res.set_facecolor(c["ax_bg"])

    shot_pos_by_id: dict = {}
    for sid, rows in (corrected_by_shot or {}).items():
        if rows:
            shot_pos_by_id[int(sid)] = float(rows[0].get("shot_pos_m", sid))
    avg = compute_layer_averages(layer_results or {}, shot_pos_by_id)

    v0 = float((avg.get("V0", {}) or {}).get("avg", 0.0) or 0.0)
    v1 = float((avg.get("V1", {}) or {}).get("avg", 0.0) or 0.0)
    v2 = float((avg.get("V2", {}) or {}).get("avg", 0.0) or 0.0)
    ti1 = float((avg.get("ti1_ms", {}) or {}).get("avg", 0.0) or 0.0)
    ti2 = float((avg.get("ti2_ms", {}) or {}).get("avg", 0.0) or 0.0)

    # Collect observed points across all shots.
    x_obs: list = []
    t_obs: list = []
    side_obs: list = []
    shot_obs: list = []
    shot_pos_sorted = sorted(
        ((int(sid), float(rows[0].get("shot_pos_m", sid)))
         for sid, rows in (corrected_by_shot or {}).items() if rows),
        key=lambda p: p[1],
    )
    shot_order = [sid for sid, _ in shot_pos_sorted]
    shot_rank = {sid: i for i, sid in enumerate(shot_order)}
    n_shots = len(shot_order)

    for sid, rows in (corrected_by_shot or {}).items():
        for r in rows:
            x_obs.append(float(r.get("true_off_m", 0.0)))
            t_obs.append(float(r.get("fb_interp_inline_ms", 0.0)))
            side_obs.append(str(r.get("side", "R")))
            shot_obs.append(int(sid))

    if not x_obs:
        ax_tx.text(0.5, 0.5, "No corrected picks", transform=ax_tx.transAxes,
                   ha="center", va="center", color=c["text"])
    else:
        xa = np.asarray(x_obs, dtype=float)
        ta = np.asarray(t_obs, dtype=float)
        sid_a = np.asarray(shot_obs, dtype=int)
        side_a = np.asarray(side_obs, dtype=object)
        order = np.argsort(xa)
        xa = xa[order]
        ta = ta[order]
        sid_a = sid_a[order]
        side_a = side_a[order]

        x_grid = np.linspace(0.0, float(np.max(xa)), 320)

        lines: list = []  # (name, t(x), color)
        if v0 > 0.0:
            m0 = 1000.0 / v0
            lines.append(("L1", m0 * x_grid, "#e63946"))
        if v1 > 0.0:
            m1 = 1000.0 / v1
            lines.append(("L2", m1 * x_grid + ti1, "#2a9d8f"))
        if v2 > 0.0:
            m2 = 1000.0 / v2
            lines.append(("L3", m2 * x_grid + ti2, "#457b9d"))

        if lines:
            mats = np.vstack([tt for _nm, tt, _cc in lines])
            env = np.min(mats, axis=0)

            # Prediction and residual assignment for observed points.
            pred: list = []
            lay_idx: list = []
            for xv in xa:
                cand: list = []
                if v0 > 0.0:
                    cand.append((0, (1000.0 / v0) * xv))
                if v1 > 0.0:
                    cand.append((1, (1000.0 / v1) * xv + ti1))
                if v2 > 0.0:
                    cand.append((2, (1000.0 / v2) * xv + ti2))
                if not cand:
                    pred.append(float("nan"))
                    lay_idx.append(-1)
                    continue
                idx, tp = min(cand, key=lambda p: p[1])
                pred.append(float(tp))
                lay_idx.append(int(idx))

            pa = np.asarray(pred, dtype=float)
            mask = np.isfinite(pa)
            if np.any(mask):
                rms_all = float(np.sqrt(np.mean((ta[mask] - pa[mask]) ** 2)))
            else:
                rms_all = 0.0

            def _obs_group_label(sid: int, side: str) -> str:
                if n_shots <= 0:
                    return "Observed"
                rnk = shot_rank.get(int(sid), 0)
                side_s = "-" if str(side).upper() == "L" else "+"
                if n_shots == 1:
                    return f"Centre{side_s}"
                if n_shots == 2:
                    return "Off-end+" if side_s == "+" else "Off-end-"
                if n_shots == 3:
                    if rnk == 1:
                        return f"Centre{side_s}"
                    return "Off-end+" if side_s == "+" else "Off-end-"

                # Optional far groups for >3 shots.
                c_mid = (n_shots - 1) / 2.0
                d = rnk - c_mid
                mag = abs(d)
                if mag <= 0.5:
                    return f"Centre{side_s}"
                if mag <= 1.5:
                    return "Off-end+" if side_s == "+" else "Off-end-"
                if mag <= 2.5:
                    return "Far+" if d > 0 else "Far-"
                return "Far++" if d > 0 else "Far--"

            grp_colors = {
                "Off-end-": "#6d597a",
                "Centre-": "#457b9d",
                "Centre+": "#2a9d8f",
                "Off-end+": "#f4a261",
                "Far-": "#8d99ae",
                "Far+": "#e9c46a",
                "Far--": "#5c677d",
                "Far++": "#b08968",
                "Observed": "#111111",
            }

            grp_xy: dict = {}
            for xv, tv, sid, sside in zip(xa, ta, sid_a, side_a):
                g = _obs_group_label(int(sid), str(sside))
                grp_xy.setdefault(g, [[], []])
                grp_xy[g][0].append(float(xv))
                grp_xy[g][1].append(float(tv))

            grp_order = ["Off-end-", "Centre-", "Centre+", "Off-end+", "Far-", "Far+", "Far--", "Far++", "Observed"]
            for g in grp_order:
                if g not in grp_xy:
                    continue
                gx, gt = grp_xy[g]
                ax_tx.plot(gx, gt, "o", ms=3.5, alpha=0.85,
                           color=grp_colors.get(g, "#111111"), label=g)

            for nm, tt, cc in lines:
                ax_tx.plot(x_grid, tt, "--", lw=1.4, color=cc, alpha=0.9, label=f"{nm} computed")
            ax_tx.plot(x_grid, env, "-", lw=2.0, color="#000000", alpha=0.85,
                       label=f"Envelope (RMS={rms_all:.3f} ms)")

            # Residuals and per-layer RMS.
            if np.any(mask):
                res = ta[mask] - pa[mask]
                x_m = xa[mask]
                idx_m = np.asarray(lay_idx, dtype=int)[mask]
                cols = ["#e63946", "#2a9d8f", "#457b9d"]
                for li in (0, 1, 2):
                    mli = (idx_m == li)
                    if not np.any(mli):
                        continue
                    rms_li = float(np.sqrt(np.mean((res[mli]) ** 2)))
                    ax_res.plot(x_m[mli], res[mli], ".", ms=5, color=cols[li], alpha=0.85,
                                label=f"Layer {li+1} RMS={rms_li:.3f} ms")
                ax_res.axhline(0.0, color="#666666", lw=1.0, ls="--")
                ax_res.text(
                    0.02, 0.02,
                    f"Global RMS: {rms_all:.3f} ms\nN: {int(np.sum(mask))}",
                    transform=ax_res.transAxes,
                    va="bottom", ha="left",
                    fontsize=9,
                    color=c["text"],
                    bbox=dict(boxstyle="round,pad=0.35",
                              facecolor=c["ax_bg"], edgecolor=c["spine"], alpha=0.90),
                )

            # Compact key:value summary in a padded lower-left box.
            h0 = float(avg.get("h1_m", 0.0) or 0.0)
            h1v = float(avg.get("h2_m", 0.0) or 0.0)
            summary_txt = (
                "Layer model\n"
                f"V0: {v0:.1f} m/s\n"
                f"V1: {v1:.1f} m/s\n"
                f"V2: {v2:.1f} m/s\n"
                f"ti1: {ti1:.2f} ms\n"
                f"ti2: {ti2:.2f} ms\n"
                f"h0: {h0:.2f} m\n"
                f"h1: {h1v:.2f} m"
            )
            ax_tx.text(
                0.02, 0.02,
                summary_txt,
                transform=ax_tx.transAxes,
                va="bottom", ha="left",
                fontsize=8,
                color=c["text"],
                bbox=dict(boxstyle="round,pad=0.4",
                          facecolor=c["ax_bg"], edgecolor=c["spine"], alpha=0.92),
            )
        else:
            ax_tx.plot(xa, ta, "o", ms=3.3, color="#111111", alpha=0.75, label="Observed")
            ax_tx.text(0.5, 0.1, "No valid average layer velocities", transform=ax_tx.transAxes,
                       ha="center", va="center", color=c["label"], fontsize=8)

    ax_tx.set_ylim(T_MAX_MS, 0.0)
    ax_tx.set_xlabel("True offset XO (m)", color=c["label"])
    ax_tx.set_ylabel("Arrival time (ms)", color=c["label"])
    ax_tx.set_title("Average-layer fit over offset", color=c["text"], fontsize=10)

    ax_res.set_xlabel("True offset XO (m)", color=c["label"])
    ax_res.set_ylabel("Residual (obs-comp) ms", color=c["label"])
    ax_res.set_title("Residuals by controlling layer", color=c["text"], fontsize=10)

    for ax in (ax_tx, ax_res):
        ax.grid(True, lw=0.3, alpha=0.3, color=c["grid"])
        ax.tick_params(colors=c["tick"])
        for sp in ax.spines.values():
            sp.set_edgecolor(c["spine"])

    h1, _ = ax_tx.get_legend_handles_labels()
    if h1:
        ax_tx.legend(fontsize=8, facecolor=c["leg_face"], edgecolor=c["leg_edge"],
                     labelcolor=c["text"], loc="best")
    h2, _ = ax_res.get_legend_handles_labels()
    if h2:
        ax_res.legend(fontsize=8, facecolor=c["leg_face"], edgecolor=c["leg_edge"],
                      labelcolor=c["text"], loc="best")

    fig.suptitle(f"Profile {profile_name} - Layer-fit RMS diagnostics", color=c["text"], fontsize=11)
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.11, top=0.90, wspace=0.20)

    out_dir = _project().ensure_dir(_project().plots_dir_for(profile_name))
    out = out_dir / f"{profile_name}_layer_fit_rms{filename_suffix}.png"
    fig.savefig(str(out), dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  layer_fit -> {_project().relative(out)}")
    return out


def export_shot_qc_plot(
    profile_name: str,
    shot_id: int,
    data: Any,
    dt_s: float,
    delay_ms: float,
    picks: dict,
    recv_positions: Any = None,
    shot_pos_m: float | None = None,
    theme: str = "light",
    filename_suffix: str = "",
    display_mode: str = "both",
    clip_factor: float = 2.0,
    wiggle_scale: float = 1.0,
    cmap: str = "gray",
) -> Path:
    """Save a per-shot QC plot with picks marked - proof, for anyone
    reviewing the results later, that a given shot's picks were placed on
    genuine arrivals and not somewhere spurious.

    Uses the SAME variable-density + wiggle rendering convention as the
    interactive `SeismicDisplay` widget in lvl_studio.py (robust per-trace
    std normalization with median clipping, positive-lobe fill, VD
    background at the 98th-percentile amplitude scale) - not a simplified
    stand-in, so the saved plot actually looks like what you saw on
    screen while picking, not a generic line plot.

    This is the plot the legacy interactive v1 picker used to save
    automatically as part of its own matplotlib window (`{profile}_shot
    {id:02d}_picks.png`); the PyQt desktop picking view and Picker V2
    auto-pick don't go through that window, so nothing was producing it
    for those workflows - this closes that gap, in the same naming
    convention, callable from anywhere picks exist for a shot.

    Parameters
    ----------
    data : Any
        (n_traces, n_samples) array for this shot.
    picks : dict
        {trace_index (1-based): pick_time_ms}.
    display_mode : str
        "wiggle", "vd", or "both" - matches the GUI's own display-mode setting.
    """
    c = theme_colors(theme)
    x = np.asarray(data, dtype=np.float64)
    n_traces, n_samples = x.shape
    t_ms = delay_ms + np.arange(n_samples) * dt_s * 1000.0
    xv = np.arange(1, n_traces + 1, dtype=float)  # trace-index x positions
    dx = 1.0  # spacing between adjacent trace positions

    fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=False)
    fig.patch.set_facecolor(c["fig_bg"])
    ax.set_facecolor(c["ax_bg"])

    mode = str(display_mode).strip().lower()

    if mode in ("vd", "both"):
        vmax = float(np.percentile(np.abs(x), 98)) if x.size else 1.0
        if vmax < 1e-12:
            vmax = 1.0
        ax.imshow(
            x.T, cmap=cmap, aspect="auto", interpolation="nearest",
            vmin=-vmax, vmax=vmax,
            extent=[xv[0] - dx / 2, xv[-1] + dx / 2, float(t_ms[-1]), float(t_ms[0])],
            alpha=0.45 if mode == "both" else 1.0, zorder=1,
        )

    if mode in ("wiggle", "both"):
        # Robust per-trace normalization: clip each trace's std to
        # [0.3, 3.0] x the profile's median std, so one dead/noisy trace
        # doesn't wash out or dwarf its neighbours - matches the GUI.
        stds = x.std(axis=1).astype(float)
        valid = stds[stds > 1e-20]
        med = float(np.median(valid)) if valid.size else 1.0
        norms = np.clip(stds, med * 0.3, med * 3.0)
        norms = np.where(norms > 1e-20, norms, med)
        scale = norms * float(clip_factor)
        defl = dx * float(wiggle_scale)
        trace_col = c["tick"]
        for i in range(n_traces):
            base_x = float(xv[i])
            tr_n = np.clip(x[i] / scale[i], -1.0, 1.0)
            wig_x = base_x + tr_n * defl
            ax.plot(wig_x, t_ms, color=trace_col, linewidth=0.5, alpha=0.85, zorder=4)
            pos = np.where(tr_n > 0.0, tr_n, 0.0)
            ax.fill_betweenx(t_ms, base_x, base_x + pos * defl, color=trace_col, alpha=0.18, zorder=3)

    if picks:
        pick_x = sorted(int(k) for k in picks.keys())
        pick_y = [float(picks[k]) for k in pick_x]
        ax.plot(pick_x, pick_y, "o-", color="#e63946", ms=4, lw=1.2, zorder=5, label="Picks")
        ax.legend(loc="upper right", fontsize=8, facecolor=c["ax_bg"], labelcolor=c["text"])

    title = f"Profile {profile_name} - Shot {shot_id}"
    if shot_pos_m is not None:
        title += f" (SP={float(shot_pos_m):.1f} m)"
    ax.set_title(title, color=c["text"], fontsize=11)
    ax.set_xlabel("Trace index", color=c["label"])
    ax.set_ylabel("Time (ms)", color=c["label"])
    ax.set_xlim(xv[0] - dx, xv[-1] + dx)
    ax.invert_yaxis()
    ax.grid(True, lw=0.3, alpha=0.3, color=c["grid"])
    ax.tick_params(colors=c["tick"])
    for sp in ax.spines.values():
        sp.set_edgecolor(c["spine"])

    out_dir = plots_directory(profile_name)
    out_path = out_dir / f"{profile_name}_shot{int(shot_id):02d}_picks{filename_suffix}.png"
    fig.tight_layout()
    fig.savefig(str(out_path), dpi=150, facecolor=c["fig_bg"])
    plt.close(fig)
    print(f"  QC picks plot -> {_project().relative(out_path)}")
    return out_path