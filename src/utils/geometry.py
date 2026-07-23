"""Geometry, profile discovery, and field-report/coordinate Excel parsing utilities."""
from __future__ import annotations

import csv
import datetime
import math
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd

from src.common.settings import GEOM_FILES
from src.io.seg2_reader import read_seg2


def load_geometry(geom_type: int) -> Any:
    """Return receiver positions (m) as a 1-D ndarray."""
    path = GEOM_FILES[geom_type]
    positions: list = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.lower().startswith("trace"):
                continue
            # Handle German decimal comma:  "2,5" -> 2.5
            parts = s.replace(",", ".").split()
            if len(parts) >= 2:
                try:
                    positions.append(float(parts[1]))
                except ValueError:
                    pass
    return np.array(positions, dtype=float)

def auto_shot_positions(recv_pos: Any) -> dict:
    """
    Compute standard refraction shot positions from the receiver geometry.

    Convention used by most land refraction surveys:
      Shot 1  :  one receiver-spacing BEFORE the first receiver
      Shot 2  :  midpoint between receivers n//2 and n//2+1
                 (labelled e.g. GP24.5 for a 48-channel spread)
      Shot 3  :  one receiver-spacing AFTER the last receiver

    Returns {1: pos_m, 2: pos_m, 3: pos_m}.
    """
    n          = len(recv_pos)
    dx_start   = float(recv_pos[1]    - recv_pos[0])    # first spacing
    dx_end     = float(recv_pos[-1]   - recv_pos[-2])   # last spacing
    mid_lo     = float(recv_pos[n // 2 - 1])
    mid_hi     = float(recv_pos[n // 2])
    return {
        1: round(float(recv_pos[0]) - dx_start, 4),
        2: round((mid_lo + mid_hi) / 2.0,        4),
        3: round(float(recv_pos[-1]) + dx_end,   4),
    }

def _excel_col_to_index(col: str) -> int:
    """Convert Excel column label (A, D, AA, ...) to zero-based index."""
    s = str(col or "").strip().upper()
    if not s or not all("A" <= ch <= "Z" for ch in s):
        raise ValueError(f"Invalid Excel column label '{col}'")
    out = 0
    for ch in s:
        out = out * 26 + (ord(ch) - ord("A") + 1)
    return out - 1

def _read_geometry_table(path: Path) -> dict:
    """Read a geometry/coordinate table regardless of file format.

    Returns {sheet_name: DataFrame} for a uniform interface: Excel files
    may have several sheets; a plain-text file (.txt/.dat/.csv) is always
    treated as a single sheet with its delimiter auto-detected (comma,
    semicolon, tab, or whitespace).

    This is what lets a manually-uploaded geometry file (any of
    .xlsx/.xls/.csv/.txt/.dat) be searched the exact same way as an
    auto-discovered "LVL*.xlsx" survey table - same profile/station/X/Y/Z
    column-detection logic in `_find_station_xyz_columns` either way.
    """
    suffix = Path(path).suffix.lower()

    if suffix in (".xlsx", ".xlsm", ".xls"):
        try:
            return pd.read_excel(path, sheet_name=None, header=None, dtype=object)
        except Exception:
            engine = "xlrd" if suffix == ".xls" else "openpyxl"
            return pd.read_excel(path, sheet_name=None, header=None, dtype=object, engine=engine)

    # Plain-text table: sniff the delimiter rather than assuming one.
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    sample = "\n".join(text.splitlines()[:20])
    delimiter = ","
    try:
        delimiter = csv.Sniffer().sniff(sample, delimiters=",;\t ").delimiter
    except Exception:
        if "\t" in sample:
            delimiter = "\t"
        elif ";" in sample:
            delimiter = ";"
        elif "," in sample:
            delimiter = ","
        else:
            delimiter = r"\s+"

    df = pd.read_csv(
        path, sep=delimiter, header=None, dtype=object, engine="python",
        skip_blank_lines=True,
    )
    return {"Sheet1": df}


def _normalize_profile_token(value: Any) -> str:
    """
    Normalize profile strings for robust matching.

    Examples:
      "LVL150" -> "150"
      "150" -> "150"
      "LVL214_A" -> "214_A"
    """
    s = str(value or "").strip().upper().replace(" ", "")
    s = s.replace("_", "").replace("-", "")
    if s.startswith("LVL"):
        s = s[3:]
    return s

def discover_field_report_excels(data_dir: Path) -> list:
    """Find likely field-report Excel files in data directory (newest first)."""
    pats = [
        "field_report*.xls",
        "field_report*.xlsx",
        "fieldreport*.xls",
        "fieldreport*.xlsx",
    ]
    seen: set[str] = set()
    out: list = []
    for pat in pats:
        for p in data_dir.glob(pat):
            key = str(p.resolve()).lower()
            if key in seen:
                continue
            seen.add(key)
            out.append(p)
    out.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return out

def discover_profile_folders(data_dir: Path) -> list:
    """Return all profile folder names under data_dir that contain SEG2 files."""
    out: list = []
    for p in sorted(data_dir.iterdir() if data_dir.exists() else []):
        if not p.is_dir():
            continue
        has_seg2 = any(p.glob("*.seg2")) or any(p.glob("*.SEG2"))
        if has_seg2:
            out.append(p.name)
    return out

def infer_geometry_from_spread_length(length_m: float) -> int | None:
    """Infer geometry type (100/200) from observed spread length in meters."""
    try:
        obs = float(length_m)
    except Exception:
        return None
    if not np.isfinite(obs) or obs <= 0.0:
        return None

    refs = {
        100: float(abs(load_geometry(100)[-1] - load_geometry(100)[0])),
        200: float(abs(load_geometry(200)[-1] - load_geometry(200)[0])),
    }
    best_geom = min(refs, key=lambda g: abs(obs - refs[g]))
    rel_err = abs(obs - refs[best_geom]) / max(refs[best_geom], 1e-6)
    if rel_err <= 0.15:
        return int(best_geom)
    return None

def infer_geometry_from_seg2_file(seg2_path: Path) -> tuple[int | None, float | None]:
    """Infer geometry from SEG2 receiver locations in one shot file."""
    try:
        _data, _dt, _ntr, _ns, _sp, _ffid, _delay, recv_locs_m = read_seg2(seg2_path)
    except Exception:
        return None, None
    if recv_locs_m is None or len(recv_locs_m) < 2:
        return None, None
    length_m = float(abs(float(np.max(recv_locs_m)) - float(np.min(recv_locs_m))))
    geom = infer_geometry_from_spread_length(length_m)
    return geom, length_m

def _find_station_xyz_columns(df: Any) -> tuple | None:
    """Detect header row and columns for profile/station/X/Y/Z in a free-form table."""
    if df is None or df.empty:
        return None

    def _norm(v: Any) -> str:
        s = str(v or "").strip().lower().replace("_", " ").replace("-", " ")
        return " ".join(s.split())

    max_scan = min(30, int(df.shape[0]))
    for ridx in range(max_scan):
        row = [_norm(v) for v in df.iloc[ridx].tolist()]

        i_profile = None
        i_station = None
        i_x = None
        i_y = None
        i_z = None

        for ci, txt in enumerate(row):
            if not txt:
                continue
            if i_profile is None and (("lvl" in txt and "number" in txt) or txt in ("lvl", "profile", "line")):
                i_profile = ci
            if i_station is None and "station" in txt:
                i_station = ci
            if i_x is None and txt in ("x", "east", "easting"):
                i_x = ci
            if i_y is None and txt in ("y", "north", "northing"):
                i_y = ci
            if i_z is None and (txt == "z" or "height" in txt or "elevation" in txt):
                i_z = ci

        if None not in (i_profile, i_station, i_x, i_y, i_z):
            return ridx, i_profile, i_station, i_x, i_y, i_z
    return None

def _to_float_or_none(v: Any) -> float | None:
    """Convert numeric-like cell values to float, else None."""
    try:
        if pd.isna(v):
            return None
    except Exception:
        pass
    s = str(v).strip().replace(",", ".")
    if s == "":
        return None
    try:
        return float(s)
    except Exception:
        return None

def load_midpoint_xyz_from_geometry_excels(profile_name: str,
                                           station_mid: float,
                                           excel_paths: list) -> tuple:
    """
    Read midpoint X/Y/Z from LVL geometry spreadsheets.

    Repeated station entries are resolved by taking the last matching row.
    Returns (x, y, z, source_file) where each value may be None.
    """
    target = _normalize_profile_token(profile_name)
    best_xyz = (None, None, None)
    best_file = None
    best_dist = float("inf")

    for path in excel_paths:
        try:
            sheets = _read_geometry_table(path)
        except Exception:
            continue
        for _sheet_name, df in sheets.items():
            cols = _find_station_xyz_columns(df)
            if cols is None:
                continue
            hdr_row, i_prof, i_sta, i_x, i_y, i_z = cols

            by_station: dict = {}
            for ridx in range(hdr_row + 1, int(df.shape[0])):
                p_raw = df.iat[ridx, i_prof] if i_prof < df.shape[1] else None
                s_raw = df.iat[ridx, i_sta] if i_sta < df.shape[1] else None
                x_raw = df.iat[ridx, i_x] if i_x < df.shape[1] else None
                y_raw = df.iat[ridx, i_y] if i_y < df.shape[1] else None
                z_raw = df.iat[ridx, i_z] if i_z < df.shape[1] else None

                if _normalize_profile_token(p_raw) != target:
                    continue
                sta = _to_float_or_none(s_raw)
                x = _to_float_or_none(x_raw)
                y = _to_float_or_none(y_raw)
                z = _to_float_or_none(z_raw)
                if sta is None or x is None or y is None or z is None:
                    continue

                # Last row wins for repeated station values.
                by_station[float(sta)] = (float(x), float(y), float(z))

            if not by_station:
                continue

            for sta, xyz in by_station.items():
                dist = abs(float(sta) - float(station_mid))
                if dist <= best_dist:
                    best_dist = dist
                    best_xyz = xyz
                    best_file = path

    return best_xyz[0], best_xyz[1], best_xyz[2], best_file

def load_profile_geometry_from_excels(profile_name: str,
                                      excel_paths: list) -> tuple:
    """
    Load profile geometry summary from LVL geometry spreadsheets.

    Returns
    -------
    (length_m, x_mid, y_mid, z_mid, source_file, station_min, station_mid, station_max, acq_datetime_de)

    Notes
    -----
    - Repeated station entries use the last row value.
    - Length is computed along sorted station order using XY polyline distance.
    - Midpoint coordinate is taken at station nearest to the midpoint station.
    """
    target = _normalize_profile_token(profile_name)
    best_by_station: dict = {}
    best_file = None
    best_acq_dt_de = None
    best_span = -1.0
    best_count = -1

    def _parse_excel_datetime(v: Any) -> datetime.datetime | None:
        if v is None:
            return None
        if isinstance(v, datetime.datetime):
            return v
        if isinstance(v, datetime.date):
            return datetime.datetime(v.year, v.month, v.day, 0, 0, 0)
        s = str(v).strip()
        if not s:
            return None
        for fmt in (
            "%d.%m.%Y %H:%M:%S", "%d.%m.%Y %H:%M",
            "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M",
            "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M",
            "%d.%m.%Y", "%d/%m/%Y", "%Y-%m-%d",
        ):
            try:
                return datetime.datetime.strptime(s, fmt)
            except Exception:
                continue
        return None

    def _extract_sheet_acq_datetime_de(df: Any, hdr_row: int) -> str | None:
        tz = ZoneInfo("Europe/Berlin")
        max_rows = min(int(df.shape[0]), max(40, hdr_row + 8))
        max_cols = int(df.shape[1])

        def _to_local_text(dt_val: datetime.datetime) -> str:
            if dt_val.tzinfo is None:
                dt_local = dt_val.replace(tzinfo=tz)
            else:
                dt_local = dt_val.astimezone(tz)
            return dt_local.strftime("%Y-%m-%d %H:%M:%S %Z")

        labels = ("datum", "date", "acquisition")
        for r in range(max_rows):
            for c in range(max_cols):
                txt = str(df.iat[r, c] if c < df.shape[1] else "").strip().lower()
                if not txt:
                    continue
                if any(lbl in txt for lbl in labels):
                    cand_vals = []
                    if c + 1 < max_cols:
                        cand_vals.append(df.iat[r, c + 1])
                    if c + 2 < max_cols:
                        cand_vals.append(df.iat[r, c + 2])
                    if r + 1 < max_rows:
                        cand_vals.append(df.iat[r + 1, c])
                    for cv in cand_vals:
                        dtv = _parse_excel_datetime(cv)
                        if dtv is not None:
                            return _to_local_text(dtv)

        for r in range(max_rows):
            for c in range(max_cols):
                dtv = _parse_excel_datetime(df.iat[r, c])
                if dtv is not None:
                    return _to_local_text(dtv)
        return None

    for path in excel_paths:
        try:
            sheets = _read_geometry_table(path)
        except Exception:
            continue

        for _sheet_name, df in sheets.items():
            cols = _find_station_xyz_columns(df)
            if cols is None:
                continue
            hdr_row, i_prof, i_sta, i_x, i_y, i_z = cols
            acq_dt_de = _extract_sheet_acq_datetime_de(df, hdr_row)

            by_station: dict = {}
            for ridx in range(hdr_row + 1, int(df.shape[0])):
                p_raw = df.iat[ridx, i_prof] if i_prof < df.shape[1] else None
                if _normalize_profile_token(p_raw) != target:
                    continue

                s_raw = df.iat[ridx, i_sta] if i_sta < df.shape[1] else None
                x_raw = df.iat[ridx, i_x] if i_x < df.shape[1] else None
                y_raw = df.iat[ridx, i_y] if i_y < df.shape[1] else None
                z_raw = df.iat[ridx, i_z] if i_z < df.shape[1] else None

                sta = _to_float_or_none(s_raw)
                x = _to_float_or_none(x_raw)
                y = _to_float_or_none(y_raw)
                z = _to_float_or_none(z_raw)
                if sta is None or x is None or y is None or z is None:
                    continue

                # Last row wins for repeated station values.
                by_station[float(sta)] = (float(x), float(y), float(z))

            if not by_station:
                continue

            st_sorted = sorted(by_station)
            span = float(st_sorted[-1] - st_sorted[0]) if len(st_sorted) >= 2 else 0.0
            cnt = int(len(st_sorted))
            if (span > best_span) or (math.isclose(span, best_span) and cnt > best_count):
                best_span = span
                best_count = cnt
                best_by_station = by_station
                best_file = path
                best_acq_dt_de = acq_dt_de

    if not best_by_station:
        return None, None, None, None, None, None, None, None, None

    stations = sorted(best_by_station)
    s_min = float(stations[0])
    s_max = float(stations[-1])
    s_mid = 0.5 * (s_min + s_max)

    s_near = min(stations, key=lambda s: abs(float(s) - s_mid))
    x_mid, y_mid, z_mid = best_by_station[float(s_near)]

    length_m = 0.0
    for i in range(1, len(stations)):
        x0, y0, _z0 = best_by_station[float(stations[i - 1])]
        x1, y1, _z1 = best_by_station[float(stations[i])]
        dx = float(x1) - float(x0)
        dy = float(y1) - float(y0)
        length_m += float(math.sqrt(dx * dx + dy * dy))

    return float(length_m), float(x_mid), float(y_mid), float(z_mid), best_file, s_min, s_mid, s_max, best_acq_dt_de

def _read_excel_table(path: Path, sheet_name: str | int | None = None) -> Any:
    """Read Excel with fallback engines for both .xls and .xlsx files."""
    sheet = (0 if sheet_name is None else sheet_name)

    # First try pandas default engine resolution.
    try:
        return pd.read_excel(path, sheet_name=sheet, header=None, dtype=object)
    except Exception:
        pass

    suffix = str(path.suffix).lower()
    last_exc: Exception | None = None

    if suffix == ".xls":
        try:
            return pd.read_excel(path, sheet_name=sheet, header=None, dtype=object, engine="xlrd")
        except Exception as exc:
            last_exc = exc
    else:
        try:
            return pd.read_excel(path, sheet_name=sheet, header=None, dtype=object, engine="openpyxl")
        except Exception as exc:
            last_exc = exc

    if last_exc is not None:
        raise last_exc
    raise ValueError(f"Could not read Excel file: {path}")

def infer_geometry_from_field_report(path: Path,
                                     profile_name: str,
                                     sheet_name: str | int | None = None) -> int | None:
    """
    Infer geometry from field report notes near profile markers.

        Rule:
            - find row containing LVL<profile>
            - inspect nearby rows (+/- 1..2)
      - if any cell contains 'short' (or 'short profile'), geometry = 100
      - otherwise return None
    """
    if not path.exists():
        return None

    df = _read_excel_table(path, sheet_name=sheet_name)
    if df is None or df.empty:
        return None

    prof_tok = _normalize_profile_token(profile_name)
    lvl_tag = f"LVL{prof_tok}"
    n_rows = int(df.shape[0])

    def _row_text(ridx: int) -> str:
        vals = [str(v) for v in df.iloc[ridx].tolist() if pd.notna(v)]
        return " ".join(vals).strip().lower()

    for ridx in range(n_rows):
        row_txt_up = _row_text(ridx).upper().replace(" ", "")
        if lvl_tag not in row_txt_up:
            continue

        for delta in (-2, -1, 1, 2):
            rr = ridx + delta
            if rr >= n_rows:
                continue
            if rr < 0:
                continue
            txt = _row_text(rr)
            if "short profile" in txt or "short" in txt:
                return 100
    return None

def load_profile_offsets_from_excel(
    path: Path,
    profile_name: str,
    ffid_by_shot: dict | None = None,
    sheet_name: str | int | None = None,
    ffid_col: str = "A",
    perp_col: str = "D",
    profile_col: str = "F",
    inline_shift_col: str | None = None,
) -> tuple:
    """
    Load per-shot perpendicular offsets (and optional inline shift) from Excel.

        Expected table convention:
            - vertically organized profile blocks
            - first row contains LVLXXX in profile column
            - continuation rows follow until next LVL row
            - shot numbering often lives in column A as 1/2/3

    Returns
    -------
    (perp_by_shot, inline_shift_by_shot)
    """
    if not path.exists():
        raise FileNotFoundError(f"Excel file not found: {path}")

    df = _read_excel_table(path, sheet_name=sheet_name)
    if df.empty:
        return {}, {}

    i_ffid = _excel_col_to_index(ffid_col)
    i_perp = _excel_col_to_index(perp_col)
    i_prof = _excel_col_to_index(profile_col)
    i_shift = _excel_col_to_index(inline_shift_col) if inline_shift_col else None
    need_cols = [i_ffid, i_perp, i_prof] + ([i_shift] if i_shift is not None else [])
    if max(need_cols) >= df.shape[1]:
        raise ValueError(
            f"Excel has {df.shape[1]} columns; requested column index {max(need_cols)+1} is out of range"
        )

    target = _normalize_profile_token(profile_name)
    ffid_to_shots: dict = {}
    for sid in sorted(ffid_by_shot.keys() if ffid_by_shot else []):
        try:
            ff = int((ffid_by_shot or {}).get(sid, 0))
            ss = int(sid)
        except Exception:
            continue
        if ff <= 0:
            continue
        ffid_to_shots.setdefault(ff, []).append(ss)

    rows_target: list = []
    n_rows = int(df.shape[0])
    lvl_starts: list = []
    target_starts: list = []
    for ridx in range(n_rows):
        if i_prof >= df.shape[1]:
            continue
        p_raw = df.iat[ridx, i_prof]
        p_txt = "" if pd.isna(p_raw) else str(p_raw).strip()
        p_up = p_txt.upper().replace(" ", "")
        if "LVL" not in p_up:
            continue
        lvl_starts.append(ridx)
        if _normalize_profile_token(p_txt) == target:
            target_starts.append(ridx)

    for start in target_starts:
        next_starts = [r for r in lvl_starts if r > start]
        end = min(next_starts) if next_starts else n_rows
        for ridx in range(start, end):
            ffid_raw = df.iat[ridx, i_ffid]
            perp_raw = df.iat[ridx, i_perp]
            shift_raw = df.iat[ridx, i_shift] if i_shift is not None else 0.0

            ffid = None
            perp = None
            shift = 0.0
            try:
                if pd.notna(ffid_raw) and str(ffid_raw).strip() != "":
                    ffid = int(float(str(ffid_raw).replace(",", ".")))
            except Exception:
                ffid = None
            try:
                if pd.notna(perp_raw) and str(perp_raw).strip() != "":
                    perp = float(str(perp_raw).replace(",", "."))
            except Exception:
                perp = None
            try:
                if pd.notna(shift_raw) and str(shift_raw).strip() != "":
                    shift = float(str(shift_raw).replace(",", "."))
            except Exception:
                shift = 0.0

            if perp is None:
                continue

            rows_target.append({"ffid": ffid, "perp": float(perp), "shift": float(shift)})

    if not rows_target:
        return {}, {}

    perp_by_shot: dict = {}
    shift_by_shot: dict = {}

    # Many field sheets use column A as shot index (1/2/3).
    shot_numbers = [int(r["ffid"]) for r in rows_target if r.get("ffid") is not None]
    if shot_numbers:
        uniq = sorted(set(shot_numbers))
        if uniq and min(uniq) >= 1 and max(uniq) <= max(10, len(rows_target) + 1):
            for row in rows_target:
                sid = row.get("ffid")
                if sid is None:
                    continue
                perp_by_shot[int(sid)] = float(row["perp"])
                shift_by_shot[int(sid)] = float(row["shift"])
            return perp_by_shot, shift_by_shot

    if ffid_to_shots:
        ffid_cursor: dict = {k: 0 for k in ffid_to_shots}
        for row in rows_target:
            ffid = row["ffid"]
            if ffid is None:
                continue
            ff_key = int(ffid)
            shot_list = ffid_to_shots.get(ff_key, [])
            if not shot_list:
                continue
            cur = int(ffid_cursor.get(ff_key, 0))
            sid = shot_list[min(cur, len(shot_list) - 1)]
            ffid_cursor[ff_key] = cur + 1
            perp_by_shot[int(sid)] = float(row["perp"])
            shift_by_shot[int(sid)] = float(row["shift"])
    else:
        for i, row in enumerate(rows_target, start=1):
            perp_by_shot[int(i)] = float(row["perp"])
            shift_by_shot[int(i)] = float(row["shift"])

    return perp_by_shot, shift_by_shot

def read_perpendicular_offsets(path: Path) -> dict:
    """
    Read a simple text file with per-shot perpendicular offsets.

    Format:
        # shot_id  perp_m
        1         0.0
        2         3.5
        3         0.0
    """
    out: dict = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            s = str(line).strip()
            if not s or s.startswith("#"):
                continue
            parts = s.replace(",", ".").split()
            if len(parts) < 2:
                continue
            try:
                sid = int(parts[0])
                pov = float(parts[1])
            except Exception:
                continue
            out[int(sid)] = float(pov)
    return out