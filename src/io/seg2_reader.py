"""SEG2 shot-record readers (traces, headers, acquisition time, coordinates)."""
from __future__ import annotations

import datetime
import warnings
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

import numpy as np
from obspy import read as _read_obspy

from src.common.settings import TRIGGER_STATIC_MS


def read_seg2(path: Path) -> tuple:
    """
    Read one SEG2 shot record.

    Returns
    -------
    data          : ndarray (n_traces x n_samples, float32)
    dt_s          : sample interval (s)
    n_traces      : int
    n_samples     : int
    shot_pos      : float or None  (SOURCE_LOCATION header, first component)
    ffid          : int  (SHOT_SEQUENCE_NUMBER or filename digits)
    delay_ms      : float  (DELAY header in ms; typically negative pre-trigger)
    recv_locs_m   : ndarray or None  (RECEIVER_LOCATION per trace, m)
    """
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        st = _read_obspy(str(path), format="SEG2")

    n_traces = len(st)
    n_samp   = max(tr.stats.npts for tr in st)
    dt_s     = float(st[0].stats.delta)

    data = np.zeros((n_traces, n_samp), dtype=np.float32)
    for i, tr in enumerate(st):
        npts = tr.stats.npts
        data[i, :npts] = tr.data.astype(np.float32)

    shot_pos: Any  = None
    ffid: int      = 0
    delay_ms: float = 0.0
    recv_locs_m: Any = None

    try:
        hdr0 = dict(st[0].stats.seg2)

        sloc = str(hdr0.get("SOURCE_LOCATION", "")).strip()
        if sloc:
            # Guard against German comma decimal
            shot_pos = float(sloc.replace(",", ".").split()[0])

        ssn = str(hdr0.get("SHOT_SEQUENCE_NUMBER", "")).strip()
        if ssn.isdigit():
            ffid = int(ssn)

        # DELAY is stored in seconds; convert to ms and include trigger static.
        delay_str = str(hdr0.get("DELAY", "0")).strip()
        delay_ms  = float(delay_str.replace(",", ".")) * 1000.0 + float(TRIGGER_STATIC_MS)

        # RECEIVER_LOCATION per trace (m along profile)
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
        ffid   = int(digits) if digits else 0

    return data, dt_s, n_traces, n_samp, shot_pos, ffid, delay_ms, recv_locs_m

def read_seg2_acquisition_time_de(path: Path,
                                  tz_name: str = "Europe/Berlin") -> str | None:
    """Return acquisition datetime formatted in German timezone from SEG2 metadata."""
    import warnings

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            st = _read_obspy(str(path), format="SEG2")
    except Exception:
        return None
    if not st:
        return None

    dt_utc = None
    try:
        dt0 = st[0].stats.starttime.datetime
        if dt0.tzinfo is None:
            dt_utc = dt0.replace(tzinfo=datetime.timezone.utc)
        else:
            dt_utc = dt0.astimezone(datetime.timezone.utc)
    except Exception:
        dt_utc = None

    if dt_utc is None:
        try:
            hdr0 = dict(st[0].stats.seg2)
        except Exception:
            hdr0 = {}
        d_raw = str(hdr0.get("ACQUISITION_DATE", "")).strip()
        t_raw = str(hdr0.get("ACQUISITION_TIME", "")).strip()
        if d_raw and t_raw:
            for fmt in ("%d/%m/%Y %H:%M:%S", "%Y-%m-%d %H:%M:%S", "%d.%m.%Y %H:%M:%S"):
                try:
                    dt_naive = datetime.datetime.strptime(f"{d_raw} {t_raw}", fmt)
                    dt_utc = dt_naive.replace(tzinfo=datetime.timezone.utc)
                    break
                except Exception:
                    continue

    if dt_utc is None:
        return None

    try:
        tz = ZoneInfo(tz_name)
    except Exception:
        tz = datetime.timezone(datetime.timedelta(hours=1), name="CET")
    dt_local = dt_utc.astimezone(tz)
    return dt_local.strftime("%Y-%m-%d %H:%M:%S %Z")

def read_seg2_mid_xyz(path: Path) -> tuple:
    """Try extracting (X,Y,Z) from SEG2 headers; return (None,None,None) if unavailable."""
    import warnings

    def _to_float(v: Any) -> float | None:
        try:
            s = str(v).strip().replace(",", ".")
            if not s:
                return None
            return float(s)
        except Exception:
            return None

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            st = _read_obspy(str(path), format="SEG2")
    except Exception:
        return None, None, None
    if not st:
        return None, None, None

    try:
        h = dict(st[0].stats.seg2)
    except Exception:
        h = {}

    key_x = next((k for k in h.keys() if str(k).upper() in ("X", "X_UTM", "EASTING", "UTM_X")), None)
    key_y = next((k for k in h.keys() if str(k).upper() in ("Y", "Y_UTM", "NORTHING", "UTM_Y")), None)
    key_z = next((k for k in h.keys() if str(k).upper() in ("Z", "ELEV", "ELEVATION", "HEIGHT", "UTM_Z")), None)
    if key_x and key_y and key_z:
        xv = _to_float(h.get(key_x))
        yv = _to_float(h.get(key_y))
        zv = _to_float(h.get(key_z))
        if xv is not None and yv is not None and zv is not None:
            return xv, yv, zv

    src_loc = str(h.get("SOURCE_LOCATION", "")).strip().replace(",", ".")
    if src_loc:
        parts = src_loc.split()
        if len(parts) >= 3:
            xv = _to_float(parts[0])
            yv = _to_float(parts[1])
            zv = _to_float(parts[2])
            if xv is not None and yv is not None and zv is not None:
                return xv, yv, zv

    return None, None, None
