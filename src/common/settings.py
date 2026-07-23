"""Shared configuration defaults for the LVL picker and refraction engine.

This module is the single source of truth for tunable parameters that were
previously module-level globals in the monolithic `lvl_refraction.py`
script. Values and names are unchanged from the original implementation.
"""
from __future__ import annotations

from src.common.paths import GEOM_TEMPLATES_DIR

# ---------------------------------------------------------------------------
# App identity - used in window titles, About dialogs, and CLI banners.
# Bump APP_VERSION when the picker engine or project format changes in a
# way users should notice.
# ---------------------------------------------------------------------------
APP_NAME: str = "LVL Studio"
APP_VERSION: str = "2.0.0"

# ---------------------------------------------------------------------------
# Geometry files:
# Map geometry type -> receiver geometry definition file.
# These are bundled app-level templates (fixed receiver spacings), not
# part of any one project - see src.common.paths for the distinction.
# ---------------------------------------------------------------------------
GEOM_FILES: dict = {
    100: GEOM_TEMPLATES_DIR / "geometry100.txt",
    200: GEOM_TEMPLATES_DIR / "geometry200.txt",
}

# ---------------------------------------------------------------------------
# Profile registry
# ---------------------------------------------------------------------------
# geom    : geometry type (100 or 200)
# line_no : legacy numeric identifier (kept for backward compatibility)
# perp_m  : perpendicular distance of shot from the geophone line (m)
# shots   : {shot_id: inline_position_m}  OR  "auto"
#           "auto" = derive from geometry:
#             shot 1 = recv[0] - first_spacing     (before spread)
#             shot 2 = (recv[n//2-1] + recv[n//2]) / 2   (midpoint)
#             shot 3 = recv[-1] + last_spacing     (after spread)
#           Manual values are used when SEG2 SOURCE_LOCATION is absent.
#           SEG2 header always takes priority if present.
# ---------------------------------------------------------------------------
PROFILES: dict = {
    "120":   {"geom": 100, "line_no": 1, "perp_m": 0.0, "shots": "auto"},
    "214_A": {"geom": 100, "line_no": 2, "perp_m": 0.0, "shots": "auto"},
    "150":   {"geom": 200, "line_no": 3, "perp_m": 0.0, "shots": "auto"},
    "415_B": {"geom": 200, "line_no": 4, "perp_m": 0.0, "shots": "auto"},
}

# Shot position source policy:
# False -> use geometry/config shot positions by default (recommended)
# True  -> allow SEG2 SOURCE_LOCATION to override when available
USE_SEG2_SHOT_POSITION: bool = False

# ---------------------------------------------------------------------------
# Ormsby bandpass  (zero-phase, frequency domain)
# ProMAX convention:  f1 - f2 - f3 - f4
#   f1 : low-cut ramp start      (Hz)
#   f2 : low-cut ramp end / passband start   (Hz)
#   f3 : passband end / high-cut ramp start  (Hz)
#   f4 : high-cut ramp end       (Hz)
# ---------------------------------------------------------------------------
BP_F1: float = 2.0    # Hz
BP_F2: float = 4.0    # Hz
BP_F3: float = 140.0  # Hz
BP_F4: float = 180.0  # Hz
BP_FFT_PAD: float = 0.25   # 25 % zero-padding (ProMAX default)
BP_REAPPLY: bool  = False  # True = apply twice (squares amplitude response)
BUTTER_ORDER: int = 4      # Butterworth order (zero-phase via sosfiltfilt)
FILTER_DEBOUNCE_MS: int = 40  # 0 = immediate live update, >0 = debounce while dragging

# ---------------------------------------------------------------------------
# Bulk time static  (ProMAX: "Bulk shift static -> Add")
# The SEG2 DELAY field (e.g. -106 ms) is read from each file and used to
# build the correct absolute time axis.  Picks are stored as absolute times
# from the trigger (t = 0).  BULK_SHIFT_MS is an ADDITIONAL correction on
# top of the DELAY  (leave at 0.0 unless a systematic offset remains).
# ---------------------------------------------------------------------------
BULK_SHIFT_MS: float = 0.0   # ms  (negative = shift to earlier times)

# ---------------------------------------------------------------------------
# Display / theme
# ---------------------------------------------------------------------------
THEME: str         = "light"   # "dark"  or  "light"
T_MAX_MS: float    = 150.0    # end of display window (ms from trigger)
T_DISPLAY_PRE_MS: float = -10.0  # ms before trigger to show (headroom above first arrivals)
CLIP_FACTOR: float = 2.0      # wiggle clip level in std-dev multiples
FILTER_ON: bool    = True     # start with filter active
DEFAULT_FILTER_MODE: str = "butter"  # "none" | "butter" | "ormsby" | "cutpass"
DISPLAY_MODE: str  = "wiggle"   # "wiggle" | "vd" | "both"
SHOW_TIMELINES: bool = True
WIGGLE_STRETCH: float = 1.0

# ---------------------------------------------------------------------------
# Gain control for display/picking
# ---------------------------------------------------------------------------
GAIN_MODE: str        = "norm"   # "none" | "norm" | "agc"
AGC_WINDOW_MS: float  = 200.0
AGC_STAT: str         = "rms"    # "mean" | "rms"

# ---------------------------------------------------------------------------
# STA/LTA auto-picker
# ---------------------------------------------------------------------------
STA_MS: float      = 3.0
LTA_MS: float      = 20.0
STALTA_TRIG: float = 3.0
AUTO_PICK_MODE: str = "hilbert_env"  # "stalta" | "maxdiff_zero" | "hilbert_env"
ZERO_X_SEARCH_DIRECTION: str = "backward"  # "backward" | "forward"
HILBERT_ONSET_PCT: float = 0.08  # Fraction of peak envelope amplitude for onset pick (0.0-1.0)
HILBERT_TARGET: str = "onset"  # "onset" | "peak"
PICK_PROCESS_ORDER: str = "F>G>X"  # "F>G>X" | "G>F>X" | "RAW>X"
MANUAL_SNAP_WIN_MS: float = 6.0  # Manual pick snap window for zero-crossing search (ms)
# Auto-pick first-arrival plausibility gate (absolute time, ms)
AUTO_USE_VELOCITY_GATE: bool = True
AUTO_VMIN_M_S: float = 120.0
AUTO_VMAX_M_S: float = 3500.0
AUTO_GATE_PAD_MS: float = 20.0

# ---------------------------------------------------------------------------
# Trigger/radio correction (ms)
# Static term is added to SEG2 delay when data are loaded so display and picks
# share the same time basis. Optional linear term is applied per trace offset.
# ---------------------------------------------------------------------------
TRIGGER_STATIC_MS: float = 6.0
TRIGGER_MS_PER_M: float = 0.0

# ---------------------------------------------------------------------------
# Refraction analysis
# ---------------------------------------------------------------------------
N_LAYERS: int          = 2      # 2 or 3
AUTO_BREAKPOINTS: bool = True   # True = AIC search;  False = MANUAL_BREAKS
# MANUAL_BREAKS: {shot_id: [receiver_index_of_crossover, ...]}
# Example:  {1: [12], 2: [14], 3: [12]}
MANUAL_BREAKS: dict = {}


def theme_colors(theme: str) -> dict:
    """Return the matplotlib colour palette for a given theme name.

    Parameters
    ----------
    theme : str
        Either ``"dark"`` or ``"light"``.

    Returns
    -------
    dict
        Colour keys used throughout picker/analysis plots.

    Notes
    -----
    Takes ``theme`` as an explicit parameter (rather than reading a mutable
    module global) so callers in any module always render with the theme
    that is actually active in the GUI, without cross-module global-state
    coupling.
    """
    if str(theme).lower() == "dark":
        return dict(
            fig_bg="#0f1117",   ax_bg="#1a1d2e",
            trace="#4a8fc1",    fill_alpha=0.20,
            pick_clr="#e63946", text="white",
            label="#aaa",       grid="#888",
            tick="#888",        spine="#333344",
            leg_face="#1e2130", leg_edge="#555",
            fit_alpha=0.90,
        )
    return dict(
        fig_bg="white",     ax_bg="#f5f5f5",
        trace="#111111",    fill_alpha=0.18,
        pick_clr="#cc0000", text="#111111",
        label="#333",       grid="#bbb",
        tick="#333",        spine="#aaaaaa",
        leg_face="#f0f0f0", leg_edge="#aaa",
        fit_alpha=0.85,
    )


def _ensure_interactive_backend() -> bool:
    """Ensure matplotlib uses an interactive backend if possible."""
    import matplotlib.pyplot as plt

    def _is_non_interactive(name: str) -> bool:
        n = str(name).strip().lower()
        non_interactive = {
            "agg", "pdf", "ps", "svg", "template", "cairo",
            "module://matplotlib_inline.backend_inline",
            "module://ipykernel.pylab.backend_inline",
        }
        return n in non_interactive

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
        print("  [WARN] Matplotlib backend is non-interactive (Agg); UI plots cannot open.")
        return False

    print(f"  [INFO] Switched matplotlib backend to interactive mode: {plt.get_backend()}")
    return True