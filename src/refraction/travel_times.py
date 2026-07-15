"""Travel-time line fitting and intercept-time depth formulas."""
from __future__ import annotations

import math
from typing import Any

import numpy as np


def fit_line(x: Any, y: Any) -> tuple:
    """
    Fit y = slope*x + intercept.
    Returns (slope [ms/m], intercept [ms], r^2).
    Velocity (m/s) = 1000 / slope  when x in m, y in ms.
    """
    if len(x) < 2:
        return float("nan"), float("nan"), float("nan")
    x_a = np.asarray(x, dtype=float)
    y_a = np.asarray(y, dtype=float)
    slope, intercept = np.polyfit(x_a, y_a, 1)
    y_hat  = slope * x_a + intercept
    ss_res = float(np.sum((y_a - y_hat) ** 2))
    ss_tot = float(np.sum((y_a - y_a.mean()) ** 2))
    r2     = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else 1.0
    return float(slope), float(intercept), float(r2)

def depth_2layer(ti_ms: float, V1: float, V2: float) -> Any:
    """
    Intercept-time depth, 2-layer model.
    h = (ti_s x V1 x V2) / (2 x sqrt(V2^2 - V1^2))  [m]
    Returns None if V2 <= V1.
    """
    if V2 <= V1 or V1 <= 0.0:
        return None
    return (ti_ms / 1000.0) * V1 * V2 / (2.0 * math.sqrt(V2 ** 2 - V1 ** 2))

def depth_3layer(ti2_ms: float, V1: float, V2: float,
                 V3: float, h1: float) -> Any:
    """
    Intercept-time depth to second refractor, 3-layer model.
    Accounts for first-layer delay time.  Returns None if geometry is invalid.
    """
    if V3 <= V2 or V2 <= V1 or V1 <= 0.0 or h1 is None or h1 <= 0.0:
        return None
    # Match the standard 3-layer intercept-time formulation used in project spreadsheets:
    # h2 = V2 * (ti2 - 2*h1*cos(arcsin(V1/V3))/V1) / (2000 * cos(arcsin(V2/V3)))
    cos_i13 = math.sqrt(max(0.0, 1.0 - (V1 / V3) ** 2))
    cos_i23 = math.sqrt(max(0.0, 1.0 - (V2 / V3) ** 2))
    if cos_i23 <= 0.0:
        return None
    layer1_delay_ms = 2.0 * h1 * cos_i13 / V1 * 1000.0
    ti2_eff = ti2_ms - layer1_delay_ms
    if ti2_eff <= 0.0:
        return None
    return (V2 * ti2_eff) / (2000.0 * cos_i23)
