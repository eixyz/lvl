"""Pure geometric/math helpers shared across the refraction engine."""
from __future__ import annotations

from typing import Any

import numpy as np


def true_offset(inline_m: Any, perp_m: float) -> Any:
    """True source-receiver distance corrected for perpendicular shot offset."""
    return np.sqrt(np.asarray(inline_m, dtype=float) ** 2 + perp_m ** 2)
