"""
Temporary compatibility wrapper for the clean first-break picker.

This file intentionally stays minimal so temp/ no longer contains legacy,
unused refraction-analysis code.

Use the clean picker module directly:
  python ../scripts/fb_picker.py <profile> [geometry]
"""

from __future__ import annotations

import sys
from pathlib import Path


# Ensure ../scripts is importable when running this temp wrapper directly.
SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from fb_picker import main  # type: ignore  # noqa: E402


if __name__ == "__main__":
    main()
