"""
File-based logging for both front-ends.

Both `lvl_studio.py` and `lvl_refraction.py` are full of `print()` calls -
that's fine for interactive use, but it means nothing survives after the
terminal closes. Rather than rewrite every print() into a logger call
(hundreds of call sites, high risk of missing one), this wraps
`sys.stdout`/`sys.stderr` in a small "tee" that writes everything to both
the original stream AND a timestamped log file. Call `setup_file_logging`
once at startup; everything printed anywhere in the process from then on
is captured automatically.
"""
from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

from src.common.paths import APP_ROOT

LOGS_DIR = APP_ROOT / "logs"


class _Tee:
    """Write to two streams at once; flush eagerly so `tail -f` works."""

    def __init__(self, primary, secondary):
        self._primary = primary
        self._secondary = secondary

    def write(self, data: str) -> int:
        self._primary.write(data)
        try:
            self._secondary.write(data)
            self._secondary.flush()
        except Exception:
            pass  # never let a logging hiccup break the app
        return len(data)

    def flush(self) -> None:
        self._primary.flush()
        try:
            self._secondary.flush()
        except Exception:
            pass

    def isatty(self) -> bool:
        return getattr(self._primary, "isatty", lambda: False)()

    def __getattr__(self, name):
        # Forward anything else (fileno, encoding, ...) to the real stream,
        # so libraries that introspect sys.stdout don't break.
        return getattr(self._primary, name)


def setup_file_logging(app_name: str) -> Path:
    """Tee stdout/stderr to `logs/<app_name>_<timestamp>.log`. Idempotent-
    ish: calling it twice in the same process just opens a second file
    (harmless, just wastes a file handle) - call it once at startup.
    """
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    log_path = LOGS_DIR / f"{app_name}_{stamp}.log"

    log_file = open(log_path, "a", encoding="utf-8", buffering=1)
    log_file.write(f"=== {app_name} session started {datetime.now().isoformat(timespec='seconds')} ===\n")

    sys.stdout = _Tee(sys.stdout, log_file)
    sys.stderr = _Tee(sys.stderr, log_file)

    return log_path


__all__ = ["setup_file_logging", "LOGS_DIR"]