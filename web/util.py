"""
Projects are filesystem folders (see src.io.project_io), not database rows
with numeric IDs. To reference one in a URL, we base64url-encode its
resolved absolute path - reversible, filesystem-remains-the-source-of-
truth, and doesn't leak the raw path into route logs any more than the
path itself already would.
"""
from __future__ import annotations

import base64
from pathlib import Path


def encode_project_id(root: Path) -> str:
    raw = str(Path(root).resolve()).encode("utf-8")
    return base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")


def decode_project_id(project_id: str) -> Path:
    padding = "=" * (-len(project_id) % 4)
    raw = base64.urlsafe_b64decode(project_id + padding)
    return Path(raw.decode("utf-8"))


__all__ = ["encode_project_id", "decode_project_id"]