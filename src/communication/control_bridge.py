from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def build_picker_command(payload: dict[str, Any], command_id: int) -> dict[str, Any]:
    msg = dict(payload)
    msg["id"] = int(command_id)
    return msg


def write_picker_command(control_file: Path, payload: dict[str, Any], command_id: int) -> dict[str, Any]:
    control_file.parent.mkdir(parents=True, exist_ok=True)
    msg = build_picker_command(payload, command_id)
    with open(control_file, "w", encoding="utf-8") as fh:
        json.dump(msg, fh)
    return msg


def read_latest_picker_command(control_file: Path | None, last_id: int) -> tuple[int, dict[str, Any] | None]:
    if control_file is None or not control_file.exists():
        return last_id, None
    try:
        with open(control_file, encoding="utf-8") as fh:
            payload = json.load(fh)
    except Exception:
        return last_id, None

    if not isinstance(payload, dict):
        return last_id, None

    try:
        cid = int(payload.get("id", 0))
    except Exception:
        return last_id, None

    if cid <= int(last_id):
        return last_id, None

    return cid, payload
