"""
SQLite-backed user accounts and project permissions for the web app.

Deliberately NOT an ORM - this is a small, internal tool and plain sqlite3
keeps the whole persistence layer readable in one file, which matters
since this needs to be maintainable by whoever inherits it.

Projects themselves are still 100% filesystem-based (see
src.io.project_io) - nothing about a project's data lives in this
database. This DB only answers "which users can open which project
folders, at what permission level" - a thin access-control layer on top
of the existing project system, not a second source of truth.
"""
from __future__ import annotations

import sqlite3
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

DB_PATH = Path(__file__).resolve().parent / "lvl_web.sqlite3"

SCHEMA = """
CREATE TABLE IF NOT EXISTS users (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    username      TEXT UNIQUE NOT NULL,
    password_hash TEXT NOT NULL,
    display_name  TEXT NOT NULL DEFAULT '',
    role          TEXT NOT NULL DEFAULT 'user',   -- 'admin' | 'user'
    created_at    TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS project_access (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    user_id       INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    project_root  TEXT NOT NULL,                  -- resolved Project.root path
    project_name  TEXT NOT NULL,
    role          TEXT NOT NULL DEFAULT 'owner',   -- 'owner' | 'editor' | 'viewer'
    granted_at    TEXT NOT NULL,
    UNIQUE(user_id, project_root)
);
"""

# Permission hierarchy: viewer can read; editor can read+write (pick/compute);
# owner can additionally share/revoke access to the project.
_ROLE_RANK = {"viewer": 0, "editor": 1, "owner": 2}


def role_at_least(role: str, minimum: str) -> bool:
    return _ROLE_RANK.get(role, -1) >= _ROLE_RANK.get(minimum, 99)


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@contextmanager
def get_db():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()


def init_db() -> None:
    with get_db() as conn:
        conn.executescript(SCHEMA)


__all__ = ["DB_PATH", "get_db", "init_db", "now_iso", "role_at_least"]