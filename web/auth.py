"""
Authentication and per-project authorization.

Sessions are signed cookies (Starlette's SessionMiddleware) holding just
the user id - no server-side session table needed. Passwords are hashed
with bcrypt via passlib. Project-level permissions are checked against
`project_access` (see db.py) on every request that touches a project, not
cached, so revoking access takes effect immediately.
"""
from __future__ import annotations

from pathlib import Path

from fastapi import Depends, HTTPException, Request, status
from passlib.context import CryptContext

from web.db import get_db, role_at_least

_pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")


def hash_password(password: str) -> str:
    return _pwd_context.hash(password)


def verify_password(password: str, password_hash: str) -> bool:
    try:
        return _pwd_context.verify(password, password_hash)
    except Exception:
        return False


def create_user(username: str, password: str, display_name: str = "", role: str = "user") -> int:
    from web.db import now_iso
    with get_db() as conn:
        cur = conn.execute(
            "INSERT INTO users (username, password_hash, display_name, role, created_at) "
            "VALUES (?, ?, ?, ?, ?)",
            (username.strip().lower(), hash_password(password), display_name, role, now_iso()),
        )
        return cur.lastrowid


def authenticate(username: str, password: str) -> dict | None:
    with get_db() as conn:
        row = conn.execute(
            "SELECT * FROM users WHERE username = ?", (username.strip().lower(),)
        ).fetchone()
    if row is None or not verify_password(password, row["password_hash"]):
        return None
    return dict(row)


def get_user_by_id(user_id: int) -> dict | None:
    with get_db() as conn:
        row = conn.execute("SELECT * FROM users WHERE id = ?", (user_id,)).fetchone()
    return dict(row) if row else None


def any_users_exist() -> bool:
    with get_db() as conn:
        row = conn.execute("SELECT COUNT(*) AS n FROM users").fetchone()
    return row["n"] > 0


# -----------------------------------------------------------------------------
# FastAPI dependencies
# -----------------------------------------------------------------------------

def get_current_user(request: Request) -> dict:
    """Require a logged-in user; raises 401 if the session has no user."""
    user_id = request.session.get("user_id")
    if not user_id:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Not logged in")
    user = get_user_by_id(user_id)
    if user is None:
        request.session.clear()
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Session invalid")
    return user


def get_current_user_optional(request: Request) -> dict | None:
    user_id = request.session.get("user_id")
    return get_user_by_id(user_id) if user_id else None


def require_admin(user: dict = Depends(get_current_user)) -> dict:
    if user["role"] != "admin":
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Admin only")
    return user


def user_role_for_project(user_id: int, project_root: Path) -> str | None:
    """The caller's role on this project ('owner'/'editor'/'viewer'), or None."""
    with get_db() as conn:
        row = conn.execute(
            "SELECT role FROM project_access WHERE user_id = ? AND project_root = ?",
            (user_id, str(Path(project_root).resolve())),
        ).fetchone()
    return row["role"] if row else None


def grant_project_access(user_id: int, project_root: Path, project_name: str, role: str = "owner") -> None:
    from web.db import now_iso
    with get_db() as conn:
        conn.execute(
            "INSERT INTO project_access (user_id, project_root, project_name, role, granted_at) "
            "VALUES (?, ?, ?, ?, ?) "
            "ON CONFLICT(user_id, project_root) DO UPDATE SET role = excluded.role",
            (user_id, str(Path(project_root).resolve()), project_name, role, now_iso()),
        )


def revoke_project_access(user_id: int, project_root: Path) -> None:
    with get_db() as conn:
        conn.execute(
            "DELETE FROM project_access WHERE user_id = ? AND project_root = ?",
            (user_id, str(Path(project_root).resolve())),
        )


def projects_for_user(user_id: int) -> list[dict]:
    with get_db() as conn:
        rows = conn.execute(
            "SELECT project_root, project_name, role FROM project_access "
            "WHERE user_id = ? ORDER BY project_name",
            (user_id,),
        ).fetchall()
    return [dict(r) for r in rows]


class RequireProjectRole:
    """FastAPI dependency factory: require at least `minimum` role on the
    project identified by the `project_root` path (query/body param, not a
    route path param - callers resolve the path first, then pass it here).

    Usage: `role: str = Depends(RequireProjectRole("editor"))` alongside a
    route that also depends on `get_current_user` and resolves the target
    project's root path from its own path parameter.
    """

    def __init__(self, minimum: str = "viewer"):
        self.minimum = minimum

    def __call__(self, project_root: Path, user: dict = Depends(get_current_user)) -> str:
        role = user_role_for_project(user["id"], project_root)
        if role is None or not role_at_least(role, self.minimum):
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail=f"No '{self.minimum}'-level access to this project.",
            )
        return role


__all__ = [
    "hash_password", "verify_password", "create_user", "authenticate",
    "get_user_by_id", "any_users_exist",
    "get_current_user", "get_current_user_optional", "require_admin",
    "user_role_for_project", "grant_project_access", "revoke_project_access",
    "projects_for_user", "RequireProjectRole",
]