from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException, Request, status
from pydantic import BaseModel

from web.auth import authenticate, create_user, get_current_user, require_admin

router = APIRouter(prefix="/auth", tags=["auth"])


class LoginRequest(BaseModel):
    username: str
    password: str


class CreateUserRequest(BaseModel):
    username: str
    password: str
    display_name: str = ""
    role: str = "user"  # 'admin' | 'user'


@router.post("/login")
def login(body: LoginRequest, request: Request):
    user = authenticate(body.username, body.password)
    if user is None:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid username or password")
    request.session["user_id"] = user["id"]
    return {"id": user["id"], "username": user["username"], "display_name": user["display_name"], "role": user["role"]}


@router.post("/logout")
def logout(request: Request):
    request.session.clear()
    return {"ok": True}


@router.get("/me")
def me(user: dict = Depends(get_current_user)):
    return {"id": user["id"], "username": user["username"], "display_name": user["display_name"], "role": user["role"]}


@router.post("/users", status_code=status.HTTP_201_CREATED)
def admin_create_user(body: CreateUserRequest, _admin: dict = Depends(require_admin)):
    """Create a new account. Admin-only - this is an internal tool, not
    open self-registration, so an admin provisions accounts for the team.
    """
    if body.role not in ("admin", "user"):
        raise HTTPException(status_code=400, detail="role must be 'admin' or 'user'")
    try:
        user_id = create_user(body.username, body.password, body.display_name, body.role)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Could not create user: {exc}") from exc
    return {"id": user_id, "username": body.username.strip().lower()}