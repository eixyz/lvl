from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from src.application.studio import Studio


app = FastAPI(
    title="LVL Studio API"
)


studio = Studio()


@app.get("/api/status")
def status():

    return {
        "application": "LVL Studio",
        "loaded": studio.project is not None
    }


@app.post("/api/project/open")
def open_project(path: str):

    studio.open_project(path)

    return {
        "status": "opened",
        "project": path
    }


app.mount(
    "/",
    StaticFiles(
        directory="apps/web/static",
        html=True
    ),
    name="web"
)