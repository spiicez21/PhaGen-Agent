from __future__ import annotations

import sys
from pathlib import Path

# Ensure repo root (parent of backend) is in sys.path so `agents` is importable.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from fastapi import FastAPI

from .config import get_settings
from .routers import jobs as jobs_router

settings = get_settings()
app = FastAPI(title="PhaGen Agentic API")
app.include_router(jobs_router.router)


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}
