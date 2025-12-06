from __future__ import annotations

import sys
from pathlib import Path
from dotenv import load_dotenv

# Ensure repo root (parent of backend) is in sys.path so `agents` is importable.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

# Load .env from repo root, overriding system environment variables
load_dotenv(ROOT / ".env", override=True)

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .config import get_settings
from .database import init_db
from .routers import jobs as jobs_router
from .routers import health as health_router
from .routers import feedback as feedback_router
from .cache import get_cache

settings = get_settings()
app = FastAPI(title="PhaGen Agentic API")
init_db()

# Initialize cache on startup
@app.on_event("startup")
async def startup_event():
    cache = get_cache()
    if cache.enabled:
        print(f"✓ Redis cache enabled: {settings.redis_url}")
    else:
        print("✓ Cache disabled, running in direct mode")

allowed_origins = {
    "http://localhost:3000",
    "http://127.0.0.1:3000",
}

app.add_middleware(
    CORSMiddleware,
    allow_origins=list(allowed_origins),
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(jobs_router.router)
app.include_router(health_router.router)
app.include_router(feedback_router.router)
