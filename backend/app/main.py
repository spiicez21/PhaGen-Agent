from __future__ import annotations

import os
import sys
from pathlib import Path
from dotenv import load_dotenv

# Ensure repo root (parent of backend) is in sys.path so `agents` is importable.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

# Load .env from repo root, overriding system environment variables
load_dotenv(ROOT / ".env", override=True)

# Configure structured logging for the entire application
from agents.logging_config import configure_logging, get_logger

configure_logging()
logger = get_logger(__name__)

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .config import get_settings
from .database import init_db
from .error_handlers import register_error_handlers
from .routers import jobs as jobs_router
from .routers import health as health_router
from .routers import feedback as feedback_router
from .routers import auth as auth_router
from .cache import get_cache

settings = get_settings()
app = FastAPI(title="PhaGen Agentic API")
init_db()
register_error_handlers(app)

# Initialize cache on startup
@app.on_event("startup")
async def startup_event():
    cache = get_cache()
    if cache.enabled:
        logger.info("redis_cache_enabled", url=settings.redis_url)
    else:
        logger.info("cache_disabled_direct_mode")

_DEFAULT_ORIGINS = "http://localhost:3000,http://127.0.0.1:3000"
allowed_origins = [
    o.strip()
    for o in os.getenv("CORS_ORIGINS", _DEFAULT_ORIGINS).split(",")
    if o.strip()
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(auth_router.router)
app.include_router(jobs_router.router)
app.include_router(health_router.router)
app.include_router(feedback_router.router)
