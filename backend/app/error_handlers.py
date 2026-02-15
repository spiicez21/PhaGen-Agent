"""Structured error handlers for the FastAPI application."""

from __future__ import annotations

from fastapi import Request
from fastapi.responses import JSONResponse

# Import agent exceptions — these may not be installed in all environments
try:
    from agents.exceptions import PhaGenError, ValidationError as AgentValidationError
except ImportError:
    PhaGenError = None  # type: ignore[assignment, misc]
    AgentValidationError = None  # type: ignore[assignment, misc]


async def phagen_error_handler(request: Request, exc: Exception) -> JSONResponse:
    """Handle PhaGenError and subclasses with structured JSON responses."""
    code = getattr(exc, "code", "INTERNAL_ERROR")
    details = getattr(exc, "details", {})
    return JSONResponse(
        status_code=500,
        content={
            "error": {
                "code": code,
                "message": str(exc),
                "details": details,
            }
        },
    )


async def validation_error_handler(request: Request, exc: Exception) -> JSONResponse:
    """Handle validation errors with 422 status."""
    code = getattr(exc, "code", "VALIDATION_ERROR")
    return JSONResponse(
        status_code=422,
        content={
            "error": {
                "code": code,
                "message": str(exc),
            }
        },
    )


def register_error_handlers(app) -> None:
    """Register structured error handlers on the FastAPI app."""
    if PhaGenError is not None:
        app.add_exception_handler(PhaGenError, phagen_error_handler)
    if AgentValidationError is not None:
        app.add_exception_handler(AgentValidationError, validation_error_handler)
