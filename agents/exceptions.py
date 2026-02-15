"""Structured exception hierarchy for PhaGen agent pipeline.

All agent-layer exceptions inherit from PhaGenError, allowing callers
to catch broadly (PhaGenError) or narrowly (e.g. LLMTimeoutError).
"""

from __future__ import annotations


class PhaGenError(Exception):
    """Base exception for all PhaGen agent errors."""

    def __init__(
        self,
        message: str,
        *,
        code: str = "PHAGEN_ERROR",
        details: dict | None = None,
    ) -> None:
        self.code = code
        self.details = details or {}
        super().__init__(message)


# -- LLM errors ---------------------------------------------------------------

class LLMError(PhaGenError):
    """LLM generation failed."""

    def __init__(self, message: str, **kw) -> None:
        super().__init__(message, code=kw.pop("code", "LLM_ERROR"), **kw)


class LLMTimeoutError(LLMError):
    """LLM call exceeded its deadline."""

    def __init__(self, message: str, **kw) -> None:
        super().__init__(message, code=kw.pop("code", "LLM_TIMEOUT"), **kw)


class LLMParseError(LLMError):
    """LLM output could not be parsed into the expected structure."""

    def __init__(self, message: str, *, raw: str = "", **kw) -> None:
        super().__init__(message, code=kw.pop("code", "LLM_PARSE_ERROR"), **kw)
        self.raw = raw


# -- Retrieval errors ----------------------------------------------------------

class RetrievalError(PhaGenError):
    """Vector store query failed."""

    def __init__(self, message: str, **kw) -> None:
        super().__init__(message, code=kw.pop("code", "RETRIEVAL_ERROR"), **kw)


# -- Tool errors ---------------------------------------------------------------

class ToolExecutionError(PhaGenError):
    """An agent tool failed during execution."""

    def __init__(self, tool_name: str, message: str, **kw) -> None:
        self.tool_name = tool_name
        super().__init__(
            f"Tool '{tool_name}' failed: {message}",
            code=kw.pop("code", "TOOL_EXECUTION_ERROR"),
            **kw,
        )


class BudgetExhaustedError(PhaGenError):
    """Agent exceeded its step, token, or API-call budget."""

    def __init__(self, message: str, **kw) -> None:
        super().__init__(message, code=kw.pop("code", "BUDGET_EXHAUSTED"), **kw)


# -- Worker / supervisor errors ------------------------------------------------

class WorkerError(PhaGenError):
    """A domain worker failed."""

    def __init__(self, worker_name: str, message: str, **kw) -> None:
        self.worker_name = worker_name
        super().__init__(
            f"Worker '{worker_name}': {message}",
            code=kw.pop("code", "WORKER_ERROR"),
            **kw,
        )


class SupervisorError(PhaGenError):
    """Supervisor-level orchestration failure."""

    def __init__(self, message: str, **kw) -> None:
        super().__init__(message, code=kw.pop("code", "SUPERVISOR_ERROR"), **kw)


# -- Validation errors ---------------------------------------------------------

class ValidationError(PhaGenError):
    """Input validation failed."""

    def __init__(self, message: str, **kw) -> None:
        super().__init__(message, code=kw.pop("code", "VALIDATION_ERROR"), **kw)


# -- Circuit breaker -----------------------------------------------------------

class CircuitOpenError(PhaGenError):
    """Circuit breaker is open — upstream service is unhealthy."""

    def __init__(self, service: str, message: str = "", **kw) -> None:
        self.service = service
        super().__init__(
            message or f"Circuit open for '{service}'; retries paused",
            code=kw.pop("code", "CIRCUIT_OPEN"),
            **kw,
        )
