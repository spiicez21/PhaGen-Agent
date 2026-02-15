"""Centralized logging configuration for the PhaGen agent pipeline.

Usage:
    from agents.logging_config import get_logger
    logger = get_logger(__name__)
    logger.info("worker_started", worker="clinical", molecule="aspirin")
"""

from __future__ import annotations

import logging
import os
import sys


_CONFIGURED = False


def configure_logging(level: str | None = None) -> None:
    """Set up structured console logging for the whole process.

    Safe to call multiple times — subsequent calls are no-ops.
    """
    global _CONFIGURED
    if _CONFIGURED:
        return
    _CONFIGURED = True

    resolved_level = (level or os.getenv("LOG_LEVEL", "INFO")).upper()
    numeric_level = getattr(logging, resolved_level, logging.INFO)

    # Use a clean format: timestamp level name message
    formatter = logging.Formatter(
        fmt="%(asctime)s %(levelname)-8s [%(name)s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(formatter)

    root = logging.getLogger()
    root.setLevel(numeric_level)
    # Remove any existing handlers to avoid duplicates
    root.handlers.clear()
    root.addHandler(handler)

    # Quiet noisy third-party loggers
    for noisy in ("httpx", "httpcore", "chromadb", "sentence_transformers", "urllib3"):
        logging.getLogger(noisy).setLevel(logging.WARNING)


def get_logger(name: str) -> logging.Logger:
    """Return a named logger, ensuring logging is configured."""
    configure_logging()
    return logging.getLogger(name)
