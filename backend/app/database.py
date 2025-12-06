from __future__ import annotations

import logging
from contextlib import contextmanager
from typing import Iterator

from sqlalchemy import create_engine
from sqlalchemy.orm import DeclarativeBase, Session, sessionmaker

from .config import get_settings


class Base(DeclarativeBase):
    """Declarative base for all ORM models."""


_settings = get_settings()
_engine = create_engine(
    _settings.database_url,
    echo=_settings.database_echo,
    future=True,
)
SessionLocal = sessionmaker(
    bind=_engine,
    autoflush=False,
    autocommit=False,
    expire_on_commit=False,
    class_=Session,
)


@contextmanager
def get_session() -> Iterator[Session]:
    session = SessionLocal()
    try:
        yield session
    finally:
        session.close()


def get_db() -> Iterator[Session]:
    """FastAPI dependency for database sessions."""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def init_db() -> None:
    """Create tables if they do not already exist."""
    from sqlalchemy.exc import SQLAlchemyError

    from . import models  # noqa: F401  ensures models register with metadata

    try:
        Base.metadata.create_all(bind=_engine)
    except SQLAlchemyError as exc:  # pragma: no cover - executed at startup
        logging.getLogger(__name__).error("Database initialization failed: %s", exc)
        raise


__all__ = ["Base", "SessionLocal", "get_session", "init_db"]
