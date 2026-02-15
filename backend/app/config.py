from __future__ import annotations

import os
import warnings
from functools import lru_cache
from pathlib import Path

from pydantic import AliasChoices, Field, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict


PROJECT_ROOT_ENV = Path(__file__).resolve().parents[2] / ".env"

_PLACEHOLDER_SECRETS = {
    "your-secret-key-change-in-production-min-32-chars-long",
    "admin",
    "phagen123",
    "changeme",
}


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_file=str(PROJECT_ROOT_ENV),
        env_file_encoding="utf-8",
        extra="allow",
        env_file_override=True,
    )

    # -- Application -------------------------------------------------------
    api_prefix: str = Field(default="/api")
    rag_top_k: int = Field(default=5)
    job_ttl_minutes: int = Field(default=60)

    # -- Database ----------------------------------------------------------
    database_url: str = Field(
        default="postgresql+psycopg://phagen:phagen@localhost:5432/phagen",
        validation_alias=AliasChoices("DATABASE_URL", "SUPABASE_URL"),
        description="SQLAlchemy connection string for the operational database.",
    )
    database_echo: bool = Field(default=False)

    # -- Object storage (S3 / MinIO) ---------------------------------------
    storage_provider: str = Field(default="s3")
    s3_endpoint_url: str = Field(default="http://localhost:9000")
    s3_region: str = Field(default="us-east-1")
    s3_access_key: str = Field(
        default="admin",
        description="S3 access key. Set via S3_ACCESS_KEY env var.",
    )
    s3_secret_key: str = Field(
        default="phagen123",
        description="S3 secret key. Set via S3_SECRET_KEY env var.",
    )
    s3_use_ssl: bool = Field(default=False)
    s3_raw_documents_bucket: str = Field(default="phagen-raw-documents")
    s3_reports_bucket: str = Field(default="phagen-report-artifacts")

    # -- Redis / Celery ----------------------------------------------------
    redis_url: str = Field(default="redis://localhost:6379/0")
    cache_enabled: bool = Field(default=True)
    celery_broker_url: str = Field(default="redis://localhost:6379/1")
    celery_result_backend: str = Field(default="redis://localhost:6379/1")

    # -- JWT Authentication ------------------------------------------------
    secret_key: str = Field(
        default="your-secret-key-change-in-production-min-32-chars-long",
        description="JWT signing secret (min 32 chars). Set via SECRET_KEY env var.",
    )
    algorithm: str = Field(default="HS256")
    access_token_expire_minutes: int = Field(default=30)

    @model_validator(mode="after")
    def _warn_placeholder_secrets(self) -> Settings:
        """Emit loud warnings when placeholder secrets are still in use."""
        env = os.getenv("PHAGEN_ENV", "development").lower()
        if env in ("production", "staging"):
            problems: list[str] = []
            if self.secret_key in _PLACEHOLDER_SECRETS:
                problems.append("SECRET_KEY is still a placeholder")
            if len(self.secret_key) < 32:
                problems.append(f"SECRET_KEY too short ({len(self.secret_key)} chars, need 32+)")
            if self.s3_access_key in _PLACEHOLDER_SECRETS:
                problems.append("S3_ACCESS_KEY is still a placeholder")
            if self.s3_secret_key in _PLACEHOLDER_SECRETS:
                problems.append("S3_SECRET_KEY is still a placeholder")
            if problems:
                raise ValueError(
                    "Insecure configuration in production: " + "; ".join(problems)
                )
        else:
            if self.secret_key in _PLACEHOLDER_SECRETS:
                warnings.warn(
                    "SECRET_KEY is a placeholder. Set a real secret via the SECRET_KEY env var.",
                    UserWarning,
                    stacklevel=2,
                )
        return self


@lru_cache
def get_settings() -> Settings:
    return Settings()
