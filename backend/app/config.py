from functools import lru_cache
from pathlib import Path
from pydantic import AliasChoices, Field
from pydantic_settings import BaseSettings, SettingsConfigDict


PROJECT_ROOT_ENV = Path(__file__).resolve().parents[2] / ".env"


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_file=str(PROJECT_ROOT_ENV),
        env_file_encoding="utf-8",
        extra="allow",
        env_file_override=True,  # .env file overrides environment variables
    )
    api_prefix: str = Field(default="/api")
    rag_top_k: int = Field(default=5)
    job_ttl_minutes: int = Field(default=60)
    database_url: str = Field(
        default="postgresql+psycopg://phagen:phagen@localhost:5432/phagen",
        validation_alias=AliasChoices("DATABASE_URL", "SUPABASE_URL"),
        description="SQLAlchemy connection string for the operational database.",
    )
    database_echo: bool = Field(default=False)
    storage_provider: str = Field(default="s3")
    s3_endpoint_url: str = Field(default="http://localhost:9000")
    s3_region: str = Field(default="us-east-1")
    s3_access_key: str = Field(default="admin")
    s3_secret_key: str = Field(default="phagen123")
    s3_use_ssl: bool = Field(default=False)
    s3_raw_documents_bucket: str = Field(default="phagen-raw-documents")
    s3_reports_bucket: str = Field(default="phagen-report-artifacts")
    redis_url: str = Field(default="redis://localhost:6379/0")
    cache_enabled: bool = Field(default=True)
    celery_broker_url: str = Field(default="redis://localhost:6379/1")
    celery_result_backend: str = Field(default="redis://localhost:6379/1")
    
    # JWT Authentication
    secret_key: str = Field(
        default="your-secret-key-change-in-production-min-32-chars-long",
        description="Secret key for JWT token generation (minimum 32 characters)"
    )
    algorithm: str = Field(default="HS256")
    access_token_expire_minutes: int = Field(default=30)


@lru_cache
def get_settings() -> Settings:
    return Settings()
