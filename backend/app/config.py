from functools import lru_cache
from pydantic import Field
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    api_prefix: str = Field(default="/api")
    rag_top_k: int = Field(default=5)
    job_ttl_minutes: int = Field(default=60)
    database_url: str = Field(
        default="postgresql+psycopg://phagen:phagen@localhost:5432/phagen",
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

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"


@lru_cache
def get_settings() -> Settings:
    return Settings()
