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

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"


@lru_cache
def get_settings() -> Settings:
    return Settings()
