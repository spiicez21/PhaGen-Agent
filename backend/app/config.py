from functools import lru_cache
from pydantic import BaseSettings, Field


class Settings(BaseSettings):
    api_prefix: str = Field(default="/api")
    rag_top_k: int = Field(default=5)
    job_ttl_minutes: int = Field(default=60)

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"


@lru_cache
def get_settings() -> Settings:
    return Settings()
