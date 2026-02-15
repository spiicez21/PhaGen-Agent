from datetime import datetime
from enum import Enum
from typing import Any, Dict, List, Literal, Optional

from pydantic import BaseModel, Field, field_validator


class Recommendation(str, Enum):
    go = "Go"
    investigate = "Investigate"
    no_go = "No-Go"


class EvidenceItem(BaseModel):
    type: str
    text: str
    url: str
    confidence: float = Field(ge=0.0, le=1.0)


class WorkerResult(BaseModel):
    summary: str
    evidence: List[EvidenceItem]
    confidence: float = Field(ge=0.0, le=1.0)
    confidence_band: Literal["low", "medium", "high"]
    metadata: Dict[str, str] = Field(default_factory=dict)


class JobStatus(str, Enum):
    pending = "PENDING"
    running = "RUNNING"
    completed = "COMPLETED"
    failed = "FAILED"


class JobCreateRequest(BaseModel):
    molecule: str = Field(..., min_length=1, max_length=200)
    synonyms: Optional[List[str]] = None
    smiles: Optional[str] = Field(default=None, max_length=500)

    @field_validator("molecule")
    @classmethod
    def _clean_molecule(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("Molecule name cannot be empty")
        return v

    @field_validator("synonyms")
    @classmethod
    def _clean_synonyms(cls, v: list[str] | None) -> list[str] | None:
        if v is None:
            return None
        if len(v) > 50:
            raise ValueError("Too many synonyms (max 50)")
        return [s.strip() for s in v if s.strip()]


class JobResponse(BaseModel):
    job_id: str
    status: JobStatus
    created_at: datetime
    updated_at: datetime
    payload: Dict[str, Any] | None = None
    recommendation: Optional[Recommendation] = None
    report_version: Optional[int] = None


# Authentication Schemas
class UserBase(BaseModel):
    email: str
    name: str


class UserCreate(UserBase):
    password: str = Field(min_length=8)


class UserLogin(BaseModel):
    email: str
    password: str


class UserResponse(UserBase):
    id: str
    is_active: bool
    created_at: datetime

    class Config:
        from_attributes = True


class Token(BaseModel):
    access_token: str
    token_type: str = "bearer"


class TokenData(BaseModel):
    email: Optional[str] = None
