from datetime import datetime
from enum import Enum
from typing import Dict, List, Optional

from pydantic import BaseModel, Field


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
    metadata: Dict[str, str] = Field(default_factory=dict)


class JobStatus(str, Enum):
    pending = "PENDING"
    running = "RUNNING"
    completed = "COMPLETED"
    failed = "FAILED"


class JobCreateRequest(BaseModel):
    molecule: str
    synonyms: Optional[List[str]] = None
    smiles: Optional[str] = None


class JobResponse(BaseModel):
    job_id: str
    status: JobStatus
    created_at: datetime
    updated_at: datetime
    payload: Dict[str, dict] | None = None
    recommendation: Optional[Recommendation] = None
