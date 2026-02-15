from __future__ import annotations

from enum import Enum
from typing import Dict, List, Optional

from pydantic import BaseModel, Field


class Recommendation(str, Enum):
    go = "Go"
    investigate = "Investigate"
    no_go = "No-Go"
    insufficient = "Insufficient"


class EvidenceItem(BaseModel):
    type: str
    text: str
    url: str
    confidence: float = Field(ge=0.0, le=1.0)


class WorkerResult(BaseModel):
    summary: str
    evidence: List[EvidenceItem]
    confidence: float = Field(ge=0.0, le=1.0)
    confidence_band: str
    metadata: Dict[str, str] = Field(default_factory=dict)
    metrics: Dict[str, float] = Field(default_factory=dict)
    alerts: List[str] = Field(default_factory=list)


class WorkerRequest(BaseModel):
    molecule: str
    synonyms: List[str] = Field(default_factory=list)
    smiles: Optional[str] = None
    top_k: int = 5
    context_tokens: int = 1200


class MasterResult(BaseModel):
    innovation_story: str
    recommendation: Recommendation
    market_score: int
    workers: Dict[str, WorkerResult]


class WorkerFailure(BaseModel):
    worker_name: str
    reason: str


class MasterRun(BaseModel):
    success: bool
    output: Optional[MasterResult] = None
    failures: List[WorkerFailure] = Field(default_factory=list)
