from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional


class Recommendation(str, Enum):
    go = "Go"
    investigate = "Investigate"
    no_go = "No-Go"
    insufficient = "Insufficient"


@dataclass
class EvidenceItem:
    type: str
    text: str
    url: str
    confidence: float


@dataclass
class WorkerResult:
    summary: str
    evidence: List[EvidenceItem]
    confidence: float
    metadata: Dict[str, str] = field(default_factory=dict)


@dataclass
class WorkerRequest:
    molecule: str
    synonyms: List[str]
    top_k: int = 5
    context_tokens: int = 1200


@dataclass
class MasterResult:
    innovation_story: str
    recommendation: Recommendation
    market_score: int
    workers: Dict[str, WorkerResult]


@dataclass
class WorkerFailure:
    worker_name: str
    reason: str


@dataclass
class MasterRun:
    success: bool
    output: Optional[MasterResult] = None
    failures: List[WorkerFailure] = field(default_factory=list)
