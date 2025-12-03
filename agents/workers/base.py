from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List

from ..models import EvidenceItem, WorkerRequest, WorkerResult
from ..retrieval import Retriever


class Worker(ABC):
    def __init__(self, name: str, retriever: Retriever) -> None:
        self.name = name
        self.retriever = retriever

    @abstractmethod
    def build_summary(self, request: WorkerRequest, passages: List[dict]) -> WorkerResult:
        raise NotImplementedError

    def run(self, request: WorkerRequest) -> WorkerResult:
        passages = self.retriever.search(
            query=request.molecule,
            source_type=self.name,
            top_k=request.top_k,
            max_tokens=request.context_tokens,
        )
        return self.build_summary(request, passages)

    def _to_evidence(self, passages: List[dict], default_type: str) -> List[EvidenceItem]:
        return [
            EvidenceItem(
                type=passage.get("source_type", default_type),
                text=passage.get("snippet", ""),
                url=passage.get("url", ""),
                confidence=max(0.4, 1 - (idx * 0.1)),
            )
            for idx, passage in enumerate(passages)
        ]
