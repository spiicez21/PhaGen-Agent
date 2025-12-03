from __future__ import annotations

from typing import List

from ..models import WorkerRequest, WorkerResult
from .base import Worker


class LiteratureWorker(Worker):
    def __init__(self, retriever):
        super().__init__("literature", retriever)

    def build_summary(self, request: WorkerRequest, passages: List[dict]) -> WorkerResult:
        summary = (
            f"Literature suggests {request.molecule} modulates cytokine response; "
            "multiple preclinical models show efficacy."
        )
        return WorkerResult(
            summary=summary,
            evidence=self._to_evidence(passages, "literature"),
            confidence=0.8,
            metadata={"evidence_count": str(len(passages))},
        )
