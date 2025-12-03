from __future__ import annotations

from typing import List

from ..models import WorkerRequest, WorkerResult
from .base import Worker


class ClinicalWorker(Worker):
    def __init__(self, retriever):
        super().__init__("clinical", retriever)

    def build_summary(self, request: WorkerRequest, passages: List[dict]) -> WorkerResult:
        trials = [
            {
                "registry_id": f"NCT{1000 + idx}",
                "phase": "Phase 2",
                "status": "Completed",
                "indication": "Mock Indication",
            }
            for idx, _ in enumerate(passages)
        ]
        summary = (
            f"Identified {len(trials)} trials mentioning {request.molecule}; "
            "most advanced signal in Phase 2."
        )
        return WorkerResult(
            summary=summary,
            evidence=self._to_evidence(passages, "clinical"),
            confidence=0.7,
            metadata={"trials": str(trials)},
        )
