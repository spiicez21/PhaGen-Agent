from __future__ import annotations

from typing import List

from ..models import WorkerRequest, WorkerResult
from .base import Worker


class PatentWorker(Worker):
    def __init__(self, retriever):
        super().__init__("patent", retriever)

    def build_summary(self, request: WorkerRequest, passages: List[dict]) -> WorkerResult:
        blockers = [
            {
                "assignee": "Example Pharma",
                "priority_date": "2018-05-01",
                "claim": f"Use of {request.molecule} for inflammatory disorders",
            }
        ]
        summary = (
            "Single broad claim detected; prior art window suggests moderate risk."
        )
        return WorkerResult(
            summary=summary,
            evidence=self._to_evidence(passages, "patent"),
            confidence=0.5,
            metadata={"possible_blockers": str(blockers)},
        )
