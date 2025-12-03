from __future__ import annotations

from typing import List

from ..models import WorkerRequest, WorkerResult
from .base import Worker


class MarketWorker(Worker):
    def __init__(self, retriever):
        super().__init__("market", retriever)

    def build_summary(self, request: WorkerRequest, passages: List[dict]) -> WorkerResult:
        score = 72
        summary = (
            f"Estimated treatable population supports a market score of {score}; "
            "competition limited to two incumbents."
        )
        return WorkerResult(
            summary=summary,
            evidence=self._to_evidence(passages, "market"),
            confidence=0.6,
            metadata={"market_score": str(score)},
        )
