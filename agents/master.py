from __future__ import annotations

from dataclasses import asdict
from typing import Dict

from .llm import LLMClient
from .models import (
    MasterResult,
    MasterRun,
    Recommendation,
    WorkerFailure,
    WorkerRequest,
    WorkerResult,
)
from .retrieval import Retriever
from .synonyms import SynonymExpander
from .workers.clinical import ClinicalWorker
from .workers.patent import PatentWorker
from .workers.literature import LiteratureWorker
from .workers.market import MarketWorker


class MasterAgent:
    def __init__(
        self,
        top_k: int = 5,
        context_tokens: int = 1200,
        llm_client: LLMClient | None = None,
        synonym_expander: SynonymExpander | None = None,
    ) -> None:
        retriever = Retriever(top_k=top_k, context_tokens=context_tokens)
        self.context_tokens = context_tokens
        self.llm = llm_client or LLMClient()
        self.synonyms = synonym_expander or SynonymExpander()
        self.workers = {
            "clinical": ClinicalWorker(retriever, self.llm),
            "patent": PatentWorker(retriever, self.llm),
            "literature": LiteratureWorker(retriever, self.llm),
            "market": MarketWorker(retriever, self.llm),
        }

    def run(
        self,
        molecule: str,
        synonyms: list[str] | None = None,
        smiles: str | None = None,
    ) -> MasterRun:
        expanded_synonyms = self.synonyms.expand(
            molecule,
            provided_synonyms=synonyms,
            smiles=smiles,
        )
        request = WorkerRequest(
            molecule=molecule,
            synonyms=expanded_synonyms,
            smiles=smiles,
            top_k=len(self.workers),
            context_tokens=self.context_tokens,
        )
        worker_outputs: Dict[str, WorkerResult] = {}
        failures: list[WorkerFailure] = []
        for name, worker in self.workers.items():
            try:
                worker_outputs[name] = worker.run(request)
            except Exception as exc:  # pragma: no cover - placeholder
                failures.append(WorkerFailure(worker_name=name, reason=str(exc)))
        if failures:
            return MasterRun(
                success=False,
                failures=failures,
            )
        market_score = int(
            float(worker_outputs["market"].metadata.get("market_score", "70"))
        )
        recommendation = self._recommend(market_score, worker_outputs)
        innovation_story = self._synthesize_story(molecule, worker_outputs)
        return MasterRun(
            success=True,
            output=MasterResult(
                innovation_story=innovation_story,
                recommendation=recommendation,
                market_score=market_score,
                workers=worker_outputs,
            ),
        )

    def _recommend(
        self, market_score: int, worker_outputs: Dict[str, WorkerResult]
    ) -> Recommendation:
        avg_conf = sum(worker.confidence for worker in worker_outputs.values()) / len(
            worker_outputs
        )
        if market_score >= 75 and avg_conf >= 0.7:
            return Recommendation.go
        if market_score >= 60:
            return Recommendation.investigate
        return Recommendation.no_go

    def _synthesize_story(
        self, molecule: str, worker_outputs: Dict[str, WorkerResult]
    ) -> str:
        clinical = worker_outputs["clinical"].summary
        patent = worker_outputs["patent"].summary
        literature = worker_outputs["literature"].summary
        market = worker_outputs["market"].summary
        return (
            f"{molecule} shows encouraging clinical signals ({clinical}). "
            f"Patent view: {patent}. Literature review: {literature}. "
            f"Market stance: {market}."
        )

    def serialize(self, result: MasterResult) -> dict:
        return asdict(result)
