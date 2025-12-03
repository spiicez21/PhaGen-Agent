from __future__ import annotations

import json
import logging
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeout
from dataclasses import asdict
from typing import Dict

from .llm import LLMClient, LLMClientError, format_structured_context
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
    DEFAULT_WORKER_TIMEOUTS = {
        "clinical": 45.0,
        "patent": 35.0,
        "literature": 35.0,
        "market": 30.0,
    }
    DEFAULT_RETRY_BUDGET = {
        "clinical": 1,
        "patent": 2,
        "literature": 1,
        "market": 1,
    }

    def __init__(
        self,
        top_k: int = 5,
        context_tokens: int = 1200,
        llm_client: LLMClient | None = None,
        synonym_expander: SynonymExpander | None = None,
        worker_timeouts: Dict[str, float] | None = None,
        worker_retry_budget: Dict[str, int] | None = None,
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
        self.worker_timeouts = {
            **self.DEFAULT_WORKER_TIMEOUTS,
            **(worker_timeouts or {}),
        }
        self.worker_retry_budget = {
            **self.DEFAULT_RETRY_BUDGET,
            **(worker_retry_budget or {}),
        }
        self._logger = logging.getLogger(__name__)

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
                worker_outputs[name] = self._run_with_retries(name, worker, request)
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
        fallback_story = self._fallback_story(molecule, worker_outputs)
        fallback_recommendation = self._recommend(market_score, worker_outputs)
        innovation_story, recommendation = self._synthesize_story_and_recommendation(
            molecule,
            worker_outputs,
            market_score,
            fallback_story,
            fallback_recommendation,
        )
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

    def _fallback_story(
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

    def _synthesize_story_and_recommendation(
        self,
        molecule: str,
        worker_outputs: Dict[str, WorkerResult],
        market_score: int,
        fallback_story: str,
        fallback_recommendation: Recommendation,
    ) -> tuple[str, Recommendation]:
        if not self.llm:
            return fallback_story, fallback_recommendation
        payload = {
            "molecule": molecule,
            "market_score": market_score,
            "workers": self._build_worker_payload(worker_outputs),
        }
        rubric = (
            "Go: compelling clinical efficacy, manageable risk, clear TAM; "
            "Investigate: mixed or emerging signals requiring validation; "
            "No-Go: weak efficacy, blocking IP/regulatory risk, or limited market."  # noqa: E501
        )
        user_prompt = (
            "You are the final reviewer who writes the Innovation Story and rubric-based recommendation.\n"
            "Context:\n"
            f"{format_structured_context(payload)}\n\n"
            "Respond as minified JSON with keys innovation_story, recommendation, rationale."
            " Recommendation must be one of Go, Investigate, No-Go per rubric: "
            f"{rubric}"
        )
        try:
            raw = self.llm.generate(
                prompt=user_prompt,
                system_prompt=
                "Synthesize cross-domain biotech diligence into a crisp Innovation Story and rubric-based recommendation.",
                temperature=0.2,
                max_tokens=400,
            )
            parsed = self._parse_master_response(raw)
        except (LLMClientError, ValueError) as exc:  # pragma: no cover - runtime
            self._logger.warning("Master synthesis fallback: %s", exc)
            return fallback_story, fallback_recommendation

        story = parsed.get("innovation_story") or fallback_story
        recommendation = self._resolve_recommendation(
            parsed.get("recommendation"),
            fallback_recommendation,
        )
        return story, recommendation

    def _build_worker_payload(
        self, worker_outputs: Dict[str, WorkerResult]
    ) -> Dict[str, dict]:
        payload: Dict[str, dict] = {}
        for name, result in worker_outputs.items():
            payload[name] = {
                "summary": result.summary,
                "confidence": round(result.confidence, 3),
                "confidence_band": result.confidence_band,
                "metadata": result.metadata,
                "evidence": [
                    {
                        "type": item.type,
                        "confidence": round(item.confidence, 3),
                        "text": (item.text or "")[:280],
                        "url": item.url,
                    }
                    for item in result.evidence[:3]
                ],
            }
        return payload

    def _parse_master_response(self, raw: str) -> dict:
        snippet = raw.strip()
        start = snippet.find("{")
        end = snippet.rfind("}")
        if start == -1 or end == -1 or end <= start:
            raise ValueError("Master synthesis did not return JSON")
        json_blob = snippet[start : end + 1]
        return json.loads(json_blob)

    def _resolve_recommendation(
        self,
        label: str | None,
        fallback: Recommendation,
    ) -> Recommendation:
        if not label:
            return fallback
        normalized = label.strip().lower()
        mapping = {
            "go": Recommendation.go,
            "investigate": Recommendation.investigate,
            "no-go": Recommendation.no_go,
            "no go": Recommendation.no_go,
        }
        return mapping.get(normalized, fallback)

    def serialize(self, result: MasterResult) -> dict:
        return asdict(result)

    # Execution helpers ------------------------------------------------
    def _run_with_retries(
        self,
        name: str,
        worker,
        request: WorkerRequest,
    ) -> WorkerResult:
        retry_budget = max(0, self.worker_retry_budget.get(name, 0))
        attempts = 0
        last_error: Exception | None = None
        while attempts <= retry_budget:
            attempts += 1
            try:
                timeout = self.worker_timeouts.get(name, 30.0)
                return self._run_with_timeout(worker, request, timeout)
            except Exception as exc:  # pragma: no cover - runtime dependent
                last_error = exc
                self._logger.warning(
                    "%s attempt %s/%s failed: %s",
                    name,
                    attempts,
                    retry_budget + 1,
                    exc,
                )
        raise RuntimeError(
            f"{name} failed after {attempts} attempts: {last_error}",
        )

    def _run_with_timeout(
        self,
        worker,
        request: WorkerRequest,
        timeout_seconds: float,
    ) -> WorkerResult:
        executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix=f"{worker.name}-worker")
        future = executor.submit(worker.run, request)
        try:
            return future.result(timeout=timeout_seconds)
        except FuturesTimeout as exc:
            future.cancel()
            raise TimeoutError(
                f"{worker.name} exceeded {timeout_seconds:.1f}s timeout",
            ) from exc
        finally:
            executor.shutdown(wait=False, cancel_futures=True)
