from __future__ import annotations

import json
import logging
import time
from abc import ABC, abstractmethod
from typing import List, Optional

from ..api_budget import API_PROVIDER_FOR_WORKER, api_budget_monitor
from ..llm import LLMClient, LLMClientError, format_structured_context
from ..models import EvidenceItem, WorkerRequest, WorkerResult
from ..retrieval import Retriever
from ..temperature import resolve_worker_temperature


logger = logging.getLogger(__name__)


class WorkerSLAMetrics:
    """SLA metrics tracker for worker performance monitoring."""
    def __init__(self, worker_name: str):
        self.worker_name = worker_name
        self.start_time = time.time()
        self.retrieval_time = 0.0
        self.llm_time = 0.0
        self.total_tokens = 0
        self.retrieval_calls = 0
        self.llm_calls = 0
        self.failures = []
    
    def record_retrieval(self, duration: float):
        self.retrieval_time += duration
        self.retrieval_calls += 1
    
    def record_llm(self, duration: float, tokens: int):
        self.llm_time += duration
        self.llm_calls += 1
        self.total_tokens += tokens
    
    def record_failure(self, stage: str, error: str):
        self.failures.append({"stage": stage, "error": str(error)[:200]})
    
    def finalize(self) -> dict:
        total_time = time.time() - self.start_time
        return {
            "worker": self.worker_name,
            "total_latency_ms": round(total_time * 1000, 2),
            "retrieval_latency_ms": round(self.retrieval_time * 1000, 2),
            "llm_latency_ms": round(self.llm_time * 1000, 2),
            "retrieval_calls": self.retrieval_calls,
            "llm_calls": self.llm_calls,
            "total_tokens": self.total_tokens,
            "failure_count": len(self.failures),
            "failures": self.failures if self.failures else None,
        }
    
    def log_summary(self):
        metrics = self.finalize()
        logger.info(
            f"[SLA] {self.worker_name}: "
            f"latency={metrics['total_latency_ms']}ms "
            f"(retrieval={metrics['retrieval_latency_ms']}ms, "
            f"llm={metrics['llm_latency_ms']}ms), "
            f"tokens={metrics['total_tokens']}, "
            f"failures={metrics['failure_count']}"
        )


class Worker(ABC):
    SOURCE_PRIORITY = {
        "clinical": 0,
        "regulatory": 1,
        "literature": 2,
        "patent": 3,
    }
    CONFIDENCE_THRESHOLDS = (
        (0.8, "high"),
        (0.6, "medium"),
        (0.0, "low"),
    )
    MIN_EVIDENCE_THRESHOLD = 2
    MIN_COVERAGE_RATIO = 0.5
    MIN_PRECISION_RATIO = 0.5

    def __init__(
        self,
        name: str,
        retriever: Retriever,
        llm_client: Optional[LLMClient] = None,
        temperature: float | None = None,
    ) -> None:
        self.name = name
        self.retriever = retriever
        self.llm = llm_client
        self.temperature = resolve_worker_temperature(name, temperature)

    @abstractmethod
    def build_summary(self, request: WorkerRequest, passages: List[dict]) -> WorkerResult:
        raise NotImplementedError

    def run(self, request: WorkerRequest) -> WorkerResult:
        sla = WorkerSLAMetrics(self.name)
        logger.info(f"      \ud83d\udd0d [{self.name}] Building search queries for: {request.molecule}")
        
        try:
            queries = self._build_query_terms(request)
            logger.info(f"      \ud83d\udd0d [{self.name}] Generated {len(queries)} search queries")

            gathered: List[dict] = []
            seen_ids: set[str] = set()
            for term in queries:
                self._record_api_budget()
                retrieval_start = time.time()
                try:
                    results = self.retriever.search(
                        query=term,
                        source_type=self.name,
                        top_k=request.top_k,
                        max_tokens=request.context_tokens,
                    )
                    sla.record_retrieval(time.time() - retrieval_start)
                except Exception as e:
                    sla.record_failure("retrieval", str(e))
                    sla.record_retrieval(time.time() - retrieval_start)
                    raise
                
                for result in results:
                    record_id = result.get("id") or f"{term}-{result.get('rank', len(gathered))}"
                    if record_id in seen_ids:
                        continue
                    seen_ids.add(record_id)
                    gathered.append(result)
                    if len(gathered) >= request.top_k:
                        break
                if len(gathered) >= request.top_k:
                    break

            if not gathered:
                logger.info(f"      \ud83d\udd04 [{self.name}] No results found, using fallback query")
                self._record_api_budget()
                retrieval_start = time.time()
                try:
                    gathered = self.retriever.search(
                        query=request.molecule,
                        source_type=self.name,
                        top_k=request.top_k,
                        max_tokens=request.context_tokens,
                    )
                    sla.record_retrieval(time.time() - retrieval_start)
                except Exception as e:
                    sla.record_failure("retrieval_fallback", str(e))
                    sla.record_retrieval(time.time() - retrieval_start)
                    raise

            logger.info(f"      \ud83d\udccb [{self.name}] Retrieved {len(gathered)} passages from knowledge base")
            ordered = sorted(gathered, key=self._source_rank)
            limited = ordered[: request.top_k]
            
            # Track LLM usage in build_summary
            self._current_sla = sla
            logger.info(f"      \ud83e\udde0 [{self.name}] Generating summary with LLM...")
            result = self.build_summary(request, limited)
            logger.info(f"      \u2705 [{self.name}] Summary complete (confidence: {result.confidence:.2f})")
            
            result.metrics = self._compute_metrics(
                gathered=gathered,
                limited=limited,
                queries=queries,
                request=request,
                result=result,
            )
            result.alerts = self._evaluate_guardrails(result.metrics)
            
            # Add SLA metrics to result
            result.metrics["sla"] = sla.finalize()
            sla.log_summary()
            
            return result
        except Exception as e:
            sla.record_failure("worker_run", str(e))
            sla.log_summary()
            raise

    def _build_query_terms(self, request: WorkerRequest) -> List[str]:
        terms: List[str] = []
        seen: set[str] = set()
        for candidate in [request.molecule, *(request.synonyms or [])]:
            normalized = (candidate or "").strip()
            if not normalized:
                continue
            lowered = normalized.lower()
            if lowered in seen:
                continue
            seen.add(lowered)
            terms.append(normalized)
        return terms or [request.molecule]

    def _to_evidence(self, passages: List[dict], default_type: str) -> List[EvidenceItem]:
        return [
            EvidenceItem(
                type=passage.get("source_type", default_type),
                text=passage.get("snippet", ""),
                url=self._format_evidence_url(passage, default_type),
                confidence=max(0.4, 1 - (idx * 0.1)),
            )
            for idx, passage in enumerate(passages)
        ]
    
    def _format_evidence_url(self, passage: dict, source_type: str) -> str:
        """Format evidence URL with descriptive placeholder if missing."""
        url = passage.get("url", "").strip()
        if url and url != "https://example.org/evidence":
            return url
        # Generate descriptive source reference
        doc_id = passage.get("id", "unknown")
        origin = passage.get("origin", "indexed")
        return f"source://{source_type}/{origin}/{doc_id}"

    # Quality metrics & guardrails ---------------------------------
    def _compute_metrics(
        self,
        gathered: List[dict],
        limited: List[dict],
        queries: List[str],
        request: WorkerRequest,
        result: WorkerResult,
    ) -> dict:
        total_retrieved = len(gathered)
        final_passages = len(limited)
        evidence_count = len(result.evidence)
        unique_sources = len(
            {
                passage.get("url")
                or passage.get("id")
                or f"{idx}:{passage.get('rank', 0)}"
                for idx, passage in enumerate(limited)
            }
        )
        high_conf_count = sum(1 for item in result.evidence if item.confidence >= 0.8)
        coverage_ratio = final_passages / max(1, request.top_k)
        precision_proxy = (
            high_conf_count / max(1, evidence_count)
            if evidence_count
            else 0.0
        )
        return {
            "query_terms": len(queries),
            "retrieved_passages": total_retrieved,
            "final_passages": final_passages,
            "evidence_count": evidence_count,
            "unique_sources": unique_sources,
            "coverage_ratio": round(coverage_ratio, 3),
            "precision_proxy": round(precision_proxy, 3),
            "high_conf_evidence": high_conf_count,
            "retriever_top_k": request.top_k,
        }

    def _evaluate_guardrails(self, metrics: dict) -> List[str]:
        alerts: List[str] = []
        evidence_count = metrics.get("evidence_count", 0)
        coverage_ratio = metrics.get("coverage_ratio", 0.0)
        precision_proxy = metrics.get("precision_proxy", 0.0)
        unique_sources = metrics.get("unique_sources", 0)
        final_passages = metrics.get("final_passages", 0)
        retrieved = metrics.get("retrieved_passages", 0)

        if retrieved == 0:
            alerts.append("[ANOMALY] No passages retrieved from index")
        elif final_passages == 0:
            alerts.append("[ANOMALY] All passages filtered out by budget")

        if evidence_count < self.MIN_EVIDENCE_THRESHOLD:
            alerts.append(
                f"Low evidence coverage ({evidence_count}/{self.MIN_EVIDENCE_THRESHOLD})"
            )
        if coverage_ratio < self.MIN_COVERAGE_RATIO:
            alerts.append(
                f"Low retrieval coverage ({coverage_ratio:.2f} < {self.MIN_COVERAGE_RATIO:.2f})"
            )
        if evidence_count and precision_proxy < self.MIN_PRECISION_RATIO:
            alerts.append(
                f"Low precision proxy ({precision_proxy:.2f} < {self.MIN_PRECISION_RATIO:.2f})"
            )
        if unique_sources <= 1 and final_passages >= 2:
            alerts.append("Evidence monoculture (single unique source)")
        return alerts

    def _record_api_budget(self) -> None:
        provider = API_PROVIDER_FOR_WORKER.get(self.name)
        if provider:
            api_budget_monitor.record(provider)

    # LLM helpers ----------------------------------------------------
    def _format_passages(self, passages: List[dict], limit: int = 5) -> str:
        formatted: List[str] = []
        for idx, passage in enumerate(passages[:limit], start=1):
            snippet = (passage.get("snippet") or "").strip().replace("\n", " ")
            if len(snippet) > 400:
                snippet = snippet[:397] + "..."
            source = passage.get("source_type") or self.name
            url = passage.get("url") or ""
            citation = f"[{idx}] ({source}) {snippet}"
            if url:
                citation += f"\nSource: {url}"
            formatted.append(citation)
        return "\n\n".join(formatted) if formatted else "No evidence retrieved."

    def _summarize_with_llm(
        self,
        request: WorkerRequest,
        instructions: str,
        passages: List[dict],
        extra_context: object | None,
        fallback: str,
    ) -> str:
        if not self.llm:
            return fallback
        context = self._format_passages(passages)
        if extra_context:
            context += "\n\nStructured context:\n" + format_structured_context(extra_context)
        user_prompt = (
            f"Molecule: {request.molecule}\n"
            f"Synonyms: {', '.join(request.synonyms) if request.synonyms else 'n/a'}\n"
            f"Objective: {instructions}\n\n"
            f"Evidence:\n{context}\n\n"
            "Write 2-3 crisp sentences, cite references like [1], and stay factual."
        )
        
        llm_start = time.time()
        try:
            response = self.llm.generate(
                prompt=user_prompt,
                system_prompt=
                "You are an expert analyst in a multi-agent molecule repurposing pipeline.",
                temperature=self.temperature,
                max_tokens=280,
            )
            # Estimate tokens (4 chars per token heuristic)
            estimated_tokens = len(user_prompt + response) // 4
            if hasattr(self, '_current_sla'):
                self._current_sla.record_llm(time.time() - llm_start, estimated_tokens)
            return response
        except LLMClientError as exc:  # pragma: no cover - runtime dependency
            if hasattr(self, '_current_sla'):
                self._current_sla.record_failure("llm_summary", str(exc))
                self._current_sla.record_llm(time.time() - llm_start, 0)
            logger.warning("LLM summary failed for %s: %s", self.name, exc)
            return fallback

    def _source_rank(self, passage: dict) -> tuple[int, int]:
        source = (passage.get("source_type") or self.name).lower()
        priority = self.SOURCE_PRIORITY.get(source, 10)
        rank = passage.get("rank") or 999
        return priority, rank

    def _calibrate_confidence(self, score: float) -> tuple[float, str]:
        clamped = max(0.0, min(1.0, score))
        for threshold, label in self.CONFIDENCE_THRESHOLDS:
            if clamped >= threshold:
                return clamped, label
        return clamped, "low"
