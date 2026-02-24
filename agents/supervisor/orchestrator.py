"""Supervisor orchestrator — coordinates workers, quality checks, and synthesis.

Replaces the monolithic MasterAgent with a clean supervisor pattern.
Supports both classic sequential execution and ReAct-enhanced workers.
"""

from __future__ import annotations

import asyncio
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeout
from typing import Dict, List, Optional

from ..exceptions import SupervisorError, WorkerError
from ..llm import LLMClient
from ..logging_config import get_logger
from ..models import (
    MasterResult,
    MasterRun,
    Recommendation,
    WorkerFailure,
    WorkerRequest,
    WorkerResult,
)
from ..retrieval import Retriever
from ..synonyms import SynonymExpander
from ..temperature import resolve_master_temperature, resolve_worker_temperatures
from ..tools.clinical import build_clinical_tools
from ..tools.literature import build_literature_tools
from ..tools.market import build_market_tools
from ..tools.patent import build_patent_tools
from ..tools.registry import ToolRegistry
from ..workers.clinical import ClinicalWorker
from ..workers.literature import LiteratureWorker
from ..workers.market import MarketWorker
from ..workers.patent import PatentWorker
from .serializer import serialize_result
from .synthesizer import (
    fallback_story,
    recommend,
    synthesize_story_and_recommendation,
)

logger = get_logger(__name__)


class SupervisorAgent:
    """Orchestrates domain workers with quality checking and retry logic."""

    DEFAULT_WORKER_TIMEOUTS = {
        "clinical": 300.0,
        "patent": 300.0,
        "literature": 300.0,
        "market": 300.0,
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
        worker_temperatures: Dict[str, float] | None = None,
        master_temperature: float | None = None,
        enable_react: bool = True,
    ) -> None:
        retriever = Retriever(top_k=top_k, context_tokens=context_tokens)
        self.context_tokens = context_tokens
        self.llm = llm_client or LLMClient()
        self.synonyms = synonym_expander or SynonymExpander()
        self.worker_temperatures = resolve_worker_temperatures(worker_temperatures)
        self.master_temperature = resolve_master_temperature(master_temperature)

        # Build tool registries per domain
        tool_registries = self._build_tool_registries(retriever) if enable_react else {}

        self.workers = {
            "clinical": ClinicalWorker(
                retriever, self.llm,
                temperature=self.worker_temperatures.get("clinical"),
                tools=tool_registries.get("clinical"),
            ),
            "patent": PatentWorker(
                retriever, self.llm,
                temperature=self.worker_temperatures.get("patent"),
                tools=tool_registries.get("patent"),
            ),
            "literature": LiteratureWorker(
                retriever, self.llm,
                temperature=self.worker_temperatures.get("literature"),
                tools=tool_registries.get("literature"),
            ),
            "market": MarketWorker(
                retriever, self.llm,
                temperature=self.worker_temperatures.get("market"),
                tools=tool_registries.get("market"),
            ),
        }
        self.worker_timeouts = {**self.DEFAULT_WORKER_TIMEOUTS, **(worker_timeouts or {})}
        self.worker_retry_budget = {**self.DEFAULT_RETRY_BUDGET, **(worker_retry_budget or {})}

    def run(
        self,
        molecule: str,
        synonyms: list[str] | None = None,
        smiles: str | None = None,
    ) -> MasterRun:
        """Execute the full supervisor pipeline."""
        logger.info("supervisor_start", extra={"molecule": molecule})

        expanded_synonyms = self.synonyms.expand(
            molecule, provided_synonyms=synonyms, smiles=smiles,
        )
        request = WorkerRequest(
            molecule=molecule,
            synonyms=expanded_synonyms,
            smiles=smiles,
            top_k=len(self.workers),
            context_tokens=self.context_tokens,
        )

        # Phase 1: Dispatch workers
        worker_outputs, failures = self._dispatch_workers(request)
        if failures:
            return MasterRun(success=False, failures=failures)

        # Phase 2: Extract market score
        market_score = int(
            float(worker_outputs["market"].metadata.get("market_score", "70"))
        )

        # Phase 3: Generate fallback (deterministic)
        fb_story = fallback_story(molecule, worker_outputs)
        fb_recommendation = recommend(market_score, worker_outputs)

        # Phase 4: LLM synthesis
        try:
            innovation_story, recommendation = synthesize_story_and_recommendation(
                self.llm,
                molecule,
                worker_outputs,
                market_score,
                fb_story,
                fb_recommendation,
                temperature=self.master_temperature,
            )
        except (SupervisorError, Exception) as exc:
            logger.warning("synthesis_fallback", extra={"error": str(exc)})
            innovation_story = fb_story
            recommendation = fb_recommendation

        logger.info(
            "supervisor_complete",
            extra={
                "molecule": molecule,
                "recommendation": recommendation.value,
                "market_score": market_score,
                "workers": len(worker_outputs),
            },
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

    def serialize(self, result: MasterResult) -> dict:
        """Serialize a MasterResult with validation and quality checks."""
        return serialize_result(result)

    # -- Worker dispatch ---------------------------------------------------

    def _dispatch_workers(
        self, request: WorkerRequest
    ) -> tuple[Dict[str, WorkerResult], List[WorkerFailure]]:
        worker_outputs: Dict[str, WorkerResult] = {}
        failures: List[WorkerFailure] = []

        for name, worker in self.workers.items():
            try:
                logger.info("worker_dispatch", extra={"worker": name})
                worker_outputs[name] = self._run_with_retries(name, worker, request)
                logger.info(
                    "worker_success",
                    extra={"worker": name, "confidence": worker_outputs[name].confidence},
                )
            except Exception as exc:
                logger.error("worker_failed", extra={"worker": name, "error": str(exc)})
                failures.append(WorkerFailure(worker_name=name, reason=str(exc)))

        return worker_outputs, failures

    def _run_with_retries(self, name: str, worker, request: WorkerRequest) -> WorkerResult:
        retry_budget = max(0, self.worker_retry_budget.get(name, 0))
        attempts = 0
        last_error: Exception | None = None
        while attempts <= retry_budget:
            attempts += 1
            try:
                timeout = self.worker_timeouts.get(name, 30.0)
                return self._run_with_timeout(worker, request, timeout)
            except Exception as exc:
                last_error = exc
                logger.warning(
                    "worker_retry",
                    extra={"worker": name, "attempt": attempts, "max": retry_budget + 1, "error": str(exc)},
                )
        raise WorkerError(name, f"Failed after {attempts} attempts: {last_error}")

    def _run_with_timeout(self, worker, request: WorkerRequest, timeout_seconds: float) -> WorkerResult:
        executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix=f"{worker.name}-worker")
        future = executor.submit(worker.run, request)
        try:
            return future.result(timeout=timeout_seconds)
        except FuturesTimeout as exc:
            future.cancel()
            raise TimeoutError(f"{worker.name} exceeded {timeout_seconds:.1f}s timeout") from exc
        finally:
            executor.shutdown(wait=False, cancel_futures=True)

    # -- Tool registry building --------------------------------------------

    def _build_tool_registries(self, retriever: Retriever) -> Dict[str, ToolRegistry]:
        """Build per-domain tool registries."""
        registries: Dict[str, ToolRegistry] = {}

        clinical_reg = ToolRegistry()
        for tool in build_clinical_tools(retriever):
            clinical_reg.register(tool, domain="clinical")
        registries["clinical"] = clinical_reg

        patent_reg = ToolRegistry()
        for tool in build_patent_tools(retriever):
            patent_reg.register(tool, domain="patent")
        registries["patent"] = patent_reg

        lit_reg = ToolRegistry()
        for tool in build_literature_tools(retriever):
            lit_reg.register(tool, domain="literature")
        registries["literature"] = lit_reg

        market_reg = ToolRegistry()
        for tool in build_market_tools(retriever):
            market_reg.register(tool, domain="market")
        registries["market"] = market_reg

        return registries
