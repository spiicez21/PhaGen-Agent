"""Backward-compatible MasterAgent — thin wrapper around SupervisorAgent.

All orchestration, quality checking, and synthesis logic has been moved to
the ``agents.supervisor`` package.  This module exists solely for import
compatibility with existing code (backend routers, CLI, tests).
"""

from __future__ import annotations

from typing import Dict

from .llm import LLMClient
from .models import MasterResult, MasterRun
from .supervisor.orchestrator import SupervisorAgent
from .synonyms import SynonymExpander


class MasterAgent:
    """Backward-compatible facade delegating to SupervisorAgent."""

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
        self._supervisor = SupervisorAgent(
            top_k=top_k,
            context_tokens=context_tokens,
            llm_client=llm_client,
            synonym_expander=synonym_expander,
            worker_timeouts=worker_timeouts,
            worker_retry_budget=worker_retry_budget,
            worker_temperatures=worker_temperatures,
            master_temperature=master_temperature,
            enable_react=enable_react,
        )

    def run(
        self,
        molecule: str,
        synonyms: list[str] | None = None,
        smiles: str | None = None,
    ) -> MasterRun:
        return self._supervisor.run(molecule, synonyms, smiles)

    def serialize(self, result: MasterResult) -> dict:
        return self._supervisor.serialize(result)
