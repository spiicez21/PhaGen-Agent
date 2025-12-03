from __future__ import annotations

import json
import logging
from abc import ABC, abstractmethod
from typing import List, Optional

from ..llm import LLMClient, LLMClientError, format_structured_context
from ..models import EvidenceItem, WorkerRequest, WorkerResult
from ..retrieval import Retriever


logger = logging.getLogger(__name__)


class Worker(ABC):
    SOURCE_PRIORITY = {
        "clinical": 0,
        "regulatory": 1,
        "literature": 2,
        "patent": 3,
    }

    def __init__(
        self,
        name: str,
        retriever: Retriever,
        llm_client: Optional[LLMClient] = None,
    ) -> None:
        self.name = name
        self.retriever = retriever
        self.llm = llm_client

    @abstractmethod
    def build_summary(self, request: WorkerRequest, passages: List[dict]) -> WorkerResult:
        raise NotImplementedError

    def run(self, request: WorkerRequest) -> WorkerResult:
        queries = [request.molecule]
        for synonym in request.synonyms:
            if synonym and synonym.lower() != request.molecule.lower():
                queries.append(synonym)

        gathered: List[dict] = []
        seen_ids: set[str] = set()
        for term in queries:
            results = self.retriever.search(
                query=term,
                source_type=self.name,
                top_k=request.top_k,
                max_tokens=request.context_tokens,
            )
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
            gathered = self.retriever.search(
                query=request.molecule,
                source_type=self.name,
                top_k=request.top_k,
                max_tokens=request.context_tokens,
            )

        ordered = sorted(gathered, key=self._source_rank)
        limited = ordered[: request.top_k]

        return self.build_summary(request, limited)

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
        try:
            return self.llm.generate(
                prompt=user_prompt,
                system_prompt=
                "You are an expert analyst in a multi-agent molecule repurposing pipeline.",
                temperature=0.15,
                max_tokens=280,
            )
        except LLMClientError as exc:  # pragma: no cover - runtime dependency
            logger.warning("LLM summary failed for %s: %s", self.name, exc)
            return fallback

    def _source_rank(self, passage: dict) -> tuple[int, int]:
        source = (passage.get("source_type") or self.name).lower()
        priority = self.SOURCE_PRIORITY.get(source, 10)
        rank = passage.get("rank") or 999
        return priority, rank
