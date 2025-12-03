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
        passages = self.retriever.search(
            query=request.molecule,
            source_type=self.name,
            top_k=request.top_k,
            max_tokens=request.context_tokens,
        )
        return self.build_summary(request, passages)

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
