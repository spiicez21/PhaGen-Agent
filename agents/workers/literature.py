from __future__ import annotations

import json
from typing import Dict, List, Optional

from ..models import EvidenceItem, WorkerRequest, WorkerResult
from .base import Worker


class LiteratureWorker(Worker):
    def __init__(self, retriever, llm=None):
        super().__init__("literature", retriever, llm)

    def build_summary(self, request: WorkerRequest, passages: List[dict]) -> WorkerResult:
        articles = self._extract_articles(passages)
        fallback_summary = self._summarize(request.molecule, articles)
        summary = self._summarize_with_llm(
            request,
            instructions=
            "Summarize mechanism-of-action, preclinical/clinical findings, and safety notes from the literature evidence.",
            passages=passages,
            extra_context={"articles": articles} if articles else None,
            fallback=fallback_summary,
        )
        confidence, confidence_band = self._calibrate_confidence(
            self._score_confidence(articles)
        )
        return WorkerResult(
            summary=summary,
            evidence=self._build_evidence(articles, passages),
            confidence=confidence,
            confidence_band=confidence_band,
            metadata=self._build_metadata(articles),
        )

    # Parsing helpers -------------------------------------------------
    def _extract_articles(self, passages: List[dict]) -> List[Dict[str, str]]:
        articles: List[Dict[str, str]] = []
        for passage in passages:
            parsed = self._parse_passage(passage)
            if parsed:
                articles.append(parsed)
        return articles

    def _parse_passage(self, passage: dict) -> Optional[Dict[str, str]]:
        snippet = (passage.get("snippet") or "").strip()
        if not snippet:
            return None
        data = self._load_json(snippet)
        if data:
            article = self._article_from_payload(data)
        else:
            article = None
        if not article:
            article = {
                "title": passage.get("title") or snippet[:80],
                "mechanism": "",
                "summary": snippet,
                "doi": "",
            }
        if passage.get("url"):
            article.setdefault("url", passage["url"])
        article.setdefault("summary", snippet)
        return article

    def _load_json(self, raw: str) -> Optional[dict]:
        try:
            return json.loads(raw)
        except (json.JSONDecodeError, TypeError):
            return None

    def _article_from_payload(self, data: dict) -> Optional[Dict[str, str]]:
        payload = data.get("payload") if "payload" in data else data
        article = payload.get("article") if isinstance(payload, dict) else None
        if not article:
            article = payload if isinstance(payload, dict) else None
        if not isinstance(article, dict):
            return None
        summary = ""
        key_findings = article.get("key_findings")
        if isinstance(key_findings, list) and key_findings:
            summary = key_findings[0]
            if len(key_findings) > 1:
                summary += f"; {key_findings[1]}"
        summary = summary or article.get("snippet", "") or article.get("text", "")
        mechanism = article.get("mechanism") or article.get("focus") or ""
        return {
            "title": article.get("title", ""),
            "mechanism": mechanism,
            "summary": summary,
            "doi": article.get("doi", ""),
            "url": article.get("url", ""),
            "journal": article.get("journal", ""),
            "year": str(article.get("year", "")),
        }

    # Reporting helpers -----------------------------------------------
    def _summarize(self, molecule: str, articles: List[Dict[str, str]]) -> str:
        if not articles:
            return f"No peer-reviewed signals parsed for {molecule}; re-run once literature corpora are indexed."
        mechanisms = [article.get("mechanism") for article in articles if article.get("mechanism")]
        mechanism_text = mechanisms[0] if mechanisms else "broad anti-fibrotic activity"
        cited = sum(1 for article in articles if article.get("doi"))
        return (
            f"Reviewed {len(articles)} publications; {molecule} shows {mechanism_text}. "
            f"{cited} include DOI-backed citations."
        )

    def _build_evidence(self, articles: List[Dict[str, str]], passages: List[dict]) -> List[EvidenceItem]:
        if articles:
            evidence: List[EvidenceItem] = []
            for idx, article in enumerate(articles):
                cite = article.get("doi") or article.get("journal") or "Peer-reviewed source"
                text = f"{article.get('title') or 'Literature evidence'} — {article.get('summary')}"
                evidence.append(
                    EvidenceItem(
                        type="Literature",
                        text=text,
                        url=article.get("url", ""),
                        confidence=max(0.45, 0.85 - idx * 0.05),
                    )
                )
            return evidence
        return self._to_evidence(passages, "literature")

    def _score_confidence(self, articles: List[Dict[str, str]]) -> float:
        if not articles:
            return 0.55
        doi_ratio = sum(1 for article in articles if article.get("doi")) / len(articles)
        return min(0.9, 0.6 + 0.2 * doi_ratio + 0.03 * min(len(articles), 5))

    def _build_metadata(self, articles: List[Dict[str, str]]) -> Dict[str, str]:
        metadata: Dict[str, str] = {
            "evidence_count": str(len(articles)),
        }
        dois = [article.get("doi") for article in articles if article.get("doi")]
        if dois:
            metadata["citations"] = json.dumps(dois, ensure_ascii=False)
        journals = sorted(
            {article.get("journal") for article in articles if article.get("journal")}
        )
        if journals:
            metadata["journals"] = ", ".join(journals)
        return {k: v for k, v in metadata.items() if v}
