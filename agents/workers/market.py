from __future__ import annotations

import json
from typing import Dict, List, Optional

from ..models import WorkerRequest, WorkerResult
from .base import Worker


class MarketWorker(Worker):
    def __init__(self, retriever):
        super().__init__("market", retriever)

    def build_summary(self, request: WorkerRequest, passages: List[dict]) -> WorkerResult:
        market_data = self._extract_market(passages)
        score = self._calculate_score(market_data)
        summary = self._summarize(request.molecule, market_data, score)
        return WorkerResult(
            summary=summary,
            evidence=self._to_evidence(passages, "market"),
            confidence=self._score_confidence(market_data),
            metadata=self._build_metadata(market_data, score),
        )

    # Parsing helpers -------------------------------------------------
    def _extract_market(self, passages: List[dict]) -> List[Dict[str, str]]:
        records: List[Dict[str, str]] = []
        for passage in passages:
            parsed = self._parse_passage(passage)
            if parsed:
                records.append(parsed)
        return records

    def _parse_passage(self, passage: dict) -> Optional[Dict[str, str]]:
        snippet = (passage.get("snippet") or "").strip()
        if not snippet:
            return None
        data = self._load_json(snippet)
        if isinstance(data, dict):
            market = data.get("market") or data
            if isinstance(market, dict):
                return {
                    "tam": self._to_str(market.get("tam")),
                    "incidence": self._to_str(market.get("incidence")),
                    "growth": self._to_str(market.get("growth")),
                    "competition": self._to_str(market.get("competition")),
                    "notes": market.get("notes", snippet),
                }
        return {
            "tam": self._extract_number(snippet, unit="$"),
            "incidence": self._extract_number(snippet, unit="pts"),
            "growth": self._extract_percentage(snippet),
            "competition": "limited" if "limited" in snippet.lower() else "",
            "notes": snippet,
        }

    def _load_json(self, raw: str) -> Optional[dict]:
        try:
            return json.loads(raw)
        except (json.JSONDecodeError, TypeError):
            return None

    def _to_str(self, value) -> str:
        if value is None:
            return ""
        return str(value)

    def _extract_number(self, text: str, unit: str) -> str:
        unit_lower = unit.lower()
        for token in text.split():
            if unit_lower in token.lower():
                return token
        return ""

    def _extract_percentage(self, text: str) -> str:
        for token in text.split():
            if token.endswith("%"):
                return token
        return ""

    # Reporting helpers -----------------------------------------------
    def _calculate_score(self, records: List[Dict[str, str]]) -> int:
        if not records:
            return 60
        score = 65
        for record in records:
            if record.get("tam"):
                score += 5
            if record.get("incidence"):
                score += 5
            if record.get("competition") and "limited" in record["competition"].lower():
                score += 5
        return min(score, 90)

    def _summarize(self, molecule: str, records: List[Dict[str, str]], score: int) -> str:
        if not records:
            return f"Insufficient market intel for {molecule}; ingest IQVIA/analog datasets to unblock."
        tam = next((r["tam"] for r in records if r.get("tam")), "$?" )
        competition = next((r["competition"] for r in records if r.get("competition")), "limited head-to-head data")
        return (
            f"Derived market score {score} for {molecule}; TAM around {tam} with {competition}."
        )

    def _score_confidence(self, records: List[Dict[str, str]]) -> float:
        if not records:
            return 0.5
        fields = ["tam", "incidence", "competition"]
        coverage = sum(1 for record in records for field in fields if record.get(field))
        max_coverage = len(records) * len(fields)
        ratio = coverage / max_coverage if max_coverage else 0
        return min(0.85, 0.55 + 0.3 * ratio)

    def _build_metadata(self, records: List[Dict[str, str]], score: int) -> Dict[str, str]:
        metadata = {"market_score": str(score)}
        tam = next((r["tam"] for r in records if r.get("tam")), "")
        if tam:
            metadata["tam"] = tam
        incidence = next((r["incidence"] for r in records if r.get("incidence")), "")
        if incidence:
            metadata["incidence"] = incidence
        competitors = next((r["competition"] for r in records if r.get("competition")), "")
        if competitors:
            metadata["competition"] = competitors
        return metadata
