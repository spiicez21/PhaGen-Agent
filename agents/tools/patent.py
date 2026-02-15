"""Patent and regulatory domain tools for the ReAct agent loop.

Tools:
  - chroma_patent_search: Retrieves patent evidence from the vector index.
  - chroma_regulatory_search: Retrieves regulatory evidence from the vector index.
  - extract_patent_claims: Extracts assignee, priority date, blocking claims from text.
  - ip_risk_score: Scores freedom-to-operate risk from patent data.
"""

from __future__ import annotations

import json
import re
from html import unescape
from typing import Optional

from pydantic import Field

from ..logging_config import get_logger
from .base import Tool, ToolInput, ToolOutput

logger = get_logger(__name__)


# -- Chroma patent retrieval -----------------------------------------------

class ChromaPatentInput(ToolInput):
    query: str = Field(description="Search query for patent evidence")
    top_k: int = Field(default=5, description="Number of results")


class ChromaPatentSearchTool(Tool):
    name = "chroma_patent_search"
    description = "Search the local patent evidence vector index for relevant filings."

    def __init__(self, retriever) -> None:
        self._retriever = retriever

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = ChromaPatentInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        try:
            results = self._retriever.search(query=inp.query, source_type="patent", top_k=inp.top_k)
            return ToolOutput(success=True, data=results, metadata={"count": len(results)})
        except Exception as exc:
            return ToolOutput(success=False, error=str(exc))

    def _input_class(self):
        return ChromaPatentInput


# -- Chroma regulatory retrieval -------------------------------------------

class ChromaRegulatoryInput(ToolInput):
    query: str = Field(description="Search query for regulatory evidence")
    top_k: int = Field(default=3, description="Number of results")


class ChromaRegulatorySearchTool(Tool):
    name = "chroma_regulatory_search"
    description = "Search the local regulatory evidence vector index for label notes and contraindications."

    def __init__(self, retriever) -> None:
        self._retriever = retriever

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = ChromaRegulatoryInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        try:
            results = self._retriever.search(query=inp.query, source_type="regulatory", top_k=inp.top_k)
            return ToolOutput(success=True, data=results, metadata={"count": len(results)})
        except Exception as exc:
            return ToolOutput(success=False, error=str(exc))

    def _input_class(self):
        return ChromaRegulatoryInput


# -- Patent claim extractor ------------------------------------------------

class PatentClaimInput(ToolInput):
    text: str = Field(description="Raw patent snippet (JSON or HTML) to parse")


class PatentClaimExtractorTool(Tool):
    name = "extract_patent_claims"
    description = "Extract assignee, priority date, and blocking claims from a patent snippet."

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = PatentClaimInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        parsed = self._parse(inp.text)
        if parsed:
            return ToolOutput(success=True, data=parsed)
        return ToolOutput(success=False, error="Could not extract patent metadata from text")

    def _parse(self, snippet: str) -> Optional[dict]:
        try:
            data = json.loads(snippet)
            if isinstance(data, dict):
                return {
                    "assignee": data.get("assignee", ""),
                    "priority_date": data.get("priority_date", ""),
                    "blocking_claim": data.get("claim", ""),
                    "title": data.get("title", ""),
                }
        except (json.JSONDecodeError, TypeError):
            pass
        return self._parse_html(snippet)

    def _parse_html(self, html: str) -> Optional[dict]:
        assignees = self._extract_meta_values(html, scheme="assignee") or [self._extract_meta(html, "DC.contributor") or ""]
        priority = self._extract_meta(html, "DC.date", scheme="dateSubmitted") or self._extract_meta(html, "DC.date") or ""
        description = self._extract_meta(html, "DC.description") or self._extract_meta(html, "description") or ""
        title = self._extract_meta(html, "DC.title") or ""
        if not any([assignees[0] if assignees else "", priority, description]):
            return None
        claim = description.split(".")[0].strip() if description else ""
        return {
            "assignee": ", ".join(a for a in assignees if a),
            "priority_date": priority,
            "blocking_claim": unescape(claim),
            "title": title,
        }

    def _extract_meta(self, html: str, name: str, scheme: str | None = None) -> Optional[str]:
        pattern = rf'<meta[^>]*name="{re.escape(name)}"[^>]*>'
        for tag in re.finditer(pattern, html, re.IGNORECASE):
            snippet = tag.group(0)
            if scheme and not re.search(rf'scheme="{re.escape(scheme)}"', snippet, re.IGNORECASE):
                continue
            content_match = re.search(r'content="([^"]+)"', snippet, re.IGNORECASE)
            if content_match:
                return unescape(content_match.group(1))
        return None

    def _extract_meta_values(self, html: str, scheme: str) -> list[str]:
        pattern = rf'<meta[^>]*scheme="{re.escape(scheme)}"[^>]*content="([^"]+)"'
        return [unescape(m) for m in re.findall(pattern, html, re.IGNORECASE)]

    def _input_class(self):
        return PatentClaimInput


# -- IP risk scorer --------------------------------------------------------

class IPRiskInput(ToolInput):
    patents: list[dict] = Field(description="List of extracted patent records")
    regulatory_notes: list[str] = Field(default_factory=list, description="Regulatory guardrail notes")


class IPRiskScorerTool(Tool):
    name = "ip_risk_score"
    description = "Score freedom-to-operate risk based on extracted patent data and regulatory notes."

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = IPRiskInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        risk_score = self._score(inp.patents, inp.regulatory_notes)
        level = "high" if risk_score >= 0.7 else "medium" if risk_score >= 0.4 else "low"
        return ToolOutput(
            success=True,
            data={
                "risk_score": round(risk_score, 3),
                "risk_level": level,
                "patent_count": len(inp.patents),
                "regulatory_flags": len(inp.regulatory_notes),
            },
        )

    def _score(self, patents: list[dict], notes: list[str]) -> float:
        score = 0.2
        if patents:
            score += 0.15 * min(len(patents), 3)
            if any(p.get("blocking_claim") for p in patents):
                score += 0.2
        if notes:
            score += 0.1 * min(len(notes), 2)
        return min(1.0, score)

    def _input_class(self):
        return IPRiskInput


def build_patent_tools(retriever) -> list[Tool]:
    """Factory: create all patent domain tools."""
    return [
        ChromaPatentSearchTool(retriever),
        ChromaRegulatorySearchTool(retriever),
        PatentClaimExtractorTool(),
        IPRiskScorerTool(),
    ]
