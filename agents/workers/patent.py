from __future__ import annotations

import json
import re
from html import unescape
from typing import Dict, List, Optional

from ..models import WorkerRequest, WorkerResult
from .base import Worker


class PatentWorker(Worker):
    def __init__(self, retriever):
        super().__init__("patent", retriever)

    def run(self, request: WorkerRequest) -> WorkerResult:
        primary = self.retriever.search(
            query=request.molecule,
            source_type="patent",
            top_k=request.top_k,
            max_tokens=request.context_tokens,
        )
        regulatory = self.retriever.search(
            query=request.molecule,
            source_type="regulatory",
            top_k=max(1, request.top_k // 2),
            max_tokens=max(200, request.context_tokens // 2),
        )
        return self.build_summary(request, primary + regulatory)

    def build_summary(self, request: WorkerRequest, passages: List[dict]) -> WorkerResult:
        patents = self._extract_patents(passages)
        guardrails = self._extract_regulatory(passages)
        summary = self._summarize(request.molecule, patents, guardrails)
        return WorkerResult(
            summary=summary,
            evidence=self._to_evidence(passages, "patent"),
            confidence=self._score_confidence(patents, guardrails),
            metadata=self._build_metadata(patents, guardrails),
        )

    # Patent parsing --------------------------------------------------
    def _extract_patents(self, passages: List[dict]) -> List[Dict[str, str]]:
        records: List[Dict[str, str]] = []
        for passage in passages:
            if passage.get("source_type") != "patent":
                continue
            snippet = passage.get("snippet", "")
            parsed = self._parse_patent_snippet(snippet)
            if not parsed:
                continue
            if passage.get("url"):
                parsed.setdefault("url", passage["url"])
            records.append(parsed)
        return records

    def _parse_patent_snippet(self, snippet: str) -> Optional[Dict[str, str]]:
        structured = self._load_json(snippet)
        if isinstance(structured, dict):
            return {
                "assignee": structured.get("assignee", ""),
                "priority_date": structured.get("priority_date", ""),
                "blocking_claim": structured.get("claim", ""),
                "title": structured.get("title", ""),
            }
        return self._parse_patent_html(snippet)

    def _load_json(self, raw: str) -> Optional[dict]:
        try:
            return json.loads(raw)
        except (TypeError, json.JSONDecodeError):
            return None

    def _parse_patent_html(self, html: str) -> Optional[Dict[str, str]]:
        assignees = self._extract_meta_values(html, scheme="assignee") or self._extract_meta(html, "DC.contributor")
        priority = self._extract_meta(html, "DC.date", scheme="dateSubmitted") or self._extract_meta(html, "DC.date")
        description = self._extract_meta(html, "DC.description") or self._extract_meta(html, "description")
        title = self._extract_meta(html, "DC.title") or self._extract_title(html)
        if not any([assignees, priority, description]):
            return None
        claim = description.split(".")[0].strip() if description else ""
        return {
            "assignee": assignees if isinstance(assignees, str) else ", ".join(assignees),
            "priority_date": priority or "",
            "blocking_claim": unescape(claim),
            "title": title or "",
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

    def _extract_meta_values(self, html: str, scheme: str) -> List[str]:
        pattern = rf'<meta[^>]*scheme="{re.escape(scheme)}"[^>]*content="([^"]+)"'
        return [unescape(match) for match in re.findall(pattern, html, re.IGNORECASE)]

    def _extract_title(self, html: str) -> str:
        match = re.search(r"<title>([^<]+)</title>", html, re.IGNORECASE)
        return unescape(match.group(1).strip()) if match else ""

    # Regulatory parsing ----------------------------------------------
    def _extract_regulatory(self, passages: List[dict]) -> List[str]:
        notes: List[str] = []
        for passage in passages:
            if passage.get("source_type") != "regulatory":
                continue
            snippet = passage.get("snippet", "")
            payload = self._load_json(snippet)
            if not payload:
                continue
            reg = payload.get("regulatory") or payload
            highlights = reg.get("label_highlights") if isinstance(reg, dict) else None
            if highlights:
                notes.extend(highlights)
        return notes

    # Reporting helpers -----------------------------------------------
    def _summarize(self, molecule: str, patents: List[Dict[str, str]], notes: List[str]) -> str:
        if not patents and not notes:
            return f"No patent or regulatory blockers detected for {molecule}; additional sources required."
        if patents:
            lead = patents[0]
            assignee = lead.get("assignee") or "unknown assignee"
            phase = f"priority {lead.get('priority_date')}" if lead.get("priority_date") else "undated"
            claim = lead.get("blocking_claim") or "broad composition claim"
            patent_summary = f"Lead filing lists {assignee} with {phase}; claim focuses on {claim.lower()}."
        else:
            patent_summary = "No relevant patent filings retrieved."
        reg_summary = "Contraindication data not retrieved."
        if notes:
            reg_summary = f"Label notes highlight {notes[0].lower()}"
        return f"{patent_summary} {reg_summary}".strip()

    def _score_confidence(self, patents: List[Dict[str, str]], notes: List[str]) -> float:
        points = 0.4
        if patents:
            coverage = sum(1 for field in ("assignee", "priority_date", "blocking_claim") if patents[0].get(field))
            points += 0.15 * coverage
        if notes:
            points += 0.1
        return min(0.9, points)

    def _build_metadata(self, patents: List[Dict[str, str]], notes: List[str]) -> Dict[str, str]:
        metadata: Dict[str, str] = {}
        if patents:
            metadata["assignees"] = ", ".join(
                sorted({patent.get("assignee", "") for patent in patents if patent.get("assignee")})
            )
            metadata["priority_dates"] = ", ".join(
                sorted({patent.get("priority_date", "") for patent in patents if patent.get("priority_date")})
            )
            claims = [patent.get("blocking_claim") for patent in patents if patent.get("blocking_claim")]
            if claims:
                metadata["blocking_claims"] = json.dumps(claims, ensure_ascii=False)
        if notes:
            metadata["contraindications"] = json.dumps(notes, ensure_ascii=False)
        return {k: v for k, v in metadata.items() if v}
