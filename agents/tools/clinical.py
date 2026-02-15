"""Clinical domain tools for the ReAct agent loop.

Tools:
  - chroma_clinical_search: Retrieves clinical evidence from the vector index.
  - pubmed_search: Searches PubMed via NCBI E-utilities.
  - clinicaltrials_search: Queries ClinicalTrials.gov v2 API.
  - parse_trial_phase: Extracts structured trial metadata from text.
"""

from __future__ import annotations

import re
from typing import Optional

import httpx
from pydantic import Field

from ..exceptions import ToolExecutionError
from ..logging_config import get_logger
from .base import Tool, ToolInput, ToolOutput
from .circuit_breaker import CircuitBreaker

logger = get_logger(__name__)

# Shared circuit breakers for external APIs
_pubmed_cb = CircuitBreaker("pubmed", failure_threshold=3, recovery_timeout=120.0)
_ctgov_cb = CircuitBreaker("clinicaltrials.gov", failure_threshold=3, recovery_timeout=120.0)


# -- Chroma retrieval tool -----------------------------------------------

class ChromaClinicalInput(ToolInput):
    query: str = Field(description="Search query for clinical evidence")
    top_k: int = Field(default=5, description="Number of results to return")


class ChromaClinicalSearchTool(Tool):
    name = "chroma_clinical_search"
    description = "Search the local clinical evidence vector index (ChromaDB) for relevant passages."

    def __init__(self, retriever) -> None:
        self._retriever = retriever

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = ChromaClinicalInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        try:
            results = self._retriever.search(
                query=inp.query,
                source_type="clinical",
                top_k=inp.top_k,
            )
            return ToolOutput(
                success=True,
                data=results,
                metadata={"count": len(results)},
            )
        except Exception as exc:
            return ToolOutput(success=False, error=str(exc))

    def _input_class(self):
        return ChromaClinicalInput


# -- PubMed search tool --------------------------------------------------

class PubMedInput(ToolInput):
    query: str = Field(description="PubMed search query (e.g. molecule name + condition)")
    max_results: int = Field(default=5, description="Maximum results to return")


class PubMedSearchTool(Tool):
    name = "pubmed_search"
    description = "Search PubMed for clinical studies and publications via NCBI E-utilities API."

    _ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    _ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = PubMedInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        try:
            results = _pubmed_cb.call(self._search, inp.query, inp.max_results)
            return ToolOutput(
                success=True,
                data=results,
                metadata={"count": len(results), "source": "pubmed"},
            )
        except Exception as exc:
            return ToolOutput(success=False, error=str(exc))

    def _search(self, query: str, max_results: int) -> list[dict]:
        search_resp = httpx.get(
            self._ESEARCH_URL,
            params={
                "db": "pubmed",
                "term": query,
                "retmax": max_results,
                "retmode": "json",
            },
            timeout=15.0,
        )
        search_resp.raise_for_status()
        ids = search_resp.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []

        summary_resp = httpx.get(
            self._ESUMMARY_URL,
            params={
                "db": "pubmed",
                "id": ",".join(ids),
                "retmode": "json",
            },
            timeout=15.0,
        )
        summary_resp.raise_for_status()
        result_data = summary_resp.json().get("result", {})

        articles = []
        for pmid in ids:
            info = result_data.get(pmid, {})
            if not isinstance(info, dict):
                continue
            articles.append({
                "pmid": pmid,
                "title": info.get("title", ""),
                "source": info.get("source", ""),
                "pubdate": info.get("pubdate", ""),
                "doi": next(
                    (aid.get("value", "") for aid in info.get("articleids", []) if aid.get("idtype") == "doi"),
                    "",
                ),
            })
        return articles

    def _input_class(self):
        return PubMedInput


# -- ClinicalTrials.gov tool ---------------------------------------------

class ClinicalTrialsInput(ToolInput):
    query: str = Field(description="Search query for clinical trials")
    max_results: int = Field(default=5, description="Maximum trials to return")


class ClinicalTrialsGovTool(Tool):
    name = "clinicaltrials_search"
    description = "Search ClinicalTrials.gov for registered clinical trials related to a molecule."

    _API_URL = "https://clinicaltrials.gov/api/v2/studies"

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = ClinicalTrialsInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        try:
            results = _ctgov_cb.call(self._search, inp.query, inp.max_results)
            return ToolOutput(
                success=True,
                data=results,
                metadata={"count": len(results), "source": "clinicaltrials.gov"},
            )
        except Exception as exc:
            return ToolOutput(success=False, error=str(exc))

    def _search(self, query: str, max_results: int) -> list[dict]:
        resp = httpx.get(
            self._API_URL,
            params={
                "query.term": query,
                "pageSize": max_results,
                "format": "json",
            },
            timeout=15.0,
        )
        resp.raise_for_status()
        studies = resp.json().get("studies", [])
        results = []
        for study in studies:
            protocol = study.get("protocolSection", {})
            ident = protocol.get("identificationModule", {})
            status_mod = protocol.get("statusModule", {})
            design = protocol.get("designModule", {})
            outcomes = protocol.get("outcomesModule", {})
            primary_outcomes = outcomes.get("primaryOutcomes", [])
            results.append({
                "nct_id": ident.get("nctId", ""),
                "title": ident.get("briefTitle", ""),
                "phase": ", ".join(design.get("phases", [])),
                "status": status_mod.get("overallStatus", ""),
                "primary_endpoint": primary_outcomes[0].get("measure", "") if primary_outcomes else "",
            })
        return results

    def _input_class(self):
        return ClinicalTrialsInput


# -- Trial phase parser tool ---------------------------------------------

class TrialParseInput(ToolInput):
    text: str = Field(description="Raw text passage to extract trial metadata from")


class TrialPhaseParseTool(Tool):
    name = "parse_trial_phase"
    description = "Extract structured trial metadata (NCT ID, phase, status, endpoint, population) from a text passage."

    _STATUS_LIST = [
        "Recruiting", "Active, not recruiting", "Completed", "Terminated",
        "Withdrawn", "Suspended", "Not yet recruiting", "Enrolling by invitation",
    ]

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = TrialParseInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        text = inp.text
        result = {
            "nct_id": self._find_nct(text),
            "phase": self._find_phase(text),
            "status": self._find_status(text),
            "primary_endpoint": self._find_endpoint(text),
            "population": self._find_population(text),
        }
        has_data = any(v for v in result.values())
        return ToolOutput(success=has_data, data=result)

    def _find_nct(self, text: str) -> str:
        match = re.search(r"NCT\d{8}", text, re.IGNORECASE)
        return match.group(0).upper() if match else ""

    def _find_phase(self, text: str) -> str:
        match = re.search(r"Phase\s+([0-4](?:[ab])?|I{1,3}V?)", text, re.IGNORECASE)
        if not match:
            return ""
        token = match.group(1)
        roman_map = {"I": "1", "II": "2", "III": "3", "IV": "4"}
        normalized = roman_map.get(token.upper(), token)
        return f"Phase {normalized}"

    def _find_status(self, text: str) -> str:
        for status in self._STATUS_LIST:
            if re.search(re.escape(status), text, re.IGNORECASE):
                return status
        return ""

    def _find_endpoint(self, text: str) -> str:
        match = re.search(
            r"primary\s+(?:endpoint|outcome)[:\-]?\s*(.+?)(?:\. |;|,|$)",
            text, re.IGNORECASE,
        )
        return match.group(1).strip() if match else ""

    def _find_population(self, text: str) -> str:
        match = re.search(r"in\s+([A-Za-z0-9 \-/]+?)\s+patients", text, re.IGNORECASE)
        return match.group(1).strip() if match else ""

    def _input_class(self):
        return TrialParseInput


def build_clinical_tools(retriever) -> list[Tool]:
    """Factory: create all clinical domain tools."""
    return [
        ChromaClinicalSearchTool(retriever),
        PubMedSearchTool(),
        ClinicalTrialsGovTool(),
        TrialPhaseParseTool(),
    ]
