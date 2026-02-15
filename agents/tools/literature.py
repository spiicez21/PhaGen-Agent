"""Literature domain tools for the ReAct agent loop.

Tools:
  - chroma_literature_search: Retrieves literature evidence from the vector index.
  - pubmed_central_search: Searches PubMed Central for full-text articles.
  - doi_resolve: Resolves a DOI to fetch article metadata.
  - citation_graph: Fetches citing/cited papers via Semantic Scholar.
"""

from __future__ import annotations

import httpx
from pydantic import Field

from ..logging_config import get_logger
from .base import Tool, ToolInput, ToolOutput
from .circuit_breaker import CircuitBreaker

logger = get_logger(__name__)

_pmc_cb = CircuitBreaker("pmc", failure_threshold=3, recovery_timeout=120.0)
_doi_cb = CircuitBreaker("doi.org", failure_threshold=3, recovery_timeout=120.0)
_s2_cb = CircuitBreaker("semantic_scholar", failure_threshold=3, recovery_timeout=120.0)


# -- Chroma literature retrieval -------------------------------------------

class ChromaLitInput(ToolInput):
    query: str = Field(description="Search query for literature evidence")
    top_k: int = Field(default=5, description="Number of results")


class ChromaLitSearchTool(Tool):
    name = "chroma_literature_search"
    description = "Search the local literature evidence vector index for peer-reviewed passages."

    def __init__(self, retriever) -> None:
        self._retriever = retriever

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = ChromaLitInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        try:
            results = self._retriever.search(query=inp.query, source_type="literature", top_k=inp.top_k)
            return ToolOutput(success=True, data=results, metadata={"count": len(results)})
        except Exception as exc:
            return ToolOutput(success=False, error=str(exc))

    def _input_class(self):
        return ChromaLitInput


# -- PubMed Central search -------------------------------------------------

class PMCInput(ToolInput):
    query: str = Field(description="Full-text search query for PubMed Central")
    max_results: int = Field(default=5, description="Maximum articles to return")


class PubMedCentralTool(Tool):
    name = "pubmed_central_search"
    description = "Search PubMed Central for full-text open-access articles."

    _ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    _ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = PMCInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        try:
            results = _pmc_cb.call(self._search, inp.query, inp.max_results)
            return ToolOutput(success=True, data=results, metadata={"count": len(results), "source": "pmc"})
        except Exception as exc:
            return ToolOutput(success=False, error=str(exc))

    def _search(self, query: str, max_results: int) -> list[dict]:
        search_resp = httpx.get(
            self._ESEARCH_URL,
            params={"db": "pmc", "term": query, "retmax": max_results, "retmode": "json"},
            timeout=15.0,
        )
        search_resp.raise_for_status()
        ids = search_resp.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        summary_resp = httpx.get(
            self._ESUMMARY_URL,
            params={"db": "pmc", "id": ",".join(ids), "retmode": "json"},
            timeout=15.0,
        )
        summary_resp.raise_for_status()
        result_data = summary_resp.json().get("result", {})
        articles = []
        for pmcid in ids:
            info = result_data.get(pmcid, {})
            if not isinstance(info, dict):
                continue
            articles.append({
                "pmcid": pmcid,
                "title": info.get("title", ""),
                "journal": info.get("source", ""),
                "pubdate": info.get("pubdate", ""),
                "doi": next(
                    (aid.get("value", "") for aid in info.get("articleids", []) if aid.get("idtype") == "doi"),
                    "",
                ),
            })
        return articles

    def _input_class(self):
        return PMCInput


# -- DOI resolver ----------------------------------------------------------

class DOIInput(ToolInput):
    doi: str = Field(description="DOI to resolve (e.g. 10.1234/example)")


class DOIResolverTool(Tool):
    name = "doi_resolve"
    description = "Resolve a DOI to fetch article metadata (title, journal, authors, year)."

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = DOIInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        try:
            result = _doi_cb.call(self._resolve, inp.doi)
            return ToolOutput(success=True, data=result)
        except Exception as exc:
            return ToolOutput(success=False, error=str(exc))

    def _resolve(self, doi: str) -> dict:
        resp = httpx.get(
            f"https://doi.org/{doi}",
            headers={"Accept": "application/citeproc+json"},
            follow_redirects=True,
            timeout=10.0,
        )
        resp.raise_for_status()
        data = resp.json()
        authors = data.get("author", [])
        author_names = [
            f"{a.get('given', '')} {a.get('family', '')}".strip()
            for a in authors[:5]
        ]
        issued = data.get("issued", {}).get("date-parts", [[]])
        year = str(issued[0][0]) if issued and issued[0] else ""
        return {
            "doi": doi,
            "title": data.get("title", ""),
            "journal": data.get("container-title", ""),
            "authors": author_names,
            "year": year,
            "type": data.get("type", ""),
        }

    def _input_class(self):
        return DOIInput


# -- Citation graph tool ---------------------------------------------------

class CitationInput(ToolInput):
    doi: str = Field(description="DOI of the paper to get citations for")
    direction: str = Field(default="citations", description="'citations' (papers citing this) or 'references' (papers this cites)")


class CitationGraphTool(Tool):
    name = "citation_graph"
    description = "Get citing or referenced papers for a DOI via Semantic Scholar API."

    _API_URL = "https://api.semanticscholar.org/graph/v1/paper"

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = CitationInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        try:
            results = _s2_cb.call(self._fetch, inp.doi, inp.direction)
            return ToolOutput(
                success=True,
                data=results,
                metadata={"count": len(results), "direction": inp.direction},
            )
        except Exception as exc:
            return ToolOutput(success=False, error=str(exc))

    def _fetch(self, doi: str, direction: str) -> list[dict]:
        endpoint = "citations" if direction == "citations" else "references"
        resp = httpx.get(
            f"{self._API_URL}/DOI:{doi}/{endpoint}",
            params={"fields": "title,year,citationCount,externalIds", "limit": 10},
            timeout=10.0,
        )
        resp.raise_for_status()
        data = resp.json().get("data", [])
        results = []
        for item in data:
            paper = item.get("citingPaper" if direction == "citations" else "citedPaper", {})
            if not paper:
                continue
            ext_ids = paper.get("externalIds", {}) or {}
            results.append({
                "title": paper.get("title", ""),
                "year": paper.get("year"),
                "citation_count": paper.get("citationCount", 0),
                "doi": ext_ids.get("DOI", ""),
            })
        return results

    def _input_class(self):
        return CitationInput


def build_literature_tools(retriever) -> list[Tool]:
    """Factory: create all literature domain tools."""
    return [
        ChromaLitSearchTool(retriever),
        PubMedCentralTool(),
        DOIResolverTool(),
        CitationGraphTool(),
    ]
