"""Market domain tools for the ReAct agent loop.

Tools:
  - chroma_market_search: Retrieves market evidence from the vector index.
  - market_data_lookup: Aggregates market data (TAM, competition, growth).
  - tam_calculator: Computes Total Addressable Market from inputs.
  - competition_analyzer: Analyzes competitive landscape.
"""

from __future__ import annotations

import json
from typing import Optional

from pydantic import Field

from ..logging_config import get_logger
from .base import Tool, ToolInput, ToolOutput

logger = get_logger(__name__)


# -- Chroma market retrieval -----------------------------------------------

class ChromaMarketInput(ToolInput):
    query: str = Field(description="Search query for market evidence")
    top_k: int = Field(default=5, description="Number of results")


class ChromaMarketSearchTool(Tool):
    name = "chroma_market_search"
    description = "Search the local market evidence vector index for TAM, competition, and pricing data."

    def __init__(self, retriever) -> None:
        self._retriever = retriever

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = ChromaMarketInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        try:
            results = self._retriever.search(query=inp.query, source_type="market", top_k=inp.top_k)
            return ToolOutput(success=True, data=results, metadata={"count": len(results)})
        except Exception as exc:
            return ToolOutput(success=False, error=str(exc))

    def _input_class(self):
        return ChromaMarketInput


# -- Market data lookup tool -----------------------------------------------

class MarketDataInput(ToolInput):
    text: str = Field(description="Raw market text snippet to extract data from")


class MarketDataTool(Tool):
    name = "market_data_lookup"
    description = "Extract structured market data (TAM, incidence, growth, competition) from a text snippet."

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = MarketDataInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        parsed = self._extract(inp.text)
        return ToolOutput(success=True, data=parsed)

    def _extract(self, text: str) -> dict:
        try:
            data = json.loads(text)
            if isinstance(data, dict):
                market = data.get("market", data)
                if isinstance(market, dict):
                    return {
                        "tam": str(market.get("tam", "")),
                        "incidence": str(market.get("incidence", "")),
                        "growth": str(market.get("growth", "")),
                        "competition": str(market.get("competition", "")),
                    }
        except (json.JSONDecodeError, TypeError):
            pass
        return {
            "tam": self._find_currency(text),
            "incidence": self._find_patients(text),
            "growth": self._find_percentage(text),
            "competition": "limited" if "limited" in text.lower() else "",
        }

    def _find_currency(self, text: str) -> str:
        for token in text.split():
            if "$" in token:
                return token
        return ""

    def _find_patients(self, text: str) -> str:
        for token in text.split():
            if "pts" in token.lower() or "patients" in token.lower():
                return token
        return ""

    def _find_percentage(self, text: str) -> str:
        for token in text.split():
            if token.endswith("%"):
                return token
        return ""

    def _input_class(self):
        return MarketDataInput


# -- TAM calculator --------------------------------------------------------

class TAMInput(ToolInput):
    incidence: int = Field(description="Number of patients (annual incidence)")
    price_per_patient: float = Field(description="Annual treatment cost per patient in USD")
    market_share_pct: float = Field(default=100.0, description="Addressable market share percentage")


class TAMCalculatorTool(Tool):
    name = "tam_calculator"
    description = "Calculate Total Addressable Market (TAM) from incidence, pricing, and market share data."

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = TAMInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        tam = inp.incidence * inp.price_per_patient * (inp.market_share_pct / 100.0)
        return ToolOutput(
            success=True,
            data={
                "tam_usd": round(tam, 2),
                "tam_display": self._format_currency(tam),
                "incidence": inp.incidence,
                "price_per_patient": inp.price_per_patient,
                "market_share_pct": inp.market_share_pct,
            },
        )

    def _format_currency(self, value: float) -> str:
        if value >= 1_000_000_000:
            return f"${value / 1_000_000_000:.1f}B"
        if value >= 1_000_000:
            return f"${value / 1_000_000:.1f}M"
        return f"${value:,.0f}"

    def _input_class(self):
        return TAMInput


# -- Competition analyzer --------------------------------------------------

class CompetitionInput(ToolInput):
    market_records: list[dict] = Field(description="List of market data records with competition info")


class CompetitionAnalyzerTool(Tool):
    name = "competition_analyzer"
    description = "Analyze competitive landscape intensity from market data records."

    def execute(self, tool_input: ToolInput) -> ToolOutput:
        inp = CompetitionInput.model_validate(tool_input.model_dump() if isinstance(tool_input, ToolInput) else tool_input)
        competitors = set()
        signals = []
        for record in inp.market_records:
            comp = record.get("competition", "")
            if comp:
                competitors.add(comp)
            notes = record.get("notes", "")
            if "limited" in notes.lower():
                signals.append("limited competition")
            elif "crowded" in notes.lower() or "competitive" in notes.lower():
                signals.append("high competition")
        intensity = "high" if len(competitors) > 3 else "medium" if len(competitors) > 1 else "low"
        return ToolOutput(
            success=True,
            data={
                "intensity": intensity,
                "competitor_count": len(competitors),
                "competitors": sorted(competitors),
                "signals": signals,
            },
        )

    def _input_class(self):
        return CompetitionInput


def build_market_tools(retriever) -> list[Tool]:
    """Factory: create all market domain tools."""
    return [
        ChromaMarketSearchTool(retriever),
        MarketDataTool(),
        TAMCalculatorTool(),
        CompetitionAnalyzerTool(),
    ]
