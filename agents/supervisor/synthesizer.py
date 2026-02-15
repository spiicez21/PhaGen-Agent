"""Story synthesis and recommendation logic.

Extracted from the original MasterAgent (master.py:344-403).
"""

from __future__ import annotations

from typing import Dict, List

from ..exceptions import LLMError, LLMParseError, SupervisorError
from ..llm import LLMClient
from ..logging_config import get_logger
from ..models import Recommendation, WorkerResult
from ..react.parser import parse_structured_response

logger = get_logger(__name__)


def recommend(
    market_score: int, worker_outputs: Dict[str, WorkerResult]
) -> Recommendation:
    """Deterministic recommendation from market score + average confidence."""
    avg_conf = sum(w.confidence for w in worker_outputs.values()) / max(1, len(worker_outputs))
    if market_score >= 75 and avg_conf >= 0.7:
        return Recommendation.go
    if market_score >= 60:
        return Recommendation.investigate
    return Recommendation.no_go


def fallback_story(
    molecule: str, worker_outputs: Dict[str, WorkerResult]
) -> str:
    """Build a lite summary when LLM synthesis is unavailable."""
    header = f"{molecule} — Lite innovation summary (LLM fallback)."
    ordered = ("clinical", "literature", "patent", "market")
    lines: List[str] = []
    for name in ordered:
        result = worker_outputs.get(name)
        if not result:
            continue
        summary = (result.summary or "").strip()
        snippet = summary[:260] + "..." if len(summary) > 260 else summary
        if not snippet.endswith((".", "!", "?")):
            snippet += "."
        lines.append(f"- {name.title()} ({result.confidence_band} confidence): {snippet}")
    for name, result in worker_outputs.items():
        if name in ordered:
            continue
        summary = (result.summary or "").strip()[:260]
        lines.append(f"- {name.title()}: {summary}")
    if not lines:
        lines.append("- No worker insights available; re-run once retrieval succeeds.")
    return "\n".join([header, *lines])


def synthesize_story_and_recommendation(
    llm: LLMClient,
    molecule: str,
    worker_outputs: Dict[str, WorkerResult],
    market_score: int,
    fallback_story_text: str,
    fallback_recommendation: Recommendation,
    temperature: float = 0.1,
) -> tuple[str, Recommendation]:
    """Use the LLM to generate the final innovation story and recommendation."""

    clinical_summary = worker_outputs["clinical"].summary if "clinical" in worker_outputs else "N/A"
    patent_summary = worker_outputs["patent"].summary if "patent" in worker_outputs else "N/A"
    literature_summary = worker_outputs["literature"].summary if "literature" in worker_outputs else "N/A"

    user_prompt = (
        f"Molecule: {molecule}\n"
        f"Market score: {market_score}/10\n"
        f"Clinical: {clinical_summary}\n"
        f"Patent: {patent_summary}\n"
        f"Literature: {literature_summary}\n\n"
        "TASK: Generate a structured mechanistic inference analysis.\n"
        "REQUIREMENTS:\n"
        "- Synthesize a story based on target prediction, pathway ranking, and MoA hypotheses.\n"
        "- Do NOT return fallback text or simple summaries.\n"
        "- Recommendation must be one of: Go, Investigate, No-Go.\n"
        "- Use the exact format below.\n\n"
        "FORMAT:\n"
        "STORY: <detailed mechanistic story or failure analysis>\n"
        "RECOMMENDATION: <Go|Investigate|No-Go>\n"
        "RATIONALE: <brief reason>"
    )

    try:
        raw = llm.generate(
            prompt=user_prompt,
            system_prompt="You are a chemical expert. Follow the output format strictly.",
            temperature=temperature,
            max_tokens=600,
        )
        parsed = parse_structured_response(raw, ["STORY", "RECOMMENDATION", "RATIONALE"])
    except LLMError as exc:
        logger.error("synthesis_llm_failed", extra={"error": str(exc)})
        raise SupervisorError(f"Synthesis LLM call failed: {exc}") from exc

    story = parsed.get("story")
    if not story:
        raise LLMParseError(
            "Generated output missing 'STORY' section",
            raw=raw[:500],
        )

    rec_label = parsed.get("recommendation")
    recommendation = _resolve_recommendation(rec_label, fallback_recommendation)
    return story, recommendation


def _resolve_recommendation(
    label: str | None, fallback: Recommendation
) -> Recommendation:
    if not label:
        return fallback
    normalized = label.strip().lower()
    mapping = {
        "go": Recommendation.go,
        "investigate": Recommendation.investigate,
        "no-go": Recommendation.no_go,
        "no go": Recommendation.no_go,
    }
    return mapping.get(normalized, fallback)
