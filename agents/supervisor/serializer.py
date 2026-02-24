"""Output serialization for the supervisor pipeline.

Extracted from MasterAgent.serialize() (master.py:485-507).
"""

from __future__ import annotations

from typing import Dict, List

from ..api_budget import api_budget_monitor, summarize_budget_status
from ..models import MasterResult, WorkerResult
from .quality_checker import (
    attach_evidence_ids,
    build_validation_summary,
    detect_story_gaps,
    link_story_claims,
)


def serialize_result(result: MasterResult) -> dict:
    """Convert a MasterResult into a fully validated, serializable dict."""
    payload = result.model_dump()
    workers: Dict[str, Dict[str, object]] = payload.get("workers", {})

    worker_evidence_map, evidence_catalog = attach_evidence_ids(workers)
    claim_links = link_story_claims(
        payload.get("innovation_story", ""),
        workers,
        worker_evidence_map,
    )
    story_checks = detect_story_gaps(
        payload.get("innovation_story", ""),
        claim_links,
        evidence_catalog,
    )

    api_budgets = api_budget_monitor.snapshot()

    payload["validation"] = build_validation_summary(claim_links)
    payload["quality"] = _build_quality_summary(
        worker_outputs=result.workers,
        story_checks=story_checks,
        api_budgets=api_budgets,
    )
    payload["workers"] = workers
    return payload


def _build_quality_summary(
    worker_outputs: Dict[str, WorkerResult],
    story_checks: dict | None = None,
    api_budgets: dict | None = None,
) -> dict:
    metrics: Dict[str, dict] = {}
    alerts: Dict[str, List[str]] = {}
    has_anomaly = False
    for name, result in (worker_outputs or {}).items():
        metrics[name] = result.metrics or {}
        if result.alerts:
            alerts[name] = result.alerts
            if any(a.lower().startswith("[anomaly]") for a in result.alerts):
                has_anomaly = True
    status = "pass"
    if alerts:
        status = "investigate" if has_anomaly else "needs_attention"
    if story_checks:
        status = _max_quality_status(status, story_checks.get("status", "pass"))
    if api_budgets:
        budget_status = summarize_budget_status(api_budgets)
        if budget_status:
            status = _max_quality_status(status, budget_status)
    summary = {
        "status": status,
        "metrics": metrics,
        "alerts": alerts,
    }
    if story_checks:
        summary["story_checks"] = story_checks
    if api_budgets:
        summary["api_budgets"] = api_budgets
    return summary


def _max_quality_status(left: str, right: str) -> str:
    order = {"pass": 0, "needs_attention": 1, "investigate": 2}
    return left if order.get(left, 0) >= order.get(right, 0) else right
