"""Quality checking — claim linking, story gap detection, hallucination flags.

Extracted from the original MasterAgent validation helpers (master.py:556-762).
"""

from __future__ import annotations

import re
from typing import Dict, List

from ..logging_config import get_logger
from ..models import WorkerResult

logger = get_logger(__name__)

STORY_SUPPORT_THRESHOLD = 0.2
STORY_STOPWORDS = {
    "the", "and", "with", "from", "that", "this", "have", "has", "will",
    "into", "than", "once", "when", "over", "after", "more", "less",
    "such", "also", "upon", "case", "cases", "data",
}

WORKER_KEYWORDS = {
    "clinical": ("clinical", "trial", "phase", "patient", "fvc", "registry"),
    "literature": ("literature", "mechanism", "paper", "doi", "study"),
    "patent": ("patent", "ip", "claim", "regulatory", "label"),
    "market": ("market", "tam", "competition", "commercial", "access", "pricing"),
}


def attach_evidence_ids(
    workers: Dict[str, Dict[str, object]],
) -> tuple[Dict[str, List[str]], Dict[str, dict]]:
    """Assign stable IDs to each piece of evidence across workers."""
    worker_map: Dict[str, List[str]] = {}
    catalog: Dict[str, dict] = {}
    for worker_name, data in workers.items():
        evidence_records = list(data.get("evidence", []) or [])
        enriched: List[dict] = []
        evidence_ids: List[str] = []
        for idx, record in enumerate(evidence_records, start=1):
            evidence_id = f"{worker_name}-{idx}"
            enriched_record = dict(record)
            enriched_record["evidence_id"] = evidence_id
            enriched.append(enriched_record)
            evidence_ids.append(evidence_id)
            catalog[evidence_id] = enriched_record
        data["evidence"] = enriched
        worker_map[worker_name] = evidence_ids
    return worker_map, catalog


def link_story_claims(
    story: str,
    workers: Dict[str, Dict[str, object]],
    worker_evidence_map: Dict[str, List[str]],
) -> List[dict]:
    """Link each sentence in the story to the most relevant worker's evidence."""
    sentences = [s.strip() for s in re.split(r"(?<=[.!?])\s+", story) if s.strip()]
    if not sentences:
        sentences = [
            str(data.get("summary", "")).strip()
            for data in workers.values()
            if data.get("summary")
        ]
    all_evidence_ids = [
        eid for ids in worker_evidence_map.values() for eid in ids
    ]
    fallback_worker = next(
        (w for w, ids in worker_evidence_map.items() if ids), "clinical"
    )
    claim_links: List[dict] = []
    for idx, sentence in enumerate(sentences, start=1):
        worker = _match_worker(sentence, fallback_worker)
        linked_ids = worker_evidence_map.get(worker) or all_evidence_ids[:1]
        status = "linked" if linked_ids else "missing"
        claim_links.append({
            "claim_id": f"claim-{idx}",
            "claim_text": sentence,
            "worker": worker,
            "evidence_ids": linked_ids,
            "status": status,
        })
    return claim_links


def detect_story_gaps(
    story: str,
    claim_links: List[dict],
    evidence_catalog: Dict[str, dict],
) -> dict:
    """Detect unsupported claims in the innovation story."""
    flagged: List[dict] = []
    total = len(claim_links)
    if not (story or "").strip() or total == 0:
        return {
            "status": "pass",
            "claims_total": total,
            "claims_flagged": 0,
            "citation_gap_ratio": 0.0,
            "flagged_claims": [],
        }
    for claim in claim_links:
        claim_text = claim.get("claim_text", "")
        evidence_ids = claim.get("evidence_ids") or []
        support_score = _claim_support_score(claim_text, evidence_ids, evidence_catalog)
        reasons: List[str] = []
        if not evidence_ids:
            reasons.append("No citations linked to claim")
        if evidence_ids and support_score < STORY_SUPPORT_THRESHOLD:
            reasons.append(f"Low lexical overlap ({support_score:.2f})")
        has_numbers = any(ch.isdigit() for ch in claim_text)
        evidence_texts = [
            (evidence_catalog[eid].get("text") or "")
            for eid in evidence_ids
            if eid in evidence_catalog
        ]
        if has_numbers and evidence_ids and evidence_texts and not any(
            any(ch.isdigit() for ch in t) for t in evidence_texts
        ):
            reasons.append("Numeric claim lacks cited numbers")
        if reasons:
            flagged.append({
                "claim_id": claim.get("claim_id"),
                "claim_text": claim_text,
                "reason": "; ".join(reasons),
                "support_score": round(support_score, 3),
                "evidence_ids": evidence_ids,
            })
    flagged_count = len(flagged)
    ratio = flagged_count / max(1, total)
    if flagged_count == 0:
        status = "pass"
    elif ratio >= 0.4 or flagged_count >= 2:
        status = "investigate"
    else:
        status = "needs_attention"
    return {
        "status": status,
        "claims_total": total,
        "claims_flagged": flagged_count,
        "citation_gap_ratio": round(ratio, 3),
        "flagged_claims": flagged,
    }


def build_validation_summary(claim_links: List[dict]) -> dict:
    total = len(claim_links)
    linked = sum(1 for c in claim_links if c.get("status") == "linked")
    status = "pass" if total and linked == total else "needs_review"
    return {
        "status": status,
        "claims_total": total,
        "claims_linked": linked,
        "claim_links": claim_links,
    }


# -- Internal helpers ------------------------------------------------------

def _match_worker(sentence: str, fallback: str) -> str:
    lowered = sentence.lower()
    for worker, keywords in WORKER_KEYWORDS.items():
        if any(kw in lowered for kw in keywords):
            return worker
    return fallback


def _claim_support_score(
    claim_text: str,
    evidence_ids: List[str],
    evidence_catalog: Dict[str, dict],
) -> float:
    claim_tokens = _tokenize_text(claim_text)
    if not claim_tokens:
        return 0.0
    evidence_tokens: set[str] = set()
    for eid in evidence_ids:
        record = evidence_catalog.get(eid) or {}
        evidence_tokens.update(_tokenize_text(record.get("text", "")))
    if not evidence_tokens:
        return 0.0
    overlap = len(claim_tokens & evidence_tokens)
    return overlap / max(1, len(claim_tokens))


def _tokenize_text(text: str) -> set[str]:
    tokens = re.findall(r"[a-z0-9]+", (text or "").lower())
    return {t for t in tokens if len(t) >= 4 and t not in STORY_STOPWORDS}
