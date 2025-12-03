from __future__ import annotations

import json
import re
from typing import Dict, List, Optional

from ..models import WorkerRequest, WorkerResult
from .base import Worker


class ClinicalWorker(Worker):
    def __init__(self, retriever, llm=None):
        super().__init__("clinical", retriever, llm)

    def build_summary(self, request: WorkerRequest, passages: List[dict]) -> WorkerResult:
        trials = self._extract_trials(passages)
        fallback_summary = self._summarize_trials(request.molecule, trials)
        summary = self._summarize_with_llm(
            request,
            instructions=
            "Highlight trial phases, statuses, endpoints, and population insights relevant to clinical readiness.",
            passages=passages,
            extra_context={"trials": trials} if trials else None,
            fallback=fallback_summary,
        )
        return WorkerResult(
            summary=summary,
            evidence=self._to_evidence(passages, "clinical"),
            confidence=self._score_confidence(trials),
            metadata=self._build_metadata(trials),
        )

    # Parsing helpers -------------------------------------------------
    def _extract_trials(self, passages: List[dict]) -> List[Dict[str, str]]:
        trials: Dict[str, Dict[str, str]] = {}
        for passage in passages:
            trial = self._parse_passage(passage)
            if not trial:
                continue
            nct_id = trial.get("nct_id") or f"trial-{len(trials) + 1}"
            merged = trials.get(nct_id, {"nct_id": nct_id})
            for key, value in trial.items():
                if value:
                    merged[key] = value
            trials[nct_id] = merged
        return list(trials.values())

    def _parse_passage(self, passage: dict) -> Optional[Dict[str, str]]:
        snippet = passage.get("snippet", "").strip()
        if not snippet:
            return None
        data = self._load_json(snippet)
        if data:
            study = data.get("study") or data
            trial = self._trial_from_study(study)
        else:
            trial = self._trial_from_text(snippet)
        if trial and passage.get("url"):
            trial.setdefault("url", passage["url"])
        return trial

    def _load_json(self, raw: str) -> Optional[dict]:
        try:
            return json.loads(raw)
        except (json.JSONDecodeError, TypeError):
            return None

    def _trial_from_study(self, study: dict) -> Optional[Dict[str, str]]:
        if not isinstance(study, dict):
            return None
        return {
            "nct_id": study.get("nct_id") or self._find_nct(str(study.get("title", ""))),
            "title": study.get("title", ""),
            "phase": study.get("phase", ""),
            "status": study.get("status", ""),
            "population": study.get("population") or study.get("indication", ""),
            "primary_endpoint": study.get("primary_outcome", ""),
            "enrollment": str(study.get("enrollment", "")) if study.get("enrollment") else "",
        }

    def _trial_from_text(self, text: str) -> Dict[str, str]:
        return {
            "nct_id": self._find_nct(text),
            "phase": self._find_phase(text),
            "status": self._find_status(text),
            "primary_endpoint": self._find_endpoint(text),
            "population": self._find_population(text),
        }

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
        statuses = [
            "Recruiting",
            "Active, not recruiting",
            "Completed",
            "Terminated",
            "Withdrawn",
            "Suspended",
            "Not yet recruiting",
            "Enrolling by invitation",
        ]
        for status in statuses:
            if re.search(re.escape(status), text, re.IGNORECASE):
                return status
        return ""

    def _find_endpoint(self, text: str) -> str:
        match = re.search(r"primary\s+(?:endpoint|outcome)[:\-]?\s*(.+?)(?:\. |;|,|$)", text, re.IGNORECASE)
        return match.group(1).strip() if match else ""

    def _find_population(self, text: str) -> str:
        match = re.search(r"in\s+([A-Za-z0-9 \-/]+?)\s+patients", text, re.IGNORECASE)
        return match.group(1).strip() if match else ""

    # Reporting helpers -----------------------------------------------
    def _summarize_trials(self, molecule: str, trials: List[Dict[str, str]]) -> str:
        if not trials:
            return f"No structured trials parsed for {molecule}; showing retrieved evidence only."
        advanced = max(trials, key=lambda trial: self._phase_rank(trial.get("phase", "")))
        completed = sum(1 for trial in trials if trial.get("status", "").lower() == "completed")
        endpoint = advanced.get("primary_endpoint") or "primary endpoint not captured"
        return (
            f"Parsed {len(trials)} trials for {molecule}; most advanced signal is {advanced.get('phase') or 'Phase ?'} "
            f"({advanced.get('nct_id') or 'NCT?'}) targeting {endpoint}. {completed} marked completed."
        )

    def _phase_rank(self, phase: str) -> int:
        order = {
            "phase 0": 0,
            "phase 1": 1,
            "phase 1b": 1,
            "phase 2": 2,
            "phase 2b": 3,
            "phase 3": 4,
            "phase 4": 5,
        }
        return order.get((phase or "").lower(), 0)

    def _score_confidence(self, trials: List[Dict[str, str]]) -> float:
        if not trials:
            return 0.5
        fields = ["phase", "status", "primary_endpoint", "population"]
        filled = sum(1 for trial in trials for field in fields if trial.get(field))
        coverage = filled / (len(trials) * len(fields)) if trials else 0
        return min(0.9, 0.55 + 0.35 * coverage)

    def _build_metadata(self, trials: List[Dict[str, str]]) -> Dict[str, str]:
        metadata: Dict[str, str] = {}
        if trials:
            metadata["trials"] = json.dumps(trials, ensure_ascii=False)
            populations = sorted({trial.get("population", "") for trial in trials if trial.get("population")})
            if populations:
                metadata["populations"] = ", ".join(populations)
        return metadata
