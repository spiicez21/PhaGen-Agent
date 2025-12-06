from __future__ import annotations

import json
import logging
import re
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeout
from dataclasses import asdict
from typing import Dict, List

from .api_budget import api_budget_monitor, summarize_budget_status
from .llm import LLMClient, LLMClientError, format_structured_context
from .temperature import resolve_master_temperature, resolve_worker_temperatures
from .models import (
    MasterResult,
    MasterRun,
    Recommendation,
    WorkerFailure,
    WorkerRequest,
    WorkerResult,
)
from .retrieval import Retriever
from .synonyms import SynonymExpander
from .workers.clinical import ClinicalWorker
from .workers.patent import PatentWorker
from .workers.literature import LiteratureWorker
from .workers.market import MarketWorker


class MasterAgent:
    DEFAULT_WORKER_TIMEOUTS = {
        "clinical": 45.0,
        "patent": 35.0,
        "literature": 35.0,
        "market": 30.0,
    }
    WORKER_KEYWORDS = {
        "clinical": ("clinical", "trial", "phase", "patient", "fvc", "registry"),
        "literature": ("literature", "mechanism", "paper", "doi", "study"),
        "patent": ("patent", "ip", "claim", "regulatory", "label"),
        "market": ("market", "tam", "competition", "commercial", "access", "pricing"),
    }
    DEFAULT_RETRY_BUDGET = {
        "clinical": 1,
        "patent": 2,
        "literature": 1,
        "market": 1,
    }
    STORY_SUPPORT_THRESHOLD = 0.2
    STORY_STOPWORDS = {
        "the",
        "and",
        "with",
        "from",
        "that",
        "this",
        "have",
        "has",
        "will",
        "into",
        "than",
        "once",
        "when",
        "over",
        "after",
        "more",
        "less",
        "such",
        "also",
        "upon",
        "into",
        "case",
        "cases",
        "data",
    }

    def __init__(
        self,
        top_k: int = 5,
        context_tokens: int = 1200,
        llm_client: LLMClient | None = None,
        synonym_expander: SynonymExpander | None = None,
        worker_timeouts: Dict[str, float] | None = None,
        worker_retry_budget: Dict[str, int] | None = None,
        worker_temperatures: Dict[str, float] | None = None,
        master_temperature: float | None = None,
    ) -> None:
        retriever = Retriever(top_k=top_k, context_tokens=context_tokens)
        self.context_tokens = context_tokens
        self.llm = llm_client or LLMClient()
        self.synonyms = synonym_expander or SynonymExpander()
        self.worker_temperatures = resolve_worker_temperatures(worker_temperatures)
        self.workers = {
            "clinical": ClinicalWorker(
                retriever,
                self.llm,
                temperature=self.worker_temperatures.get("clinical"),
            ),
            "patent": PatentWorker(
                retriever,
                self.llm,
                temperature=self.worker_temperatures.get("patent"),
            ),
            "literature": LiteratureWorker(
                retriever,
                self.llm,
                temperature=self.worker_temperatures.get("literature"),
            ),
            "market": MarketWorker(
                retriever,
                self.llm,
                temperature=self.worker_temperatures.get("market"),
            ),
        }
        self.worker_timeouts = {
            **self.DEFAULT_WORKER_TIMEOUTS,
            **(worker_timeouts or {}),
        }
        self.worker_retry_budget = {
            **self.DEFAULT_RETRY_BUDGET,
            **(worker_retry_budget or {}),
        }
        self.master_temperature = resolve_master_temperature(master_temperature)
        self._logger = logging.getLogger(__name__)

    def run(
        self,
        molecule: str,
        synonyms: list[str] | None = None,
        smiles: str | None = None,
    ) -> MasterRun:
        self._logger.info(f"🚀 Starting agent analysis for molecule: {molecule}")
        self._logger.info(f"   Provided synonyms: {synonyms or 'None'}")
        self._logger.info(f"   SMILES: {smiles or 'Not provided'}")
        
        expanded_synonyms = self.synonyms.expand(
            molecule,
            provided_synonyms=synonyms,
            smiles=smiles,
        )
        self._logger.info(f"   Expanded synonyms: {expanded_synonyms}")
        request = WorkerRequest(
            molecule=molecule,
            synonyms=expanded_synonyms,
            smiles=smiles,
            top_k=len(self.workers),
            context_tokens=self.context_tokens,
        )
        worker_outputs: Dict[str, WorkerResult] = {}
        failures: list[WorkerFailure] = []
        self._logger.info(f"\n📋 Executing {len(self.workers)} workers: {list(self.workers.keys())}")
        
        for name, worker in self.workers.items():
            try:
                self._logger.info(f"   ⚙️  Running {name} worker...")
                worker_outputs[name] = self._run_with_retries(name, worker, request)
                self._logger.info(f"   ✅ {name} worker completed (confidence: {worker_outputs[name].confidence:.2f})")
            except Exception as exc:  # pragma: no cover - placeholder
                self._logger.error(f"   ❌ {name} worker failed: {str(exc)}")
                failures.append(WorkerFailure(worker_name=name, reason=str(exc)))
        if failures:
            return MasterRun(
                success=False,
                failures=failures,
            )
        market_score = int(
            float(worker_outputs["market"].metadata.get("market_score", "70"))
        )
        self._logger.info(f"\n🎯 Market score: {market_score}")
        
        self._logger.info("\n🔄 Generating fallback story and recommendation...")
        fallback_story = self._fallback_story(molecule, worker_outputs)
        fallback_recommendation = self._recommend(market_score, worker_outputs)
        
        self._logger.info("\n🧠 Synthesizing final innovation story and recommendation...")
        innovation_story, recommendation = self._synthesize_story_and_recommendation(
            molecule,
            worker_outputs,
            market_score,
            fallback_story,
            fallback_recommendation,
        )
        self._logger.info(f"   Final recommendation: {recommendation}")
        self._logger.info(f"\n✅ Agent analysis completed successfully for {molecule}")
        self._logger.info(f"   Story length: {len(innovation_story)} chars")
        self._logger.info(f"   Workers completed: {len(worker_outputs)}/{len(self.workers)}\n")
        
        return MasterRun(
            success=True,
            output=MasterResult(
                innovation_story=innovation_story,
                recommendation=recommendation,
                market_score=market_score,
                workers=worker_outputs,
            ),
        )

    def _recommend(
        self, market_score: int, worker_outputs: Dict[str, WorkerResult]
    ) -> Recommendation:
        avg_conf = sum(worker.confidence for worker in worker_outputs.values()) / len(
            worker_outputs
        )
        if market_score >= 75 and avg_conf >= 0.7:
            return Recommendation.go
        if market_score >= 60:
            return Recommendation.investigate
        return Recommendation.no_go

    def _fallback_story(
        self, molecule: str, worker_outputs: Dict[str, WorkerResult]
    ) -> str:
        header = f"{molecule} — Lite innovation summary (LLM fallback)."
        lines = self._build_lite_story_lines(worker_outputs)
        return "\n".join([header, *lines])

    def _build_lite_story_lines(self, worker_outputs: Dict[str, WorkerResult]) -> List[str]:
        ordered = ("clinical", "literature", "patent", "market")
        lines: List[str] = []
        for name in ordered:
            result = worker_outputs.get(name)
            if not result:
                continue
            highlight = self._summarize_worker_highlight(name, result)
            lines.append(
                f"- {name.title()} ({result.confidence_band} confidence): {highlight}"
            )
        for name, result in worker_outputs.items():
            if name in ordered:
                continue
            highlight = self._summarize_worker_highlight(name, result)
            lines.append(f"- {name.title()}: {highlight}")
        if not lines:
            lines.append("- No worker insights available; re-run once retrieval succeeds.")
        return lines

    def _summarize_worker_highlight(
        self,
        worker_name: str,
        result: WorkerResult,
    ) -> str:
        summary = (result.summary or "").strip()
        metadata_hint = self._metadata_hint(worker_name, result.metadata)
        parts = [part for part in (summary, metadata_hint) if part]
        if not parts:
            parts = [self._evidence_hint(result)]
        text = " | ".join(parts)
        return self._ensure_sentence(self._truncate(text, 260))

    def _metadata_hint(self, worker_name: str, metadata: Dict[str, str]) -> str:
        if not metadata:
            return ""
        hints: List[str] = []
        if worker_name == "clinical":
            trials_hint = self._describe_trial_count(metadata.get("trials"))
            if trials_hint:
                hints.append(trials_hint)
            populations = metadata.get("populations")
            if populations:
                hints.append(f"populations: {self._truncate(populations, 60)}")
        elif worker_name == "patent":
            assignees = metadata.get("assignees")
            if assignees:
                hints.append(f"assignees: {self._truncate(assignees, 60)}")
            priority = metadata.get("priority_dates")
            if priority:
                hints.append(f"priority: {self._truncate(priority, 40)}")
            contra = self._first_json_entry(metadata.get("contraindications"))
            if contra:
                hints.append(f"contra: {self._truncate(contra, 60)}")
        elif worker_name == "literature":
            count = metadata.get("evidence_count")
            if count:
                hints.append(f"{count} studies reviewed")
            journals = metadata.get("journals")
            if journals:
                hints.append(f"journals: {self._truncate(journals, 60)}")
        elif worker_name == "market":
            score = metadata.get("market_score")
            if score:
                hints.append(f"score {score}")
            tam = metadata.get("tam")
            if tam:
                hints.append(f"TAM {self._truncate(tam, 30)}")
            comp = metadata.get("competition")
            if comp:
                hints.append(f"competition: {self._truncate(comp, 40)}")

        if not hints:
            for key, value in metadata.items():
                if value:
                    hints.append(f"{key}: {self._truncate(str(value), 60)}")
                    break
        return "; ".join(hints)

    def _describe_trial_count(self, raw: str | None) -> str:
        if not raw:
            return ""
        try:
            data = json.loads(raw)
        except (json.JSONDecodeError, TypeError):
            return ""
        if isinstance(data, list):
            count = len(data)
            if count:
                return f"{count} trial{'s' if count != 1 else ''} parsed"
        return ""

    def _first_json_entry(self, raw: str | None) -> str:
        if not raw:
            return ""
        try:
            data = json.loads(raw)
        except (json.JSONDecodeError, TypeError):
            return ""
        if isinstance(data, list) and data:
            return str(data[0])
        if isinstance(data, dict):
            for key, value in data.items():
                if value:
                    return f"{key}: {value}"
        return str(data)

    def _evidence_hint(self, result: WorkerResult) -> str:
        if result.evidence:
            snippet = (result.evidence[0].text or "").strip()
            if snippet:
                return self._truncate(snippet, 200)
        return "No structured insight captured yet"

    def _truncate(self, value: str, limit: int = 160) -> str:
        if len(value) <= limit:
            return value
        return value[: limit - 3].rstrip() + "..."

    def _ensure_sentence(self, text: str) -> str:
        cleaned = text.strip()
        if not cleaned:
            return "No structured insight captured yet."
        if cleaned[-1] in ".!?":
            return cleaned
        return f"{cleaned}."

    def _synthesize_story_and_recommendation(
        self,
        molecule: str,
        worker_outputs: Dict[str, WorkerResult],
        market_score: int,
        fallback_story: str,
        fallback_recommendation: Recommendation,
    ) -> tuple[str, Recommendation]:
        if not self.llm:
            return fallback_story, fallback_recommendation
        payload = {
            "molecule": molecule,
            "market_score": market_score,
            "workers": self._build_worker_payload(worker_outputs),
        }
        rubric = (
            "Go: compelling clinical efficacy, manageable risk, clear TAM; "
            "Investigate: mixed or emerging signals requiring validation; "
            "No-Go: weak efficacy, blocking IP/regulatory risk, or limited market."  # noqa: E501
        )
        # Use simpler prompt optimized for small models
        user_prompt = (
            f"Molecule: {molecule}\n"
            f"Market score: {market_score}/10\n"
            f"Clinical: {worker_outputs.get('clinical', WorkerResult()).summary[:150]}\n"
            f"Patent: {worker_outputs.get('patent', WorkerResult()).summary[:150]}\n\n"
            "Output valid JSON only (no markdown):\n"
            "{\"innovation_story\": \"2-3 sentence story\", \"recommendation\": \"Go\", \"rationale\": \"brief reason\"}\n\n"
            f"Recommendation must be: Go, Investigate, or No-Go based on: {rubric}"
        )
        try:
            raw = self.llm.generate(
                prompt=user_prompt,
                system_prompt="Return only valid JSON. No explanations or markdown.",
                temperature=0.1,  # Lower temp for more reliable JSON
                max_tokens=300,
            )
            parsed = self._parse_master_response(raw)
        except (LLMClientError, ValueError) as exc:  # pragma: no cover - runtime
            self._logger.debug("Master synthesis JSON parse failed: %s, using fallback", exc)
            return fallback_story, fallback_recommendation

        story = parsed.get("innovation_story") or fallback_story
        recommendation = self._resolve_recommendation(
            parsed.get("recommendation"),
            fallback_recommendation,
        )
        return story, recommendation

    def _build_worker_payload(
        self, worker_outputs: Dict[str, WorkerResult]
    ) -> Dict[str, dict]:
        payload: Dict[str, dict] = {}
        for name, result in worker_outputs.items():
            payload[name] = {
                "summary": result.summary,
                "confidence": round(result.confidence, 3),
                "confidence_band": result.confidence_band,
                "metadata": result.metadata,
                "evidence": [
                    {
                        "type": item.type,
                        "confidence": round(item.confidence, 3),
                        "text": (item.text or "")[:280],
                        "url": item.url,
                    }
                    for item in result.evidence[:3]
                ],
                "metrics": result.metrics,
            }
            if result.alerts:
                payload[name]["alerts"] = result.alerts
        return payload

    def _parse_master_response(self, raw: str) -> dict:
        snippet = raw.strip()
        # Try direct parse first
        try:
            return json.loads(snippet)
        except json.JSONDecodeError:
            pass
        # Try extracting JSON from markdown code blocks
        if "```json" in snippet:
            start = snippet.find("```json") + 7
            end = snippet.find("```", start)
            if end > start:
                snippet = snippet[start:end].strip()
                try:
                    return json.loads(snippet)
                except json.JSONDecodeError:
                    pass
        # Try finding first { to last }
        start = snippet.find("{")
        end = snippet.rfind("}")
        if start != -1 and end != -1 and end > start:
            json_blob = snippet[start : end + 1]
            try:
                return json.loads(json_blob)
            except json.JSONDecodeError as e:
                self._logger.debug(f"JSON parse failed: {e}, raw response: {raw[:200]}")
        raise ValueError("Master synthesis did not return valid JSON")

    def _resolve_recommendation(
        self,
        label: str | None,
        fallback: Recommendation,
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

    def serialize(self, result: MasterResult) -> dict:
        payload = asdict(result)
        workers: Dict[str, Dict[str, object]] = payload.get("workers", {})
        worker_evidence_map, evidence_catalog = self._attach_evidence_ids(workers)
        claim_links = self._link_story_claims(
            payload.get("innovation_story", ""),
            workers,
            worker_evidence_map,
        )
        story_checks = self._detect_story_gaps(
            payload.get("innovation_story", ""),
            claim_links,
            evidence_catalog,
        )
        api_budgets = api_budget_monitor.snapshot()
        payload["validation"] = self._build_validation_summary(claim_links)
        payload["quality"] = self._build_quality_summary(
            worker_outputs=result.workers,
            story_checks=story_checks,
            api_budgets=api_budgets,
        )
        payload["workers"] = workers
        return payload

    # Execution helpers ------------------------------------------------
    def _run_with_retries(
        self,
        name: str,
        worker,
        request: WorkerRequest,
    ) -> WorkerResult:
        retry_budget = max(0, self.worker_retry_budget.get(name, 0))
        attempts = 0
        last_error: Exception | None = None
        while attempts <= retry_budget:
            attempts += 1
            try:
                timeout = self.worker_timeouts.get(name, 30.0)
                return self._run_with_timeout(worker, request, timeout)
            except Exception as exc:  # pragma: no cover - runtime dependent
                last_error = exc
                self._logger.warning(
                    "%s attempt %s/%s failed: %s",
                    name,
                    attempts,
                    retry_budget + 1,
                    exc,
                )
        raise RuntimeError(
            f"{name} failed after {attempts} attempts: {last_error}",
        )

    def _run_with_timeout(
        self,
        worker,
        request: WorkerRequest,
        timeout_seconds: float,
    ) -> WorkerResult:
        executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix=f"{worker.name}-worker")
        future = executor.submit(worker.run, request)
        try:
            return future.result(timeout=timeout_seconds)
        except FuturesTimeout as exc:
            future.cancel()
            raise TimeoutError(
                f"{worker.name} exceeded {timeout_seconds:.1f}s timeout",
            ) from exc
        finally:
            executor.shutdown(wait=False, cancel_futures=True)

    # Validation helpers ----------------------------------------------
    def _attach_evidence_ids(
        self, workers: Dict[str, Dict[str, object]]
    ) -> tuple[Dict[str, List[str]], Dict[str, dict]]:
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

    def _link_story_claims(
        self,
        story: str,
        workers: Dict[str, Dict[str, object]],
        worker_evidence_map: Dict[str, List[str]],
    ) -> List[dict]:
        sentences = [
            sentence.strip()
            for sentence in re.split(r"(?<=[.!?])\s+", story)
            if sentence.strip()
        ]
        if not sentences:
            sentences = [
                str(data.get("summary", "")).strip()
                for data in workers.values()
                if data.get("summary")
            ]
        all_evidence_ids = [
            evidence_id
            for evidence_list in worker_evidence_map.values()
            for evidence_id in evidence_list
        ]
        fallback_worker = next(
            (worker for worker, ids in worker_evidence_map.items() if ids),
            "clinical",
        )
        claim_links: List[dict] = []
        for idx, sentence in enumerate(sentences, start=1):
            worker = self._match_worker(sentence, fallback_worker)
            linked_ids = worker_evidence_map.get(worker) or all_evidence_ids[:1]
            status = "linked" if linked_ids else "missing"
            claim_links.append(
                {
                    "claim_id": f"claim-{idx}",
                    "claim_text": sentence,
                    "worker": worker,
                    "evidence_ids": linked_ids,
                    "status": status,
                }
            )
        return claim_links

    def _match_worker(self, sentence: str, fallback_worker: str) -> str:
        lowered = sentence.lower()
        for worker, keywords in self.WORKER_KEYWORDS.items():
            if any(keyword in lowered for keyword in keywords):
                return worker
        return fallback_worker

    def _build_validation_summary(self, claim_links: List[dict]) -> dict:
        total = len(claim_links)
        linked = sum(1 for claim in claim_links if claim.get("status") == "linked")
        status = "pass" if total and linked == total else "needs_review"
        return {
            "status": status,
            "claims_total": total,
            "claims_linked": linked,
            "claim_links": claim_links,
        }

    def _build_quality_summary(
        self,
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
                if any(alert.lower().startswith("[anomaly]") for alert in result.alerts):
                    has_anomaly = True
        status = "pass"
        if alerts:
            status = "investigate" if has_anomaly else "needs_attention"
        if story_checks:
            status = self._max_quality_status(status, story_checks.get("status", "pass"))
        if api_budgets:
            budget_status = summarize_budget_status(api_budgets)
            if budget_status:
                status = self._max_quality_status(status, budget_status)
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

    def _detect_story_gaps(
        self,
        story: str,
        claim_links: List[dict],
        evidence_catalog: Dict[str, dict],
    ) -> dict:
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
            support_score = self._claim_support_score(claim_text, evidence_ids, evidence_catalog)
            reasons: List[str] = []
            if not evidence_ids:
                reasons.append("No citations linked to claim")
            if evidence_ids and support_score < self.STORY_SUPPORT_THRESHOLD:
                reasons.append(f"Low lexical overlap ({support_score:.2f})")
            has_numbers = any(ch.isdigit() for ch in claim_text)
            evidence_texts = [
                (evidence_catalog[e_id].get("text") or "")
                for e_id in evidence_ids
                if e_id in evidence_catalog
            ]
            if has_numbers and evidence_ids and evidence_texts and not any(
                any(ch.isdigit() for ch in text) for text in evidence_texts
            ):
                reasons.append("Numeric claim lacks cited numbers")
            if reasons:
                flagged.append(
                    {
                        "claim_id": claim.get("claim_id"),
                        "claim_text": claim_text,
                        "reason": "; ".join(reasons),
                        "support_score": round(support_score, 3),
                        "evidence_ids": evidence_ids,
                    }
                )
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

    def _claim_support_score(
        self,
        claim_text: str,
        evidence_ids: List[str],
        evidence_catalog: Dict[str, dict],
    ) -> float:
        claim_tokens = self._tokenize_text(claim_text)
        if not claim_tokens:
            return 0.0
        evidence_tokens: set[str] = set()
        for evidence_id in evidence_ids:
            record = evidence_catalog.get(evidence_id) or {}
            evidence_tokens.update(self._tokenize_text(record.get("text", "")))
        if not evidence_tokens:
            return 0.0
        overlap = len(claim_tokens & evidence_tokens)
        return overlap / max(1, len(claim_tokens))

    def _tokenize_text(self, text: str) -> set[str]:
        tokens = re.findall(r"[a-z0-9]+", (text or "").lower())
        return {
            token
            for token in tokens
            if len(token) >= 4 and token not in self.STORY_STOPWORDS
        }

    def _max_quality_status(self, left: str, right: str) -> str:
        order = {"pass": 0, "needs_attention": 1, "investigate": 2}
        left_rank = order.get(left, 0)
        right_rank = order.get(right, 0)
        return left if left_rank >= right_rank else right
