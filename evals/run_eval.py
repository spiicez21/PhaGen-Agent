#!/usr/bin/env python
"""PhaGen eval harness.

Runs deterministic checks against fixture payloads or the live API to ensure
recommendations, guardrails, and worker coverage stay within expected bounds.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Sequence

DEFAULT_API_BASE = "http://localhost:8000/api"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run PhaGen evaluation suite")
    parser.add_argument(
        "--cases-file",
        type=Path,
        default=Path("evals/cases.json"),
        help="Path to cases.json (defaults to evals/cases.json)",
    )
    parser.add_argument(
        "--mode",
        choices=("fixtures", "live"),
        default="fixtures",
        help="fixtures = validate stored payloads, live = hit running API",
    )
    parser.add_argument(
        "--api-base",
        default=DEFAULT_API_BASE,
        help="API base URL when mode=live (default: %(default)s)",
    )
    parser.add_argument(
        "--case",
        action="append",
        dest="case_ids",
        help="Run a specific case id (can be repeated)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=180.0,
        help="Maximum seconds to wait for a live job to finish",
    )
    parser.add_argument(
        "--poll-interval",
        type=float,
        default=3.0,
        help="Seconds between polling /api/jobs/{id} in live mode",
    )
    return parser.parse_args()


def load_cases(cases_path: Path) -> List[dict]:
    with cases_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, list):  # pragma: no cover - sanity guard
        raise ValueError("cases.json must contain a list")
    return data


def load_fixture(base_dir: Path, relative_path: str) -> dict:
    fixture_path = (base_dir / relative_path).resolve()
    with fixture_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if "payload" in data:
        return data["payload"]
    return data


def run_live_case(
    api_base: str,
    case: dict,
    timeout: float,
    poll_interval: float,
) -> dict:
    request_body: Dict[str, Any] = {"molecule": case.get("molecule")}
    if case.get("synonyms"):
        request_body["synonyms"] = case["synonyms"]
    if case.get("smiles"):
        request_body["smiles"] = case["smiles"]

    job = _http_json(f"{api_base.rstrip('/')}/jobs", method="POST", payload=request_body)
    job_id = job.get("job_id")
    if not job_id:
        raise RuntimeError(f"API did not return a job_id for {case['id']}")

    deadline = time.time() + timeout
    while time.time() <= deadline:
        response = _http_json(f"{api_base.rstrip('/')}/jobs/{job_id}", method="GET")
        status = response.get("status")
        if status == "FAILED":
            raise RuntimeError(f"Job {job_id} failed: {response.get('payload')}")
        if status == "COMPLETED" and response.get("payload"):
            return response["payload"]
        time.sleep(poll_interval)
    raise TimeoutError(f"Job {job_id} did not complete within {timeout}s")


def _http_json(url: str, method: str = "GET", payload: dict | None = None) -> dict:
    data: bytes | None = None
    headers = {"Accept": "application/json"}
    if payload is not None:
        headers["Content-Type"] = "application/json"
        data = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(url, data=data, headers=headers, method=method.upper())
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            body = response.read().decode("utf-8")
            return json.loads(body)
    except urllib.error.HTTPError as exc:  # pragma: no cover - network guard
        detail = exc.read().decode("utf-8", errors="ignore")
        raise RuntimeError(f"HTTP {exc.code} for {url}: {detail}") from exc


def evaluate_case(case: dict, payload: dict) -> List[str]:
    expectations: dict = case.get("expectations", {})
    failures: List[str] = []

    expected_rec = expectations.get("recommendation")
    if expected_rec and payload.get("recommendation") != expected_rec:
        failures.append(
            f"Recommendation mismatch (expected {expected_rec}, got {payload.get('recommendation')})"
        )

    score_rules = expectations.get("market_score")
    if score_rules:
        score = payload.get("market_score") or 0
        minimum = score_rules.get("min")
        maximum = score_rules.get("max")
        if minimum is not None and score < minimum:
            failures.append(f"Market score {score} below minimum {minimum}")
        if maximum is not None and score > maximum:
            failures.append(f"Market score {score} above maximum {maximum}")

    quality_expected = expectations.get("quality_status")
    quality_status = ((payload.get("quality") or {}).get("status") or "").lower()
    if quality_expected:
        normalized_targets = {value.lower() for value in quality_expected}
        if quality_status not in normalized_targets:
            failures.append(
                f"Quality status {quality_status or 'n/a'} not in {sorted(normalized_targets)}"
            )

    workers: dict = payload.get("workers") or {}
    required_workers: Sequence[str] = expectations.get("min_workers") or []
    for worker in required_workers:
        if worker not in workers:
            failures.append(f"Missing worker output for {worker}")

    evidence_expectations: dict = expectations.get("min_evidence") or {}
    for worker, min_evidence in evidence_expectations.items():
        evidence = ((workers.get(worker) or {}).get("evidence") or [])
        if len(evidence) < int(min_evidence):
            failures.append(
                f"Worker {worker} only has {len(evidence)} evidence items (expected >= {min_evidence})"
            )

    story_min_words = expectations.get("story_min_words")
    if story_min_words:
        story_words = len((payload.get("innovation_story") or "").split())
        if story_words < story_min_words:
            failures.append(
                f"Innovation story only {story_words} words (expected >= {story_min_words})"
            )

    return failures


def main() -> None:
    args = parse_args()
    cases = load_cases(args.cases_file)
    selected = set(args.case_ids or [])
    base_dir = args.cases_file.parent

    errors = False
    for case in cases:
        if selected and case.get("id") not in selected:
            continue
        print(f"→ Evaluating {case.get('id')} ({case.get('molecule')}) ...", end=" ")
        try:
            if args.mode == "fixtures":
                payload = load_fixture(base_dir, case["fixture"])
            else:
                payload = run_live_case(
                    args.api_base,
                    case,
                    timeout=args.timeout,
                    poll_interval=args.poll_interval,
                )
            failures = evaluate_case(case, payload)
            if failures:
                errors = True
                print("❌")
                for failure in failures:
                    print(f"   - {failure}")
            else:
                print("✅")
        except Exception as exc:  # pragma: no cover - CLI guard
            errors = True
            print("❌")
            print(f"   - Error: {exc}")

    if errors:
        sys.exit(1)
    print("All evaluations passed.")


if __name__ == "__main__":
    main()
