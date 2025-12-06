from __future__ import annotations

import time
from collections import deque
from dataclasses import dataclass
from datetime import datetime, timezone
from threading import Lock
from typing import Deque, Dict, Tuple


@dataclass(frozen=True)
class ApiBudgetConfig:
    label: str
    per_minute: int
    per_day: int
    warn_ratio: float = 0.8


class ApiBudgetMonitor:
    """Tracks rolling API usage and remaining budget for external providers."""

    def __init__(self, configs: Dict[str, ApiBudgetConfig]):
        self._configs = configs
        self._lock = Lock()
        self._history: Dict[str, Dict[str, Deque[Tuple[float, int]]]] = {
            name: {
                "minute": deque(),
                "day": deque(),
            }
            for name in configs
        }
        self._totals: Dict[str, int] = {name: 0 for name in configs}
        self._last_call: Dict[str, float | None] = {name: None for name in configs}

    def record(self, api_name: str, cost: int = 1) -> None:
        if api_name not in self._configs or cost <= 0:
            return
        now = time.time()
        with self._lock:
            trackers = self._history.setdefault(api_name, {"minute": deque(), "day": deque()})
            self._purge(trackers["minute"], window_seconds=60, now=now)
            self._purge(trackers["day"], window_seconds=86_400, now=now)
            trackers["minute"].append((now, cost))
            trackers["day"].append((now, cost))
            self._totals[api_name] = self._totals.get(api_name, 0) + cost
            self._last_call[api_name] = now

    def snapshot(self) -> Dict[str, dict]:
        now = time.time()
        summary: Dict[str, dict] = {}
        with self._lock:
            for name, cfg in self._configs.items():
                trackers = self._history.setdefault(name, {"minute": deque(), "day": deque()})
                self._purge(trackers["minute"], 60, now)
                self._purge(trackers["day"], 86_400, now)
                minute_used = self._count(trackers["minute"])
                day_used = self._count(trackers["day"])
                status = self._resolve_status(cfg, minute_used, day_used)
                summary[name] = {
                    "label": cfg.label,
                    "per_minute_limit": cfg.per_minute,
                    "per_day_limit": cfg.per_day,
                    "minute_used": minute_used,
                    "day_used": day_used,
                    "minute_remaining": max(cfg.per_minute - minute_used, 0),
                    "day_remaining": max(cfg.per_day - day_used, 0),
                    "status": status,
                    "last_call_at": self._iso_timestamp(self._last_call.get(name)),
                    "reset_in_seconds": {
                        "minute": self._seconds_until_reset(trackers["minute"], 60, now),
                        "day": self._seconds_until_reset(trackers["day"], 86_400, now),
                    },
                    "lifetime_calls": self._totals.get(name, 0),
                }
        return summary

    def _resolve_status(
        self,
        cfg: ApiBudgetConfig,
        minute_used: int,
        day_used: int,
    ) -> str:
        if minute_used >= cfg.per_minute or day_used >= cfg.per_day:
            return "exceeded"
        warn_minute = cfg.per_minute * cfg.warn_ratio
        warn_day = cfg.per_day * cfg.warn_ratio
        if minute_used >= warn_minute or day_used >= warn_day:
            return "warning"
        return "ok"

    @staticmethod
    def _count(records: Deque[Tuple[float, int]]) -> int:
        return sum(entry[1] for entry in records)

    @staticmethod
    def _seconds_until_reset(
        records: Deque[Tuple[float, int]],
        window_seconds: int,
        now: float,
    ) -> int:
        if not records:
            return 0
        oldest_timestamp = records[0][0]
        remaining = int(max(window_seconds - (now - oldest_timestamp), 0))
        return remaining

    @staticmethod
    def _iso_timestamp(epoch: float | None) -> str | None:
        if epoch is None:
            return None
        return datetime.fromtimestamp(epoch, tz=timezone.utc).isoformat()

    @staticmethod
    def _purge(records: Deque[Tuple[float, int]], window_seconds: int, now: float) -> None:
        while records and (now - records[0][0]) >= window_seconds:
            records.popleft()


DEFAULT_API_BUDGETS: Dict[str, ApiBudgetConfig] = {
    "ncbi": ApiBudgetConfig(
        label="NCBI E-utilities",
        per_minute=180,  # 3 requests/sec guideline
        per_day=10_000,
        warn_ratio=0.75,
    ),
    "openfda": ApiBudgetConfig(
        label="OpenFDA",
        per_minute=240,  # documented public limit
        per_day=5_000,
        warn_ratio=0.75,
    ),
}

API_PROVIDER_FOR_WORKER = {
    "clinical": "ncbi",
    "literature": "ncbi",
    "regulatory": "openfda",
}

api_budget_monitor = ApiBudgetMonitor(DEFAULT_API_BUDGETS)


def summarize_budget_status(snapshot: Dict[str, dict]) -> str | None:
    worst: str | None = None
    for data in snapshot.values():
        status = data.get("status")
        if status == "exceeded":
            return "investigate"
        if status == "warning":
            worst = "needs_attention"
    return worst
