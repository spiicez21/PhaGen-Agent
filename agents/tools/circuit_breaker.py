"""Three-state circuit breaker for external API calls.

States: CLOSED (normal) -> OPEN (failing) -> HALF_OPEN (testing recovery).
"""

from __future__ import annotations

import time
import threading
from typing import Callable, TypeVar

from ..exceptions import CircuitOpenError
from ..logging_config import get_logger

logger = get_logger(__name__)

T = TypeVar("T")


class CircuitBreaker:
    """Wraps a callable with circuit-breaker protection."""

    def __init__(
        self,
        service: str,
        failure_threshold: int = 5,
        recovery_timeout: float = 60.0,
    ) -> None:
        self.service = service
        self._failure_threshold = failure_threshold
        self._recovery_timeout = recovery_timeout
        self._failure_count = 0
        self._state = "closed"
        self._last_failure_time = 0.0
        self._lock = threading.Lock()

    @property
    def state(self) -> str:
        with self._lock:
            if self._state == "open":
                if time.monotonic() - self._last_failure_time > self._recovery_timeout:
                    self._state = "half_open"
            return self._state

    def call(self, fn: Callable[..., T], *args, **kwargs) -> T:
        """Execute *fn* through the circuit breaker."""
        current_state = self.state

        if current_state == "open":
            raise CircuitOpenError(
                self.service,
                f"Circuit open for '{self.service}'; retry after {self._recovery_timeout}s",
            )

        try:
            result = fn(*args, **kwargs)
        except Exception as exc:
            self._record_failure()
            raise
        else:
            self._record_success()
            return result

    def _record_failure(self) -> None:
        with self._lock:
            self._failure_count += 1
            self._last_failure_time = time.monotonic()
            if self._failure_count >= self._failure_threshold:
                logger.warning(
                    "circuit_opened",
                    extra={"service": self.service, "failures": self._failure_count},
                )
                self._state = "open"

    def _record_success(self) -> None:
        with self._lock:
            if self._state == "half_open":
                logger.info(
                    "circuit_closed",
                    extra={"service": self.service},
                )
            self._failure_count = 0
            self._state = "closed"

    def reset(self) -> None:
        """Manually reset the circuit breaker."""
        with self._lock:
            self._failure_count = 0
            self._state = "closed"
            self._last_failure_time = 0.0
