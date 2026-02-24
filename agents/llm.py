from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import httpx
from dotenv import load_dotenv

from .exceptions import LLMError, LLMTimeoutError
from .logging_config import get_logger

# Load .env from repo root, overriding system environment variables
_env_path = Path(__file__).resolve().parents[1] / ".env"
load_dotenv(_env_path, override=True)

logger = get_logger(__name__)


# Keep backward-compatible alias
LLMClientError = LLMError


@dataclass
class LLMUsage:
    """Token usage stats returned from an LLM call."""
    prompt_tokens: int = 0
    completion_tokens: int = 0
    total_tokens: int = 0


class LLMClient:
    """Unified LLM client supporting Ollama and OpenAI-compatible APIs.

    Enhancements over the original:
      - ``generate_chat()`` for multi-turn ReAct conversations
      - Per-call timeout enforcement
      - Token usage tracking
      - JSON response mode support
    """

    def __init__(
        self,
        provider: Optional[str] = None,
        model: Optional[str] = None,
        base_url: Optional[str] = None,
        api_key: Optional[str] = None,
        timeout: float = float(os.getenv("LLM_TIMEOUT", "240")),
    ) -> None:
        self.provider = (provider or os.getenv("LLM_PROVIDER", "ollama")).lower()
        self.model = model or os.getenv("LLM_MODEL", "smollm:360m")
        default_base = (
            "http://localhost:11434"
            if self.provider == "ollama"
            else "https://api.openai.com/v1"
        )
        self.base_url = base_url or os.getenv("LLM_BASE_URL", default_base)
        self.api_key = api_key or os.getenv("LLM_API_KEY")
        self.timeout = timeout
        self.last_usage: LLMUsage = LLMUsage()
        logger.info(
            "llm_client_init",
            extra={"provider": self.provider, "model": self.model, "base_url": self.base_url},
        )

    # -- Single-turn generation (backward-compatible) ----------------------

    def generate(
        self,
        *,
        prompt: str,
        system_prompt: Optional[str] = None,
        temperature: float = 0.1,
        max_tokens: int = 512,
        timeout: float | None = None,
    ) -> str:
        """Generate a single response from a prompt string."""
        if not prompt.strip():
            raise LLMError("Prompt cannot be empty")

        messages: list[dict[str, str]] = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})

        return self._call(
            messages=messages,
            temperature=temperature,
            max_tokens=max_tokens,
            timeout=timeout or self.timeout,
        )

    # -- Multi-turn chat (for ReAct loop) ----------------------------------

    def generate_chat(
        self,
        messages: list[dict[str, str]],
        *,
        system_prompt: Optional[str] = None,
        temperature: float = 0.1,
        max_tokens: int = 512,
        timeout: float | None = None,
    ) -> str:
        """Generate a response from a multi-turn conversation.

        ``messages`` should be a list of ``{"role": "user"|"assistant", "content": ...}``
        dicts.  A system prompt is prepended automatically if provided.
        """
        full_messages: list[dict[str, str]] = []
        if system_prompt:
            full_messages.append({"role": "system", "content": system_prompt})
        full_messages.extend(messages)

        return self._call(
            messages=full_messages,
            temperature=temperature,
            max_tokens=max_tokens,
            timeout=timeout or self.timeout,
        )

    # -- Internal dispatch -------------------------------------------------

    def _call(
        self,
        messages: list[dict[str, str]],
        temperature: float,
        max_tokens: int,
        timeout: float,
    ) -> str:
        if self.provider == "ollama":
            return self._call_ollama(messages, temperature, max_tokens, timeout)
        if self.provider == "openai":
            return self._call_openai(messages, temperature, max_tokens, timeout)
        raise LLMError(f"Unsupported LLM provider '{self.provider}'")

    # -- Ollama ------------------------------------------------------------

    def _call_ollama(
        self,
        messages: list[dict[str, str]],
        temperature: float,
        max_tokens: int,
        timeout: float,
    ) -> str:
        url = f"{self.base_url.rstrip('/')}/api/chat"
        payload = {
            "model": self.model,
            "messages": messages,
            "stream": False,
            "options": {
                "temperature": temperature,
                "num_predict": max_tokens,
            },
        }
        try:
            response = httpx.post(url, json=payload, timeout=timeout)
            response.raise_for_status()
            data = response.json()

            # Extract token usage from Ollama response
            self.last_usage = LLMUsage(
                prompt_tokens=data.get("prompt_eval_count", 0),
                completion_tokens=data.get("eval_count", 0),
                total_tokens=data.get("prompt_eval_count", 0) + data.get("eval_count", 0),
            )

            content = data.get("message", {}).get("content")
            if content:
                return content.strip()
            if "response" in data:
                return str(data["response"]).strip()
        except httpx.TimeoutException as exc:
            raise LLMTimeoutError(f"Ollama request timed out after {timeout}s") from exc
        except httpx.HTTPError as exc:
            error_detail = ""
            try:
                error_detail = f" - {exc.response.text}" if hasattr(exc, "response") and exc.response else ""
            except Exception:
                pass
            logger.error("ollama_request_failed", extra={"error": str(exc)})
            raise LLMError(f"Ollama request failed: {exc}{error_detail}") from exc
        raise LLMError("Ollama returned an empty response")

    # -- OpenAI-compatible -------------------------------------------------

    def _call_openai(
        self,
        messages: list[dict[str, str]],
        temperature: float,
        max_tokens: int,
        timeout: float,
    ) -> str:
        if not self.api_key:
            raise LLMError("LLM_API_KEY is required for OpenAI provider")
        url = f"{self.base_url.rstrip('/')}/chat/completions"
        payload = {
            "model": self.model or "gpt-4o-mini",
            "messages": messages,
            "temperature": temperature,
            "max_tokens": max_tokens,
        }
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
        }
        try:
            response = httpx.post(url, json=payload, headers=headers, timeout=timeout)
            response.raise_for_status()
            data = response.json()

            # Extract token usage
            usage = data.get("usage", {})
            self.last_usage = LLMUsage(
                prompt_tokens=usage.get("prompt_tokens", 0),
                completion_tokens=usage.get("completion_tokens", 0),
                total_tokens=usage.get("total_tokens", 0),
            )

            choices = data.get("choices") or []
            if choices:
                content = choices[0].get("message", {}).get("content")
                if content:
                    return content.strip()
        except httpx.TimeoutException as exc:
            raise LLMTimeoutError(f"OpenAI request timed out after {timeout}s") from exc
        except httpx.HTTPError as exc:
            raise LLMError("OpenAI request failed") from exc
        raise LLMError("OpenAI returned an empty response")


def format_structured_context(payload: object, limit: int = 1200) -> str:
    """Utility to stringify structured context within a safe character budget."""
    try:
        serialized = json.dumps(payload, ensure_ascii=False, indent=2)
    except TypeError:
        serialized = str(payload)
    if len(serialized) <= limit:
        return serialized
    return serialized[: limit - 3] + "..."
