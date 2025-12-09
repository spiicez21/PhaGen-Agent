from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Optional

import httpx
from dotenv import load_dotenv

# Load .env from repo root, overriding system environment variables
_env_path = Path(__file__).resolve().parents[1] / ".env"
load_dotenv(_env_path, override=True)

logger = logging.getLogger(__name__)


class LLMClientError(RuntimeError):
    """Raised when the LLM runtime cannot fulfill a request."""


class LLMClient:
    def __init__(
        self,
        provider: Optional[str] = None,
        model: Optional[str] = None,
        base_url: Optional[str] = None,
        api_key: Optional[str] = None,
        timeout: float = 60.0,
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
        logger.info(f"LLMClient initialized: provider={self.provider}, model={self.model}, base_url={self.base_url}")

    def generate(
        self,
        *,
        prompt: str,
        system_prompt: Optional[str] = None,
        temperature: float = 0.1,
        max_tokens: int = 512,
    ) -> str:
        if not prompt.strip():
            raise LLMClientError("Prompt cannot be empty")

        if self.provider == "ollama":
            return self._generate_with_ollama(
                prompt=prompt,
                system_prompt=system_prompt,
                temperature=temperature,
                max_tokens=max_tokens,
            )
        if self.provider == "openai":
            return self._generate_with_openai(
                prompt=prompt,
                system_prompt=system_prompt,
                temperature=temperature,
                max_tokens=max_tokens,
            )
        raise LLMClientError(f"Unsupported LLM provider '{self.provider}'")

    # Provider implementations ---------------------------------------
    def _generate_with_ollama(
        self,
        *,
        prompt: str,
        system_prompt: Optional[str],
        temperature: float,
        max_tokens: int,
    ) -> str:
        url = f"{self.base_url.rstrip('/')}/api/chat"
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})
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
            logger.debug(f"Sending Ollama request to {url} with model={self.model}")
            response = httpx.post(url, json=payload, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            message = data.get("message", {}).get("content")
            if message:
                return message.strip()
            if "response" in data:
                return str(data["response"]).strip()
        except httpx.HTTPError as exc:
            error_detail = ""
            try:
                error_detail = f" - {exc.response.text}" if hasattr(exc, 'response') else ""
            except:
                pass
            logger.error(f"Ollama request failed: {exc}{error_detail}")
            raise LLMClientError(f"Ollama request failed: {exc}{error_detail}") from exc
        raise LLMClientError("Ollama returned an empty response")

    def _generate_with_openai(
        self,
        *,
        prompt: str,
        system_prompt: Optional[str],
        temperature: float,
        max_tokens: int,
    ) -> str:
        if not self.api_key:
            raise LLMClientError("LLM_API_KEY is required for OpenAI provider")
        url = f"{self.base_url.rstrip('/')}/chat/completions"
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})
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
            response = httpx.post(
                url, json=payload, headers=headers, timeout=self.timeout
            )
            response.raise_for_status()
            data = response.json()
            choices = data.get("choices") or []
            if choices:
                content = choices[0].get("message", {}).get("content")
                if content:
                    return content.strip()
        except httpx.HTTPError as exc:
            raise LLMClientError("OpenAI request failed") from exc
        raise LLMClientError("OpenAI returned an empty response")


def format_structured_context(payload: object, limit: int = 1200) -> str:
    """Utility to stringify structured context within a safe character budget."""
    try:
        serialized = json.dumps(payload, ensure_ascii=False, indent=2)
    except TypeError:
        serialized = str(payload)
    if len(serialized) <= limit:
        return serialized
    return serialized[: limit - 3] + "..."
