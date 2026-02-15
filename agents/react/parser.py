"""ReAct output parser — extracts Thought/Action/Action Input/Final Answer
from LLM responses using multiple strategies for robustness."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from typing import Optional


@dataclass
class ParsedReActResponse:
    thought: str | None
    action: str | None
    action_input: dict | None
    final_answer: str | None
    raw: str


def parse_react_response(raw: str) -> ParsedReActResponse:
    """Parse LLM output in ReAct format.

    Tries multiple strategies:
      1. Regex extraction for Thought/Action/Action Input/Final Answer blocks
      2. JSON-embedded action input fallback
      3. Treats entire output as final answer if nothing matches
    """
    text = raw.strip()

    # Check for Final Answer first
    final_match = re.search(
        r"Final\s+Answer\s*:\s*(.*)",
        text,
        re.DOTALL | re.IGNORECASE,
    )
    if final_match:
        thought = _extract_thought(text)
        return ParsedReActResponse(
            thought=thought,
            action=None,
            action_input=None,
            final_answer=final_match.group(1).strip(),
            raw=raw,
        )

    # Extract Thought
    thought = _extract_thought(text)

    # Extract Action
    action_match = re.search(
        r"Action\s*:\s*(\S+)",
        text,
        re.IGNORECASE,
    )
    action = action_match.group(1).strip() if action_match else None

    # Extract Action Input
    action_input = _extract_action_input(text)

    # If we have an action, return it
    if action:
        return ParsedReActResponse(
            thought=thought,
            action=action,
            action_input=action_input or {},
            final_answer=None,
            raw=raw,
        )

    # Nothing matched — treat as implicit final answer
    return ParsedReActResponse(
        thought=thought,
        action=None,
        action_input=None,
        final_answer=text if text else None,
        raw=raw,
    )


def _extract_thought(text: str) -> str | None:
    match = re.search(
        r"Thought\s*:\s*(.*?)(?=\n(?:Action|Final\s+Answer)\s*:|$)",
        text,
        re.DOTALL | re.IGNORECASE,
    )
    return match.group(1).strip() if match else None


def _extract_action_input(text: str) -> dict | None:
    """Extract the Action Input value, handling both JSON and simple text."""
    match = re.search(
        r"Action\s+Input\s*:\s*(.*?)(?=\n(?:Thought|Action|Observation|Final\s+Answer)\s*:|$)",
        text,
        re.DOTALL | re.IGNORECASE,
    )
    if not match:
        return None

    raw_input = match.group(1).strip()

    # Try JSON parsing
    try:
        parsed = json.loads(raw_input)
        if isinstance(parsed, dict):
            return parsed
    except (json.JSONDecodeError, TypeError):
        pass

    # Try extracting embedded JSON
    json_match = re.search(r"(\{[\s\S]*\})", raw_input)
    if json_match:
        try:
            parsed = json.loads(json_match.group(1))
            if isinstance(parsed, dict):
                return parsed
        except (json.JSONDecodeError, TypeError):
            pass

    # Fallback: treat as a single "query" parameter
    if raw_input:
        return {"query": raw_input}

    return None


def parse_structured_response(raw: str, expected_keys: list[str]) -> dict:
    """Parse a KEY: VALUE formatted LLM response (used for synthesis).

    Example input:
        STORY: some text here
        RECOMMENDATION: Go
        RATIONALE: because of X

    Returns dict with lowercased keys.
    """
    text = raw.strip().replace("\r\n", "\n")
    data: dict = {}

    for key in expected_keys:
        # Build pattern that matches KEY: value until the next KEY: or end
        other_keys = [k for k in expected_keys if k != key]
        if other_keys:
            lookahead = "|".join(re.escape(k) for k in other_keys)
            pattern = rf"{re.escape(key)}\s*:\s*(.*?)(?=\n(?:{lookahead})\s*:|$)"
        else:
            pattern = rf"{re.escape(key)}\s*:\s*(.*?)$"

        match = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
        if match:
            data[key.lower()] = match.group(1).strip()

    # Fallback: try JSON
    if not data:
        try:
            data = json.loads(text)
        except (json.JSONDecodeError, TypeError):
            pass

        json_match = re.search(r"(\{[\s\S]*\})", text)
        if json_match and not data:
            try:
                data = json.loads(json_match.group(1))
            except (json.JSONDecodeError, TypeError):
                pass

    return data
