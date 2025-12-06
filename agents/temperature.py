from __future__ import annotations

import os
from typing import Dict

DEFAULT_MASTER_TEMPERATURE = 0.2
DEFAULT_WORKER_TEMPERATURES: Dict[str, float] = {
    "clinical": 0.12,
    "literature": 0.18,
    "patent": 0.1,
    "market": 0.22,
}


def _clamp_temperature(value: float) -> float:
    return max(0.0, min(1.0, value))


def _parse_temperature(raw: str | float | None) -> float | None:
    if raw is None:
        return None
    if isinstance(raw, (int, float)):
        return _clamp_temperature(float(raw))
    try:
        return _clamp_temperature(float(raw))
    except (ValueError, TypeError):
        return None


def resolve_worker_temperature(worker_name: str, override: float | None = None) -> float:
    direct = _parse_temperature(override)
    if direct is not None:
        return direct
    env_key = f"LLM_TEMP_{worker_name.upper()}"
    env_value = _parse_temperature(os.getenv(env_key))
    if env_value is not None:
        return env_value
    default_env = _parse_temperature(os.getenv("LLM_TEMP_DEFAULT"))
    if default_env is not None:
        return default_env
    return DEFAULT_WORKER_TEMPERATURES.get(worker_name, 0.15)


def resolve_worker_temperatures(overrides: Dict[str, float] | None = None) -> Dict[str, float]:
    resolved: Dict[str, float] = {}
    names = set(DEFAULT_WORKER_TEMPERATURES.keys())
    if overrides:
        names.update(overrides.keys())
    for name in names:
        override_value = overrides.get(name) if overrides else None
        resolved[name] = resolve_worker_temperature(name, override_value)
    return resolved


def resolve_master_temperature(override: float | None = None) -> float:
    direct = _parse_temperature(override)
    if direct is not None:
        return direct
    env_value = _parse_temperature(os.getenv("LLM_TEMP_MASTER"))
    if env_value is not None:
        return env_value
    default_env = _parse_temperature(os.getenv("LLM_TEMP_DEFAULT"))
    if default_env is not None:
        return default_env
    return DEFAULT_MASTER_TEMPERATURE
