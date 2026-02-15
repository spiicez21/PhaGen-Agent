"""Base abstractions for agent tools.

Every capability an agent can invoke (API calls, parsers, calculators)
is modeled as a Tool with typed input/output and a prompt-friendly schema.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

from pydantic import BaseModel, Field


class ToolInput(BaseModel):
    """Base class for all tool inputs.  Each concrete tool subclasses this."""

    class Config:
        extra = "forbid"


class ToolOutput(BaseModel):
    """Standardized tool output returned to the ReAct loop."""

    success: bool
    data: Any = None
    error: str | None = None
    metadata: dict = Field(default_factory=dict)

    def __str__(self) -> str:
        if self.success:
            # Truncate large data for inclusion in ReAct observation
            text = str(self.data)
            if len(text) > 2000:
                text = text[:2000] + "... [truncated]"
            return text
        return f"ERROR: {self.error}"


class Tool(ABC):
    """Abstract tool that agents can invoke via the ReAct loop."""

    name: str
    description: str  # Injected into the ReAct system prompt

    @abstractmethod
    def execute(self, tool_input: ToolInput) -> ToolOutput:
        """Run the tool and return a structured output."""
        ...

    def input_schema(self) -> dict:
        """Return the JSON schema for this tool's input model."""
        input_cls = self._input_class()
        return input_cls.model_json_schema() if input_cls else {}

    def _input_class(self) -> type[ToolInput] | None:
        """Override to return the concrete ToolInput subclass for this tool."""
        return None

    def schema_for_prompt(self) -> str:
        """Render a human-readable description for injection into ReAct prompts."""
        schema = self.input_schema()
        props = schema.get("properties", {})
        required = set(schema.get("required", []))
        lines = [f"  {self.name}: {self.description}"]
        if props:
            lines.append("    Parameters:")
            for param_name, param_info in props.items():
                req_marker = " (required)" if param_name in required else ""
                desc = param_info.get("description", param_info.get("type", ""))
                lines.append(f"      - {param_name}{req_marker}: {desc}")
        return "\n".join(lines)
