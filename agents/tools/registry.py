"""Tool registry for organizing and discovering agent tools."""

from __future__ import annotations

from .base import Tool


class ToolRegistry:
    """Registry of available tools, optionally organized by domain."""

    def __init__(self) -> None:
        self._tools: dict[str, Tool] = {}
        self._domain_map: dict[str, list[str]] = {}

    def register(self, tool: Tool, domain: str | None = None) -> None:
        """Register a tool, optionally tagging it with a domain."""
        self._tools[tool.name] = tool
        if domain:
            self._domain_map.setdefault(domain, []).append(tool.name)

    def get(self, name: str) -> Tool:
        """Retrieve a tool by name. Raises KeyError if not found."""
        try:
            return self._tools[name]
        except KeyError:
            available = ", ".join(sorted(self._tools)) or "(none)"
            raise KeyError(f"Tool '{name}' not found. Available: {available}")

    def list_tools(self, domain: str | None = None) -> list[Tool]:
        """List all tools, optionally filtered by domain."""
        if domain is None:
            return list(self._tools.values())
        names = self._domain_map.get(domain, [])
        return [self._tools[n] for n in names if n in self._tools]

    def has(self, name: str) -> bool:
        return name in self._tools

    def render_tool_descriptions(self, domain: str | None = None) -> str:
        """Render all tool schemas for injection into a ReAct prompt."""
        tools = self.list_tools(domain)
        if not tools:
            return "No tools available."
        return "\n\n".join(tool.schema_for_prompt() for tool in tools)
