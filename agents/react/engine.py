"""ReAct (Reason + Act) reasoning loop engine.

Each agent worker gets a ReActEngine that orchestrates a think→act→observe
loop until the agent produces a Final Answer or exhausts its step budget.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Optional

from ..exceptions import BudgetExhaustedError, LLMError, ToolExecutionError
from ..logging_config import get_logger
from ..tools.base import ToolOutput
from ..tools.registry import ToolRegistry
from .parser import ParsedReActResponse, parse_react_response

logger = get_logger(__name__)


@dataclass
class ReActStep:
    """One iteration of the think→act→observe loop."""
    thought: str
    action: str
    action_input: dict
    observation: str
    duration_ms: float = 0.0


@dataclass
class ReActTrace:
    """Full trace of a ReAct execution."""
    steps: list[ReActStep] = field(default_factory=list)
    final_answer: str = ""
    total_tokens: int = 0
    total_tool_calls: int = 0
    total_duration_ms: float = 0.0

    def summary(self) -> str:
        lines = [f"ReAct trace: {len(self.steps)} steps, {self.total_tool_calls} tool calls"]
        for i, step in enumerate(self.steps, 1):
            lines.append(f"  Step {i}: {step.action}({list(step.action_input.keys())}) → {step.observation[:80]}...")
        return "\n".join(lines)


REACT_SYSTEM_PROMPT = """You are {role}. You analyze evidence methodically using available tools.

Available tools:
{tool_descriptions}

For EACH reasoning step, use this exact format:
Thought: <what you need to find out next>
Action: <tool_name>
Action Input: {{"param1": "value1", "param2": "value2"}}

After a tool runs, you will see:
Observation: <tool result>

When you have gathered enough evidence to answer, respond with:
Thought: I have sufficient evidence to provide a comprehensive answer.
Final Answer: <your detailed, structured answer>

Important rules:
- You have at most {max_steps} tool calls. Use them wisely.
- Always start with the vector index search tool to check existing evidence.
- Use external API tools only when local evidence is insufficient.
- Cite specific evidence in your Final Answer.
- Stay factual — do not hallucinate data."""


class ReActEngine:
    """Executes a ReAct reasoning loop with a tool registry."""

    def __init__(
        self,
        llm_client,
        tools: ToolRegistry,
        max_steps: int = 6,
        max_tokens_per_turn: int = 512,
    ) -> None:
        self.llm = llm_client
        self.tools = tools
        self.max_steps = max_steps
        self.max_tokens_per_turn = max_tokens_per_turn

    def run(
        self,
        role: str,
        query: str,
        context: str = "",
        temperature: float = 0.1,
    ) -> ReActTrace:
        """Execute the ReAct loop until Final Answer or budget exhaustion."""
        start_time = time.monotonic()
        trace = ReActTrace()

        system_prompt = REACT_SYSTEM_PROMPT.format(
            role=role,
            tool_descriptions=self.tools.render_tool_descriptions(),
            max_steps=self.max_steps,
        )

        # Build initial conversation
        conversation: list[dict[str, str]] = []
        if context:
            conversation.append({
                "role": "user",
                "content": f"Context:\n{context}\n\nQuery: {query}",
            })
        else:
            conversation.append({"role": "user", "content": query})

        for step_idx in range(self.max_steps):
            step_start = time.monotonic()

            # Generate next ReAct step
            try:
                response_text = self.llm.generate_chat(
                    messages=conversation,
                    system_prompt=system_prompt,
                    temperature=temperature,
                    max_tokens=self.max_tokens_per_turn,
                )
            except Exception as exc:
                logger.warning("react_llm_failed", extra={"step": step_idx, "error": str(exc)})
                # Force synthesis from what we have
                trace.final_answer = self._emergency_synthesis(trace)
                break

            # Estimate tokens
            trace.total_tokens += len(response_text) // 4

            # Parse the response
            parsed = parse_react_response(response_text)

            # Check for Final Answer
            if parsed.final_answer:
                trace.final_answer = parsed.final_answer
                break

            # Execute the tool
            if not parsed.action:
                # LLM didn't follow format — treat as final answer
                trace.final_answer = response_text
                break

            tool_output = self._execute_tool(parsed.action, parsed.action_input or {})
            observation = str(tool_output)

            step_duration = (time.monotonic() - step_start) * 1000
            step = ReActStep(
                thought=parsed.thought or "",
                action=parsed.action,
                action_input=parsed.action_input or {},
                observation=observation,
                duration_ms=step_duration,
            )
            trace.steps.append(step)
            trace.total_tool_calls += 1

            logger.info(
                "react_step",
                extra={
                    "step": step_idx + 1,
                    "action": parsed.action,
                    "success": tool_output.success,
                    "duration_ms": round(step_duration, 1),
                },
            )

            # Add to conversation for next turn
            conversation.append({"role": "assistant", "content": response_text})
            conversation.append({"role": "user", "content": f"Observation: {observation}"})

        else:
            # Budget exhausted — force a final answer
            logger.warning("react_budget_exhausted", extra={"steps": self.max_steps})
            trace.final_answer = self._force_final_answer(trace, conversation, system_prompt, temperature)

        trace.total_duration_ms = (time.monotonic() - start_time) * 1000
        return trace

    def _execute_tool(self, tool_name: str, tool_input: dict) -> ToolOutput:
        """Look up and execute a tool by name."""
        if not self.tools.has(tool_name):
            return ToolOutput(
                success=False,
                error=f"Unknown tool '{tool_name}'. Available: {', '.join(t.name for t in self.tools.list_tools())}",
            )
        tool = self.tools.get(tool_name)
        try:
            # Build the proper ToolInput from the dict
            input_cls = tool._input_class()
            if input_cls:
                typed_input = input_cls.model_validate(tool_input)
            else:
                from ..tools.base import ToolInput as BaseInput
                typed_input = BaseInput()
            return tool.execute(typed_input)
        except Exception as exc:
            return ToolOutput(success=False, error=f"Tool '{tool_name}' error: {str(exc)[:200]}")

    def _force_final_answer(
        self,
        trace: ReActTrace,
        conversation: list[dict],
        system_prompt: str,
        temperature: float,
    ) -> str:
        """Force the LLM to produce a Final Answer given observations so far."""
        force_prompt = (
            "You have used all your tool calls. Based on the evidence gathered so far, "
            "you MUST now provide your Final Answer.\n\n"
            "Final Answer:"
        )
        conversation.append({"role": "user", "content": force_prompt})
        try:
            response = self.llm.generate_chat(
                messages=conversation,
                system_prompt=system_prompt,
                temperature=temperature,
                max_tokens=self.max_tokens_per_turn,
            )
            parsed = parse_react_response(response)
            return parsed.final_answer or response
        except Exception:
            return self._emergency_synthesis(trace)

    def _emergency_synthesis(self, trace: ReActTrace) -> str:
        """Build a minimal answer from observations when the LLM fails."""
        if not trace.steps:
            return "Unable to gather evidence. No tool calls succeeded."
        observations = [
            f"- {step.action}: {step.observation[:200]}"
            for step in trace.steps
            if step.observation and "ERROR" not in step.observation
        ]
        if observations:
            return "Based on available evidence:\n" + "\n".join(observations)
        return "Evidence gathering encountered errors. Manual review required."
