from __future__ import annotations

from pathlib import Path
from typing import List


class Retriever:
    """Thin facade for semantic + keyword retrieval (mock for now)."""

    def __init__(self, index_path: Path | None = None, top_k: int = 5) -> None:
        self.index_path = index_path or Path("indexes/mock-index.json")
        self.top_k = top_k

    def search(self, query: str, source_type: str | None = None) -> List[dict]:
        # Placeholder deterministic snippets; replace with FAISS/Chroma binding later.
        base = {
            "query": query,
            "source_type": source_type or "generic",
            "url": "https://example.org/evidence",
            "snippet": f"Mock passage referencing {query}",
        }
        return [base | {"rank": r + 1} for r in range(self.top_k)]
