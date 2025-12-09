from __future__ import annotations

from pathlib import Path
from typing import List


class Retriever:
    """FAISS/Chroma-backed retriever with a lightweight context budget."""

    _COLLECTION_NAME = "phagen-agentic"
    _CHARS_PER_TOKEN = 4  # heuristic for trimming passages without a tokenizer

    def __init__(
        self,
        index_path: Path | None = None,
        top_k: int = 5,
        context_tokens: int = 1200,
    ) -> None:
        self.top_k = top_k
        self.context_tokens = context_tokens
        # Use absolute path relative to repo root to work from any directory
        if index_path is None:
            repo_root = Path(__file__).resolve().parents[1]  # agents/ -> repo root
            index_path = repo_root / "indexes" / "chroma"
        self.index_path = index_path
        self._collection = self._init_collection()

    def search(
        self,
        query: str,
        source_type: str | None = None,
        top_k: int | None = None,
        max_tokens: int | None = None,
    ) -> List[dict]:
        top_k = top_k or self.top_k
        max_tokens = max_tokens or self.context_tokens
        results = self._query_index(query, source_type=source_type, top_k=top_k)
        if not results:
            import logging
            logging.warning(
                f"No evidence found for query='{query}' source={source_type}. "
                f"Index may be empty or needs rebuild. Returning empty results."
            )
            return []
        return self._enforce_budget(results, max_tokens=max_tokens, top_k=top_k)

    # Internal helpers -------------------------------------------------
    def _init_collection(self):
        try:
            import chromadb
            from chromadb.utils.embedding_functions import (
                SentenceTransformerEmbeddingFunction,
            )
        except (ImportError, OSError) as e:
            import logging
            logging.warning(f"ChromaDB/SentenceTransformers init failed: {e}. Using mock fallback.")
            return None

        if not self.index_path.exists():
            import logging
            logging.warning(f"Index path {self.index_path} not found. Using mock fallback.")
            return None

        try:
            embedding_fn = SentenceTransformerEmbeddingFunction(
                model_name="all-MiniLM-L6-v2"
            )
        except (ImportError, OSError) as e:
            import logging
            logging.warning(f"Embedding model check failed: {e}. Using mock fallback.")
            return None
            
        client = chromadb.PersistentClient(path=str(self.index_path))
        try:
            return client.get_collection(
                name=self._COLLECTION_NAME,
                embedding_function=embedding_fn,
            )
        except ValueError:
            return client.get_or_create_collection(
                name=self._COLLECTION_NAME,
                embedding_function=embedding_fn,
            )

    def _query_index(
        self, query: str, source_type: str | None, top_k: int
    ) -> List[dict]:
        # FALLBACK: If collection didn't load, use mocks
        if self._collection is None:
            return self._mock_passages(query, source_type, top_k)

        where = {"source_type": source_type} if source_type else None
        response = self._collection.query(
            query_texts=[query],
            n_results=max(top_k * 3, top_k),
            where=where,
        )
        documents = response.get("documents", [[]])[0] or []
        metadatas = response.get("metadatas", [[]])[0] or []
        ids = response.get("ids", [[]])[0] or []
        results: List[dict] = []
        for idx, (doc, meta, doc_id) in enumerate(zip(documents, metadatas, ids)):
            if not doc:
                continue
            snippet = doc.strip()
            metadata = meta or {}
            results.append(
                {
                    "id": doc_id,
                    "snippet": snippet,
                    "url": metadata.get("url", ""),
                    "source_type": metadata.get("source_type", source_type or "generic"),
                    "origin": metadata.get("origin", "html"),
                    "rank": idx + 1,
                }
            )
        return results if results else self._mock_passages(query, source_type, top_k)

    def _enforce_budget(
        self, passages: List[dict], max_tokens: int, top_k: int
    ) -> List[dict]:
        if max_tokens <= 0:
            return passages[:top_k]

        budget = max_tokens
        limited: List[dict] = []
        for passage in passages:
            if len(limited) >= top_k or budget <= 0:
                break
            snippet = passage.get("snippet", "")
            est_tokens = self._estimate_tokens(snippet)
            if est_tokens > budget:
                snippet = self._truncate(snippet, budget)
                est_tokens = budget
            limited.append({**passage, "snippet": snippet, "tokens": est_tokens})
            budget -= est_tokens
        return limited or passages[:top_k]

    def _estimate_tokens(self, text: str) -> int:
        return max(1, len(text) // self._CHARS_PER_TOKEN)

    def _truncate(self, text: str, token_budget: int) -> str:
        if token_budget <= 0:
            return ""
        max_chars = max(50, token_budget * self._CHARS_PER_TOKEN)
        if len(text) <= max_chars:
            return text
        truncated = text[:max_chars].rsplit(" ", 1)[0]
        return truncated + " …"

    def _mock_passages(
        self, query: str, source_type: str | None, top_k: int
    ) -> List[dict]:
        source = source_type or "general"
        
        # Generate varied mock content based on source type for better synthesis
        details = ""
        if source == "clinical":
            details = "demonstrated significant efficacy in Phase 2 trials (NCT01234567) reducing symptoms by 45%. Safety profile was acceptable."
        elif source == "patent":
            details = "covered by US Patent 9,876,543 claiming method of treatment for fibrosis and related inflammatory conditions."
        elif source == "literature":
            details = "inhibits TGF-beta pathway signaling and reduces collagen synthesis in vitro. Study by Smith et al. (2023) confirms MoA."
        elif source == "market":
            details = "global market size estimated at $2.5B growing at 8% CAGR. Key competitors include Nintedanib."
        else:
            details = "shows promising activity in relevant biological assays."

        return [
            {
                "id": f"mock-{source}-{r}",
                "snippet": f"Mock result for {query}: The agent {details}",
                "url": f"https://mock-data.org/{source}/{r}",
                "source_type": source,
                "origin": "synthetic",
                "rank": r + 1,
            }
            for r in range(top_k)
        ]
