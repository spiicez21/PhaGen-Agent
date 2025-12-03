"""Build a local Chroma (FAISS-compatible) index from crawler outputs.

Usage:
    python indexes/build_index.py

Requirements:
    pip install chromadb sentence-transformers
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Iterable, List

try:
    import chromadb
    from chromadb.utils.embedding_functions import (
        SentenceTransformerEmbeddingFunction,
    )
except ImportError as exc:  # pragma: no cover - runtime guard
    raise SystemExit(
        "chromadb and sentence-transformers are required. Run 'pip install chromadb sentence-transformers'"
    ) from exc

DATASET_DIR = Path(__file__).resolve().parents[1] / "crawler" / "storage" / "datasets" / "default"
PERSIST_DIR = Path(__file__).resolve().parent / "chroma"
COLLECTION_NAME = "phagen-agentic"


def iter_records(directory: Path) -> Iterable[dict]:
    if not directory.exists():
        raise SystemExit(f"Dataset directory {directory} does not exist. Run the crawler first.")

    for file in sorted(directory.glob("*.json")):
        with file.open("r", encoding="utf-8") as handle:
            yield json.load(handle)


def record_to_document(record: dict) -> tuple[str, dict, str]:
    text = ""
    if record.get("via") == "api" and record.get("payload"):
        text = json.dumps(record["payload"], ensure_ascii=False)
    else:
        text = record.get("snippet") or record.get("text") or ""
    metadata = {
        "source_type": record.get("source_type", "unknown"),
        "origin": record.get("via", "html"),
        "url": record.get("url", ""),
    }
    doc_id = record.get("id") or record.get("url") or str(hash(text))
    return doc_id, metadata, text


def build_index(records: Iterable[dict]) -> None:
    embedding_fn = SentenceTransformerEmbeddingFunction(model_name="all-MiniLM-L6-v2")
    client = chromadb.PersistentClient(path=str(PERSIST_DIR))

    try:
        client.delete_collection(COLLECTION_NAME)
    except ValueError:
        pass

    collection = client.get_or_create_collection(
        name=COLLECTION_NAME,
        embedding_function=embedding_fn,
    )

    documents: List[str] = []
    metadatas: List[dict] = []
    ids: List[str] = []
    for record in records:
        doc_id, metadata, text = record_to_document(record)
        if not text:
            continue
        documents.append(text)
        metadatas.append(metadata)
        ids.append(doc_id)

    if not documents:
        raise SystemExit("No documents found to index.")

    collection.add(documents=documents, metadatas=metadatas, ids=ids)
    print(f"Indexed {len(documents)} documents into {PERSIST_DIR}")


if __name__ == "__main__":  # pragma: no cover
    build_index(iter_records(DATASET_DIR))
