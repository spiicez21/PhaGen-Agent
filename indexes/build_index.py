"""Build a local Chroma (FAISS-compatible) index from crawler outputs.

Adds optional snapshotting so daily/monthly builds stay reproducible.

Usage examples:
    python indexes/build_index.py                             # build + daily snapshot (default)
    python indexes/build_index.py --cadence monthly           # build + monthly snapshot name
    python indexes/build_index.py --snapshot-name demo-run    # custom snapshot folder name
    python indexes/build_index.py --no-snapshot               # build only (legacy behavior)

Requirements:
    pip install chromadb sentence-transformers
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
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
SNAPSHOT_ROOT = Path(__file__).resolve().parent / "chroma_snapshots"
COLLECTION_NAME = "phagen-agentic"
MANIFEST_NAME = "manifest.json"


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Chroma index with optional snapshots")
    parser.add_argument(
        "--dataset",
        type=Path,
        default=DATASET_DIR,
        help="Path to crawler dataset directory (default: crawler/storage/datasets/default)",
    )
    parser.add_argument(
        "--cadence",
        choices=["daily", "monthly"],
        default="daily",
        help="Snapshot cadence label used when auto-naming snapshots (default: daily)",
    )
    parser.add_argument(
        "--snapshot-name",
        type=str,
        help="Explicit snapshot folder name (overrides cadence-based naming)",
    )
    parser.add_argument(
        "--no-snapshot",
        action="store_true",
        help="Skip copying the built index into the snapshots directory",
    )
    parser.add_argument(
        "--keep-daily",
        type=int,
        default=14,
        help="How many daily snapshots to retain (default: 14 days)",
    )
    parser.add_argument(
        "--keep-monthly",
        type=int,
        default=12,
        help="How many monthly snapshots to retain (default: 12 months)",
    )
    parser.add_argument(
        "--no-retention",
        action="store_true",
        help="Disable automatic snapshot retention cleanup",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing snapshot that has the same resolved name",
    )
    return parser.parse_args()


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


def build_index(records: Iterable[dict], persist_dir: Path) -> int:
    persist_dir.mkdir(parents=True, exist_ok=True)
    embedding_fn = SentenceTransformerEmbeddingFunction(model_name="all-MiniLM-L6-v2")
    client = chromadb.PersistentClient(path=str(persist_dir))

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
    print(f"Indexed {len(documents)} documents into {persist_dir}")
    return len(documents)


def sanitize_snapshot_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]", "-", name.strip())
    return cleaned or "snapshot"


def default_snapshot_name(cadence: str) -> str:
    now = utc_now()
    if cadence == "monthly":
        return f"monthly-{now.strftime('%Y-%m')}"
    return f"daily-{now.strftime('%Y-%m-%d')}"


def snapshot_dataset(persist_dir: Path, snapshot_name: str, force: bool) -> Path:
    SNAPSHOT_ROOT.mkdir(parents=True, exist_ok=True)
    snapshot_dir = SNAPSHOT_ROOT / snapshot_name
    if snapshot_dir.exists():
        if force:
            shutil.rmtree(snapshot_dir)
        else:
            raise SystemExit(
                f"Snapshot '{snapshot_name}' already exists. Use --force to overwrite or pick a new name."
            )
    shutil.copytree(persist_dir, snapshot_dir)
    return snapshot_dir


def write_manifest(snapshot_dir: Path, manifest: dict) -> None:
    manifest_path = snapshot_dir / MANIFEST_NAME
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def dataset_digest(dataset_dir: Path) -> str:
    hasher = hashlib.sha256()
    for file in sorted(dataset_dir.glob("*.json")):
        hasher.update(file.name.encode("utf-8"))
        with file.open("rb") as handle:
            while chunk := handle.read(8192):
                hasher.update(chunk)
    return hasher.hexdigest()


def current_git_commit() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (subprocess.SubprocessError, FileNotFoundError):  # pragma: no cover - env dependent
        return None
    return result.stdout.strip() or None


def cleanup_snapshots(keep_daily: int, keep_monthly: int) -> list[str]:
    removed: list[str] = []
    today = utc_now().date()

    if SNAPSHOT_ROOT.exists():
        for snapshot_dir in SNAPSHOT_ROOT.iterdir():
            if not snapshot_dir.is_dir():
                continue
            name = snapshot_dir.name
            if keep_daily >= 0 and name.startswith("daily-"):
                try:
                    parsed_date = datetime.strptime(name[len("daily-") :], "%Y-%m-%d").date()
                except ValueError:
                    parsed_date = None
                if parsed_date and (today - parsed_date).days > keep_daily:
                    shutil.rmtree(snapshot_dir)
                    removed.append(name)
                    continue
            if keep_monthly >= 0 and name.startswith("monthly-"):
                try:
                    parsed_month = datetime.strptime(name[len("monthly-") :], "%Y-%m").date()
                except ValueError:
                    parsed_month = None
                if parsed_month:
                    month_delta = (today.year - parsed_month.year) * 12 + (today.month - parsed_month.month)
                    if month_delta > keep_monthly:
                        shutil.rmtree(snapshot_dir)
                        removed.append(name)
    return removed


def main() -> None:
    args = parse_args()
    dataset_dir = args.dataset.resolve()
    record_count = build_index(iter_records(dataset_dir), PERSIST_DIR)

    if args.no_snapshot:
        return

    resolved_name = sanitize_snapshot_name(
        args.snapshot_name or default_snapshot_name(args.cadence)
    )
    snapshot_dir = snapshot_dataset(PERSIST_DIR, resolved_name, args.force)
    digest = dataset_digest(dataset_dir)
    manifest = {
        "snapshot_name": resolved_name,
        "cadence": args.cadence,
        "created_at": utc_now().isoformat().replace("+00:00", "Z"),
        "records_indexed": record_count,
        "dataset_dir": str(dataset_dir),
        "dataset_hash": digest,
        "git_commit": current_git_commit(),
    }
    write_manifest(snapshot_dir, manifest)
    print(f"Snapshot '{resolved_name}' saved under {snapshot_dir}")

    if not args.no_retention:
        removed = cleanup_snapshots(args.keep_daily, args.keep_monthly)
        if removed:
            print("Removed expired snapshots:", ", ".join(sorted(removed)))


if __name__ == "__main__":  # pragma: no cover
    main()
