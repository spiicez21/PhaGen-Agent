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
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import chromadb
    from chromadb.utils.embedding_functions import (
        SentenceTransformerEmbeddingFunction,
    )
except ImportError as exc:  # pragma: no cover - runtime guard
    raise SystemExit(
        "chromadb and sentence-transformers are required. Run 'pip install chromadb sentence-transformers'"
    ) from exc

BASE_DIR = Path(__file__).resolve().parent
DATASET_DIR = Path(__file__).resolve().parents[1] / "crawler" / "storage" / "datasets" / "default"
PERSIST_DIR = BASE_DIR / "chroma"
SNAPSHOT_ROOT = BASE_DIR / "chroma_snapshots"
EMBED_CACHE_PATH = BASE_DIR / ".embedding_cache.json"
STRUCTURE_INPUT = BASE_DIR / "data" / "normalized_smiles.jsonl"
STRUCTURE_IMAGES_DIR = BASE_DIR / "data" / "structures" / "images"
STRUCTURE_METADATA_DIR = BASE_DIR / "data" / "structures" / "metadata"
STRUCTURE_MANIFEST = BASE_DIR / "data" / "structures" / "structures.manifest.json"
COLLECTION_NAME = "phagen-agentic"
MANIFEST_NAME = "manifest.json"
EMBED_MODEL_NAME = "all-MiniLM-L6-v2"
DEDUP_SOURCE_TYPES = {"clinical", "literature"}
SOURCE_PRIORITY = {"clinical": 2, "literature": 1}


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
    parser.add_argument(
        "--cache-path",
        type=Path,
        default=EMBED_CACHE_PATH,
        help="Path to embedding cache JSON file (default: indexes/.embedding_cache.json)",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Disable embedding cache reuse (always re-embed)",
    )
    parser.add_argument(
        "--no-structures",
        action="store_true",
        help="Skip RDKit structure rendering/catalog refresh",
    )
    parser.add_argument(
        "--structure-records",
        type=Path,
        default=STRUCTURE_INPUT,
        help="Canonical SMILES JSONL feeding the structure catalog (default: indexes/data/normalized_smiles.jsonl)",
    )
    parser.add_argument(
        "--structure-output-dir",
        type=Path,
        default=STRUCTURE_IMAGES_DIR,
        help="Directory for rendered SVG assets",
    )
    parser.add_argument(
        "--structure-metadata-dir",
        type=Path,
        default=STRUCTURE_METADATA_DIR,
        help="Directory for per-structure metadata JSON files",
    )
    parser.add_argument(
        "--structure-manifest",
        type=Path,
        default=STRUCTURE_MANIFEST,
        help="Manifest JSON describing all rendered structures",
    )
    parser.add_argument(
        "--structure-width",
        type=int,
        default=420,
        help="Width of rendered SVGs (default: 420)",
    )
    parser.add_argument(
        "--structure-height",
        type=int,
        default=320,
        help="Height of rendered SVGs (default: 320)",
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


def compute_doc_hash(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


@dataclass
class PreparedRecord:
    doc_id: str
    metadata: dict
    text: str


def ensure_list(vector) -> List[float]:  # type: ignore[override]
    if isinstance(vector, list):
        return vector
    try:
        return vector.tolist()  # type: ignore[attr-defined]
    except AttributeError:
        return list(vector)


def prepare_records(records: Iterable[dict]) -> List[PreparedRecord]:
    prepared: List[PreparedRecord] = []
    for record in records:
        doc_id, metadata, text = record_to_document(record)
        if not text:
            continue
        metadata = dict(metadata)
        metadata.setdefault("source_type", record.get("source_type", metadata.get("source_type", "unknown")))
        metadata["doc_hash"] = metadata.get("doc_hash") or compute_doc_hash(text)
        prepared.append(PreparedRecord(doc_id=doc_id, metadata=metadata, text=text))
    return prepared


def apply_dedup(records: Sequence[PreparedRecord]) -> Tuple[List[PreparedRecord], dict]:
    keep_flags = [True] * len(records)
    seen: Dict[str, Tuple[int, PreparedRecord]] = {}
    skipped = 0
    replaced = 0

    for idx, rec in enumerate(records):
        source_type = rec.metadata.get("source_type", "unknown")
        if source_type not in DEDUP_SOURCE_TYPES:
            continue
        dedup_key = rec.metadata.get("doc_hash") or compute_doc_hash(rec.text)
        existing = seen.get(dedup_key)
        if existing is None:
            seen[dedup_key] = (idx, rec)
            continue

        existing_idx, existing_rec = existing
        existing_priority = SOURCE_PRIORITY.get(existing_rec.metadata.get("source_type", "unknown"), 0)
        new_priority = SOURCE_PRIORITY.get(source_type, 0)

        if new_priority > existing_priority:
            keep_flags[existing_idx] = False
            existing_rec.metadata.setdefault("dedup", {})
            existing_rec.metadata["dedup"].update(
                {
                    "duplicate_of": rec.doc_id,
                    "reason": "replaced-by-higher-priority",
                }
            )
            rec.metadata.setdefault("dedup", {})
            rec.metadata["dedup"].update(
                {
                    "preferred_over": existing_rec.doc_id,
                    "reason": "higher-priority-source",
                }
            )
            seen[dedup_key] = (idx, rec)
            replaced += 1
        else:
            keep_flags[idx] = False
            rec.metadata.setdefault("dedup", {})
            rec.metadata["dedup"].update(
                {
                    "duplicate_of": existing_rec.doc_id,
                    "reason": "lower-priority-source",
                }
            )
            skipped += 1

    filtered = [rec for rec, keep in zip(records, keep_flags) if keep]
    stats = {
        "evaluated": len(records),
        "kept": len(filtered),
        "duplicates_skipped": skipped,
        "duplicates_replaced": replaced,
        "source_types": sorted(DEDUP_SOURCE_TYPES),
    }
    return filtered, stats


class EmbeddingCache:
    def __init__(self, path: Path, model_name: str) -> None:
        self.path = path
        self.model_name = model_name
        self.entries: Dict[str, Dict[str, List[float]]] = {}
        self.loaded = False
        self.reused = 0
        self.created = 0
        self._load()

    def _load(self) -> None:
        if not self.path.exists():
            return
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return
        if data.get("model") != self.model_name:
            return
        entries = data.get("entries") or {}
        self.entries = {
            key: {"embedding": list(value.get("embedding", []))}
            for key, value in entries.items()
            if isinstance(value, dict) and "embedding" in value
        }
        self.loaded = True

    def lookup(self, text: str, key: Optional[str] = None) -> Tuple[str, Optional[List[float]]]:
        key = key or compute_doc_hash(text)
        entry = self.entries.get(key)
        if entry:
            self.reused += 1
            return key, list(entry["embedding"])
        return key, None

    def store(self, key: str, embedding: List[float]) -> None:
        self.entries[key] = {"embedding": list(embedding)}
        self.created += 1

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "model": self.model_name,
            "entries": self.entries,
            "updated_at": utc_now().isoformat().replace("+00:00", "Z"),
        }
        self.path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")


def maybe_init_cache(args: argparse.Namespace) -> Optional[EmbeddingCache]:
    if args.no_cache:
        return None
    return EmbeddingCache(args.cache_path, EMBED_MODEL_NAME)


def build_index(
    records: Sequence[PreparedRecord],
    persist_dir: Path,
    embedding_fn: SentenceTransformerEmbeddingFunction,
    cache: Optional[EmbeddingCache] = None,
) -> tuple[int, int, int]:
    persist_dir.mkdir(parents=True, exist_ok=True)
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
    embeddings: List[List[float]] = []
    for record in records:
        doc_id = record.doc_id
        metadata = dict(record.metadata)
        text = record.text

        cached_key = metadata.get("doc_hash")
        cached_embedding: Optional[List[float]] = None
        if cache:
            cached_key, cached_embedding = cache.lookup(text, cached_key)

        if cached_embedding is None:
            vector = embedding_fn([text])[0]
            embedding_values = ensure_list(vector)
            if cache:
                cache.store(cached_key or compute_doc_hash(text), embedding_values)
        else:
            embedding_values = cached_embedding

        documents.append(text)
        metadata["doc_hash"] = cached_key or metadata.get("doc_hash") or compute_doc_hash(text)
        metadatas.append(metadata)
        ids.append(doc_id)
        embeddings.append(embedding_values)

    if not documents:
        raise SystemExit("No documents found to index.")

    collection.add(
        documents=documents,
        metadatas=metadatas,
        ids=ids,
        embeddings=embeddings,
    )
    reused = cache.reused if cache else 0
    created = cache.created if cache else len(documents)
    print(
        f"Indexed {len(documents)} documents into {persist_dir} (embeddings reused: {reused}, newly encoded: {created})"
    )
    return len(documents), reused, created


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


def maybe_render_structures(args: argparse.Namespace) -> dict | None:
    if args.no_structures:
        return None

    records_path = (args.structure_records or STRUCTURE_INPUT).resolve()
    if not records_path.exists():
        print(f"Structure records file {records_path} not found; skipping structure catalog refresh.")
        return None

    try:
        from structure_renderer import StructureRenderConfig, render_structure_catalog
    except Exception as exc:  # noqa: BLE001
        print(f"Skipping structure rendering because structure_renderer could not load: {exc}")
        return None

    config = StructureRenderConfig(
        input_path=records_path,
        output_dir=(args.structure_output_dir or STRUCTURE_IMAGES_DIR).resolve(),
        metadata_dir=(args.structure_metadata_dir or STRUCTURE_METADATA_DIR).resolve(),
        manifest_path=(args.structure_manifest or STRUCTURE_MANIFEST).resolve(),
        width=args.structure_width,
        height=args.structure_height,
    )
    try:
        summary = render_structure_catalog(config)
    except Exception as exc:  # noqa: BLE001
        print(f"Structure rendering failed: {exc}")
        return None

    print(
        f"Rendered {summary['rendered']} structure(s) into {config.output_dir} (manifest {summary['manifest']})."
    )
    return summary


def main() -> None:
    args = parse_args()
    dataset_dir = args.dataset.resolve()
    embedding_fn = SentenceTransformerEmbeddingFunction(model_name=EMBED_MODEL_NAME)
    cache = maybe_init_cache(args)
    prepared_records = prepare_records(iter_records(dataset_dir))
    deduped_records, dedup_stats = apply_dedup(prepared_records)
    record_count, reused, created = build_index(
        deduped_records,
        PERSIST_DIR,
        embedding_fn,
        cache=cache,
    )

    if cache and not args.no_cache:
        cache.save()

    structure_summary = maybe_render_structures(args)

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
        "embeddings_reused": reused,
        "embeddings_encoded": created,
        "dedup": dedup_stats,
        "dataset_dir": str(dataset_dir),
        "dataset_hash": digest,
        "git_commit": current_git_commit(),
        "embedding_model": EMBED_MODEL_NAME,
    }
    if structure_summary:
        manifest["structures"] = structure_summary
    write_manifest(snapshot_dir, manifest)
    print(f"Snapshot '{resolved_name}' saved under {snapshot_dir}")

    if not args.no_retention:
        removed = cleanup_snapshots(args.keep_daily, args.keep_monthly)
        if removed:
            print("Removed expired snapshots:", ", ".join(sorted(removed)))


if __name__ == "__main__":  # pragma: no cover
    main()
