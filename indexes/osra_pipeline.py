"""Run OSRA on approved structure diagrams and emit SMILES records.

Usage:
    python indexes/osra_pipeline.py \
        --input indexes/data/sample_diagrams.jsonl \
        --output indexes/data/osra_results.jsonl

Each input line must be JSON with at least `image_path` and `allow_osra`. Records are
skipped unless `allow_osra` is true so sites with restrictive terms of service are
never processed accidentally.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PureWindowsPath
from typing import Iterable, List, Optional

from rdkit import Chem
from rdkit.Chem import inchi

DEFAULT_BASE_DIR = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = Path(__file__).resolve().parent / "data" / "sample_diagrams.jsonl"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "data" / "osra_results.jsonl"
DEFAULT_MANIFEST = Path(__file__).resolve().parent / "data" / "manifests" / "osra_results.manifest.json"
STDOUT_LIMIT = 2000
STDERR_LIMIT = 2000


@dataclass
class DiagramRecord:
    raw: dict
    image_path: Path


@dataclass
class OsraResult:
    payload: dict
    error: Optional[str] = None


class OsraNotAvailable(RuntimeError):
    """Raised when the OSRA binary is missing."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert diagrams to SMILES via OSRA")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="JSONL file with diagram metadata")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Destination JSONL for OSRA results")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST, help="Summary report JSON")
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE_DIR, help="Resolve relative image paths from here")
    parser.add_argument("--osra-binary", default="osra", help="Path to OSRA executable (default: osra in PATH)")
    parser.add_argument(
        "--extra-osra-args",
        nargs="*",
        default=[],
        help="Additional flags to pass to OSRA before the image path (e.g., --resolution 300)",
    )
    parser.add_argument("--timeout", type=int, default=45, help="Per-image timeout in seconds")
    parser.add_argument(
        "--skip-missing",
        action="store_true",
        help="Skip records whose image files are missing instead of treating them as failures",
    )
    parser.add_argument(
        "--image-path-mode",
        choices=["windows", "wsl"],
        default="windows",
        help="Format image paths for native Windows binaries (default) or WSL-mounted paths (/mnt/d/...)",
    )
    return parser.parse_args()


def load_records(path: Path, base_dir: Path) -> Iterable[DiagramRecord]:
    if not path.exists():
        raise SystemExit(f"Input file {path} does not exist")
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                raw = json.loads(line)
            except json.JSONDecodeError as exc:  # pragma: no cover - CLI guard
                raise SystemExit(f"Invalid JSON at line {line_number}: {exc}") from exc

            image_path_value = raw.get("image_path")
            if not image_path_value:
                raise SystemExit(f"Record at line {line_number} is missing image_path")

            image_path = Path(image_path_value)
            if not image_path.is_absolute():
                image_path = (base_dir / image_path).resolve()

            yield DiagramRecord(raw=raw, image_path=image_path)


def ensure_osra_available(binary: str) -> str:
    candidate_path = Path(binary)
    if candidate_path.exists():
        return str(candidate_path)

    resolved = shutil.which(binary)
    if resolved:
        return resolved

    raise OsraNotAvailable("OSRA binary not found. Install OSRA and/or set --osra-binary to its location.")


def compute_sha256(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def truncate(value: str, limit: int) -> str:
    if len(value) <= limit:
        return value
    return value[: limit - 3] + "..."


def extract_smiles(stdout: str) -> str:
    tokens: List[str] = []
    for line in stdout.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        for piece in stripped.replace("\t", " ").split():
            tokens.append(piece.strip())

    for token in tokens:
        mol = Chem.MolFromSmiles(token)
        if mol:
            try:
                Chem.SanitizeMol(mol)
            except Exception:  # pragma: no cover - RDKit sanitization rarely fails for valid SMILES
                continue
            return Chem.MolToSmiles(mol, canonical=True)
    raise ValueError("OSRA output did not contain a valid SMILES string")


def format_image_path_for_osra(path: Path, mode: str) -> str:
    if mode == "windows":
        return str(path)

    win_path = PureWindowsPath(path)
    drive = win_path.drive.rstrip(":")
    if not drive:
        raise ValueError("WSL path mode requires an absolute Windows path with a drive letter")
    relative_parts = win_path.parts[1:]
    relative = "/".join(part.strip("\\/") for part in relative_parts if part not in {"\\", "/"})
    relative = relative.replace("\\", "/")
    return f"/mnt/{drive.lower()}/{relative}"


def process_record(record: DiagramRecord, args: argparse.Namespace) -> OsraResult:
    metadata = record.raw
    if not metadata.get("allow_osra"):
        return OsraResult(
            payload={
                "diagram_id": metadata.get("diagram_id"),
                "image_path": str(record.image_path),
                "skipped": True,
                "reason": "allow_osra flag not set",
            }
        )

    if not record.image_path.exists():
        if args.skip_missing:
            return OsraResult(
                payload={
                    "diagram_id": metadata.get("diagram_id"),
                    "image_path": str(record.image_path),
                    "skipped": True,
                    "reason": "image missing",
                }
            )
        raise FileNotFoundError(f"Image file not found: {record.image_path}")

    osra_image_path = format_image_path_for_osra(record.image_path, args.image_path_mode)
    command = [args.osra_binary, "-o", "smi", *args.extra_osra_args, osra_image_path]
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=args.timeout,
            check=False,
        )
    except FileNotFoundError as exc:  # pragma: no cover - guard for race conditions
        raise OsraNotAvailable("OSRA binary not found during execution") from exc

    stdout = completed.stdout.strip()
    stderr = completed.stderr.strip()
    if completed.returncode != 0:
        raise RuntimeError(f"OSRA exited with {completed.returncode}: {stderr or stdout}")

    canonical_smiles = extract_smiles(stdout)
    mol = Chem.MolFromSmiles(canonical_smiles)
    canonical_inchi = inchi.MolToInchi(mol) if mol else None
    inchikey_value = inchi.MolToInchiKey(mol) if mol else None

    synonyms = [value for value in [metadata.get("molecule_hint"), canonical_smiles, canonical_inchi, inchikey_value] if value]
    deduped_synonyms = sorted({alias for alias in synonyms if alias})

    payload = {
        "diagram_id": metadata.get("diagram_id") or record.image_path.stem,
        "image_path": str(record.image_path),
        "osra_image_path": osra_image_path,
        "image_sha256": compute_sha256(record.image_path),
        "source_url": metadata.get("source_url"),
        "molecule_hint": metadata.get("molecule_hint"),
        "notes": metadata.get("notes"),
        "allow_osra": True,
        "osra": {
            "command": command,
            "stdout": truncate(stdout, STDOUT_LIMIT),
            "stderr": truncate(stderr, STDERR_LIMIT),
            "exit_code": completed.returncode,
        },
        "smiles": canonical_smiles,
        "inchi": canonical_inchi,
        "canonical_smiles": canonical_smiles,
        "canonical_inchi": canonical_inchi,
        "inchikey": inchikey_value,
        "synonyms": deduped_synonyms,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    return OsraResult(payload=payload)


def write_results(path: Path, results: List[OsraResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for result in results:
            handle.write(json.dumps(result.payload, ensure_ascii=False) + "\n")


def write_manifest(
    path: Path,
    processed: int,
    converted: int,
    skipped: int,
    skipped_records: List[dict],
    failures: List[dict],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "processed": processed,
        "converted": converted,
        "skipped": skipped,
        "skipped_records": skipped_records,
        "failed": len(failures),
        "failures": failures,
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    resolved_binary = ensure_osra_available(args.osra_binary)
    args.osra_binary = resolved_binary

    results: List[OsraResult] = []
    failures: List[dict] = []
    skipped = 0
    skipped_records: List[dict] = []

    for record in load_records(args.input, args.base_dir):
        try:
            outcome = process_record(record, args)
        except Exception as exc:  # pragma: no cover - CLI guard
            failures.append({"diagram_id": record.raw.get("diagram_id"), "image_path": str(record.image_path), "error": str(exc)})
            continue

        if outcome.payload.get("skipped"):
            skipped += 1
            skipped_records.append({
                "diagram_id": outcome.payload.get("diagram_id"),
                "image_path": outcome.payload.get("image_path"),
                "reason": outcome.payload.get("reason"),
            })
        else:
            results.append(outcome)

    write_results(args.output, results)
    write_manifest(
        args.manifest,
        processed=len(results) + skipped + len(failures),
        converted=len(results),
        skipped=skipped,
        skipped_records=skipped_records,
        failures=failures,
    )

    if failures:
        print(f"⚠️  {len(failures)} diagram(s) failed — see {args.manifest} for details.")
    else:
        print(f"✅ Converted {len(results)} diagram(s); {skipped} skipped per policy.")


if __name__ == "__main__":  # pragma: no cover
    main()
