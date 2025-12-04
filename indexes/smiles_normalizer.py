"""Canonicalize SMILES/InChI records ahead of indexing and synonym expansion.

Usage:
    python indexes/smiles_normalizer.py \
        --input indexes/data/sample_smiles.jsonl \
        --output indexes/data/normalized_smiles.jsonl

The input file must contain newline-delimited JSON objects with the following
fields:
    {
        "name": "Pirfenidone",
        "smiles": "O=C1NC(=O)N(C)C(=O)N1C",
        "inchi": "optional InChI string",
        "synonyms": ["alias-1", "alias-2"]
    }

The script emits a matching JSONL file with canonical SMILES, InChI, and
InChIKey plus the deduplicated synonym set and provenance metadata.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional

from rdkit import Chem
from rdkit.Chem import inchi

DEFAULT_INPUT = Path(__file__).resolve().parent / "data" / "sample_smiles.jsonl"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "data" / "normalized_smiles.jsonl"
DEFAULT_MANIFEST = Path(__file__).resolve().parent / "data" / "normalized_smiles.manifest.json"


@dataclass
class NormalizationResult:
    payload: dict
    error: Optional[str] = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Canonicalize SMILES/InChI via RDKit")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Path to JSONL file with raw SMILES records")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Destination JSONL file with canonicalized entries",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="Optional path to write normalization summary metadata",
    )
    parser.add_argument(
        "--append",
        action="store_true",
        help="Append to the output file instead of overwriting it",
    )
    return parser.parse_args()


def load_records(path: Path) -> Iterable[dict]:
    if not path.exists():
        raise SystemExit(f"Input file {path} does not exist")
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                yield json.loads(line)
            except json.JSONDecodeError as exc:
                raise SystemExit(f"Invalid JSON at line {line_number}: {exc}") from exc


def to_mol(record: dict):
    smiles = record.get("smiles") or ""
    inchi_value = record.get("inchi") or ""
    mol = None
    source = None

    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        source = "smiles"
    if mol is None and inchi_value:
        mol = Chem.MolFromInchi(inchi_value)
        source = "inchi"

    if mol is None:
        raise ValueError("Unable to parse SMILES or InChI")

    Chem.SanitizeMol(mol)
    return mol, source


def normalize_record(record: dict) -> NormalizationResult:
    mol, source = to_mol(record)
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    canonical_inchi = inchi.MolToInchi(mol)
    inchikey_value = inchi.MolToInchiKey(mol)

    synonyms: List[str] = []
    for value in record.get("synonyms", []) or []:
        if value:
            synonyms.append(value.strip())
    for alias in (canonical_smiles, canonical_inchi, inchikey_value):
        if alias:
            synonyms.append(alias)
    if record.get("name"):
        synonyms.append(record["name"].strip())

    deduped_synonyms = sorted({alias for alias in synonyms if alias})

    payload = {
        "name": record.get("name", ""),
        "source": source,
        "input": {
            "smiles": record.get("smiles"),
            "inchi": record.get("inchi"),
        },
        "canonical_smiles": canonical_smiles,
        "canonical_inchi": canonical_inchi,
        "inchikey": inchikey_value,
        "synonyms": deduped_synonyms,
    }
    return NormalizationResult(payload=payload)


def write_results(path: Path, results: List[NormalizationResult], append: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if append else "w"
    with path.open(mode, encoding="utf-8") as handle:
        for result in results:
            handle.write(json.dumps(result.payload, ensure_ascii=False) + "\n")


def write_manifest(path: Path, processed: int, failures: List[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "processed": processed,
        "normalized": processed - len(failures),
        "failed": len(failures),
        "failures": failures,
    }
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    normalized: List[NormalizationResult] = []
    failures: List[dict] = []

    for record in load_records(args.input):
        try:
            normalized.append(normalize_record(record))
        except Exception as exc:  # pragma: no cover - simple CLI guard
            failures.append({"record": record, "error": str(exc)})

    if failures:
        print(f"⚠️  {len(failures)} record(s) failed normalization. See manifest for details.")
    else:
        print("✅ All records normalized successfully.")

    write_results(args.output, normalized, append=args.append)
    write_manifest(args.manifest, processed=len(normalized) + len(failures), failures=failures)


if __name__ == "__main__":  # pragma: no cover
    main()
