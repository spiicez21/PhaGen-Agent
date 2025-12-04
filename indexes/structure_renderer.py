"""Shared helpers for rendering RDKit structure assets for the index pipeline."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

try:  # pragma: no cover - heavy optional dependency
    from rdkit import Chem
    from rdkit.Chem import AllChem, inchi
    from rdkit.Chem.Draw import rdMolDraw2D
except Exception as exc:  # noqa: BLE001
    Chem = None  # type: ignore[assignment]
    AllChem = None  # type: ignore[assignment]
    rdMolDraw2D = None  # type: ignore[assignment]
    inchi = None  # type: ignore[assignment]
    _RDKIT_IMPORT_ERROR: Exception | None = exc
else:
    _RDKIT_IMPORT_ERROR = None


@dataclass
class StructureRenderConfig:
    input_path: Path
    output_dir: Path
    metadata_dir: Path
    manifest_path: Path
    width: int = 420
    height: int = 320


@dataclass
class StructureAsset:
    image_id: str
    svg_path: Path
    metadata_path: Path
    canonical_smiles: str
    canonical_inchi: str
    inchikey: str
    name: str
    generated_at: str
    source_type: str
    source_ref: str
    license_info: str


def render_structure_catalog(config: StructureRenderConfig) -> Dict[str, object]:
    """Render SVG + metadata assets for every SMILES entry in the input JSONL."""

    if Chem is None or AllChem is None or rdMolDraw2D is None or inchi is None:
        raise RuntimeError("RDKit is required to render structures") from _RDKIT_IMPORT_ERROR

    records = list(_load_records(config.input_path))
    config.output_dir.mkdir(parents=True, exist_ok=True)
    config.metadata_dir.mkdir(parents=True, exist_ok=True)
    assets: List[StructureAsset] = []
    failures: List[Dict[str, str]] = []
    used_ids: Dict[str, int] = {}

    for record in records:
        try:
            asset = _render_record(record, config, used_ids)
        except Exception as exc:  # noqa: BLE001
            failures.append({
                "name": str(record.get("name") or record.get("input", {}).get("smiles") or "unknown"),
                "error": str(exc),
            })
            continue
        if asset:
            assets.append(asset)

    manifest_payload = _build_manifest(config, assets, failures, len(records))
    config.manifest_path.parent.mkdir(parents=True, exist_ok=True)
    config.manifest_path.write_text(
        json.dumps(manifest_payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    return {
        "processed": len(records),
        "rendered": len(assets),
        "failures": len(failures),
        "manifest": str(config.manifest_path),
    }


def _load_records(path: Path) -> Iterable[dict]:
    if not path.exists():
        raise SystemExit(f"Structure input file {path} does not exist")
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text:
                continue
            try:
                yield json.loads(text)
            except json.JSONDecodeError as exc:  # pragma: no cover - CLI guard
                raise SystemExit(f"Invalid JSON at line {line_number}: {exc}") from exc


def _render_record(record: dict, config: StructureRenderConfig, used_ids: Dict[str, int]) -> StructureAsset | None:
    smiles = (
        record.get("canonical_smiles")
        or record.get("canonicalSmiles")
        or record.get("smiles")
        or record.get("input", {}).get("smiles")
    )
    if not smiles:
        raise ValueError("Record missing canonical_smiles")

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("RDKit could not parse SMILES")

    AllChem.Compute2DCoords(molecule)
    drawer = rdMolDraw2D.MolDraw2DSVG(config.width, config.height)
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
    canonical_inchi = record.get("canonical_inchi") or inchi.MolToInchi(molecule)
    inchikey = record.get("inchikey") or inchi.MolToInchiKey(molecule)

    image_id = _derive_image_id(record, canonical_smiles, inchikey, used_ids)
    svg_path = config.output_dir / f"{image_id}.svg"
    metadata_path = config.metadata_dir / f"{image_id}.json"

    svg_path.write_text(svg, encoding="utf-8")
    generated_at = _utcnow()
    source_type = record.get("source_type") or record.get("source") or "rdkit"
    source_ref = record.get("source_ref") or record.get("source_reference") or canonical_smiles
    license_info = record.get("license") or "internal-use-only"
    metadata = {
        "image_id": image_id,
        "asset_path": str(svg_path),
        "source_type": source_type,
        "source_reference": source_ref,
        "smiles": canonical_smiles,
        "inchi": canonical_inchi or "",
        "inchikey": inchikey or "",
        "name": record.get("name", ""),
        "synonyms": record.get("synonyms", []),
        "generated_at": generated_at,
        "license": license_info,
        "origin": "indexes.render_structures",
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    return StructureAsset(
        image_id=image_id,
        svg_path=svg_path,
        metadata_path=metadata_path,
        canonical_smiles=canonical_smiles,
        canonical_inchi=canonical_inchi or "",
        inchikey=inchikey or "",
        name=record.get("name", ""),
        generated_at=generated_at,
        source_type=source_type,
        source_ref=source_ref,
        license_info=license_info,
    )


def _derive_image_id(record: dict, smiles: str, inchikey: str | None, used_ids: Dict[str, int]) -> str:
    bases = [record.get("name"), inchikey, smiles]
    for base in bases:
        candidate = _sanitize_identifier(str(base or ""))
        if candidate:
            break
    else:
        candidate = "structure"

    count = used_ids.get(candidate, 0)
    used_ids[candidate] = count + 1
    if count:
        candidate = f"{candidate}-{count+1}"
    return candidate


def _sanitize_identifier(value: str) -> str:
    allowed = [ch.lower() for ch in value if ch.isalnum() or ch in {"-", "_"}]
    if not allowed:
        return ""
    return "".join(allowed)[:64]


def _build_manifest(
    config: StructureRenderConfig,
    assets: List[StructureAsset],
    failures: List[Dict[str, str]],
    processed: int,
) -> Dict[str, object]:
    manifest_dir = config.manifest_path.parent
    entries: List[Dict[str, object]] = []
    for asset in assets:
        entries.append(
            {
                "image_id": asset.image_id,
                "name": asset.name,
                "canonical_smiles": asset.canonical_smiles,
                "canonical_inchi": asset.canonical_inchi,
                "inchikey": asset.inchikey,
                "svg_path": _relativize(asset.svg_path, manifest_dir),
                "metadata_path": _relativize(asset.metadata_path, manifest_dir),
                "generated_at": asset.generated_at,
                "source_type": asset.source_type,
                "source_ref": asset.source_ref,
                "license": asset.license_info,
            }
        )

    return {
        "generated_at": _utcnow(),
        "input": str(config.input_path),
        "output_dir": str(config.output_dir),
        "metadata_dir": str(config.metadata_dir),
        "entries": entries,
        "failures": failures,
        "processed": processed,
        "rendered": len(assets),
    }


def _relativize(path: Path, base: Path) -> str:
    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


__all__ = ["StructureRenderConfig", "render_structure_catalog"]
