from __future__ import annotations

import os
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Literal

try:  # pragma: no cover - optional heavy dependency
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
except Exception as exc:  # noqa: BLE001
    Chem = None  # type: ignore[assignment]
    AllChem = None  # type: ignore[assignment]
    rdMolDraw2D = None  # type: ignore[assignment]
    _RDKIT_IMPORT_ERROR: Exception | None = exc
else:
    _RDKIT_IMPORT_ERROR = None


_DEFAULT_ASSET_ROOT = Path(__file__).resolve().parent / "report_assets"
REPORT_ASSET_ROOT = Path(os.getenv("REPORT_ASSETS_DIR", str(_DEFAULT_ASSET_ROOT)))
REPORT_IMAGES_ROOT = REPORT_ASSET_ROOT / "reports" / "images"
STRUCTURE_ASSET_DIR = REPORT_IMAGES_ROOT / "structures"
STRUCTURE_METADATA_DIR = REPORT_IMAGES_ROOT / "metadata"


@dataclass
class StructureRenderResult:
    svg: str
    path: Path
    metadata_path: Path
    metadata: Dict[str, str]


def render_structure_svg(
    smiles: str,
    *,
    output_dir: Path | None = None,
    metadata_dir: Path | None = None,
    filename_hint: str,
    width: int = 420,
    height: int = 320,
    source_type: Literal["smiles", "inchi", "pubchem"] = "smiles",
    source_reference: str | None = None,
    inchi: str | None = None,
) -> StructureRenderResult:
    if not smiles:
        raise ValueError("SMILES string is required to render a structure.")
    if Chem is None or rdMolDraw2D is None or AllChem is None:
        raise RuntimeError(
            "RDKit is not available. Install rdkit-pypi to enable structure rendering."
        ) from _RDKIT_IMPORT_ERROR

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("RDKit could not parse the provided SMILES string.")

    AllChem.Compute2DCoords(molecule)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    target_dir = output_dir or STRUCTURE_ASSET_DIR
    target_dir.mkdir(parents=True, exist_ok=True)
    safe_stem = _sanitize_filename(filename_hint)
    asset_path = target_dir / f"{safe_stem}.svg"
    asset_path.write_text(svg, encoding="utf-8")
    meta_dir = metadata_dir or STRUCTURE_METADATA_DIR
    meta_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = meta_dir / f"{safe_stem}.json"
    generated_at = datetime.utcnow().replace(tzinfo=timezone.utc).isoformat()
    source_ref = source_reference or smiles
    metadata = {
        "image_id": safe_stem,
        "asset_path": str(asset_path),
        "source_type": source_type,
        "source_reference": source_ref,
        "smiles": smiles,
        "inchi": inchi or "",
        "generated_at": generated_at,
    }
    metadata_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    return StructureRenderResult(
        svg=svg,
        path=asset_path,
        metadata_path=metadata_path,
        metadata=metadata,
    )


def build_structure_payload(
    *,
    smiles: str,
    job_id: str,
    molecule_label: str,
    output_dir: Path | None = None,
    metadata_dir: Path | None = None,
    inchi: str | None = None,
    source_type: Literal["smiles", "inchi", "pubchem"] = "smiles",
    source_reference: str | None = None,
) -> Dict[str, str]:
    filename_hint = f"{molecule_label}-{job_id[:8]}"
    fallback_image_id = _sanitize_filename(filename_hint)
    try:
        result = render_structure_svg(
            smiles,
            output_dir=output_dir,
            metadata_dir=metadata_dir,
            filename_hint=filename_hint,
            inchi=inchi,
            source_type=source_type,
            source_reference=source_reference,
        )
    except Exception as exc:  # noqa: BLE001
        return {
            "svg": "",
            "path": "",
            "error": str(exc),
            "smiles": smiles,
            "metadata_path": "",
            "source_type": source_type,
            "source_reference": source_reference or smiles,
            "inchi": inchi or "",
            "generated_at": "",
            "image_id": fallback_image_id,
        }

    return {
        "svg": result.svg,
        "path": str(result.path),
        "error": "",
        "smiles": smiles,
        "metadata_path": str(result.metadata_path),
        "source_type": result.metadata.get("source_type", source_type),
        "source_reference": result.metadata.get("source_reference", source_reference or smiles),
        "inchi": inchi or "",
        "generated_at": result.metadata.get("generated_at", ""),
        "image_id": result.metadata.get("image_id", fallback_image_id),
    }


def _sanitize_filename(filename: str) -> str:
    cleaned = [ch for ch in filename if ch.isalnum() or ch in ("-", "_")]
    if not cleaned:
        return "structure"
    return "".join(cleaned)


__all__ = [
    "StructureRenderResult",
    "render_structure_svg",
    "build_structure_payload",
    "STRUCTURE_ASSET_DIR",
    "STRUCTURE_METADATA_DIR",
]
