from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict

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
STRUCTURE_ASSET_DIR = REPORT_ASSET_ROOT / "structures"


@dataclass
class StructureRenderResult:
    svg: str
    path: Path


def render_structure_svg(
    smiles: str,
    *,
    output_dir: Path | None = None,
    filename_hint: str,
    width: int = 420,
    height: int = 320,
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

    return StructureRenderResult(svg=svg, path=asset_path)


def build_structure_payload(
    *,
    smiles: str,
    job_id: str,
    molecule_label: str,
    output_dir: Path | None = None,
) -> Dict[str, str]:
    filename_hint = f"{molecule_label}-{job_id[:8]}"
    try:
        result = render_structure_svg(
            smiles,
            output_dir=output_dir,
            filename_hint=filename_hint,
        )
    except Exception as exc:  # noqa: BLE001
        return {
            "svg": "",
            "path": "",
            "error": str(exc),
            "smiles": smiles,
        }

    return {
        "svg": result.svg,
        "path": str(result.path),
        "error": "",
        "smiles": smiles,
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
]
