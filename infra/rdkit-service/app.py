from __future__ import annotations

import base64
import io
import os
from typing import Literal

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
except ImportError as exc:  # pragma: no cover - runtime guard
    raise RuntimeError("RDKit must be installed in the rdkit-service container") from exc

app = FastAPI(title="rdkit-service", version=os.getenv("RDKIT_SERVICE_VERSION", "0.1.0"))


class RenderPayload(BaseModel):
    smiles: str = Field(..., description="Input SMILES string")
    format: Literal["svg", "png"] = Field("svg", description="Image output format")
    size: int = Field(400, ge=64, le=1024, description="Square image size in pixels")


def _mol_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    Chem.SanitizeMol(mol)
    return mol


@app.get("/health")
def healthcheck() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/render")
def render(payload: RenderPayload) -> dict[str, str]:
    mol = _mol_from_smiles(payload.smiles)
    size = (payload.size, payload.size)

    if payload.format == "svg":
        svg = Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=size, useSVG=True)
        return {"format": "svg", "data": str(svg)}

    image = Draw.MolToImage(mol, size=size)
    buffer = io.BytesIO()
    image.save(buffer, format="PNG")
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return {"format": "png", "encoding": "base64", "data": encoded}
