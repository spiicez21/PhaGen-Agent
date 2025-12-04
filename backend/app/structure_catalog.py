from __future__ import annotations

import json
import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST_PATH = REPO_ROOT / "indexes" / "data" / "structures" / "structures.manifest.json"
ENV_CATALOG_PATH = "STRUCTURE_CATALOG_PATH"


@dataclass(frozen=True)
class StructureCatalogEntry:
    image_id: str
    svg_path: Path
    metadata_path: Path
    canonical_smiles: str | None
    canonical_inchi: str | None
    inchikey: str | None
    source_type: str | None
    source_ref: str | None
    license: str | None
    generated_at: str | None


class StructureCatalog:
    def __init__(self, manifest_path: Path) -> None:
        self.manifest_path = manifest_path
        self._entries_by_smiles: Dict[str, StructureCatalogEntry] = {}
        self._entries_by_inchikey: Dict[str, StructureCatalogEntry] = {}
        self._entries_by_image: Dict[str, StructureCatalogEntry] = {}
        self._last_mtime: float | None = None
        self.generated_at: str | None = None

    def lookup(
        self,
        *,
        smiles: str | None = None,
        inchikey: str | None = None,
        image_id: str | None = None,
    ) -> StructureCatalogEntry | None:
        self._refresh_if_needed()
        if smiles:
            entry = self._entries_by_smiles.get(smiles.strip())
            if entry:
                return entry
        if inchikey:
            entry = self._entries_by_inchikey.get(inchikey.strip().upper())
            if entry:
                return entry
        if image_id:
            entry = self._entries_by_image.get(image_id.strip())
            if entry:
                return entry
        return None

    # Internal helpers -------------------------------------------------
    def _refresh_if_needed(self) -> None:
        try:
            stat = self.manifest_path.stat()
        except FileNotFoundError:
            self._entries_by_smiles.clear()
            self._entries_by_inchikey.clear()
            self._entries_by_image.clear()
            self.generated_at = None
            self._last_mtime = None
            return

        if self._last_mtime and stat.st_mtime == self._last_mtime:
            return

        self._load_manifest()
        self._last_mtime = stat.st_mtime

    def _load_manifest(self) -> None:
        with self.manifest_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)

        base_dir = self.manifest_path.parent
        self._entries_by_smiles.clear()
        self._entries_by_inchikey.clear()
        self._entries_by_image.clear()
        for raw_entry in payload.get("entries", []):
            self._register_entry(base_dir, raw_entry)
        self.generated_at = payload.get("generated_at")

    def _register_entry(self, base_dir: Path, raw_entry: dict) -> None:
        svg_path = self._resolve_path(base_dir, raw_entry.get("svg_path"))
        metadata_path = self._resolve_path(base_dir, raw_entry.get("metadata_path"))
        inchikey = raw_entry.get("inchikey") or None
        if inchikey:
            inchikey = inchikey.strip().upper()
        entry = StructureCatalogEntry(
            image_id=raw_entry.get("image_id", ""),
            svg_path=svg_path,
            metadata_path=metadata_path,
            canonical_smiles=self._clean(raw_entry.get("canonical_smiles")),
            canonical_inchi=self._clean(raw_entry.get("canonical_inchi")),
            inchikey=inchikey,
            source_type=self._clean(raw_entry.get("source_type")),
            source_ref=self._clean(raw_entry.get("source_ref")),
            license=self._clean(raw_entry.get("license")),
            generated_at=self._clean(raw_entry.get("generated_at")),
        )

        if entry.canonical_smiles:
            self._entries_by_smiles[entry.canonical_smiles] = entry
        if entry.inchikey:
            self._entries_by_inchikey[entry.inchikey] = entry
        if entry.image_id:
            self._entries_by_image[entry.image_id] = entry

    @staticmethod
    def _resolve_path(base_dir: Path, raw_path: Optional[str]) -> Path:
        if not raw_path:
            return base_dir
        candidate = Path(raw_path)
        if not candidate.is_absolute():
            candidate = base_dir / candidate
        return candidate.resolve()

    @staticmethod
    def _clean(value: Optional[str]) -> Optional[str]:
        if value is None:
            return None
        text = str(value).strip()
        return text or None


def _resolve_catalog_path() -> Path:
    override = os.getenv(ENV_CATALOG_PATH)
    if override:
        return Path(override).expanduser().resolve()
    return DEFAULT_MANIFEST_PATH


@lru_cache(maxsize=None)
def _catalog_for_path(manifest_path: str) -> StructureCatalog:
    return StructureCatalog(Path(manifest_path))


def get_structure_catalog() -> StructureCatalog:
    path = _resolve_catalog_path()
    return _catalog_for_path(str(path))


__all__ = ["StructureCatalog", "StructureCatalogEntry", "get_structure_catalog"]
