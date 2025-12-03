from __future__ import annotations

import json
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

_DEFAULT_DATA = [
    {
        "canonical": "pirfenidone",
        "smiles": "O=C1NC(=O)N(C)C(C)=C1C1=CC=CC=C1",
        "synonyms": [
            "Esbriet",
            "Pirespa",
            "5-methyl-1-phenyl-2(1H)-pyridone",
            "PFD",
        ],
    },
    {
        "canonical": "nintedanib",
        "smiles": "CC1=CC(=CC=C1OC)NC(=O)C2=NC(=NC(=N2)N)N(C)C",
        "synonyms": [
            "Ofev",
            "BIBF 1120",
            "Vargatef",
        ],
    },
]


class SynonymExpander:
    """Lightweight synonym expander using heuristics + a small built-in catalog."""

    def __init__(
        self,
        extra_catalog: Optional[Iterable[Dict[str, str | List[str]]]] = None,
    ) -> None:
        self.catalog = list(_DEFAULT_DATA)
        if extra_catalog:
            self.catalog.extend(list(extra_catalog))
        self._name_index = {
            self._key(entry["canonical"]): entry for entry in self.catalog
        }
        self._smiles_index = {
            entry.get("smiles", "").strip(): entry
            for entry in self.catalog
            if entry.get("smiles")
        }

    def expand(
        self,
        molecule: str,
        *,
        provided_synonyms: Optional[Iterable[str]] = None,
        smiles: Optional[str] = None,
    ) -> List[str]:
        candidates: OrderedDict[str, str] = OrderedDict()
        self._add_candidate(candidates, molecule)
        if provided_synonyms:
            for synonym in provided_synonyms:
                self._add_candidate(candidates, synonym)
        if smiles:
            self._add_candidate(candidates, smiles.strip())
        matched = self._match_entry(molecule, smiles)
        if matched:
            self._add_candidate(candidates, matched.get("canonical", ""))
            for synonym in matched.get("synonyms", []) or []:
                self._add_candidate(candidates, synonym)
        # Heuristic variants (case, punctuation stripping)
        base_variants = self._heuristic_variants(molecule)
        for variant in base_variants:
            self._add_candidate(candidates, variant)
        return [value for value in candidates.values() if value]

    # Internal helpers -------------------------------------------------
    def _match_entry(self, molecule: str, smiles: Optional[str]):
        if smiles and smiles.strip() in self._smiles_index:
            return self._smiles_index[smiles.strip()]
        key = self._key(molecule)
        return self._name_index.get(key)

    def _add_candidate(self, store: "OrderedDict[str, str]", value: Optional[str]) -> None:
        if not value:
            return
        normalized = self._key(value)
        if not normalized:
            return
        if normalized not in store:
            store[normalized] = value.strip()

    def _heuristic_variants(self, molecule: str) -> List[str]:
        token = molecule.strip()
        if not token:
            return []
        variants = {
            token.lower(),
            token.upper(),
            token.title(),
            token.replace(" ", ""),
            token.replace("-", ""),
            token.replace("-", " "),
        }
        return [variant for variant in variants if variant and variant != token]

    def _key(self, value: str) -> str:
        return value.strip().lower().replace(" ", "").replace("-", "")


def load_catalog_from_file(path: Path) -> List[Dict[str, str | List[str]]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)
