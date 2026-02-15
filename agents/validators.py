"""Input validation utilities for molecule names and SMILES strings."""

from __future__ import annotations

import re

from .exceptions import ValidationError

# Allow alphanumeric, hyphens, spaces, parentheses, brackets, commas, periods
# Max 200 chars — enough for long IUPAC names.
_MOLECULE_NAME_RE = re.compile(r"^[A-Za-z0-9\-\s\(\)\[\],\.'\u00C0-\u024F]{1,200}$")

# SMILES charset (covers organic subset + stereochemistry + charges)
_SMILES_RE = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)\\/=#%\$\.\:~]{1,500}$")


def validate_molecule_name(name: str) -> str:
    """Return a cleaned molecule name or raise ValidationError."""
    name = name.strip()
    if not name:
        raise ValidationError("Molecule name cannot be empty", code="EMPTY_MOLECULE")
    if len(name) > 200:
        raise ValidationError(
            f"Molecule name too long ({len(name)} chars, max 200)",
            code="MOLECULE_NAME_TOO_LONG",
        )
    if not _MOLECULE_NAME_RE.match(name):
        raise ValidationError(
            f"Invalid molecule name: '{name[:50]}...'",
            code="INVALID_MOLECULE_NAME",
        )
    return name


def validate_smiles(smiles: str | None) -> str | None:
    """Return a cleaned SMILES string or raise ValidationError."""
    if smiles is None:
        return None
    smiles = smiles.strip()
    if not smiles:
        return None
    if len(smiles) > 500:
        raise ValidationError(
            "SMILES string too long (max 500 chars)",
            code="SMILES_TOO_LONG",
        )
    if not _SMILES_RE.match(smiles):
        raise ValidationError("Invalid SMILES string", code="INVALID_SMILES")
    return smiles


def validate_synonyms(synonyms: list[str] | None, *, max_count: int = 50) -> list[str]:
    """Validate and deduplicate a list of synonym strings."""
    if not synonyms:
        return []
    if len(synonyms) > max_count:
        raise ValidationError(
            f"Too many synonyms ({len(synonyms)}, max {max_count})",
            code="TOO_MANY_SYNONYMS",
        )
    seen: set[str] = set()
    result: list[str] = []
    for syn in synonyms:
        cleaned = syn.strip()
        if not cleaned:
            continue
        lower = cleaned.lower()
        if lower in seen:
            continue
        seen.add(lower)
        result.append(cleaned)
    return result
