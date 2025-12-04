"""CLI wrapper for rendering RDKit structure assets from canonical SMILES records."""

from __future__ import annotations

import argparse
from pathlib import Path

from structure_renderer import StructureRenderConfig, render_structure_catalog

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = BASE_DIR / "data" / "normalized_smiles.jsonl"
DEFAULT_OUTPUT_DIR = BASE_DIR / "data" / "structures" / "images"
DEFAULT_METADATA_DIR = BASE_DIR / "data" / "structures" / "metadata"
DEFAULT_MANIFEST = BASE_DIR / "data" / "structures" / "structures.manifest.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render RDKit structure SVG assets for downstream pipelines")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="JSONL file with canonical SMILES entries")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="Where to store SVG assets")
    parser.add_argument(
        "--metadata-dir",
        type=Path,
        default=DEFAULT_METADATA_DIR,
        help="Directory for per-structure metadata JSON files",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="Catalog manifest JSON file (records every rendered structure)",
    )
    parser.add_argument("--width", type=int, default=420, help="SVG width in pixels (default: 420)")
    parser.add_argument("--height", type=int, default=320, help="SVG height in pixels (default: 320)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = StructureRenderConfig(
        input_path=args.input.resolve(),
        output_dir=args.output_dir.resolve(),
        metadata_dir=args.metadata_dir.resolve(),
        manifest_path=args.manifest.resolve(),
        width=args.width,
        height=args.height,
    )
    summary = render_structure_catalog(config)
    print(
        f"Rendered {summary['rendered']} structure(s) from {summary['processed']} record(s). "
        f"Manifest: {summary['manifest']}"
    )


if __name__ == "__main__":  # pragma: no cover
    main()
