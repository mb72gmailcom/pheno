"""CLI for model score column statistics."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pheno.model_stats import compute_model_stats


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pheno-model-stats",
        description="Summarize the distribution of each column in a model score table.",
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Parent of chrN directories of model output",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for chrN/columns.by_row.tsv and chrN/columns.by_variant.tsv",
    )
    parser.add_argument(
        "--otari",
        action="store_true",
        help="Read chrN/{start}_{end}/variant_effects_comprehensive.tsv",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.otari:
        print("error: set --otari", file=sys.stderr)
        return 1
    try:
        written = compute_model_stats(args.input_dir, args.output_dir, otari=True)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote column statistics for {written} chromosomes to {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
