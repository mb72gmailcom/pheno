"""CLI for per-person Otari burden."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pheno.burden import compute_burden, load_people
from pheno.mapping import load_column_map


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pheno-asd-burden",
        description="Compute per-person Otari burden for variants in short-format TSVs.",
    )
    parser.add_argument("--family-file", required=True, type=Path, help="Family TSV with a person id and an ASD column")
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Parent of chrN directories of {prefix}_{start}_{end}.tsv files",
    )
    parser.add_argument(
        "--otari-dir",
        required=True,
        type=Path,
        help="Parent of {prefix}/chrN/{start}_{end}/variant_effects_comprehensive.tsv",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Genome-wide JSON file. Each chromosome is also written to {parent}/chrN/{name}",
    )
    parser.add_argument("--file-pattern", required=True, help="Filename prefix, e.g. inherited or denovo")
    parser.add_argument(
        "--column-map",
        default=None,
        type=Path,
        help="JSON map of family-file column names and ASD value aliases",
    )
    parser.add_argument(
        "--otari-columns",
        default="max_effect",
        help="Comma-separated Otari score columns. Default: max_effect",
    )
    parser.add_argument("--threshold", required=True, type=float, help="Damaging when the collapsed score is greater than this value")
    parser.add_argument(
        "--abs",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Take the absolute value of each transcript score before collapsing (default: true)",
    )
    transcripts = parser.add_mutually_exclusive_group()
    transcripts.add_argument("--max-transcripts", action="store_true", help="Collapse transcripts by maximum (default)")
    transcripts.add_argument("--mean-transcripts", action="store_true", help="Collapse transcripts by mean")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    columns = _parse_columns(args.otari_columns)
    transcripts = "mean" if args.mean_transcripts else "max"
    try:
        column_map = load_column_map(args.column_map)
        people = load_people(args.family_file, column_map)
        compute_burden(
            args.input_dir,
            args.otari_dir,
            args.file_pattern,
            people,
            columns,
            threshold=args.threshold,
            use_abs=args.abs,
            transcripts=transcripts,
            output=args.output,
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote burden to {args.output}")
    return 0


def _parse_columns(text: str) -> list[str]:
    columns: list[str] = []
    for part in text.split(","):
        name = part.strip()
        if name and name not in columns:
            columns.append(name)
    if not columns:
        raise ValueError("at least one Otari column is required")
    return columns


if __name__ == "__main__":
    sys.exit(main())
