"""CLI for the within-family ASD dosage contrast."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pheno.gene_burden import load_analysis_groups
from pheno.mapping import load_column_map
from pheno.transmission import compute_transmission


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pheno-asd-transmission",
        description=(
            "For each variant in an ASD and unaffected sibship, write the "
            "number of informative families, the dosage z-score, and the "
            "exact two-sided permutation p-value."
        ),
    )
    parser.add_argument(
        "--family-file",
        required=True,
        type=Path,
        help="Family TSV with a person id and an ASD column",
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Parent of chrN, chrX, and chrY directories of TRIO_CALLS TSVs",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output TSV with one row per variant",
    )
    parser.add_argument(
        "--file-pattern",
        required=True,
        help="Filename prefix, e.g. inherited or denovo",
    )
    parser.add_argument(
        "--column-map",
        default=None,
        type=Path,
        help="JSON map of family-file column names and ASD value aliases",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        column_map = load_column_map(args.column_map)
        people, _n_asd, _n_siblings = load_analysis_groups(args.family_file, column_map)
        written = compute_transmission(args.input_dir, args.file_pattern, people, args.output)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote {written} variants to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
