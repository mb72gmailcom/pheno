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
        default=None,
        type=Path,
        help="Single TSV with one row per variant. Mutually exclusive with --output-dir",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        type=Path,
        help=(
            "Directory of chrN/{file-pattern}.tsv files, ordered by position. "
            "Mutually exclusive with --output"
        ),
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
    if (args.output is None) == (args.output_dir is None):
        print("error: set exactly one of --output and --output-dir", file=sys.stderr)
        return 1
    try:
        column_map = load_column_map(args.column_map)
        people, _n_asd, _n_siblings = load_analysis_groups(args.family_file, column_map)
        written = compute_transmission(
            args.input_dir,
            args.file_pattern,
            people,
            args.output,
            output_dir=args.output_dir,
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    destination = args.output if args.output is not None else args.output_dir
    print(f"Wrote {written} variants to {destination}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
