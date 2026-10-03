"""CLI for per-patient variant counts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pheno.mapping import load_column_map
from pheno.patient_counts import count_patient_variants, write_patient_counts
from pheno.split import load_asd_status


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pheno-asd-counts",
        description="Count variants per ASD patient and per unaffected patient.",
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
        help="Parent of chrN, chrX, and chrY directories of short-format variant TSVs",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output JSON file with asd and unaffected counts",
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
        status = load_asd_status(args.family_file, column_map)
        counts = count_patient_variants(args.input_dir, args.file_pattern, status)
        write_patient_counts(args.output, counts)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote patient counts to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
