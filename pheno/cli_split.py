"""CLI for splitting variant TSVs by ASD and unaffected carriers."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pheno.mapping import load_column_map
from pheno.split import load_asd_status, split_directory


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pheno-asd-split",
        description=(
            "Split variant TSVs into ASD-only, unaffected-only, and both-carrier files."
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
        help="Parent of chrN, chrX, and chrY directories of short-format variant TSVs",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Output directory. Files are {prefix}_{asd,unaffected,both}_*.tsv",
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
        split_directory(args.input_dir, args.output_dir, args.file_pattern, status)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote split TSVs to {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
