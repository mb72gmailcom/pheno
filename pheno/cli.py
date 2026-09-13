"""CLI for per-chromosome ASD association counts."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from pheno.counts import count_directory
from pheno.mapping import load_column_map
from pheno.pedigree import load_cohort
from pheno.stats import add_stats

logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pheno-asd",
        description=(
            "Count inherited-variant carriers among trio children and test ASD "
            "association (McNemar exact primary; Fisher extra columns)."
        ),
    )
    parser.add_argument("--family-file", required=True, help="Family TSV (person, parents, sex, asd)")
    parser.add_argument(
        "--idir",
        required=True,
        help="Directory of variant TSVs for a single chromosome",
    )
    parser.add_argument(
        "--file-pattern",
        required=True,
        help="Filename prefix, e.g. inherited → inherited.tsv and inherited_*.tsv",
    )
    parser.add_argument(
        "--column-map",
        default=None,
        help="JSON map of family-file column names and value aliases",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV path (one file for this chromosome)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )
    column_map = load_column_map(args.column_map)
    cohort = load_cohort(args.family_file, column_map)
    counts = count_directory(args.idir, args.file_pattern, cohort)
    result = add_stats(counts)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False)
    logger.info("wrote %s variants to %s", len(result), output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
