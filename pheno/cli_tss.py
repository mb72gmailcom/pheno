"""CLI for distances from variants to transcript start sites."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pheno.tss import write_tss_distances


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pheno-tss",
        description=(
            "Measure variants from an inherited or denovo run against "
            "protein-coding transcript start sites."
        ),
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Parent of chrN directories of variant TSVs",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory of chrN/tss.tsv files",
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Filename prefix, such as inherited or denovo. Reads chrN/{prefix}*.tsv",
    )
    parser.add_argument(
        "--gtf",
        required=True,
        type=Path,
        help="GENCODE GTF or GTF.gz annotation",
    )
    parser.add_argument(
        "--window",
        default=1000,
        type=int,
        help="Bases on each side of a transcript TSS (default 1000)",
    )
    parser.add_argument(
        "--indels",
        action="store_true",
        help="Keep only alleles whose reference and alternative lengths differ",
    )
    parser.add_argument(
        "--chr",
        default=None,
        help="Process only this chromosome, such as 1 or X. Without it, every chromosome directory is processed.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        written = write_tss_distances(
            args.input_dir,
            args.output_dir,
            args.prefix,
            args.gtf,
            window=args.window,
            indels_only=args.indels,
            chrom=args.chr,
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote {written} rows to {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
