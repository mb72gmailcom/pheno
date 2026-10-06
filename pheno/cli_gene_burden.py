"""CLI for per-gene Otari burden."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pheno.gene_burden import compute_gene_burden, load_analysis_groups
from pheno.mapping import load_column_map


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pheno-asd-burden-gene",
        description="Compute per-gene Otari burden for ASD children and unaffected siblings.",
    )
    parser.add_argument("--family-file", required=True, type=Path, help="Family TSV with person, family, parents, and ASD columns")
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
        help="Parent of {otari-prefix}/chrN/{start}_{end}/variant_effects_comprehensive.tsv",
    )
    parser.add_argument(
        "--annotation",
        required=True,
        type=Path,
        help="Otari clean gene table used to assign variants to genes",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for chrN/genes.json, chrN/gene_summary.json, and gene_summary.json",
    )
    parser.add_argument(
        "--file-pattern",
        required=True,
        help="Input TSV prefix, e.g. inherited or inherited_asd",
    )
    parser.add_argument(
        "--otari-prefix",
        default=None,
        help="Otari directory name under --otari-dir. Default: --file-pattern",
    )
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
        people, n_asd, n_siblings = load_analysis_groups(args.family_file, column_map)
        compute_gene_burden(
            args.input_dir,
            args.otari_dir,
            args.annotation,
            args.file_pattern,
            people,
            columns,
            threshold=args.threshold,
            use_abs=args.abs,
            transcripts=transcripts,
            n_asd=n_asd,
            n_unaffected_sibling=n_siblings,
            output_dir=args.output_dir,
            otari_prefix=args.otari_prefix,
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote gene burden to {args.output_dir}")
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
