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
        description="Compute per-gene burden for ASD children and unaffected siblings.",
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
        type=Path,
        help="Parent of chrN/{start}_{end}/variant_effects_comprehensive.tsv",
    )
    parser.add_argument(
        "--score-dir",
        type=Path,
        help="Parent of chrN/{score-prefix}_{start}_{end}.tsv score files",
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
        help="Directory for chrN/genes.tsv.gz, chrN/gene_summary.json, and gene_summary.json",
    )
    parser.add_argument(
        "--file-pattern",
        required=True,
        help="Input TSV prefix, e.g. inherited or inherited_asd",
    )
    parser.add_argument(
        "--score-prefix",
        default=None,
        help="Filename prefix under --score-dir/chrN. Required with --score-dir",
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
    parser.add_argument(
        "--score-columns",
        default=None,
        help="Comma-separated columns in the score TSV. Required with --score-dir",
    )
    parser.add_argument("--threshold", required=True, type=float, help="Damaging when the score is greater than this value")
    parser.add_argument(
        "--window",
        type=int,
        default=2000,
        help="Base pairs past each end of the gene. The gene body is included. Default: 2000",
    )
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
    transcripts = "mean" if args.mean_transcripts else "max"
    try:
        if (args.otari_dir is None) == (args.score_dir is None):
            raise ValueError("set exactly one of --otari-dir and --score-dir")
        if args.score_dir is not None:
            if not args.score_prefix:
                raise ValueError("--score-prefix is required with --score-dir")
            if not args.score_columns:
                raise ValueError("--score-columns is required with --score-dir")
            columns = _parse_columns(args.score_columns, "score")
        else:
            columns = _parse_columns(args.otari_columns, "Otari")
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
            score_dir=args.score_dir,
            score_prefix=args.score_prefix,
            window=args.window,
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote gene burden to {args.output_dir}")
    return 0


def _parse_columns(text: str, kind: str) -> list[str]:
    columns: list[str] = []
    for part in text.split(","):
        name = part.strip()
        if name and name not in columns:
            columns.append(name)
    if not columns:
        raise ValueError(f"at least one {kind} column is required")
    return columns


if __name__ == "__main__":
    sys.exit(main())
