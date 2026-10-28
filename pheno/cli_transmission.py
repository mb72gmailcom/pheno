"""CLI for the within-family ASD dosage contrast."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pheno.gene_burden import load_analysis_groups
from pheno.mapping import load_column_map
from pheno.transmission import compute_transmission, gene_table_path
from pheno.transmission_scores import gene_scores_path


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
    parser.add_argument(
        "--genes",
        default=None,
        type=Path,
        help="Annotation table. Also write a gene-level dosage contrast",
    )
    parser.add_argument(
        "--window",
        default=None,
        type=int,
        help="Bases past both ends of a gene when --genes is set (default 2000)",
    )
    parser.add_argument(
        "--fcap",
        default=None,
        type=int,
        help=(
            "Exact p-value when the number of informative families is at most this value. "
            "Above it, use the two-sided normal tail of z"
        ),
    )
    parser.add_argument(
        "--otari-dir",
        default=None,
        type=Path,
        help=(
            "Parent of chrN/{start}_{end}/variant_effects_comprehensive.tsv. "
            "Weight each gene by the signed tissue scores. Mutually exclusive with --seqweaver"
        ),
    )
    parser.add_argument(
        "--seqweaver",
        default=None,
        type=Path,
        help=(
            "Parent of chrN directories of seqweaver_human diffs. "
            "Weight each gene by those columns. Mutually exclusive with --otari-dir"
        ),
    )
    parser.add_argument(
        "--column-names",
        default=None,
        type=Path,
        help="Seqweaver column names, one per line, in matrix order. Requires --seqweaver",
    )
    parser.add_argument(
        "--chr",
        default=None,
        help="Process only this chromosome, such as 1 or X. Required with --debug-gene.",
    )
    parser.add_argument(
        "--debug-gene",
        default=None,
        help="Print how this gene id is assigned and whether its row is written. Requires --chr. A version suffix is ignored.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if (args.output is None) == (args.output_dir is None):
        print("error: set exactly one of --output and --output-dir", file=sys.stderr)
        return 1
    if args.window is not None and args.genes is None:
        print("error: --window requires --genes", file=sys.stderr)
        return 1
    try:
        column_map = load_column_map(args.column_map)
        people, _n_asd, _n_siblings = load_analysis_groups(args.family_file, column_map)
        written, gene_written = compute_transmission(
            args.input_dir,
            args.file_pattern,
            people,
            args.output,
            output_dir=args.output_dir,
            genes=args.genes,
            window=2000 if args.window is None else args.window,
            fcap=args.fcap,
            otari_dir=args.otari_dir,
            seqweaver_dir=args.seqweaver,
            column_names=args.column_names,
            debug_gene=args.debug_gene,
            chrom=args.chr,
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    destination = args.output if args.output is not None else args.output_dir
    if args.otari_dir is not None or args.seqweaver is not None:
        score_destination = gene_scores_path(args.output) if args.output is not None else args.output_dir
        print(f"Wrote {written} variants to {destination} and {gene_written} gene scores to {score_destination}")
    elif args.genes is None:
        print(f"Wrote {written} variants to {destination}")
    else:
        gene_destination = gene_table_path(args.output) if args.output is not None else args.output_dir
        print(f"Wrote {written} variants to {destination} and {gene_written} genes to {gene_destination}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
