"""CLI for gene-level ASD association (gene + upstream/downstream flanks)."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from pheno.aggregate import aggregate_genes
from pheno.gtf import DEFAULT_GENE_TYPES
from pheno.mapping import load_column_map
from pheno.pedigree import load_cohort
from pheno.stats import add_stats

logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pheno-asd-gene",
        description=(
            "Union rare-variant carriers in each gene plus strand-aware flanks, "
            "then test ASD association (McNemar exact primary; Fisher extra columns)."
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
        "--gtf",
        required=True,
        help="GENCODE GTF for this chromosome (gene features)",
    )
    parser.add_argument(
        "--column-map",
        default=None,
        help="JSON map of family-file column names and value aliases",
    )
    parser.add_argument(
        "--upstream",
        type=int,
        default=5000,
        help="Bases upstream of TSS (strand-aware, default 5000)",
    )
    parser.add_argument(
        "--downstream",
        type=int,
        default=1000,
        help="Bases downstream of TES (strand-aware, default 1000)",
    )
    parser.add_argument(
        "--gene-types",
        default="protein_coding",
        help="Comma-separated GTF gene_type values (default protein_coding)",
    )
    parser.add_argument(
        "--all-gene-types",
        action="store_true",
        help="Use every GTF gene feature (overrides --gene-types)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV path (one file for this chromosome, one row per gene)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )
    if args.upstream < 0 or args.downstream < 0:
        raise ValueError("upstream and downstream must be >= 0")
    gene_types: tuple[str, ...] | None
    if args.all_gene_types:
        gene_types = None
    else:
        gene_types = tuple(token.strip() for token in args.gene_types.split(",") if token.strip())
        if not gene_types:
            gene_types = DEFAULT_GENE_TYPES

    column_map = load_column_map(args.column_map)
    cohort = load_cohort(args.family_file, column_map)
    counts = aggregate_genes(
        args.idir,
        args.file_pattern,
        args.gtf,
        cohort,
        upstream=args.upstream,
        downstream=args.downstream,
        gene_types=gene_types,
    )
    if counts.empty:
        logger.info("no genes with overlapping variants; writing empty table")
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        counts.to_csv(output, sep="\t", index=False)
        return 0

    result = add_stats(counts)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False)
    logger.info("wrote %s genes to %s", len(result), output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
