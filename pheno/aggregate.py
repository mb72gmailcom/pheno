"""Union carriers of all rare variants overlapping a gene ± flanks."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from pheno.counts import (
    COUNT_COLUMNS,
    _normalize_variant_columns,
    counts_from_patient_lists,
    discover_files,
)
from pheno.gtf import load_genes, normalize_chrom, overlap_variants
from pheno.pedigree import Cohort

logger = logging.getLogger(__name__)


def _read_variants(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    frame = _normalize_variant_columns(frame)
    frame = frame.reset_index(drop=True)
    frame["vid"] = range(len(frame))
    frame["variant_key"] = (
        frame["CHROM"].map(normalize_chrom)
        + "_"
        + frame["POS"].astype(str)
        + "_"
        + frame["REF"].astype(str)
        + "_"
        + frame["ALT"].astype(str)
    )
    return frame


def _union_patients(hits: pd.DataFrame, variants: pd.DataFrame) -> pd.DataFrame:
    merged = hits.merge(variants[["vid", "variant_key", "PATIENTS"]], on="vid", how="inner")
    n_variants = (
        merged[["gene_id", "variant_key"]]
        .drop_duplicates()
        .groupby("gene_id", sort=False)
        .size()
        .rename("n_variants")
    )
    exploded = merged.assign(person=merged["PATIENTS"].fillna("").astype(str).str.split(";"))
    exploded = exploded.explode("person", ignore_index=True)
    exploded["person"] = exploded["person"].str.strip()
    exploded = exploded[exploded["person"].ne("")]
    exploded = exploded.drop_duplicates(["gene_id", "person"])
    patients = exploded.groupby("gene_id", sort=False)["person"].agg(";".join).rename("PATIENTS")
    return pd.concat([n_variants, patients], axis=1).reset_index()


def aggregate_directory(
    idir: str | Path,
    file_pattern: str,
    genes: pd.DataFrame,
    cohort: Cohort,
) -> pd.DataFrame:
    pair_frames: list[pd.DataFrame] = []
    files = discover_files(idir, file_pattern)
    logger.info("aggregating %s variant files from %s", len(files), idir)
    for path in files:
        logger.info("  %s", path.name)
        variants = _read_variants(path)
        hits = overlap_variants(variants, genes)
        if hits.empty:
            continue
        pair_frames.append(_union_patients(hits, variants))

    if not pair_frames:
        logger.info("no variants overlapped any gene interval")
        return pd.DataFrame(
            columns=[
                "CHROM",
                "gene_id",
                "gene_name",
                "gene_type",
                "strand",
                "gene_start",
                "gene_end",
                "region_start",
                "region_end",
                "n_variants",
                "key",
                *COUNT_COLUMNS,
            ]
        )

    combined = pd.concat(pair_frames, ignore_index=True)
    n_variants = combined.groupby("gene_id", sort=False)["n_variants"].sum()
    exploded = combined.assign(person=combined["PATIENTS"].fillna("").astype(str).str.split(";"))
    exploded = exploded.explode("person", ignore_index=True)
    exploded["person"] = exploded["person"].str.strip()
    exploded = exploded[exploded["person"].ne("")]
    exploded = exploded.drop_duplicates(["gene_id", "person"])
    patients = exploded.groupby("gene_id", sort=False)["person"].agg(";".join)

    units = genes.loc[genes["gene_id"].isin(n_variants.index)].copy()
    units = units.drop_duplicates(subset="gene_id", keep="first")
    units["n_variants"] = units["gene_id"].map(n_variants).fillna(0).astype("int64")
    units["PATIENTS"] = units["gene_id"].map(patients).fillna("")
    units["key"] = units["gene_id"]
    units["vid"] = range(len(units))
    counts = counts_from_patient_lists(units, units["PATIENTS"], cohort)
    return counts.drop(columns=["PATIENTS", "vid"], errors="ignore")


def aggregate_genes(
    idir: str | Path,
    file_pattern: str,
    gtf_path: str | Path,
    cohort: Cohort,
    upstream: int = 5000,
    downstream: int = 1000,
    gene_types: tuple[str, ...] | None = ("protein_coding",),
) -> pd.DataFrame:
    genes = load_genes(gtf_path, upstream=upstream, downstream=downstream, gene_types=gene_types)
    return aggregate_directory(idir, file_pattern, genes, cohort)
