"""Load GENCODE gene intervals and expand strand-aware flanks."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

DEFAULT_GENE_TYPES = ("protein_coding",)


def normalize_chrom(value: str) -> str:
    token = str(value).strip()
    if token.lower().startswith("chr"):
        return "chr" + token[3:]
    return "chr" + token


def _parse_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in raw.strip().split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, value = part.split(None, 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def _open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def expand_flanks(start: int, end: int, strand: str, upstream: int, downstream: int) -> tuple[int, int]:
    if strand == "-":
        region_start = start - downstream
        region_end = end + upstream
    else:
        region_start = start - upstream
        region_end = end + downstream
    return max(1, region_start), max(max(1, region_start), region_end)


def load_genes(
    gtf_path: str | Path,
    upstream: int = 5000,
    downstream: int = 1000,
    gene_types: tuple[str, ...] | None = DEFAULT_GENE_TYPES,
) -> pd.DataFrame:
    """Return one row per GTF gene feature, with flank-expanded coordinates."""
    path = Path(gtf_path)
    rows: list[dict[str, object]] = []
    with _open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = _parse_attributes(fields[8])
            gene_type = attrs.get("gene_type") or attrs.get("gene_biotype") or ""
            if gene_types is not None and gene_type not in gene_types:
                continue
            gene_id = attrs.get("gene_id", "")
            if not gene_id:
                continue
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            region_start, region_end = expand_flanks(start, end, strand, upstream, downstream)
            rows.append(
                {
                    "CHROM": normalize_chrom(fields[0]),
                    "gene_id": gene_id,
                    "gene_name": attrs.get("gene_name", ""),
                    "gene_type": gene_type,
                    "strand": strand,
                    "gene_start": start,
                    "gene_end": end,
                    "region_start": region_start,
                    "region_end": region_end,
                }
            )
    genes = pd.DataFrame(rows)
    if genes.empty:
        raise ValueError(f"no gene features loaded from {path}")
    n_before = len(genes)
    genes = genes.drop_duplicates(subset="gene_id", keep="first")
    if len(genes) < n_before:
        logger.warning("dropped %s duplicate gene_id rows in GTF", n_before - len(genes))
    logger.info(
        "loaded %s genes from %s (upstream=%s, downstream=%s)",
        len(genes),
        path,
        upstream,
        downstream,
    )
    return genes.reset_index(drop=True)


def overlap_variants(variants: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """Return variant/gene pairs where POS lies in a gene's flank-expanded interval."""
    if variants.empty or genes.empty:
        return pd.DataFrame(columns=["vid", "gene_id"])

    work = variants[["vid", "CHROM", "POS"]].copy()
    work["CHROM"] = work["CHROM"].map(normalize_chrom)
    work["POS"] = work["POS"].astype(int)

    hits: list[tuple[int, str]] = []
    for chrom, var_chr in work.groupby("CHROM", sort=False):
        gene_chr = genes.loc[genes["CHROM"] == chrom]
        if gene_chr.empty:
            continue
        events: list[tuple[int, int, int]] = []
        gene_ids = gene_chr["gene_id"].to_numpy()
        for i, (start, end) in enumerate(
            zip(gene_chr["region_start"].to_numpy(), gene_chr["region_end"].to_numpy())
        ):
            events.append((int(start), 0, i))
            events.append((int(end), 2, i))
        vids = var_chr["vid"].to_numpy()
        for vid, pos in zip(vids, var_chr["POS"].to_numpy()):
            events.append((int(pos), 1, int(vid)))
        events.sort()
        active: set[int] = set()
        for _coord, kind, idx in events:
            if kind == 0:
                active.add(idx)
            elif kind == 2:
                active.discard(idx)
            else:
                for gene_i in active:
                    hits.append((idx, str(gene_ids[gene_i])))

    if not hits:
        return pd.DataFrame(columns=["vid", "gene_id"])
    return pd.DataFrame(hits, columns=["vid", "gene_id"]).drop_duplicates()
