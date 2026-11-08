"""Distance from variants to annotated protein-coding splice sites."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from pheno.gtf import _open_text, _parse_attributes
from pheno.tss import (
    Allele,
    _distance,
    _overlaps,
    _read_variants,
    _select_dirs,
    _write_table,
    canonical_chrom,
    chromosome_name,
)

_OUTPUT_NAME = "splice.tsv"
_FIELDS = (
    "chrom",
    "pos",
    "ref",
    "alt",
    "patients",
    "transcript_id",
    "gene_id",
    "gene_name",
    "strand",
    "site_type",
    "site",
    "distance",
    "covers_site",
)


@dataclass(frozen=True)
class SpliceSite:
    transcript_id: str
    gene_id: str
    gene_name: str
    strand: str
    site_type: str
    site: int


@dataclass
class _TranscriptExons:
    gene_id: str = ""
    gene_name: str = ""
    strand: str = ""
    tx_type: str = ""
    saw_transcript: bool = False
    exons: list[tuple[int, int]] = field(default_factory=list)


def load_splice_sites(gtf_path: Path) -> dict[str, list[SpliceSite]]:
    """Load donor and acceptor bases for protein-coding transcripts.

    Each intron contributes two sites. The donor is the first intron base in
    the direction of transcription, and the acceptor is the last. A one-exon
    transcript has no splice site. ``transcript_type`` is preferred over
    ``gene_type``.
    """
    if not gtf_path.is_file():
        raise ValueError(f"annotation file not found: {gtf_path}")
    grouped = _collect_exons(gtf_path)
    sites: dict[str, list[SpliceSite]] = {}
    for (chrom, transcript_id), transcript in grouped.items():
        if transcript.tx_type != "protein_coding" or transcript.strand not in {"+", "-"}:
            continue
        if not transcript.gene_id:
            continue
        for site_type, site in _junctions(transcript.exons, transcript.strand):
            sites.setdefault(chrom, []).append(
                SpliceSite(
                    transcript_id=transcript_id,
                    gene_id=transcript.gene_id,
                    gene_name=transcript.gene_name,
                    strand=transcript.strand,
                    site_type=site_type,
                    site=site,
                )
            )
    if not sites:
        raise ValueError(f"no protein_coding splice sites loaded from {gtf_path}")
    return sites


def write_splice_distances(
    input_dir: Path,
    output_dir: Path,
    prefix: str,
    gtf_path: Path,
    *,
    window: int = 1000,
    indels_only: bool = False,
    chrom: str | None = None,
) -> int:
    """Write ``{output_dir}/chrN/splice.tsv`` for each processed chromosome.

    A row is one variant and one annotated splice site within ``window`` bases
    on either side. Distance is signed in the direction of transcription:
    negative is upstream of that site, and zero means the variant interval
    contains it. ``covers_site`` is 1 only when a deleted or substituted
    reference base is the splice-site base.
    """
    if window < 0:
        raise ValueError("--window must be at least 0")
    if not prefix:
        raise ValueError("--prefix is empty")
    if not input_dir.is_dir():
        raise ValueError(f"input directory not found: {input_dir}")
    only = None if chrom is None else chromosome_name(chrom)
    sites = load_splice_sites(gtf_path)
    targets = _select_dirs(input_dir, only)
    output_dir.mkdir(parents=True, exist_ok=True)
    total = 0
    processed = 0
    for chrom_dir in targets:
        name = chrom_dir.name
        files = sorted(
            (path for path in chrom_dir.glob(f"{prefix}*.tsv") if path.is_file()),
            key=lambda path: path.name,
        )
        if not files:
            if only is not None:
                raise ValueError(f"no {prefix}*.tsv files in {chrom_dir}")
            continue
        alleles: list[Allele] = []
        for path in files:
            alleles.extend(_read_variants(path, name, indels_only))
        rows = _match(alleles, sites.get(name, []), window, name)
        _write_table(output_dir / name / _OUTPUT_NAME, rows, _FIELDS)
        total += len(rows)
        processed += 1
        print(f"finished {name}: {len(rows)} rows", flush=True)
    if processed == 0:
        raise ValueError(f"no {prefix}*.tsv files under {input_dir}")
    return total


def _collect_exons(gtf_path: Path) -> dict[tuple[str, str], _TranscriptExons]:
    grouped: dict[tuple[str, str], _TranscriptExons] = {}
    with _open_text(gtf_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"transcript", "exon"}:
                continue
            chrom = canonical_chrom(fields[0])
            if chrom is None:
                continue
            attrs = _parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id", "")
            if not transcript_id:
                continue
            record = grouped.setdefault((chrom, transcript_id), _TranscriptExons())
            tx_type = (
                attrs.get("transcript_type")
                or attrs.get("transcript_biotype")
                or attrs.get("gene_type")
                or attrs.get("gene_biotype")
                or ""
            )
            if fields[2] == "transcript":
                record.saw_transcript = True
                record.tx_type = tx_type
                record.strand = fields[6]
                record.gene_id = attrs.get("gene_id", "") or record.gene_id
                record.gene_name = attrs.get("gene_name", "") or record.gene_name
                continue
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError as exc:
                raise ValueError(f"{gtf_path}:{line_number} has a non-numeric coordinate") from exc
            if (start, end) not in record.exons:
                record.exons.append((start, end))
            if record.saw_transcript:
                continue
            record.tx_type = tx_type or record.tx_type
            record.strand = fields[6] or record.strand
            record.gene_id = attrs.get("gene_id", "") or record.gene_id
            record.gene_name = attrs.get("gene_name", "") or record.gene_name
    return grouped


def _junctions(exons: list[tuple[int, int]], strand: str) -> list[tuple[str, int]]:
    ordered = sorted(exons)
    sites: list[tuple[str, int]] = []
    for (_left_start, left_end), (right_start, _right_end) in zip(ordered, ordered[1:]):
        intron_start = left_end + 1
        intron_end = right_start - 1
        if intron_end < intron_start:
            continue
        if strand == "+":
            sites.append(("donor", intron_start))
            sites.append(("acceptor", intron_end))
        else:
            sites.append(("donor", intron_end))
            sites.append(("acceptor", intron_start))
    return sites


def _match(
    alleles: list[Allele],
    sites: list[SpliceSite],
    window: int,
    chrom: str,
) -> list[list[object]]:
    if not alleles or not sites:
        return []
    bucket = max(window, 1)
    index: dict[int, list[int]] = {}
    for i, site in enumerate(sites):
        index.setdefault(site.site // bucket, []).append(i)
    rows: list[list[object]] = []
    for allele in alleles:
        first = (allele.start - window) // bucket
        last = (allele.end + window) // bucket
        for key in range(first, last + 1):
            for i in index.get(key, ()):
                site = sites[i]
                if not _overlaps(allele.start, allele.end, site.site, window):
                    continue
                covers = int(allele.covers_reference and allele.start <= site.site <= allele.end)
                rows.append(
                    [
                        chrom,
                        allele.pos,
                        allele.ref,
                        allele.alt,
                        allele.patients,
                        site.transcript_id,
                        site.gene_id,
                        site.gene_name,
                        site.strand,
                        site.site_type,
                        site.site,
                        _distance(site.strand, site.site, allele.start, allele.end),
                        covers,
                    ]
                )
    rows.sort(key=lambda row: (int(row[1]), str(row[2]), str(row[3]), str(row[5]), int(row[10]), str(row[9])))
    return rows
