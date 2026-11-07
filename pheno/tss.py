"""Distance from variants to protein-coding transcript start sites."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from pheno.gtf import _open_text, _parse_attributes

_CHROM_DIR = re.compile(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
_OUTPUT_NAME = "tss.tsv"
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
    "tss",
    "distance",
    "covers_tss",
)
_PATIENT_COLUMNS = frozenset({"PATIENTS", "TRIO_CALLS"})


@dataclass(frozen=True)
class Transcript:
    transcript_id: str
    gene_id: str
    gene_name: str
    strand: str
    tss: int


@dataclass(frozen=True)
class Allele:
    pos: int
    ref: str
    alt: str
    patients: str
    start: int
    end: int
    covers_reference: bool


def canonical_chrom(value: str) -> str | None:
    """Return chr1–chr22, chrX, or chrY. Other contigs return None."""
    token = str(value).strip()
    if not token:
        return None
    suffix = token[3:] if token.lower().startswith("chr") else token
    if suffix.upper() in {"X", "Y"}:
        suffix = suffix.upper()
    name = f"chr{suffix}"
    if _CHROM_DIR.fullmatch(name):
        return name
    return None


def chromosome_name(value: str) -> str:
    name = canonical_chrom(value)
    if name is None:
        raise ValueError(f"--chr must be a chromosome number, X, or Y, not {value}")
    return name


def load_transcripts(gtf_path: Path) -> dict[str, list[Transcript]]:
    """Load protein-coding transcript TSS positions from a GTF or GTF.gz file.

    The TSS is the transcript start on ``+`` and the transcript end on ``-``.
    ``transcript_type`` is preferred over ``gene_type``, so a non-coding
    transcript of a protein-coding gene is skipped.
    """
    if not gtf_path.is_file():
        raise ValueError(f"annotation file not found: {gtf_path}")
    grouped: dict[str, list[Transcript]] = {}
    seen: set[tuple[str, str]] = set()
    with _open_text(gtf_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            chrom = canonical_chrom(fields[0])
            if chrom is None:
                continue
            attrs = _parse_attributes(fields[8])
            tx_type = (
                attrs.get("transcript_type")
                or attrs.get("transcript_biotype")
                or attrs.get("gene_type")
                or attrs.get("gene_biotype")
                or ""
            )
            if tx_type != "protein_coding":
                continue
            transcript_id = attrs.get("transcript_id", "")
            gene_id = attrs.get("gene_id", "")
            if not transcript_id or not gene_id:
                continue
            strand = fields[6]
            if strand not in {"+", "-"}:
                continue
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError as exc:
                raise ValueError(f"{gtf_path}:{line_number} has a non-numeric coordinate") from exc
            tss = start if strand == "+" else end
            key = (chrom, transcript_id)
            if key in seen:
                continue
            seen.add(key)
            grouped.setdefault(chrom, []).append(
                Transcript(
                    transcript_id=transcript_id,
                    gene_id=gene_id,
                    gene_name=attrs.get("gene_name", ""),
                    strand=strand,
                    tss=tss,
                )
            )
    if not grouped:
        raise ValueError(f"no protein_coding transcripts loaded from {gtf_path}")
    return grouped


def allele_interval(pos: int, ref: str, alt: str) -> tuple[int, int, bool]:
    """Return the reference interval used for distance, and whether it can cover a base.

    Equal-length alleles use the written reference span. A pure deletion omits
    the anchor base. A pure insertion is the reference base immediately
    downstream of the inserted sequence and does not cover a reference base.
    """
    if len(ref) == len(alt):
        return pos, pos + len(ref) - 1, True
    if len(ref) > len(alt):
        if ref.startswith(alt):
            return pos + len(alt), pos + len(ref) - 1, True
        if ref.endswith(alt):
            return pos, pos + (len(ref) - len(alt)) - 1, True
        return pos, pos + len(ref) - 1, True
    if alt.startswith(ref):
        site = pos + len(ref)
        return site, site, False
    if alt.endswith(ref):
        return pos, pos, False
    return pos, pos + len(ref) - 1, True


def write_tss_distances(
    input_dir: Path,
    output_dir: Path,
    prefix: str,
    gtf_path: Path,
    *,
    window: int = 1000,
    indels_only: bool = False,
    chrom: str | None = None,
) -> int:
    """Write ``{output_dir}/chrN/tss.tsv`` for each processed chromosome.

    A row is one variant and one transcript whose interval falls within
    ``window`` bases on either side of the TSS. Distance is signed in the
    direction of transcription: negative is upstream, and zero means the
    interval contains the TSS. ``covers_tss`` is 1 only when a deleted or
    substituted reference base is the TSS.
    """
    if window < 0:
        raise ValueError("--window must be at least 0")
    if not prefix:
        raise ValueError("--prefix is empty")
    if not input_dir.is_dir():
        raise ValueError(f"input directory not found: {input_dir}")
    only = None if chrom is None else chromosome_name(chrom)
    transcripts = load_transcripts(gtf_path)
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
        rows = _match(alleles, transcripts.get(name, []), window, name)
        _write_table(output_dir / name / _OUTPUT_NAME, rows)
        total += len(rows)
        processed += 1
        print(f"finished {name}: {len(rows)} rows", flush=True)
    if processed == 0:
        raise ValueError(f"no {prefix}*.tsv files under {input_dir}")
    return total


def _select_dirs(input_dir: Path, only: str | None) -> list[Path]:
    dirs = _chrom_dirs(input_dir)
    if only is None:
        if not dirs:
            raise ValueError(f"no chromosome directories under {input_dir}")
        return dirs
    match = [path for path in dirs if path.name == only]
    if not match:
        raise ValueError(f"chromosome directory not found: {input_dir / only}")
    return match


def _chrom_dirs(input_dir: Path) -> list[Path]:
    dirs = [path for path in input_dir.iterdir() if path.is_dir() and _CHROM_DIR.fullmatch(path.name)]

    def sort_key(path: Path) -> tuple[int, int]:
        suffix = path.name[3:]
        if suffix == "X":
            return (1, 0)
        if suffix == "Y":
            return (2, 0)
        return (0, int(suffix))

    return sorted(dirs, key=sort_key)


def _read_variants(path: Path, chrom: str, indels_only: bool) -> list[Allele]:
    alleles: list[Allele] = []
    header_seen = False
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                if header_seen:
                    continue
                header_seen = True
                fields = line[1:].split("\t")
                if len(fields) < 6 or fields[5] not in _PATIENT_COLUMNS:
                    raise ValueError(f"{path} requires a PATIENTS or TRIO_CALLS column")
                continue
            if not header_seen:
                raise ValueError(f"{path} is missing a header")
            fields = line.split("\t")
            if len(fields) < 6:
                raise ValueError(f"{path}:{line_number} expected at least 6 columns")
            row_chrom = canonical_chrom(fields[0])
            if row_chrom != chrom:
                raise ValueError(f"{path}:{line_number} chromosome {fields[0]} is not {chrom}")
            try:
                pos = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number} invalid POS {fields[1]!r}") from exc
            if pos < 1:
                raise ValueError(f"{path}:{line_number} invalid POS {fields[1]!r}")
            ref = fields[3]
            alt = fields[4]
            if _skip_allele(ref, alt, indels_only):
                continue
            start, end, covers_reference = allele_interval(pos, ref, alt)
            alleles.append(
                Allele(
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    patients=fields[5],
                    start=start,
                    end=end,
                    covers_reference=covers_reference,
                )
            )
    if not header_seen:
        raise ValueError(f"{path} is empty")
    return alleles


def _skip_allele(ref: str, alt: str, indels_only: bool) -> bool:
    if not ref or not alt or "," in ref or "," in alt:
        return True
    if _is_symbolic(ref) or _is_symbolic(alt):
        return True
    return indels_only and len(ref) == len(alt)


def _is_symbolic(allele: str) -> bool:
    if allele in {".", "*"}:
        return True
    return "<" in allele or ">" in allele or "[" in allele or "]" in allele


def _match(
    alleles: list[Allele],
    transcripts: list[Transcript],
    window: int,
    chrom: str,
) -> list[list[object]]:
    if not alleles or not transcripts:
        return []
    bucket = max(window, 1)
    index: dict[int, list[int]] = {}
    for i, transcript in enumerate(transcripts):
        index.setdefault(transcript.tss // bucket, []).append(i)
    rows: list[list[object]] = []
    for allele in alleles:
        first = (allele.start - window) // bucket
        last = (allele.end + window) // bucket
        for key in range(first, last + 1):
            for i in index.get(key, ()):
                transcript = transcripts[i]
                if not _overlaps(allele.start, allele.end, transcript.tss, window):
                    continue
                covers = int(allele.covers_reference and allele.start <= transcript.tss <= allele.end)
                rows.append(
                    [
                        chrom,
                        allele.pos,
                        allele.ref,
                        allele.alt,
                        allele.patients,
                        transcript.transcript_id,
                        transcript.gene_id,
                        transcript.gene_name,
                        transcript.strand,
                        transcript.tss,
                        _distance(transcript.strand, transcript.tss, allele.start, allele.end),
                        covers,
                    ]
                )
    rows.sort(key=lambda row: (int(row[1]), str(row[2]), str(row[3]), str(row[5]), int(row[9])))
    return rows


def _overlaps(start: int, end: int, tss: int, window: int) -> bool:
    left = max(1, tss - window)
    right = tss + window
    return start <= right and end >= left


def _distance(strand: str, tss: int, start: int, end: int) -> int:
    if start <= tss <= end:
        return 0
    if strand == "-":
        if end < tss:
            return tss - end
        return tss - start
    if end < tss:
        return end - tss
    return start - tss


def _write_table(path: Path, rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(_FIELDS)
        writer.writerows(rows)
    temporary.replace(path)
