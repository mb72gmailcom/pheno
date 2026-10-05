"""Per-gene Otari burden for ASD children and their unaffected siblings."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

from pheno.burden import (
    _add_score,
    _chrom_dirs,
    _collapse,
    _empty_column,
    _iter_variant_rows,
    _merge_column,
    _shard_files,
    _variant_id,
    write_burden,
)
from pheno.mapping import ColumnMap

_OTARI_NAME = "variant_effects_comprehensive.tsv"
_GENE_MAP_NAME = "interpretability_analysis.tsv"
_GENES_NAME = "genes.json"
_SUMMARY_NAME = "gene_summary.json"
_WINDOW = 2000


def load_analysis_groups(
    family_file: Path,
    column_map: ColumnMap,
) -> tuple[dict[str, tuple[str, str]], int, int]:
    """Return ASD children and unaffected siblings, plus the two group sizes.

    A child has both parents recorded. An unaffected sibling is an unaffected
    child who shares a family id with an ASD child. Parents are omitted.
    Children who carry nothing in the variant TSVs are still included.
    """
    if not family_file.is_file():
        raise ValueError(f"family file not found: {family_file}")
    with family_file.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        required = [column_map.person, column_map.family, column_map.mother, column_map.father, column_map.asd]
        if column_map.mother_sequenced:
            required.append(column_map.mother_sequenced)
        if column_map.father_sequenced:
            required.append(column_map.father_sequenced)
        missing = [name for name in required if name not in fieldnames]
        if missing:
            raise ValueError(f"family file missing columns: {missing}")
        children: dict[str, tuple[str, str]] = {}
        for row in reader:
            person = (row.get(column_map.person) or "").strip()
            if not person or person in children:
                continue
            if not _is_child(row, column_map):
                continue
            kind = column_map.asd_status(row.get(column_map.asd) or "")
            if kind is None:
                continue
            status = "asd" if kind == "asd" else "unaffected"
            family_id = (row.get(column_map.family) or "").strip()
            children[person] = (status, family_id)
    asd_families = {family_id for status, family_id in children.values() if status == "asd"}
    people = {
        person: info
        for person, info in children.items()
        if info[0] == "asd" or (info[0] == "unaffected" and info[1] in asd_families)
    }
    if not any(status == "asd" for status, _family_id in people.values()):
        raise ValueError("no ASD children in the family file")
    n_asd = sum(status == "asd" for status, _family_id in people.values())
    n_siblings = sum(status == "unaffected" for status, _family_id in people.values())
    return people, n_asd, n_siblings


def compute_gene_burden(
    input_dir: Path,
    otari_dir: Path,
    annotation: Path,
    prefix: str,
    people: dict[str, tuple[str, str]],
    columns: list[str],
    *,
    threshold: float,
    use_abs: bool,
    transcripts: str,
    n_asd: int,
    n_unaffected_sibling: int,
    output_dir: Path | None = None,
    otari_prefix: str | None = None,
) -> tuple[dict[str, object], dict[str, object]]:
    """Per-gene burden for ``people``, then the cohort summary.

    A variant is assigned to every gene whose body, or the 2000 bp around
    either end, contains it. Transcript scores from the Otari shard are
    collapsed within each of those genes. ``genes.json`` keeps only people
    who carry a variant in that gene. ``gene_summary.json`` divides by every
    ASD child and every unaffected sibling, counting non-carriers as zero.
    """
    if transcripts not in ("max", "mean"):
        raise ValueError("transcripts must be 'max' or 'mean'")
    if not input_dir.is_dir():
        raise ValueError(f"input directory not found: {input_dir}")
    if not otari_dir.is_dir():
        raise ValueError(f"Otari directory not found: {otari_dir}")
    if not prefix:
        raise ValueError("file pattern must not be empty")
    if otari_prefix is None:
        otari_prefix = prefix
    if not otari_prefix:
        raise ValueError("Otari prefix must not be empty")
    if not columns:
        raise ValueError("at least one Otari column is required")

    gene_index = load_gene_index(annotation)
    payloads: list[dict[str, object]] = []
    for chrom_dir in _chrom_dirs(input_dir):
        genes: dict[str, dict[str, dict[str, object]]] = {}
        for path, start, end in _shard_files(chrom_dir, prefix):
            rows = list(_iter_variant_rows(path))
            needed = {_variant_id(chrom, pos, ref, alt) for chrom, pos, ref, alt, _patients in rows}
            otari_path = otari_dir / otari_prefix / chrom_dir.name / f"{start}_{end}"
            scores = _load_gene_scores(otari_path, needed, columns)
            assigned = _assign_rows(rows, gene_index)
            for (chrom, pos, ref, alt, carriers), gene_names in zip(rows, assigned):
                variant_id = _variant_id(chrom, pos, ref, alt)
                for gene_name in gene_names:
                    _add_gene_variant(
                        genes,
                        people,
                        columns,
                        gene_name,
                        scores.get((variant_id, gene_name)),
                        carriers,
                        threshold=threshold,
                        use_abs=use_abs,
                        transcripts=transcripts,
                    )
        payload = _genes_payload(
            genes,
            threshold=threshold,
            use_abs=use_abs,
            transcripts=transcripts,
            columns=columns,
        )
        if output_dir is not None:
            write_burden(output_dir / chrom_dir.name / _GENES_NAME, payload)
        payloads.append(payload)
        print(f"finished processing {chrom_dir.name}", flush=True)

    merged = _merge_gene_payloads(
        payloads,
        columns,
        threshold=threshold,
        use_abs=use_abs,
        transcripts=transcripts,
    )
    summary = _summarize(
        merged,
        columns,
        n_asd=n_asd,
        n_unaffected_sibling=n_unaffected_sibling,
    )
    if output_dir is not None:
        write_burden(output_dir / _GENES_NAME, merged)
        write_burden(output_dir / _SUMMARY_NAME, summary)
    return merged, summary


def load_gene_index(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    if not path.is_file():
        raise ValueError(f"annotation not found: {path}")
    opener = gzip.open if path.suffix == ".gz" else open
    indexed: dict[str, list[tuple[int, int, str]]] = {}
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        required = ("feature", "chr", "start", "end", "strand", "name")
        absent = [name for name in required if name not in fieldnames]
        if absent:
            raise ValueError(f"{path} is missing columns: {absent}")
        for line_number, row in enumerate(reader, start=2):
            if (row.get("feature") or "").strip() != "gene":
                continue
            name = (row.get("name") or "").strip()
            chrom = (row.get("chr") or "").strip()
            if not name or not chrom:
                raise ValueError(f"{path}:{line_number} gene row is missing chr or name")
            try:
                start = int(row["start"])
                end = int(row["end"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{line_number} has a non-integer coordinate") from exc
            strand = (row.get("strand") or "").strip()
            for key in _chrom_aliases(chrom):
                bucket = indexed.setdefault(key, [])
                for lo, hi in _gene_intervals(start, end, strand):
                    bucket.append((lo, hi, name))
    if not indexed:
        raise ValueError(f"no gene rows in {path}")
    for intervals in indexed.values():
        intervals.sort()
    return indexed


def _is_child(row: dict[str, str], column_map: ColumnMap) -> bool:
    mother = row.get(column_map.mother) or ""
    father = row.get(column_map.father) or ""
    if column_map.is_bad_id(mother) or column_map.is_bad_id(father):
        return False
    if column_map.mother_sequenced and not column_map.is_sequenced(row.get(column_map.mother_sequenced) or ""):
        return False
    if column_map.father_sequenced and not column_map.is_sequenced(row.get(column_map.father_sequenced) or ""):
        return False
    return True


def _chrom_aliases(chrom: str) -> list[str]:
    aliases = [chrom]
    other = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    if other not in aliases:
        aliases.append(other)
    return aliases


def _gene_intervals(start: int, end: int, strand: str) -> list[tuple[int, int]]:
    tss = start if strand == "+" else end
    tel = end if strand == "+" else start
    spans = [(tss - _WINDOW, tss + _WINDOW), (tel - _WINDOW, tel + _WINDOW)]
    if start <= end:
        spans.append((start, end))
    spans.sort()
    merged: list[list[int]] = []
    for lo, hi in spans:
        if hi < lo:
            continue
        if merged and lo <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [(lo, hi) for lo, hi in merged]


def _assign_rows(
    rows: list[tuple[str, str, str, str, list[str]]],
    gene_index: dict[str, list[tuple[int, int, str]]],
) -> list[list[str]]:
    positions = [int(pos) for _chrom, pos, _ref, _alt, _carriers in rows]
    assigned: list[list[str]] = [[] for _row in rows]
    by_chrom: dict[str, list[int]] = {}
    for index, (chrom, _pos, _ref, _alt, _carriers) in enumerate(rows):
        by_chrom.setdefault(chrom, []).append(index)
    for chrom, indexes in by_chrom.items():
        intervals = _intervals_for(chrom, gene_index)
        if not intervals:
            continue
        names = _genes_at([positions[index] for index in indexes], intervals)
        for index, gene_names in zip(indexes, names):
            assigned[index] = gene_names
    return assigned


def _intervals_for(chrom: str, gene_index: dict[str, list[tuple[int, int, str]]]) -> list[tuple[int, int, str]]:
    for key in _chrom_aliases(chrom):
        intervals = gene_index.get(key)
        if intervals:
            return intervals
    return []


def _genes_at(positions: list[int], intervals: list[tuple[int, int, str]]) -> list[list[str]]:
    order = sorted(range(len(positions)), key=lambda index: positions[index])
    found: list[list[str]] = [[] for _position in positions]
    active: list[tuple[int, int, str]] = []
    cursor = 0
    for index in order:
        pos = positions[index]
        while cursor < len(intervals) and intervals[cursor][0] <= pos:
            active.append(intervals[cursor])
            cursor += 1
        active = [interval for interval in active if interval[1] >= pos]
        names: list[str] = []
        seen: set[str] = set()
        for _lo, _hi, name in active:
            if name in seen:
                continue
            seen.add(name)
            names.append(name)
        found[index] = names
    return found


def _load_gene_scores(
    otari_dir: Path,
    needed: set[str],
    columns: list[str],
) -> dict[tuple[str, str], dict[str, list[float]]]:
    comprehensive = otari_dir / _OTARI_NAME
    gene_map = otari_dir / _GENE_MAP_NAME
    scores: dict[tuple[str, str], dict[str, list[float]]] = {}
    if not comprehensive.is_file() or not gene_map.is_file():
        return scores
    transcript_gene = _load_transcript_genes(gene_map, needed)
    with comprehensive.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing = [name for name in ("variant_id", "transcript_id", *columns) if name not in fieldnames]
        if missing:
            raise ValueError(f"{comprehensive} is missing columns: {missing}")
        for row in reader:
            variant_id = (row.get("variant_id") or "").strip()
            if variant_id not in needed:
                continue
            transcript_id = (row.get("transcript_id") or "").strip()
            gene_name = transcript_gene.get((variant_id, transcript_id))
            if not gene_name:
                continue
            bucket = scores.setdefault((variant_id, gene_name), {column: [] for column in columns})
            for column in columns:
                try:
                    bucket[column].append(float(row[column]))
                except ValueError as exc:
                    raise ValueError(
                        f"{comprehensive} has a non-numeric {column} for {variant_id}"
                    ) from exc
    return scores


def _load_transcript_genes(path: Path, needed: set[str]) -> dict[tuple[str, str], str]:
    mapping: dict[tuple[str, str], str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing = [name for name in ("variant_id", "gene_id", "transcript_id") if name not in fieldnames]
        if missing:
            raise ValueError(f"{path} is missing columns: {missing}")
        for row in reader:
            variant_id = (row.get("variant_id") or "").strip()
            if variant_id not in needed:
                continue
            transcript_id = (row.get("transcript_id") or "").strip()
            gene_name = (row.get("gene_id") or "").strip()
            if transcript_id and gene_name:
                mapping[(variant_id, transcript_id)] = gene_name
    return mapping


def _add_gene_variant(
    genes: dict[str, dict[str, dict[str, object]]],
    people: dict[str, tuple[str, str]],
    columns: list[str],
    gene_name: str,
    transcript_scores: dict[str, list[float]] | None,
    carriers: list[str],
    *,
    threshold: float,
    use_abs: bool,
    transcripts: str,
) -> None:
    collapsed: dict[str, float] | None = None
    if transcript_scores:
        collapsed = {
            column: _collapse(transcript_scores[column], use_abs=use_abs, transcripts=transcripts)
            for column in columns
        }
    bucket = genes.setdefault(gene_name, {})
    for person in dict.fromkeys(carriers):
        info = people.get(person)
        if info is None:
            continue
        record = bucket.get(person)
        if record is None:
            status, family_id = info
            record = {
                "status": status,
                "family_id": family_id,
                "n_variants": 0,
                "n_unscored": 0,
            }
            for column in columns:
                record[column] = _empty_column()
            bucket[person] = record
        record["n_variants"] = int(record["n_variants"]) + 1
        if collapsed is None:
            record["n_unscored"] = int(record["n_unscored"]) + 1
            continue
        for column, score in collapsed.items():
            _add_score(record[column], score, threshold)


def _genes_payload(
    genes: dict[str, dict[str, dict[str, object]]],
    *,
    threshold: float,
    use_abs: bool,
    transcripts: str,
    columns: list[str],
) -> dict[str, object]:
    return {
        "threshold": threshold,
        "abs": use_abs,
        "transcripts": transcripts,
        "columns": list(columns),
        "genes": {
            gene: {person: genes[gene][person] for person in sorted(genes[gene])}
            for gene in sorted(genes)
        },
    }


def _merge_gene_payloads(
    payloads: list[dict[str, object]],
    columns: list[str],
    *,
    threshold: float,
    use_abs: bool,
    transcripts: str,
) -> dict[str, object]:
    genes: dict[str, dict[str, dict[str, object]]] = {}
    for payload in payloads:
        chrom_genes = payload["genes"]
        if not isinstance(chrom_genes, dict):
            raise ValueError("gene burden payload is missing genes")
        for gene, people in chrom_genes.items():
            if not isinstance(people, dict):
                raise ValueError(f"gene burden record for {gene} is not an object")
            bucket = genes.setdefault(str(gene), {})
            for person, record in people.items():
                if not isinstance(record, dict):
                    raise ValueError(f"gene burden record for {person} is not an object")
                current = bucket.get(person)
                if current is None:
                    copied = {
                        "status": record["status"],
                        "family_id": record["family_id"],
                        "n_variants": int(record["n_variants"]),
                        "n_unscored": int(record["n_unscored"]),
                    }
                    for column in columns:
                        column_record = record[column]
                        if not isinstance(column_record, dict):
                            raise ValueError(f"burden column {column} for {person} is not an object")
                        copied[column] = dict(column_record)
                    bucket[str(person)] = copied
                    continue
                current["n_variants"] = int(current["n_variants"]) + int(record["n_variants"])
                current["n_unscored"] = int(current["n_unscored"]) + int(record["n_unscored"])
                for column in columns:
                    column_record = record[column]
                    if not isinstance(column_record, dict):
                        raise ValueError(f"burden column {column} for {person} is not an object")
                    _merge_column(current[column], column_record)
    return _genes_payload(
        genes,
        threshold=threshold,
        use_abs=use_abs,
        transcripts=transcripts,
        columns=columns,
    )


def _summarize(
    payload: dict[str, object],
    columns: list[str],
    *,
    n_asd: int,
    n_unaffected_sibling: int,
) -> dict[str, object]:
    chrom_genes = payload["genes"]
    if not isinstance(chrom_genes, dict):
        raise ValueError("gene burden payload is missing genes")
    genes: dict[str, object] = {}
    for gene, carriers in chrom_genes.items():
        if not isinstance(carriers, dict):
            raise ValueError(f"gene burden record for {gene} is not an object")
        genes[str(gene)] = {
            "asd": _group_summary(carriers, "asd", n_asd, columns),
            "unaffected_sibling": _group_summary(carriers, "unaffected", n_unaffected_sibling, columns),
        }
    return {
        "threshold": payload["threshold"],
        "abs": payload["abs"],
        "transcripts": payload["transcripts"],
        "columns": payload["columns"],
        "n_asd": n_asd,
        "n_unaffected_sibling": n_unaffected_sibling,
        "genes": genes,
    }


def _group_summary(
    carriers: dict[str, object],
    status: str,
    n_people: int,
    columns: list[str],
) -> dict[str, object]:
    totals = {column: {"sum_effect": 0.0, "n_scored": 0, "n_damaging": 0} for column in columns}
    n_carriers = 0
    for record in carriers.values():
        if not isinstance(record, dict) or record.get("status") != status:
            continue
        n_carriers += 1
        for column in columns:
            stats = record[column]
            if not isinstance(stats, dict):
                raise ValueError(f"burden column {column} is not an object")
            bucket = totals[column]
            bucket["sum_effect"] += float(stats["sum_effect"])
            bucket["n_scored"] += int(stats["n_scored"])
            bucket["n_damaging"] += int(stats["n_damaging"])
    summary: dict[str, object] = {"n_people": n_people, "n_carriers": n_carriers}
    for column in columns:
        bucket = totals[column]
        scored = bucket["n_scored"]
        summary[column] = {
            "mean_sum_effect": (bucket["sum_effect"] / n_people) if n_people else None,
            "n_scored": scored,
            "n_damaging": bucket["n_damaging"],
            "fraction_damaging": (bucket["n_damaging"] / scored) if scored else None,
        }
    return summary
