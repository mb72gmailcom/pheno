"""Per-person Otari burden for variants carried in short-format TSVs."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from pheno.mapping import ColumnMap

_CHROM_DIR = re.compile(r"^chr(\d+|X|Y)$")
_OTARI_NAME = "variant_effects_comprehensive.tsv"
_BURDEN_NAME = "burden.json"


def compute_burden(
    input_dir: Path,
    otari_dir: Path,
    prefix: str,
    people: dict[str, tuple[str, str]],
    columns: list[str],
    *,
    threshold: float,
    use_abs: bool,
    transcripts: str,
    output_dir: Path | None = None,
    otari_prefix: str | None = None,
) -> dict[str, object]:
    """Burden for people who appear in ``prefix`` TSVs and have a known status.

    ``people`` maps a person id to ``(status, family_id)``. ``status`` is
    ``asd`` or ``unaffected``. Transcript scores are made absolute when
    ``use_abs`` is set, then collapsed with ``max`` or ``mean``.

    ``prefix`` selects ``{prefix}_{start}_{end}.tsv``. ``otari_prefix`` selects
    ``{otari-dir}/{otari_prefix}/chrN/{start}_{end}/``. When ``otari_prefix``
    is omitted, it is the same as ``prefix``.

    Each chromosome is scored on its own. When ``output_dir`` is set, that
    chromosome payload is written to ``{output_dir}/{chrom}/burden.json``
    before the next chromosome starts. The returned payload, also written to
    ``{output_dir}/burden.json``, is the sum of those chromosome files.
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

    payloads: list[dict[str, object]] = []
    for chrom_dir in _chrom_dirs(input_dir):
        patients: dict[str, dict[str, object]] = {}
        for path, start, end in _shard_files(chrom_dir, prefix):
            rows = list(_iter_variant_rows(path))
            needed = {_variant_id(chrom, pos, ref, alt) for chrom, pos, ref, alt, _patients in rows}
            otari_path = (
                otari_dir / otari_prefix / chrom_dir.name / f"{start}_{end}" / _OTARI_NAME
            )
            scores = _load_otari(otari_path, needed, columns)
            for chrom, pos, ref, alt, carriers in rows:
                _add_variant(
                    patients,
                    people,
                    columns,
                    scores.get(_variant_id(chrom, pos, ref, alt)),
                    carriers,
                    threshold=threshold,
                    use_abs=use_abs,
                    transcripts=transcripts,
                )
        payload = _payload(
            patients,
            threshold=threshold,
            use_abs=use_abs,
            transcripts=transcripts,
            columns=columns,
        )
        if output_dir is not None:
            write_burden(output_dir / chrom_dir.name / _BURDEN_NAME, payload)
        payloads.append(payload)
        print(f"finished processing {chrom_dir.name}", flush=True)

    total = _sum_burdens(
        payloads,
        columns,
        threshold=threshold,
        use_abs=use_abs,
        transcripts=transcripts,
    )
    if output_dir is not None:
        write_burden(output_dir / _BURDEN_NAME, total)
    return total


def load_people(family_file: Path, column_map: ColumnMap) -> dict[str, tuple[str, str]]:
    """Map person id to ``(status, family_id)`` for people with known ASD status."""
    if not family_file.is_file():
        raise ValueError(f"family file not found: {family_file}")
    with family_file.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing = [
            name
            for name in (column_map.person, column_map.family, column_map.asd)
            if name not in fieldnames
        ]
        if missing:
            raise ValueError(f"family file missing columns: {missing}")
        people: dict[str, tuple[str, str]] = {}
        for row in reader:
            person = (row.get(column_map.person) or "").strip()
            if not person or person in people:
                continue
            kind = column_map.asd_status(row.get(column_map.asd) or "")
            if kind is None:
                continue
            status = "asd" if kind == "asd" else "unaffected"
            family_id = (row.get(column_map.family) or "").strip()
            people[person] = (status, family_id)
    if not people:
        raise ValueError("no people with known ASD status")
    return people


def write_burden(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _payload(
    patients: dict[str, dict[str, object]],
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
        "patients": {person: patients[person] for person in sorted(patients)},
    }


def _sum_burdens(
    payloads: list[dict[str, object]],
    columns: list[str],
    *,
    threshold: float,
    use_abs: bool,
    transcripts: str,
) -> dict[str, object]:
    """Add per-chromosome burdens. ``fraction_damaging`` is recomputed."""
    patients: dict[str, dict[str, object]] = {}
    for payload in payloads:
        chrom_patients = payload["patients"]
        if not isinstance(chrom_patients, dict):
            raise ValueError("burden payload is missing patients")
        for person, record in chrom_patients.items():
            if not isinstance(record, dict):
                raise ValueError(f"burden record for {person} is not an object")
            total = patients.get(person)
            if total is None:
                total = {
                    "status": record["status"],
                    "family_id": record["family_id"],
                    "n_variants": 0,
                    "n_unscored": 0,
                }
                for column in columns:
                    total[column] = _empty_column()
                patients[person] = total
            total["n_variants"] = int(total["n_variants"]) + int(record["n_variants"])
            total["n_unscored"] = int(total["n_unscored"]) + int(record["n_unscored"])
            for column in columns:
                column_record = record[column]
                if not isinstance(column_record, dict):
                    raise ValueError(f"burden column {column} for {person} is not an object")
                _merge_column(total[column], column_record)
    return _payload(
        patients,
        threshold=threshold,
        use_abs=use_abs,
        transcripts=transcripts,
        columns=columns,
    )


def _merge_column(dest: dict[str, object], src: dict[str, object]) -> None:
    dest["n_scored"] = int(dest["n_scored"]) + int(src["n_scored"])
    dest["n_damaging"] = int(dest["n_damaging"]) + int(src["n_damaging"])
    dest["sum_effect"] = float(dest["sum_effect"]) + float(src["sum_effect"])
    dest["sum_damaging"] = float(dest["sum_damaging"]) + float(src["sum_damaging"])
    src_max = src["max_score"]
    if src_max is not None:
        current = dest["max_score"]
        if current is None or float(src_max) > float(current):
            dest["max_score"] = float(src_max)
    scored = int(dest["n_scored"])
    dest["fraction_damaging"] = int(dest["n_damaging"]) / scored if scored else None


def _add_variant(
    patients: dict[str, dict[str, object]],
    people: dict[str, tuple[str, str]],
    columns: list[str],
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
    for person in dict.fromkeys(carriers):
        info = people.get(person)
        if info is None:
            continue
        record = patients.get(person)
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
            patients[person] = record
        record["n_variants"] = int(record["n_variants"]) + 1
        if collapsed is None:
            record["n_unscored"] = int(record["n_unscored"]) + 1
            continue
        for column, score in collapsed.items():
            _add_score(record[column], score, threshold)


def _add_score(stats: dict[str, object], score: float, threshold: float) -> None:
    stats["n_scored"] = int(stats["n_scored"]) + 1
    stats["sum_effect"] = float(stats["sum_effect"]) + score
    current = stats["max_score"]
    if current is None or score > float(current):
        stats["max_score"] = score
    if score > threshold:
        stats["n_damaging"] = int(stats["n_damaging"]) + 1
        stats["sum_damaging"] = float(stats["sum_damaging"]) + score
    scored = int(stats["n_scored"])
    stats["fraction_damaging"] = int(stats["n_damaging"]) / scored if scored else None


def _empty_column() -> dict[str, object]:
    return {
        "n_scored": 0,
        "n_damaging": 0,
        "fraction_damaging": None,
        "sum_effect": 0.0,
        "sum_damaging": 0.0,
        "max_score": None,
    }


def _collapse(values: list[float], *, use_abs: bool, transcripts: str) -> float:
    scores = [abs(value) for value in values] if use_abs else list(values)
    if transcripts == "max":
        return max(scores)
    return sum(scores) / len(scores)


def _variant_id(chrom: str, pos: str, ref: str, alt: str) -> str:
    contig = chrom[3:] if chrom.startswith("chr") else chrom
    return f"{contig}_{pos}_{ref}_{alt}_hg38"


def _load_otari(
    path: Path,
    needed: set[str],
    columns: list[str],
) -> dict[str, dict[str, list[float]]]:
    scores: dict[str, dict[str, list[float]]] = {}
    if not path.is_file():
        return scores
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing = [name for name in ("variant_id", *columns) if name not in fieldnames]
        if missing:
            raise ValueError(f"{path} is missing columns: {missing}")
        for row in reader:
            variant_id = (row.get("variant_id") or "").strip()
            if variant_id not in needed:
                continue
            bucket = scores.setdefault(variant_id, {column: [] for column in columns})
            for column in columns:
                try:
                    bucket[column].append(float(row[column]))
                except ValueError as exc:
                    raise ValueError(
                        f"{path} has a non-numeric {column} for {variant_id}"
                    ) from exc
    return scores


def _iter_variant_rows(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                if line.startswith("#") and line.split("\t")[-1] == "TRIO_CALLS":
                    raise ValueError(f"{path} requires input in short format")
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                raise ValueError(f"{path} expected 6 tab-separated columns")
            if "=" in fields[5]:
                raise ValueError(f"{path} requires input in short format")
            patients = [patient for patient in fields[5].split(";") if patient]
            yield fields[0], fields[1], fields[3], fields[4], patients


def _chrom_dirs(input_dir: Path) -> list[Path]:
    dirs = [
        path
        for path in input_dir.iterdir()
        if path.is_dir() and _CHROM_DIR.fullmatch(path.name)
    ]

    def sort_key(path: Path) -> tuple[int, int | str]:
        suffix = path.name[3:]
        if suffix == "X":
            return (1, "X")
        if suffix == "Y":
            return (2, "Y")
        return (0, int(suffix))

    return sorted(dirs, key=sort_key)


def _shard_files(chrom_dir: Path, prefix: str) -> list[tuple[Path, str, str]]:
    pattern = re.compile(rf"^{re.escape(prefix)}_(\d+)_(\d+)\.tsv$")
    found: list[tuple[Path, str, str]] = []
    for path in chrom_dir.iterdir():
        match = pattern.fullmatch(path.name)
        if path.is_file() and match:
            found.append((path, match.group(1), match.group(2)))
    return sorted(found, key=lambda item: item[0].name)
