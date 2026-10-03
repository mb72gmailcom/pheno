"""Split variant TSVs into ASD-only, unaffected-only, and both-carrier files."""

from __future__ import annotations

import csv
import re
from pathlib import Path

from pheno.mapping import ColumnMap

_CHROM_DIR = re.compile(r"^chr(\d+|X|Y)$")
_SKIP_SUFFIXES = ("_asd", "_unaffected", "_both")
TSV_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"
_LABELS = ("asd", "unaffected", "both")


def split_directory(
    input_dir: Path,
    output_dir: Path,
    prefix: str,
    status: dict[str, str],
) -> None:
    """Write one TSV per class under ``output_dir/chrN/``.

    ``status`` maps a person id to ``asd`` or ``ctrl``. People missing from
    the map are dropped before the row is classified.
    """
    if not input_dir.is_dir():
        raise ValueError(f"input directory not found: {input_dir}")
    if not prefix:
        raise ValueError("file pattern must not be empty")

    for chrom_dir in _chrom_dirs(input_dir):
        for path in _variant_files(chrom_dir, prefix):
            groups = _classify_file(path, status)
            for label in _LABELS:
                rows = groups[label]
                if not rows:
                    continue
                destination = (
                    output_dir / chrom_dir.name / _output_name(path.stem, prefix, label)
                )
                destination.parent.mkdir(parents=True, exist_ok=True)
                with destination.open("w", encoding="utf-8") as handle:
                    handle.write(TSV_HEADER)
                    handle.writelines(rows)


def load_asd_status(family_file: Path, column_map: ColumnMap) -> dict[str, str]:
    """Map person id to ``asd`` or ``ctrl`` from the family file."""
    if not family_file.is_file():
        raise ValueError(f"family file not found: {family_file}")
    with family_file.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing = [
            name
            for name in (column_map.person, column_map.asd)
            if name not in fieldnames
        ]
        if missing:
            raise ValueError(f"family file missing columns: {missing}")
        status: dict[str, str] = {}
        for row in reader:
            person = (row.get(column_map.person) or "").strip()
            if not person or person in status:
                continue
            kind = column_map.asd_status(row.get(column_map.asd) or "")
            if kind is None:
                continue
            status[person] = kind
    if not status:
        raise ValueError("no people with known ASD status")
    return status


def _classify_file(path: Path, status: dict[str, str]) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {label: [] for label in _LABELS}
    for columns, patients in _iter_short_rows(path):
        asd = [person for person in patients if status.get(person) == "asd"]
        unaffected = [person for person in patients if status.get(person) == "ctrl"]
        if asd and unaffected:
            kept = [person for person in patients if status.get(person) in ("asd", "ctrl")]
            label = "both"
        elif asd:
            kept = asd
            label = "asd"
        elif unaffected:
            kept = unaffected
            label = "unaffected"
        else:
            continue
        site = "\t".join(columns)
        groups[label].append(f"{site}\t{';'.join(kept)}\n")
    return groups


def _iter_short_rows(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                if line.split("\t")[-1] == "TRIO_CALLS":
                    raise ValueError(f"{path} requires input in short format")
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                raise ValueError(f"{path} expected 6 tab-separated columns")
            patients_field = fields[5]
            if "=" in patients_field:
                raise ValueError(f"{path} requires input in short format")
            patients = [patient for patient in patients_field.split(";") if patient]
            yield fields[:5], patients


def _output_name(stem: str, prefix: str, label: str) -> str:
    return f"{prefix}_{label}{stem.removeprefix(prefix)}.tsv"


def _chrom_dirs(input_dir: Path) -> list[Path]:
    dirs = [
        path
        for path in input_dir.iterdir()
        if path.is_dir() and _CHROM_DIR.fullmatch(path.name)
    ]
    return sorted(dirs, key=lambda path: _chrom_sort_key(path.name))


def _chrom_sort_key(name: str) -> tuple[int, int | str]:
    suffix = name[3:]
    if suffix == "X":
        return (1, "X")
    if suffix == "Y":
        return (2, "Y")
    return (0, int(suffix))


def _variant_files(chrom_dir: Path, prefix: str) -> list[Path]:
    files: list[Path] = []
    for path in chrom_dir.glob(f"{prefix}*.tsv"):
        if not path.is_file() or not path.stem.startswith(prefix):
            continue
        rest = path.stem.removeprefix(prefix)
        if rest.startswith(_SKIP_SUFFIXES):
            continue
        files.append(path)
    return sorted(files, key=lambda path: path.name)
