"""Count variants per ASD and unaffected patient."""

from __future__ import annotations

import json
from pathlib import Path

from pheno.split import _chrom_dirs, _iter_short_rows, _variant_files


def count_patient_variants(
    input_dir: Path,
    prefix: str,
    status: dict[str, str],
) -> dict[str, dict[str, int]]:
    """Count rows of ``prefix`` in which each person appears.

    Every person in ``status`` is included. People with no rows have count 0.
    ``status`` maps a person id to ``asd`` or ``ctrl``.
    """
    if not input_dir.is_dir():
        raise ValueError(f"input directory not found: {input_dir}")
    if not prefix:
        raise ValueError("file pattern must not be empty")

    asd = {person: 0 for person, kind in status.items() if kind == "asd"}
    unaffected = {person: 0 for person, kind in status.items() if kind == "ctrl"}
    buckets = {person: asd for person in asd}
    buckets.update({person: unaffected for person in unaffected})

    for chrom_dir in _chrom_dirs(input_dir):
        for path in _variant_files(chrom_dir, prefix):
            for _columns, patients in _iter_short_rows(path):
                for person in set(patients):
                    bucket = buckets.get(person)
                    if bucket is not None:
                        bucket[person] += 1

    return {
        "asd": dict(sorted(asd.items())),
        "unaffected": dict(sorted(unaffected.items())),
    }


def write_patient_counts(path: Path, counts: dict[str, dict[str, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(counts, handle, indent=2, sort_keys=True)
        handle.write("\n")
