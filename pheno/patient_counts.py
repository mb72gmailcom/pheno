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

    Only people who appear in those TSVs and have a known ASD status are
    included. ``status`` maps a person id to ``asd`` or ``ctrl``.
    """
    if not input_dir.is_dir():
        raise ValueError(f"input directory not found: {input_dir}")
    if not prefix:
        raise ValueError("file pattern must not be empty")

    asd: dict[str, int] = {}
    unaffected: dict[str, int] = {}
    buckets = {"asd": asd, "ctrl": unaffected}

    for chrom_dir in _chrom_dirs(input_dir):
        for path in _variant_files(chrom_dir, prefix):
            for _columns, patients in _iter_short_rows(path):
                for person in set(patients):
                    kind = status.get(person)
                    bucket = buckets.get(kind) if kind is not None else None
                    if bucket is None:
                        continue
                    bucket[person] = bucket.get(person, 0) + 1

    return {
        "asd": dict(sorted(asd.items())),
        "unaffected": dict(sorted(unaffected.items())),
    }


def write_patient_counts(path: Path, counts: dict[str, dict[str, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(counts, handle, indent=2, sort_keys=True)
        handle.write("\n")
