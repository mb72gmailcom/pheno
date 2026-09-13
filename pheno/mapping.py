"""Column-name and value-alias mapping for family files."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


DEFAULT_COLUMNS = {
    "person": "person",
    "family": "family",
    "mother": "mother",
    "father": "father",
    "mother_sequenced": None,
    "father_sequenced": None,
    "sex": "sex",
    "asd": "asd",
}

DEFAULT_VALUES = {
    "bad_ids": ["0", "false", "-", ""],
    "sequenced_true": ["true", "1", "yes"],
    "sex_male": ["male", "1"],
    "sex_female": ["female", "2"],
    "asd_affected": ["true", "2"],
    "asd_unaffected": ["false", "1"],
}

REQUIRED_COLUMNS = ("person", "family", "mother", "father", "sex", "asd")
ABSENT_TOKENS = {"", "-", "none", "null"}


def _absent_column(value: Any) -> bool:
    if value is None:
        return True
    return str(value).strip().lower() in ABSENT_TOKENS


def _as_str_list(values: Any, fallback: list[str]) -> list[str]:
    if values is None:
        return list(fallback)
    if isinstance(values, str):
        return [values]
    return [str(v) for v in values]


def _norm_set(values: list[str]) -> set[str]:
    return {str(v).strip().lower() for v in values}


@dataclass(frozen=True)
class ColumnMap:
    person: str
    family: str
    mother: str
    father: str
    sex: str
    asd: str
    mother_sequenced: str | None = None
    father_sequenced: str | None = None
    bad_ids: set[str] = field(default_factory=set)
    sequenced_true: set[str] = field(default_factory=set)
    sex_male: set[str] = field(default_factory=set)
    sex_female: set[str] = field(default_factory=set)
    asd_affected: set[str] = field(default_factory=set)
    asd_unaffected: set[str] = field(default_factory=set)

    def is_bad_id(self, value: str) -> bool:
        stripped = value.strip()
        if stripped == "":
            return True
        return stripped.lower() in self.bad_ids

    def is_sequenced(self, value: str) -> bool:
        return value.strip().lower() in self.sequenced_true

    def asd_status(self, value: str) -> str | None:
        token = value.strip().lower()
        if token in self.asd_affected:
            return "asd"
        if token in self.asd_unaffected:
            return "ctrl"
        return None

    def sex_status(self, value: str) -> str | None:
        token = value.strip().lower()
        if token in self.sex_male:
            return "male"
        if token in self.sex_female:
            return "female"
        return None


def load_column_map(path: str | Path | None) -> ColumnMap:
    payload: dict[str, Any] = {}
    if path is not None:
        with Path(path).open() as handle:
            payload = json.load(handle)

    columns = {**DEFAULT_COLUMNS, **payload.get("columns", {})}
    values = {**DEFAULT_VALUES, **payload.get("values", {})}

    resolved: dict[str, str | None] = {}
    for name in DEFAULT_COLUMNS:
        raw = columns.get(name, DEFAULT_COLUMNS[name])
        resolved[name] = None if _absent_column(raw) else str(raw).strip()

    missing = [name for name in REQUIRED_COLUMNS if not resolved[name]]
    if missing:
        raise ValueError(f"column map is missing required columns: {missing}")

    return ColumnMap(
        person=resolved["person"],  # type: ignore[arg-type]
        family=resolved["family"],  # type: ignore[arg-type]
        mother=resolved["mother"],  # type: ignore[arg-type]
        father=resolved["father"],  # type: ignore[arg-type]
        sex=resolved["sex"],  # type: ignore[arg-type]
        asd=resolved["asd"],  # type: ignore[arg-type]
        mother_sequenced=resolved["mother_sequenced"],
        father_sequenced=resolved["father_sequenced"],
        bad_ids=_norm_set(_as_str_list(values.get("bad_ids"), DEFAULT_VALUES["bad_ids"])),
        sequenced_true=_norm_set(
            _as_str_list(values.get("sequenced_true"), DEFAULT_VALUES["sequenced_true"])
        ),
        sex_male=_norm_set(_as_str_list(values.get("sex_male"), DEFAULT_VALUES["sex_male"])),
        sex_female=_norm_set(_as_str_list(values.get("sex_female"), DEFAULT_VALUES["sex_female"])),
        asd_affected=_norm_set(
            _as_str_list(values.get("asd_affected"), DEFAULT_VALUES["asd_affected"])
        ),
        asd_unaffected=_norm_set(
            _as_str_list(values.get("asd_unaffected"), DEFAULT_VALUES["asd_unaffected"])
        ),
    )
