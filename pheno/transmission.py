"""Within-family ASD dosage contrast for TRIO_CALLS variants."""

from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path

_CHROM_DIR = re.compile(r"^chr(\d+|X|Y)$")
_DOSAGE = {
    "0/0": 0,
    "0|0": 0,
    "0/1": 1,
    "1/0": 1,
    "0|1": 1,
    "1|0": 1,
    "1/1": 2,
    "1|1": 2,
}
_FIELDS = [
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "n_families",
    "n_informative",
    "z",
    "p_exact",
]


def compute_transmission(
    input_dir: Path,
    prefix: str,
    people: dict[str, tuple[str, str]],
    output: Path | None = None,
    *,
    output_dir: Path | None = None,
) -> int:
    """Write one row for each variant seen in an ASD and unaffected sibship.

    ``output`` is one TSV. ``output_dir`` writes ``chrN/{prefix}.tsv``, with
    the rows of each chromosome ordered by position. Exactly one of the two
    paths is set. A child contributes only when that child is listed on the
    variant row. Dosage is 0, 1, or 2 for a biallelic genotype. Multiallelic
    variants are skipped. ``z`` is blank when no family is informative, and
    ``p_exact`` is then 1.
    """
    if (output is None) == (output_dir is None):
        raise ValueError("set exactly one of --output and --output-dir")
    if not input_dir.is_dir():
        raise ValueError(f"input directory not found: {input_dir}")
    if not prefix:
        raise ValueError("file pattern must not be empty")
    files = _discover(input_dir, prefix)
    for path in files:
        _require_trio_header(path)

    handle = None
    writer: csv.DictWriter | None = None
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        handle = output.open("w", encoding="utf-8", newline="")
        writer = csv.DictWriter(handle, fieldnames=_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
    else:
        assert output_dir is not None
        output_dir.mkdir(parents=True, exist_ok=True)

    written = 0
    seen: set[tuple[str, str, str, str]] = set()
    pending: list[dict[str, object]] = []
    current_chrom = ""
    try:
        for path in files:
            chrom_name = path.parent.name
            if chrom_name != current_chrom:
                written += _flush_chromosome(current_chrom, pending, prefix, writer, output_dir)
                if current_chrom:
                    print(f"finished processing {current_chrom}", flush=True)
                current_chrom = chrom_name
                seen.clear()
            for row in _iter_rows(path):
                if row is None:
                    continue
                chrom, pos, variant_id, ref, alt, calls = row
                key = (chrom, pos, ref, alt)
                if key in seen:
                    raise ValueError(f"duplicate variant {chrom} {pos} {ref} {alt}")
                seen.add(key)
                result = _variant_result(calls, people)
                if result is None:
                    continue
                n_families, n_informative, z, p_exact = result
                pending.append(
                    {
                        "CHROM": chrom,
                        "POS": pos,
                        "ID": variant_id,
                        "REF": ref,
                        "ALT": alt,
                        "n_families": n_families,
                        "n_informative": n_informative,
                        "z": "" if z is None else format(z, ".12g"),
                        "p_exact": format(p_exact, ".12g"),
                    }
                )
        written += _flush_chromosome(current_chrom, pending, prefix, writer, output_dir)
        if current_chrom:
            print(f"finished processing {current_chrom}", flush=True)
    finally:
        if handle is not None:
            handle.close()
    return written


def _flush_chromosome(
    chrom: str,
    pending: list[dict[str, object]],
    prefix: str,
    writer: csv.DictWriter | None,
    output_dir: Path | None,
) -> int:
    if not chrom:
        return 0
    pending.sort(key=lambda row: (int(str(row["POS"])), str(row["REF"]), str(row["ALT"]), str(row["ID"])))
    count = len(pending)
    if output_dir is not None:
        if pending:
            path = output_dir / chrom / f"{prefix}.tsv"
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("w", encoding="utf-8", newline="") as handle:
                chrom_writer = csv.DictWriter(
                    handle, fieldnames=_FIELDS, delimiter="\t", lineterminator="\n"
                )
                chrom_writer.writeheader()
                chrom_writer.writerows(pending)
    elif writer is not None:
        writer.writerows(pending)
    pending.clear()
    return count


def _variant_result(
    calls: list[tuple[str, int]],
    people: dict[str, tuple[str, str]],
) -> tuple[int, int, float | None, float] | None:
    grouped: dict[str, list[tuple[bool, int]]] = defaultdict(list)
    seen_people: set[str] = set()
    for person, dosage in calls:
        info = people.get(person)
        if info is None or person in seen_people:
            continue
        seen_people.add(person)
        status, family_id = info
        grouped[family_id].append((status == "asd", dosage))

    families = [family for records in grouped.values() if (family := _family_contrast(records)) is not None]
    if not families:
        return None
    informative = [family for family in families if family[2] is not None]
    total = sum((family[0] for family in families), Fraction(0))
    variance = sum((family[1] for family in families), Fraction(0))
    z = None if variance == 0 else float(total) / math.sqrt(float(variance))
    p_exact = 1.0 if not informative else float(_exact_tail(informative, total))
    return len(families), len(informative), z, p_exact


def _family_contrast(
    records: list[tuple[bool, int]],
) -> tuple[Fraction, Fraction, dict[Fraction, int] | None] | None:
    """Return observed difference, null variance, and the permutation weights."""
    asd = [dosage for is_asd, dosage in records if is_asd]
    unaffected = [dosage for is_asd, dosage in records if not is_asd]
    if not asd or not unaffected:
        return None
    dosages = [dosage for _is_asd, dosage in records]
    n_asd = len(asd)
    n_unaffected = len(unaffected)
    n = n_asd + n_unaffected
    observed = Fraction(sum(asd), n_asd) - Fraction(sum(unaffected), n_unaffected)
    mean = Fraction(sum(dosages), n)
    squares = sum((Fraction(dosage) - mean) ** 2 for dosage in dosages)
    variance = Fraction(n) * squares / (n_asd * n_unaffected * (n - 1))
    if squares == 0:
        return observed, variance, None
    weights: dict[Fraction, int] = {}
    for chosen in combinations(range(n), n_asd):
        selected = set(chosen)
        asd_sum = sum(dosage for index, dosage in enumerate(dosages) if index in selected)
        difference = Fraction(asd_sum, n_asd) - Fraction(sum(dosages) - asd_sum, n_unaffected)
        weights[difference] = weights.get(difference, 0) + 1
    return observed, variance, weights


def _exact_tail(
    families: list[tuple[Fraction, Fraction, dict[Fraction, int] | None]],
    observed: Fraction,
) -> Fraction:
    grouped: dict[tuple[tuple[Fraction, int], ...], int] = {}
    for _difference, _variance, weights in families:
        if weights is None:
            continue
        key = tuple(sorted(weights.items()))
        grouped[key] = grouped.get(key, 0) + 1
    distribution: dict[Fraction, int] = {Fraction(0): 1}
    for key, copies in grouped.items():
        distribution = _convolve(distribution, _scale(dict(key), copies))
    threshold = abs(observed)
    total = sum(distribution.values())
    tail = sum(weight for value, weight in distribution.items() if abs(value) >= threshold)
    return Fraction(tail, total)


def _convolve(left: dict[Fraction, int], right: dict[Fraction, int]) -> dict[Fraction, int]:
    combined: dict[Fraction, int] = {}
    for left_value, left_weight in left.items():
        for right_value, right_weight in right.items():
            total = left_value + right_value
            combined[total] = combined.get(total, 0) + left_weight * right_weight
    return _reduce(combined)


def _scale(weights: dict[Fraction, int], copies: int) -> dict[Fraction, int]:
    result: dict[Fraction, int] = {Fraction(0): 1}
    base = dict(weights)
    remaining = copies
    while remaining:
        if remaining & 1:
            result = _convolve(result, base)
        remaining >>= 1
        if remaining:
            base = _convolve(base, base)
    return result


def _reduce(weights: dict[Fraction, int]) -> dict[Fraction, int]:
    divisor = 0
    for weight in weights.values():
        divisor = math.gcd(divisor, weight)
    if divisor <= 1:
        return weights
    return {value: weight // divisor for value, weight in weights.items()}


def _discover(input_dir: Path, prefix: str) -> list[Path]:
    found: list[Path] = []
    for chrom_dir in _chrom_dirs(input_dir):
        matches = [path for path in chrom_dir.glob(f"{prefix}*.tsv") if path.is_file()]
        found.extend(sorted(matches, key=lambda path: path.name))
    return found


def _chrom_dirs(input_dir: Path) -> list[Path]:
    dirs = [path for path in input_dir.iterdir() if path.is_dir() and _CHROM_DIR.fullmatch(path.name)]

    def sort_key(path: Path) -> tuple[int, int | str]:
        suffix = path.name[3:]
        if suffix == "X":
            return (1, "X")
        if suffix == "Y":
            return (2, "Y")
        return (0, int(suffix))

    return sorted(dirs, key=sort_key)


def _require_trio_header(path: Path) -> None:
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if not line.startswith("#"):
                raise ValueError(f"{path} is missing a TRIO_CALLS header")
            if line.split("\t")[-1] != "TRIO_CALLS":
                raise ValueError(f"{path} requires a TRIO_CALLS column")
            return
    raise ValueError(f"{path} is empty")


def _iter_rows(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                raise ValueError(f"{path}:{line_number} expected 6 tab-separated columns")
            if "," in fields[3] or "," in fields[4]:
                yield None
                continue
            try:
                calls = _parse_calls(fields[5])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number} {exc}") from exc
            if calls is None:
                yield None
                continue
            yield fields[0], fields[1], fields[2], fields[3], fields[4], calls


def _parse_calls(payload: str) -> list[tuple[str, int]] | None:
    calls: list[tuple[str, int]] = []
    for part in payload.split(";"):
        if not part:
            continue
        if "=" not in part:
            raise ValueError(f"TRIO_CALLS entry {part!r} is missing '='")
        person, genotypes = part.split("=", 1)
        if not person:
            raise ValueError("TRIO_CALLS entry is missing a child id")
        fields = genotypes.split("|")
        if len(fields) < 3:
            raise ValueError(f"TRIO_CALLS entry {part!r} has no child genotype")
        dosage = _DOSAGE.get(fields[-2])
        if dosage is None:
            if fields[-2] in {".", "./.", ".|."}:
                continue
            return None
        calls.append((person, dosage))
    return calls
