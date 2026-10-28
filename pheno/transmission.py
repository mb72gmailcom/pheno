"""Within-family ASD dosage contrast for TRIO_CALLS variants."""

from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path

from scipy.stats import norm

from pheno.gene_burden import _genes_at, _intervals_for, load_gene_index
from pheno.transmission_scores import (
    ScoreSource,
    gene_score_rows,
    gene_scores_path,
    load_chromosome_scores,
    prepare_score_source,
    score_fieldnames,
)

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
_GENE_FIELDS = [
    "CHROM",
    "gene",
    "gene_name",
    "gene_type",
    "n_variants",
    "n_families",
    "n_informative",
    "z",
    "p_exact",
]
_DEFAULT_WINDOW = 2000


def gene_table_path(output: Path) -> Path:
    """Path of the gene table written beside a single variant ``output`` file."""
    suffix = output.suffix or ".tsv"
    return output.with_name(f"{output.stem}.genes{suffix}")


def compute_transmission(
    input_dir: Path,
    prefix: str,
    people: dict[str, tuple[str, str]],
    output: Path | None = None,
    *,
    output_dir: Path | None = None,
    genes: Path | None = None,
    window: int = _DEFAULT_WINDOW,
    fcap: int | None = None,
    otari_dir: Path | None = None,
    seqweaver_dir: Path | None = None,
    column_names: Path | None = None,
) -> tuple[int, int]:
    """Write one row for each variant seen in an ASD and unaffected sibship.

    ``output`` is one TSV. ``output_dir`` writes ``chrN/{prefix}.tsv``, with
    the rows of each chromosome ordered by position. Exactly one of the two
    paths is set. A sibling in a touched family who is absent from the row is
    dosage 0. Dosage is 0, 1, or 2 for a biallelic genotype. Multiallelic
    variants are skipped. ``z`` is blank when no family is informative, and
    ``p_exact`` is then 1. When ``genes`` is set, each gene sums those dosages
    and is tested with one label shuffle per family. ``fcap`` keeps the exact
    p-value for that many informative families or fewer, and uses the
    two-sided normal tail of ``z`` above it. ``otari_dir`` or ``seqweaver_dir``
    replaces that gene table with a wide score-weighted table. Its ``p`` is
    always the normal tail.
    """
    if (output is None) == (output_dir is None):
        raise ValueError("set exactly one of --output and --output-dir")
    if fcap is not None and fcap < 0:
        raise ValueError("--fcap must be at least 0")
    score_source = prepare_score_source(
        otari_dir,
        seqweaver_dir,
        column_names,
        genes=genes,
    )
    if not input_dir.is_dir():
        raise ValueError(f"input directory not found: {input_dir}")
    if not prefix:
        raise ValueError("file pattern must not be empty")
    files = _discover(input_dir, prefix)
    for path in files:
        _require_trio_header(path)
    gene_index = gene_names = gene_types = None
    if genes is not None:
        gene_index, gene_names, gene_types = load_gene_index(genes, window)

    handle = None
    gene_handle = None
    score_handle = None
    writer: csv.DictWriter | None = None
    gene_writer: csv.DictWriter | None = None
    score_writer: csv.DictWriter | None = None
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        handle = output.open("w", encoding="utf-8", newline="")
        writer = csv.DictWriter(handle, fieldnames=_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        if gene_index is not None and score_source is None:
            gene_handle = gene_table_path(output).open("w", encoding="utf-8", newline="")
            gene_writer = csv.DictWriter(
                gene_handle, fieldnames=_GENE_FIELDS, delimiter="\t", lineterminator="\n"
            )
            gene_writer.writeheader()
        if score_source is not None:
            score_handle = gene_scores_path(output).open("w", encoding="utf-8", newline="")
            score_writer = csv.DictWriter(
                score_handle,
                fieldnames=score_fieldnames(score_source.columns),
                delimiter="\t",
                lineterminator="\n",
            )
            score_writer.writeheader()
    else:
        assert output_dir is not None
        output_dir.mkdir(parents=True, exist_ok=True)

    written = 0
    gene_written = 0
    seen: set[tuple[str, str, str, str]] = set()
    pending: list[dict[str, object]] = []
    gene_records: list[tuple[int, str, str, str, dict[str, int]]] = []
    current_chrom = ""
    members = _members_by_family(people)
    try:
        for path in files:
            chrom_name = path.parent.name
            if chrom_name != current_chrom:
                variants, gene_rows = _flush_chromosome(
                    current_chrom,
                    pending,
                    gene_records,
                    prefix,
                    writer,
                    gene_writer,
                    output_dir,
                    gene_index,
                    gene_names or {},
                    gene_types or {},
                    people,
                    members,
                    fcap,
                    score_source,
                    score_writer,
                )
                written += variants
                gene_written += gene_rows
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
                listed = _cohort_dosages(calls, people)
                if gene_index is not None and listed:
                    gene_records.append((int(pos), pos, ref, alt, listed))
                result = _contrast(listed, people, members, fcap)
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
        variants, gene_rows = _flush_chromosome(
            current_chrom,
            pending,
            gene_records,
            prefix,
            writer,
            gene_writer,
            output_dir,
            gene_index,
            gene_names or {},
            gene_types or {},
            people,
            members,
            fcap,
            score_source,
            score_writer,
        )
        written += variants
        gene_written += gene_rows
        if current_chrom:
            print(f"finished processing {current_chrom}", flush=True)
    finally:
        if handle is not None:
            handle.close()
        if gene_handle is not None:
            gene_handle.close()
        if score_handle is not None:
            score_handle.close()
    return written, gene_written


def _flush_chromosome(
    chrom: str,
    pending: list[dict[str, object]],
    gene_records: list[tuple[int, str, str, str, dict[str, int]]],
    prefix: str,
    writer: csv.DictWriter | None,
    gene_writer: csv.DictWriter | None,
    output_dir: Path | None,
    gene_index: dict[str, list[tuple[int, int, str]]] | None,
    gene_names: dict[str, str],
    gene_types: dict[str, str],
    people: dict[str, tuple[str, str]],
    members: dict[str, list[tuple[str, bool]]],
    fcap: int | None,
    score_source: ScoreSource | None = None,
    score_writer: csv.DictWriter | None = None,
) -> tuple[int, int]:
    if not chrom:
        return 0, 0
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
    gene_count = 0
    if gene_index is not None and score_source is not None:
        loaded = load_chromosome_scores(score_source, chrom, gene_records)
        gene_rows = gene_score_rows(
            chrom,
            gene_records,
            gene_index,
            gene_names,
            gene_types,
            people,
            members,
            score_source,
            loaded,
        )
        gene_count = len(gene_rows)
        fields = score_fieldnames(score_source.columns)
        if output_dir is not None:
            if gene_rows:
                path = output_dir / chrom / f"{prefix}.gene_scores.tsv"
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("w", encoding="utf-8", newline="") as handle:
                    chrom_writer = csv.DictWriter(
                        handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                    )
                    chrom_writer.writeheader()
                    chrom_writer.writerows(gene_rows)
        elif score_writer is not None:
            score_writer.writerows(gene_rows)
    elif gene_index is not None:
        gene_rows = _gene_rows(
            chrom, gene_records, gene_index, gene_names, gene_types, people, members, fcap
        )
        gene_count = len(gene_rows)
        if output_dir is not None:
            if gene_rows:
                path = output_dir / chrom / f"{prefix}.genes.tsv"
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("w", encoding="utf-8", newline="") as handle:
                    chrom_writer = csv.DictWriter(
                        handle, fieldnames=_GENE_FIELDS, delimiter="\t", lineterminator="\n"
                    )
                    chrom_writer.writeheader()
                    chrom_writer.writerows(gene_rows)
        elif gene_writer is not None:
            gene_writer.writerows(gene_rows)
    gene_records.clear()
    return count, gene_count


def _members_by_family(
    people: dict[str, tuple[str, str]],
) -> dict[str, list[tuple[str, bool]]]:
    members: dict[str, list[tuple[str, bool]]] = defaultdict(list)
    for person, (status, family_id) in people.items():
        members[family_id].append((person, status == "asd"))
    return members


def _cohort_dosages(
    calls: list[tuple[str, int]],
    people: dict[str, tuple[str, str]],
) -> dict[str, int]:
    listed: dict[str, int] = {}
    for person, dosage in calls:
        if person not in people or person in listed:
            continue
        listed[person] = dosage
    return listed


def _contrast(
    listed: dict[str, int],
    people: dict[str, tuple[str, str]],
    members: dict[str, list[tuple[str, bool]]],
    fcap: int | None = None,
) -> tuple[int, int, float | None, float] | None:
    if not listed:
        return None
    families = []
    for family_id in {people[person][1] for person in listed}:
        records = [(is_asd, listed.get(person, 0)) for person, is_asd in members[family_id]]
        family = _family_contrast(records)
        if family is not None:
            families.append(family)
    if not families:
        return None
    informative = [family for family in families if family[2] is not None]
    total = sum((family[0] for family in families), Fraction(0))
    variance = sum((family[1] for family in families), Fraction(0))
    z = None if variance == 0 else float(total) / math.sqrt(float(variance))
    if not informative or z is None:
        p_exact = 1.0
    elif fcap is not None and len(informative) > fcap:
        p_exact = float(2.0 * norm.sf(abs(z)))
    else:
        p_exact = float(_exact_tail(informative, total))
    return len(families), len(informative), z, p_exact


def _gene_rows(
    chrom: str,
    records: list[tuple[int, str, str, str, dict[str, int]]],
    gene_index: dict[str, list[tuple[int, int, str]]],
    gene_names: dict[str, str],
    gene_types: dict[str, str],
    people: dict[str, tuple[str, str]],
    members: dict[str, list[tuple[str, bool]]],
    fcap: int | None,
) -> list[dict[str, object]]:
    if not records:
        return []
    records.sort(key=lambda item: item[0])
    intervals = _intervals_for(chrom, gene_index)
    if not intervals:
        return []
    assigned = _genes_at([pos for pos, _text, _ref, _alt, _dosages in records], intervals)
    sums: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    counts: dict[str, int] = defaultdict(int)
    first_pos: dict[str, int] = {}
    for (pos, _text, _ref, _alt, dosages), gene_ids in zip(records, assigned):
        for gene_id in gene_ids:
            counts[gene_id] += 1
            if gene_id not in first_pos or pos < first_pos[gene_id]:
                first_pos[gene_id] = pos
            for person, dosage in dosages.items():
                sums[gene_id][person] += dosage
    rows: list[dict[str, object]] = []
    for gene_id in sorted(counts, key=lambda gene: (first_pos[gene], gene)):
        result = _contrast(sums[gene_id], people, members, fcap)
        if result is None:
            continue
        n_families, n_informative, z, p_exact = result
        rows.append(
            {
                "CHROM": chrom,
                "gene": gene_id,
                "gene_name": gene_names.get(gene_id, ""),
                "gene_type": gene_types.get(gene_id, ""),
                "n_variants": counts[gene_id],
                "n_families": n_families,
                "n_informative": n_informative,
                "z": "" if z is None else format(z, ".12g"),
                "p_exact": format(p_exact, ".12g"),
            }
        )
    return rows


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
