"""Within-family ASD dosage contrast for TRIO_CALLS variants."""

from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import TextIO

from scipy.stats import norm

from pheno.burden import _variant_id
from pheno.gene_burden import _chrom_aliases, _genes_at, _intervals_for, load_gene_index
from pheno.transmission_scores import (
    GeneDebug,
    GeneScores,
    ScoreSource,
    _load_seqweaver,
    format_carriers,
    gene_score_rows,
    gene_scores_path,
    load_chromosome_scores,
    load_otari_shard,
    paired_otari_shards,
    paired_seqweaver_shards,
    prepare_score_source,
    report_long_genes,
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
    "gene_id",
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
    debug_gene: str | None = None,
    chrom: str | None = None,
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
    replaces that gene table with a wide score-weighted table. Its per-score
    ``p`` is always the normal tail. That table also repeats ``n_informative``,
    ``z``, and ``p_exact`` with every variant weight set to 1, and those follow
    ``fcap``. ``debug_gene`` prints how that gene id is assigned and requires
    ``chrom``, so only that chromosome is read. A version suffix is ignored.
    ``chrom`` alone, such as ``1`` or ``chrX``, also limits the run to one
    chromosome.
    """
    if (output is None) == (output_dir is None):
        raise ValueError("set exactly one of --output and --output-dir")
    if fcap is not None and fcap < 0:
        raise ValueError("--fcap must be at least 0")
    if debug_gene and not chrom:
        raise ValueError("--debug-gene requires --chr")
    selected = _normalize_chrom(chrom) if chrom else None
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
    score_plan = None
    files: list[Path] = []
    if score_source is not None and score_source.kind == "otari":
        score_plan = paired_otari_shards(input_dir, score_source.root, prefix, selected)
    elif score_source is not None and score_source.kind == "seqweaver":
        score_plan = paired_seqweaver_shards(input_dir, score_source.root, prefix, selected)
    else:
        files = _discover(input_dir, prefix, selected)
    if selected is not None:
        if score_plan is not None:
            score_plan = [item for item in score_plan if item[0] == selected]
            if not score_plan:
                raise ValueError(f"no {prefix}_START_END.tsv files for {selected} under {input_dir}")
        else:
            files = [path for path in files if path.parent.name == selected]
            if not files:
                raise ValueError(f"no {prefix} files for {selected} under {input_dir}")
    if score_plan is not None:
        for _chrom, shards in score_plan:
            for sites_path, _shard, _start, _end in shards:
                _require_trio_header(sites_path)
    else:
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
    members = _members_by_family(people)
    debug = GeneDebug(debug_gene)
    try:
        if score_plan is not None:
            assert score_source is not None
            assert gene_index is not None
            written, gene_written = _process_score_shards(
                score_plan,
                prefix,
                writer,
                output_dir,
                gene_index,
                gene_names or {},
                gene_types or {},
                people,
                members,
                fcap,
                score_source,
                score_writer,
                debug,
            )
        else:
            written, gene_written = _process_files(
                files,
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
                debug,
            )
    finally:
        if handle is not None:
            handle.close()
        if gene_handle is not None:
            gene_handle.close()
        if score_handle is not None:
            score_handle.close()
    return written, gene_written


def _process_files(
    files: list[Path],
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
    score_source: ScoreSource | None,
    score_writer: csv.DictWriter | None,
    debug: GeneDebug | None = None,
) -> tuple[int, int]:
    written = 0
    gene_written = 0
    seen: set[tuple[str, str, str, str]] = set()
    pending: list[dict[str, object]] = []
    gene_records: list[tuple[int, str, str, str, dict[str, int]]] = []
    current_chrom = ""
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
                gene_names,
                gene_types,
                people,
                members,
                fcap,
                score_source,
                score_writer,
                debug,
            )
            written += variants
            gene_written += gene_rows
            if current_chrom:
                print(f"finished processing {current_chrom}", flush=True)
            current_chrom = chrom_name
            seen.clear()
        _read_site_file(path, seen, pending, gene_records, people, members, fcap, gene_index is not None)
    variants, gene_rows = _flush_chromosome(
        current_chrom,
        pending,
        gene_records,
        prefix,
        writer,
        gene_writer,
        output_dir,
        gene_index,
        gene_names,
        gene_types,
        people,
        members,
        fcap,
        score_source,
        score_writer,
        debug,
    )
    written += variants
    gene_written += gene_rows
    if current_chrom:
        print(f"finished processing {current_chrom}", flush=True)
    return written, gene_written


def _process_score_shards(
    plan: list[tuple[str, list[tuple[Path, Path, int, int]]]],
    prefix: str,
    writer: csv.DictWriter | None,
    output_dir: Path | None,
    gene_index: dict[str, list[tuple[int, int, str]]],
    gene_names: dict[str, str],
    gene_types: dict[str, str],
    people: dict[str, tuple[str, str]],
    members: dict[str, list[tuple[str, bool]]],
    fcap: int | None,
    score_source: ScoreSource,
    score_writer: csv.DictWriter | None,
    debug: GeneDebug | None = None,
) -> tuple[int, int]:
    written = 0
    gene_written = 0
    fields = score_fieldnames(score_source.columns)
    open_scores: dict[str, tuple[TextIO, csv.DictWriter]] = {}
    debug = debug or GeneDebug(None)
    try:
        for chrom, shards in plan:
            pending: list[dict[str, object]] = []
            seen: set[tuple[str, str, str, str]] = set()
            records: list[tuple[int, str, str, str, dict[str, int]]] = []
            intervals = _intervals_for(chrom, gene_index)
            shard_bounds = [(start, end) for _sites, _shard, start, end in shards]
            report_long_genes(chrom, intervals, shard_bounds, gene_names)
            debug.begin_chromosome(chrom, intervals, shard_bounds, gene_names)
            running = GeneScores(
                chrom,
                intervals,
                gene_names,
                gene_types,
                people,
                members,
                score_source,
                fcap,
                debug,
            )
            for index, (sites_path, shard, start, end) in enumerate(shards):
                records.clear()
                _read_site_file(sites_path, seen, pending, records, people, members, fcap, True)
                assigned = 0
                if records:
                    assigned = running.add(
                        records,
                        _load_shard_scores(score_source, chrom, shard, records),
                        f"{start}_{end}",
                    )
                next_start = shards[index + 1][2] if index + 1 < len(shards) else None
                n_rows = _append_score_rows(
                    chrom,
                    running.emit(None if next_start is None else next_start),
                    fields,
                    prefix,
                    output_dir,
                    score_writer,
                    open_scores,
                )
                gene_written += n_rows
                if debug.gene_id and debug.gene_id in running.counts and next_start is not None:
                    debug.still_open(f"{start}_{end}", running.reach.get(debug.gene_id, -1), next_start)
                print(
                    f"{chrom} shard {start}_{end}: variants {len(records)}, assigned {assigned}, "
                    f"gene rows {n_rows}, still open {len(running.counts)}",
                    flush=True,
                )
            debug.end_chromosome(chrom)
            _close_score_file(open_scores, chrom)
            variants, _gene_rows = _flush_chromosome(
                chrom,
                pending,
                [],
                prefix,
                writer,
                None,
                output_dir,
                None,
                {},
                {},
                people,
                members,
                fcap,
                None,
                None,
            )
            written += variants
            print(f"finished processing {chrom}", flush=True)
    finally:
        for handle, _score_writer in open_scores.values():
            handle.close()
    return written, gene_written


def _load_shard_scores(
    score_source: ScoreSource,
    chrom: str,
    shard: Path,
    records: list[tuple[int, str, str, str, dict[str, int]]],
):
    if score_source.kind == "otari":
        needed = {_variant_id(chrom, pos, ref, alt) for _pos, pos, ref, alt, _dosages in records}
        return load_otari_shard(shard, score_source.columns, needed)
    needed_labels: set[tuple[str, str, str, str]] = set()
    for _pos, pos, ref, alt, _dosages in records:
        for alias in _chrom_aliases(chrom):
            needed_labels.add((alias, pos, ref, alt))
    return _load_seqweaver(shard, score_source.columns, needed_labels)


def _close_score_file(open_scores: dict[str, tuple[TextIO, csv.DictWriter]], chrom: str) -> None:
    """Close this chromosome's gene-score file before the next chromosome starts."""
    opened = open_scores.pop(chrom, None)
    if opened is None:
        return
    handle, _writer = opened
    handle.close()


def _append_score_rows(
    chrom: str,
    rows: list[dict[str, object]],
    fields: list[str],
    prefix: str,
    output_dir: Path | None,
    score_writer: csv.DictWriter | None,
    open_scores: dict[str, tuple[TextIO, csv.DictWriter]],
) -> int:
    if not rows:
        return 0
    if output_dir is not None:
        if chrom not in open_scores:
            path = output_dir / chrom / f"{prefix}.gene_scores.tsv"
            path.parent.mkdir(parents=True, exist_ok=True)
            handle = path.open("w", encoding="utf-8", newline="")
            chrom_writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            chrom_writer.writeheader()
            open_scores[chrom] = (handle, chrom_writer)
        handle, chrom_writer = open_scores[chrom]
        chrom_writer.writerows(rows)
        handle.flush()
    elif score_writer is not None:
        score_writer.writerows(rows)
    return len(rows)


def _read_site_file(
    path: Path,
    seen: set[tuple[str, str, str, str]],
    pending: list[dict[str, object]],
    gene_records: list[tuple[int, str, str, str, dict[str, int]]],
    people: dict[str, tuple[str, str]],
    members: dict[str, list[tuple[str, bool]]],
    fcap: int | None,
    keep_genes: bool,
) -> None:
    for row in _iter_rows(path):
        if row is None:
            continue
        chrom, pos, variant_id, ref, alt, calls = row
        key = (chrom, pos, ref, alt)
        if key in seen:
            raise ValueError(f"duplicate variant {chrom} {pos} {ref} {alt}")
        seen.add(key)
        listed = _cohort_dosages(calls, people)
        if keep_genes and listed:
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
    debug: GeneDebug | None = None,
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
            fcap,
            debug,
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
            chrom,
            gene_records,
            gene_index,
            gene_names,
            gene_types,
            people,
            members,
            fcap,
            debug,
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
    debug: GeneDebug | None = None,
) -> list[dict[str, object]]:
    debug = debug or GeneDebug(None)
    intervals = _intervals_for(chrom, gene_index)
    debug.begin_chromosome(chrom, intervals, [], gene_names)
    if not records or not intervals:
        debug.end_chromosome(chrom)
        return []
    records.sort(key=lambda item: item[0])
    assigned = _genes_at([pos for pos, _text, _ref, _alt, _dosages in records], intervals)
    sums: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    counts: dict[str, int] = defaultdict(int)
    first_pos: dict[str, int] = {}
    for (pos, pos_text, ref, alt, dosages), gene_ids in zip(records, assigned):
        for gene_id in gene_ids:
            counts[gene_id] += 1
            if gene_id not in first_pos or pos < first_pos[gene_id]:
                first_pos[gene_id] = pos
            for person, dosage in dosages.items():
                sums[gene_id][person] += dosage
            if debug.matches(gene_id):
                debug.variant(chrom, pos_text, ref, alt, dosages)
    rows: list[dict[str, object]] = []
    for gene_id in sorted(counts, key=lambda gene: (first_pos[gene], gene)):
        result = _contrast(sums[gene_id], people, members, fcap)
        if result is None:
            if debug.matches(gene_id):
                debug.dropped(format_carriers(sums[gene_id], people, members))
            continue
        n_families, n_informative, z, p_exact = result
        row = {
            "CHROM": chrom,
            "gene_id": gene_id,
            "gene_name": gene_names.get(gene_id, ""),
            "gene_type": gene_types.get(gene_id, ""),
            "n_variants": counts[gene_id],
            "n_families": n_families,
            "n_informative": n_informative,
            "z": "" if z is None else format(z, ".12g"),
            "p_exact": format(p_exact, ".12g"),
        }
        if debug.matches(gene_id):
            debug.wrote_row(row)
        rows.append(row)
    debug.end_chromosome(chrom)
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


def _normalize_chrom(chrom: str) -> str:
    text = chrom.strip()
    suffix = text[3:] if text.lower().startswith("chr") else text
    if suffix.upper() in {"X", "Y"}:
        suffix = suffix.upper()
    name = f"chr{suffix}"
    if not _CHROM_DIR.fullmatch(name):
        raise ValueError(f"--chr must be a chromosome number, X, or Y, not {chrom}")
    return name


def _discover(input_dir: Path, prefix: str, only: str | None = None) -> list[Path]:
    found: list[Path] = []
    for chrom_dir in _chrom_dirs(input_dir):
        if only is not None and chrom_dir.name != only:
            continue
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
