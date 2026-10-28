"""Gene dosage contrasts weighted by Otari or Seqweaver predictions."""

from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import norm

from pheno.burden import _variant_id
from pheno.gene_burden import _chrom_aliases, _genes_at, _intervals_for, _unversioned_gene_id

_CHROM_DIR = re.compile(r"^chr(\d+|X|Y)$")
_SHARD_NAME = re.compile(r"^(\d+)_(\d+)$")
_SEQWEAVER_SHARD = re.compile(r"^(?:variants_)?(\d+)_(\d+)$")
_OTARI_TABLE = "variant_effects_comprehensive.tsv"
_OTARI_GENES = "interpretability_analysis.tsv"
_SKIPPED_OTARI = {"variant_id", "transcript_id", "max_effect", "mean_effect"}
_LABEL_COLUMNS = ("chrom", "pos", "ref", "alt")
_FIXED_FIELDS = [
    "CHROM",
    "gene_id",
    "gene_name",
    "gene_type",
    "n_variants",
    "n_families",
    "n_scored",
    "n_informative",
    "z",
    "p_exact",
]
_FLUSH = 100_000


class ScoreSource:
    """One model, and the score columns written for every gene."""

    def __init__(self, kind: str, root: Path, columns: list[str]) -> None:
        self.kind = kind
        self.root = root
        self.columns = columns


def gene_scores_path(output: Path) -> Path:
    """Path of the weighted gene table written beside a single variant ``output`` file."""
    suffix = output.suffix or ".tsv"
    return output.with_name(f"{output.stem}.gene_scores{suffix}")


class GeneDebug:
    """Print how one gene id is assigned and whether its row is written."""

    def __init__(self, gene_id: str | None) -> None:
        self.gene_id = _unversioned_gene_id(gene_id) if gene_id else ""
        self.on_chromosome = False
        self.assigned = False
        self.wrote = False

    def matches(self, gene_id: str) -> bool:
        return bool(self.gene_id) and gene_id == self.gene_id

    def line(self, message: str) -> None:
        if self.gene_id:
            print(f"debug {self.gene_id}: {message}", flush=True)

    def begin_chromosome(
        self,
        chrom: str,
        intervals: list[tuple[int, int, str]],
        shards: list[tuple[int, int]],
        gene_names: dict[str, str],
    ) -> None:
        self.on_chromosome = False
        self.assigned = False
        self.wrote = False
        if not self.gene_id:
            return
        lo: int | None = None
        hi: int | None = None
        for start, end, gene_id in intervals:
            if gene_id != self.gene_id:
                continue
            self.on_chromosome = True
            lo = start if lo is None else min(lo, start)
            hi = end if hi is None else max(hi, end)
        if not self.on_chromosome or lo is None or hi is None:
            self.line(f"{chrom}: not in the annotation")
            return
        name = gene_names.get(self.gene_id, "")
        label = f"{self.gene_id} {name}".rstrip()
        if not shards:
            self.line(f"{chrom}: {label} interval {lo}-{hi}")
            return
        covered = [(start, end) for start, end in shards if start <= hi and end >= lo]
        ranges = ", ".join(f"{start}_{end}" for start, end in covered) or "none"
        self.line(f"{chrom}: {label} interval {lo}-{hi}; shards {ranges}")

    def variant(
        self,
        chrom: str,
        pos: str,
        ref: str,
        alt: str,
        dosages: dict[str, int],
        *,
        shard: str = "",
        scored: bool | None = None,
    ) -> None:
        self.assigned = True
        where = f" shard {shard}" if shard else ""
        doses = ",".join(f"{person}={dosage}" for person, dosage in sorted(dosages.items()))
        score = ""
        if scored is True:
            score = "; score found"
        elif scored is False:
            score = "; no score"
        self.line(f"variant {chrom} {pos} {ref} {alt}{where}; dosages {doses}{score}")

    def still_open(self, shard: str, reach: int, next_start: int) -> None:
        self.line(f"stays open after shard {shard}; reach {reach}; next shard starts {next_start}")

    def dropped(self, detail: str) -> None:
        self.line(f"dropping; {detail}")

    def wrote_row(self, row: dict[str, object]) -> None:
        self.wrote = True
        scored = row.get("n_scored", "")
        scored_text = f" n_scored={scored}" if "n_scored" in row else ""
        self.line(
            "writing"
            f" n_variants={row['n_variants']}"
            f" n_families={row['n_families']}"
            f"{scored_text}"
            f" n_informative={row['n_informative']}"
            f" z={row['z']}"
            f" p_exact={row['p_exact']}"
        )

    def end_chromosome(self, chrom: str) -> None:
        if not self.gene_id or not self.on_chromosome:
            return
        if not self.assigned:
            self.line(f"{chrom}: no variant assigned")


def format_carriers(
    carriers: set[str] | dict[str, object],
    people: dict[str, tuple[str, str]],
    members: dict[str, list[tuple[str, bool]]],
) -> str:
    """Why a gene with these carriers has no ASD and unaffected sibship."""
    parts: list[str] = []
    people_ids = carriers if isinstance(carriers, set) else carriers.keys()
    for person in sorted(people_ids):
        status, family_id = people[person]
        group = members[family_id]
        n_asd = sum(1 for _person, is_asd in group if is_asd)
        n_unaffected = sum(1 for _person, is_asd in group if not is_asd)
        parts.append(f"{person} {status} family {family_id} asd={n_asd} unaffected={n_unaffected}")
    detail = "; ".join(parts) if parts else "no carriers"
    return f"no carrier family has both an ASD child and an unaffected sibling; {detail}"


def score_fieldnames(columns: list[str]) -> list[str]:
    names = list(_FIXED_FIELDS)
    for column in columns:
        names.append(f"{column}_n_informative")
        names.append(f"{column}_z")
        names.append(f"{column}_p")
    if len(names) != len(set(names)):
        raise ValueError("score column names collide after adding suffixes")
    return names


def prepare_score_source(
    otari_dir: Path | None,
    seqweaver_dir: Path | None,
    column_names: Path | None,
    *,
    genes: Path | None,
) -> ScoreSource | None:
    """Validate the model flags and return the score columns, or None for genotype only."""
    if otari_dir is not None and seqweaver_dir is not None:
        raise ValueError("set only one of --otari-dir and --seqweaver")
    if column_names is not None and seqweaver_dir is None:
        raise ValueError("--column-names requires --seqweaver")
    if seqweaver_dir is not None and column_names is None:
        raise ValueError("--column-names is required with --seqweaver")
    if otari_dir is None and seqweaver_dir is None:
        return None
    if genes is None:
        raise ValueError("--genes is required with --otari-dir or --seqweaver")
    if otari_dir is not None:
        if not otari_dir.is_dir():
            raise ValueError(f"otari directory not found: {otari_dir}")
        columns = _otari_columns(otari_dir)
        source = ScoreSource("otari", otari_dir, columns)
    else:
        assert seqweaver_dir is not None
        if not seqweaver_dir.is_dir():
            raise ValueError(f"seqweaver directory not found: {seqweaver_dir}")
        source = ScoreSource("seqweaver", seqweaver_dir, _load_column_names(column_names))
    score_fieldnames(source.columns)
    return source


def paired_otari_shards(
    input_dir: Path,
    otari_dir: Path,
    prefix: str,
    only: str | None = None,
) -> list[tuple[str, list[tuple[Path, Path, int, int]]]]:
    """Pair site files with Otari shard directories that use the same start and end.

    Each entry is ``(sites file, Otari shard directory, start, end)`` in position order.
    """
    return _pair_shards(input_dir, prefix, otari_dir, _otari_shard_dirs, "otari", only)


def paired_seqweaver_shards(
    input_dir: Path,
    seqweaver_dir: Path,
    prefix: str,
    only: str | None = None,
) -> list[tuple[str, list[tuple[Path, Path, int, int]]]]:
    """Pair site files with Seqweaver directories named ``variants_START_END`` or ``START_END``."""
    return _pair_shards(input_dir, prefix, seqweaver_dir, _seqweaver_shard_dirs, "seqweaver", only)


def _pair_shards(
    input_dir: Path,
    prefix: str,
    model_dir: Path,
    discover,
    label: str,
    only: str | None = None,
) -> list[tuple[str, list[tuple[Path, Path, int, int]]]]:
    pattern = re.compile(rf"^{re.escape(prefix)}_(\d+)_(\d+)\.tsv$")
    plan: list[tuple[str, list[tuple[Path, Path, int, int]]]] = []
    found_sites = False
    chrom_dirs = _chrom_dirs(input_dir)
    if only is not None:
        chrom_dirs = [path for path in chrom_dirs if path.name == only]
    for chrom_dir in chrom_dirs:
        loose = [path for path in chrom_dir.glob(f"{prefix}*.tsv") if path.is_file()]
        if not loose:
            continue
        found_sites = True
        sites: dict[tuple[int, int], Path] = {}
        for path in loose:
            match = pattern.fullmatch(path.name)
            if match is None:
                raise ValueError(f"{path.name} does not match {prefix}_START_END.tsv")
            key = (int(match.group(1)), int(match.group(2)))
            if key in sites:
                raise ValueError(f"duplicate shard {chrom_dir.name} {key[0]}_{key[1]}")
            sites[key] = path
        model = discover(model_dir / chrom_dir.name)
        if set(sites) != set(model):
            raise ValueError(
                f"{chrom_dir.name} shard splits differ: "
                f"sites has {_shard_text(sorted(sites))}, "
                f"{label} has {_shard_text(sorted(model))}"
            )
        plan.append((chrom_dir.name, [(sites[key], model[key], key[0], key[1]) for key in sorted(sites)]))
    if not found_sites:
        if only is not None:
            return []
        raise ValueError(f"no {prefix}_START_END.tsv files under {input_dir}")
    return plan


def _otari_shard_dirs(chrom_dir: Path) -> dict[tuple[int, int], Path]:
    found: dict[tuple[int, int], Path] = {}
    if not chrom_dir.is_dir():
        return found
    for path in chrom_dir.iterdir():
        match = _SHARD_NAME.fullmatch(path.name)
        if path.is_dir() and match is not None and (path / _OTARI_TABLE).is_file():
            found[(int(match.group(1)), int(match.group(2)))] = path
    return found


def _seqweaver_shard_dirs(chrom_dir: Path) -> dict[tuple[int, int], Path]:
    found: dict[tuple[int, int], Path] = {}
    if not chrom_dir.is_dir():
        return found
    for path in chrom_dir.iterdir():
        match = _SEQWEAVER_SHARD.fullmatch(path.name)
        if not path.is_dir() or match is None or not _has_seqweaver_matrix(path):
            continue
        key = (int(match.group(1)), int(match.group(2)))
        if key in found:
            raise ValueError(f"duplicate seqweaver shard {chrom_dir.name} {key[0]}_{key[1]}")
        found[key] = path
    return found


def _has_seqweaver_matrix(shard: Path) -> bool:
    return any(
        path.is_file() and path.parent.name == "seqweaver_human"
        for path in shard.rglob("annot.*diffs.h5")
    )


def _shard_text(shards: list[tuple[int, int]]) -> str:
    labels = [f"{start}_{end}" for start, end in shards]
    if len(labels) > 6:
        return ", ".join(labels[:6]) + f", ... ({len(labels)} shards)"
    if not labels:
        return "none"
    return ", ".join(labels)


def load_chromosome_scores(
    source: ScoreSource,
    chrom: str,
    records: list[tuple[int, str, str, str, dict[str, int]]],
) -> dict[tuple[str, str], np.ndarray] | dict[tuple[str, str, str, str], np.ndarray]:
    """Load one weight vector per needed variant, or per variant and gene for Otari."""
    if source.kind == "otari":
        needed = {_variant_id(chrom, pos, ref, alt) for _pos, pos, ref, alt, _dosages in records}
        return _load_otari(source.root / chrom, source.columns, needed)
    needed: set[tuple[str, str, str, str]] = set()
    for _pos, pos, ref, alt, _dosages in records:
        for alias in _chrom_aliases(chrom):
            needed.add((alias, pos, ref, alt))
    return _load_seqweaver(source.root / chrom, source.columns, needed)


def score_vector(
    loaded: dict[tuple[str, str], np.ndarray] | dict[tuple[str, str, str, str], np.ndarray],
    source: ScoreSource,
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    gene_id: str,
) -> np.ndarray | None:
    if source.kind == "otari":
        otari = loaded
        assert isinstance(otari, dict)
        return otari.get((_variant_id(chrom, pos, ref, alt), gene_id))
    seqweaver = loaded
    for alias in _chrom_aliases(chrom):
        found = seqweaver.get((alias, pos, ref, alt))
        if found is not None:
            return found
    return None


def load_otari_shard(
    shard: Path,
    columns: list[str],
    needed: set[str],
) -> dict[tuple[str, str], np.ndarray]:
    """Collapse transcripts in one ``start_end`` Otari directory."""
    if not needed:
        return {}
    table = shard / _OTARI_TABLE
    gene_map = shard / _OTARI_GENES
    if not gene_map.is_file():
        raise ValueError(f"{table} is missing {gene_map.name}")
    grouped = _read_otari_shard(table, gene_map, columns, needed)
    return {key: _collapse_transcripts(rows) for key, rows in grouped.items()}


def gene_shard_covers(
    intervals: list[tuple[int, int, str]],
    shards: list[tuple[int, int]],
) -> list[tuple[str, int, int, list[tuple[int, int]]]]:
    """How many shards a gene stays open across.

    The span is the gene's lowest start through its highest end, which already
    includes the window. A gene stays open from the first shard that overlaps
    that span through the last shard that still starts at or before the end.
    """
    bounds: dict[str, list[int]] = {}
    for lo, hi, gene_id in intervals:
        current = bounds.get(gene_id)
        if current is None:
            bounds[gene_id] = [lo, hi]
        else:
            current[0] = min(current[0], lo)
            current[1] = max(current[1], hi)
    ordered = sorted(shards)
    found: list[tuple[str, int, int, list[tuple[int, int]]]] = []
    for gene_id, (lo, hi) in bounds.items():
        covered = [(start, end) for start, end in ordered if start <= hi and end >= lo]
        if not covered:
            continue
        first = covered[0][0]
        held = [(start, end) for start, end in ordered if first <= start <= hi]
        found.append((gene_id, lo, hi, held))
    found.sort(key=lambda item: (-len(item[3]), item[1], item[0]))
    return found


def report_long_genes(
    chrom: str,
    intervals: list[tuple[int, int, str]],
    shards: list[tuple[int, int]],
    gene_names: dict[str, str],
    minimum: int = 3,
) -> int:
    """Print genes that stay open across ``minimum`` or more shards. Returns that count."""
    covers = gene_shard_covers(intervals, shards)
    long = [item for item in covers if len(item[3]) >= minimum]
    if not long:
        if covers:
            longest = max(len(item[3]) for item in covers)
            print(
                f"{chrom}: no gene spans {minimum} or more shards (longest span is {longest})",
                flush=True,
            )
        else:
            print(f"{chrom}: no gene spans {minimum} or more shards", flush=True)
        return 0
    for gene_id, lo, hi, held in long:
        name = gene_names.get(gene_id, "")
        label = f"{gene_id} {name}".rstrip()
        ranges = ", ".join(f"{start}_{end}" for start, end in held)
        print(
            f"{chrom}: {label} ({lo}-{hi}) spans {len(held)} shards and stays open across {ranges}",
            flush=True,
        )
    noun = "gene spans" if len(long) == 1 else "genes span"
    print(f"{chrom}: {len(long)} {noun} {minimum} or more shards", flush=True)
    return len(long)


class GeneScores:
    """Running per-gene totals. Genes are emitted once no later shard can add to them."""

    def __init__(
        self,
        chrom: str,
        intervals: list[tuple[int, int, str]],
        gene_names: dict[str, str],
        gene_types: dict[str, str],
        people: dict[str, tuple[str, str]],
        members: dict[str, list[tuple[str, bool]]],
        source: ScoreSource,
        fcap: int | None = None,
        debug: GeneDebug | None = None,
    ) -> None:
        self.chrom = chrom
        self.intervals = intervals
        self.gene_names = gene_names
        self.gene_types = gene_types
        self.people = people
        self.members = members
        self.source = source
        self.fcap = fcap
        self.debug = debug or GeneDebug(None)
        self.reach: dict[str, int] = {}
        for _lo, hi, gene_id in intervals:
            current = self.reach.get(gene_id)
            if current is None or hi > current:
                self.reach[gene_id] = hi
        self.counts: dict[str, int] = defaultdict(int)
        self.scored: dict[str, int] = defaultdict(int)
        self.carriers: dict[str, set[str]] = defaultdict(set)
        self.weights: dict[str, dict[str, np.ndarray]] = {}
        self.dosages: dict[str, dict[str, int]] = {}
        self.first_pos: dict[str, int] = {}

    def add(
        self,
        records: list[tuple[int, str, str, str, dict[str, int]]],
        loaded: dict[tuple[str, str], np.ndarray] | dict[tuple[str, str, str, str], np.ndarray],
        shard: str = "",
    ) -> None:
        if not records or not self.intervals:
            return
        records.sort(key=lambda item: (item[0], item[2], item[3]))
        assigned = _genes_at([pos for pos, _text, _ref, _alt, _dosages in records], self.intervals)
        width = len(self.source.columns)
        for (pos, pos_text, ref, alt, dosages), gene_ids in zip(records, assigned):
            for gene_id in gene_ids:
                self.counts[gene_id] += 1
                if gene_id not in self.first_pos or pos < self.first_pos[gene_id]:
                    self.first_pos[gene_id] = pos
                self.carriers[gene_id].update(dosages)
                dose_bucket = self.dosages.setdefault(gene_id, {})
                for person, dosage in dosages.items():
                    dose_bucket[person] = dose_bucket.get(person, 0) + dosage
                values = score_vector(loaded, self.source, self.chrom, pos_text, ref, alt, gene_id)
                if self.debug.matches(gene_id):
                    self.debug.variant(
                        self.chrom,
                        pos_text,
                        ref,
                        alt,
                        dosages,
                        shard=shard,
                        scored=values is not None,
                    )
                if values is None:
                    continue
                self.scored[gene_id] += 1
                bucket = self.weights.setdefault(gene_id, {})
                for person, dosage in dosages.items():
                    current = bucket.get(person)
                    if current is None:
                        current = np.zeros(width, dtype=np.float64)
                        bucket[person] = current
                    current += dosage * values

    def emit(self, before: int | None) -> list[dict[str, object]]:
        """Write genes that end before ``before``. ``None`` writes every gene still open."""
        ready = [
            gene_id
            for gene_id in self.counts
            if before is None or self.reach.get(gene_id, -1) < before
        ]
        ready.sort(key=lambda gene_id: (self.first_pos[gene_id], gene_id))
        rows: list[dict[str, object]] = []
        for gene_id in ready:
            row = self._row(gene_id)
            if self.debug.matches(gene_id):
                if row is None:
                    self.debug.dropped(format_carriers(self.carriers[gene_id], self.people, self.members))
                else:
                    self.debug.wrote_row(row)
            if row is not None:
                rows.append(row)
            self._drop(gene_id)
        return rows

    def _row(self, gene_id: str) -> dict[str, object] | None:
        summary = _summarize(
            self.carriers[gene_id],
            self.weights.get(gene_id, {}),
            self.people,
            self.members,
            len(self.source.columns),
        )
        if summary is None:
            return None
        n_families, informative, z_values = summary
        row: dict[str, object] = {
            "CHROM": self.chrom,
            "gene_id": gene_id,
            "gene_name": self.gene_names.get(gene_id, ""),
            "gene_type": self.gene_types.get(gene_id, ""),
            "n_variants": self.counts[gene_id],
            "n_families": n_families,
            "n_scored": self.scored[gene_id],
        }
        from pheno.transmission import _contrast

        plain = _contrast(self.dosages.get(gene_id, {}), self.people, self.members, self.fcap)
        if plain is None:
            row["n_informative"] = 0
            row["z"] = ""
            row["p_exact"] = format(1.0, ".12g")
        else:
            _families, n_informative, z_value, p_exact = plain
            row["n_informative"] = n_informative
            row["z"] = "" if z_value is None else format(z_value, ".12g")
            row["p_exact"] = format(p_exact, ".12g")
        for index, column in enumerate(self.source.columns):
            count = int(informative[index])
            row[f"{column}_n_informative"] = count
            if count == 0:
                row[f"{column}_z"] = ""
                row[f"{column}_p"] = format(1.0, ".12g")
            else:
                z_value = float(z_values[index])
                row[f"{column}_z"] = format(z_value, ".12g")
                row[f"{column}_p"] = format(float(2.0 * norm.sf(abs(z_value))), ".12g")
        return row

    def _drop(self, gene_id: str) -> None:
        del self.counts[gene_id]
        self.scored.pop(gene_id, None)
        self.carriers.pop(gene_id, None)
        self.weights.pop(gene_id, None)
        self.dosages.pop(gene_id, None)
        self.first_pos.pop(gene_id, None)


def gene_score_rows(
    chrom: str,
    records: list[tuple[int, str, str, str, dict[str, int]]],
    gene_index: dict[str, list[tuple[int, int, str]]],
    gene_names: dict[str, str],
    gene_types: dict[str, str],
    people: dict[str, tuple[str, str]],
    members: dict[str, list[tuple[str, bool]]],
    source: ScoreSource,
    loaded: dict[tuple[str, str], np.ndarray] | dict[tuple[str, str, str, str], np.ndarray],
    fcap: int | None = None,
    debug: GeneDebug | None = None,
) -> list[dict[str, object]]:
    """One wide row per gene. ``p`` is the two-sided normal tail of ``z``."""
    debug = debug or GeneDebug(None)
    intervals = _intervals_for(chrom, gene_index)
    debug.begin_chromosome(chrom, intervals, [], gene_names)
    if not records or not intervals:
        debug.end_chromosome(chrom)
        return []
    running = GeneScores(
        chrom,
        intervals,
        gene_names,
        gene_types,
        people,
        members,
        source,
        fcap,
        debug,
    )
    running.add(records, loaded)
    rows = running.emit(None)
    debug.end_chromosome(chrom)
    return rows


def _summarize(
    carriers: set[str],
    person_weights: dict[str, np.ndarray],
    people: dict[str, tuple[str, str]],
    members: dict[str, list[tuple[str, bool]]],
    width: int,
) -> tuple[int, np.ndarray, np.ndarray] | None:
    family_ids: list[str] = []
    seen: set[str] = set()
    for person in carriers:
        family_id = people[person][1]
        if family_id in seen:
            continue
        seen.add(family_id)
        group = members[family_id]
        if any(is_asd for _person, is_asd in group) and any(not is_asd for _person, is_asd in group):
            family_ids.append(family_id)
    if not family_ids:
        return None
    total = np.zeros(width, dtype=np.float64)
    variance = np.zeros(width, dtype=np.float64)
    informative = np.zeros(width, dtype=np.int64)
    for family_id in family_ids:
        group = members[family_id]
        matrix = np.zeros((len(group), width), dtype=np.float64)
        asd_rows: list[int] = []
        unaffected_rows: list[int] = []
        for index, (person, is_asd) in enumerate(group):
            found = person_weights.get(person)
            if found is not None:
                matrix[index] = found
            if is_asd:
                asd_rows.append(index)
            else:
                unaffected_rows.append(index)
        n_asd = len(asd_rows)
        n_unaffected = len(unaffected_rows)
        n = n_asd + n_unaffected
        observed = matrix[asd_rows].mean(axis=0) - matrix[unaffected_rows].mean(axis=0)
        mean = matrix.mean(axis=0)
        squares = ((matrix - mean) ** 2).sum(axis=0)
        total += observed
        variance += n * squares / (n_asd * n_unaffected * (n - 1))
        informative += squares > 0
    z_values = np.zeros(width, dtype=np.float64)
    np.divide(total, np.sqrt(variance), out=z_values, where=informative > 0)
    return len(family_ids), informative, z_values


def _otari_columns(root: Path) -> list[str]:
    tables = _otari_tables(root)
    if not tables:
        raise ValueError(f"no {_OTARI_TABLE} files under {root}")
    with tables[0].open(encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
    if "variant_id" not in header or "transcript_id" not in header:
        raise ValueError(f"{tables[0]} is missing columns: variant_id, transcript_id")
    columns = [name for name in header if name not in _SKIPPED_OTARI]
    if not columns:
        raise ValueError(f"{tables[0]} has no score columns")
    return columns


def _otari_tables(root: Path) -> list[Path]:
    tables: list[Path] = []
    if not root.is_dir():
        return tables
    for chrom_dir in _chrom_dirs(root):
        for shard in sorted((path for path in chrom_dir.iterdir() if path.is_dir()), key=lambda item: item.name):
            table = shard / _OTARI_TABLE
            if table.is_file():
                tables.append(table)
    return tables


def _load_otari(
    chrom_dir: Path,
    columns: list[str],
    needed: set[str],
) -> dict[tuple[str, str], np.ndarray]:
    grouped: dict[tuple[str, str], list[np.ndarray]] = defaultdict(list)
    if not chrom_dir.is_dir() or not needed:
        return {}
    for shard in sorted((path for path in chrom_dir.iterdir() if path.is_dir()), key=lambda item: item.name):
        table = shard / _OTARI_TABLE
        if not table.is_file():
            continue
        gene_map = shard / _OTARI_GENES
        if not gene_map.is_file():
            raise ValueError(f"{table} is missing {gene_map.name}")
        for key, rows in _read_otari_shard(table, gene_map, columns, needed).items():
            grouped[key].extend(rows)
    return {key: _collapse_transcripts(rows) for key, rows in grouped.items()}


def _read_otari_shard(
    table: Path,
    gene_map: Path,
    columns: list[str],
    needed: set[str],
) -> dict[tuple[str, str], list[np.ndarray]]:
    grouped: dict[tuple[str, str], list[np.ndarray]] = defaultdict(list)
    mapping = _transcript_genes(gene_map, needed)
    with table.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        file_columns = [name for name in fieldnames if name not in _SKIPPED_OTARI]
        if file_columns != columns:
            raise ValueError(f"{table} score columns do not match")
        for row in reader:
            variant_id = (row.get("variant_id") or "").strip()
            if variant_id not in needed:
                continue
            transcript_id = (row.get("transcript_id") or "").strip()
            gene_id = mapping.get((variant_id, transcript_id))
            if not gene_id:
                continue
            parsed: list[float] = []
            complete = True
            for column in columns:
                cell = (row.get(column) or "").strip()
                if not cell:
                    complete = False
                    continue
                try:
                    value = float(cell)
                except ValueError as exc:
                    raise ValueError(f"{table} has a non-numeric {column} for {variant_id}") from exc
                if not math.isfinite(value):
                    complete = False
                parsed.append(value)
            if complete and len(parsed) == len(columns):
                grouped[(variant_id, gene_id)].append(np.asarray(parsed, dtype=np.float64))
    return grouped


def _collapse_transcripts(rows: list[np.ndarray]) -> np.ndarray:
    chosen = np.empty(rows[0].shape[0], dtype=np.float64)
    for index in range(chosen.shape[0]):
        best = float(rows[0][index])
        for row in rows[1:]:
            value = float(row[index])
            if abs(value) > abs(best) or (abs(value) == abs(best) and value > best):
                best = value
        chosen[index] = best
    return chosen


def _transcript_genes(path: Path, needed: set[str]) -> dict[tuple[str, str], str]:
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
            gene_id = _unversioned_gene_id(row.get("gene_id") or "")
            key = (variant_id, transcript_id)
            if transcript_id and gene_id and key not in mapping:
                mapping[key] = gene_id
    return mapping


def _load_seqweaver(
    chrom_dir: Path,
    columns: list[str],
    needed: set[tuple[str, str, str, str]],
) -> dict[tuple[str, str, str, str], np.ndarray]:
    found: dict[tuple[str, str, str, str], np.ndarray] = {}
    seen: set[tuple[str, str, str, str]] = set()
    if not chrom_dir.is_dir() or not needed:
        return found
    try:
        import h5py
    except ImportError as exc:
        raise ValueError("h5py is required for --seqweaver") from exc
    pairs: list[tuple[Path, Path]] = []
    for path in chrom_dir.rglob("annot.*diffs.h5"):
        if path.is_file() and path.parent.name == "seqweaver_human":
            pairs.append((path, _row_labels_path(path)))
    for diffs, labels in sorted(pairs, key=lambda item: item[0].as_posix()):
        keys = _label_keys(labels)
        with h5py.File(diffs, "r") as handle:
            if "data" not in handle:
                raise ValueError(f"{diffs} is missing data")
            data = handle["data"]
            if data.ndim != 2:
                raise ValueError(f"{diffs} data is not a matrix")
            if data.shape[1] != len(columns):
                raise ValueError(
                    f"{diffs} has {data.shape[1]} columns and the column names file has {len(columns)}"
                )
            n_rows = int(data.shape[0])
            if len(keys) != n_rows:
                raise ValueError(f"{labels} has {len(keys)} label rows for {n_rows} variants")
            for start in range(0, n_rows, _FLUSH):
                block = np.asarray(data[start : start + _FLUSH], dtype=np.float64)
                for offset, row in enumerate(block):
                    key = keys[start + offset]
                    if key not in needed or key in seen:
                        continue
                    seen.add(key)
                    if not np.isfinite(row).all():
                        continue
                    found[key] = np.array(row, copy=True)
    return found


def _row_labels_path(diffs: Path) -> Path:
    name = diffs.name
    if not name.endswith("_diffs.h5"):
        raise ValueError(f"{diffs} is not an annot.*diffs.h5 file")
    labels = diffs.with_name(name[: -len("_diffs.h5")] + "_row_labels.txt")
    if not labels.is_file():
        raise ValueError(f"{diffs} is missing {labels.name}")
    return labels


def _label_keys(path: Path) -> list[tuple[str, str, str, str]]:
    """Variant keys from a headered label table. A header-only file has no variants."""
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"{path} is empty")
    if "\t" in lines[0]:
        rows = list(csv.reader(lines, delimiter="\t"))
    else:
        rows = [line.split() for line in lines]
    header = rows[0]
    missing = [name for name in _LABEL_COLUMNS if name not in header]
    if missing:
        raise ValueError(f"{path} is missing columns: {missing}")
    index = {name: header.index(name) for name in _LABEL_COLUMNS}
    keys: list[tuple[str, str, str, str]] = []
    width = max(index.values())
    for line_number, row in enumerate(rows[1:], start=2):
        if len(row) <= width:
            raise ValueError(f"{path}:{line_number} is missing chrom, pos, ref, or alt")
        fields = tuple(row[index[name]].strip() for name in _LABEL_COLUMNS)
        if not all(fields):
            raise ValueError(f"{path}:{line_number} is missing chrom, pos, ref, or alt")
        keys.append((fields[0], fields[1], fields[2], fields[3]))
    return keys


def _load_column_names(path: Path | None) -> list[str]:
    if path is None or not path.is_file():
        raise ValueError(f"column names file not found: {path}")
    names = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not names:
        raise ValueError(f"{path} has no column names")
    if len(names) != len(set(names)):
        raise ValueError(f"{path} has duplicate column names")
    return names


def _chrom_dirs(root: Path) -> list[Path]:
    dirs = [path for path in root.iterdir() if path.is_dir() and _CHROM_DIR.fullmatch(path.name)]

    def sort_key(path: Path) -> tuple[int, int | str]:
        suffix = path.name[3:]
        if suffix == "X":
            return (1, "X")
        if suffix == "Y":
            return (2, "Y")
        return (0, int(suffix))

    return sorted(dirs, key=sort_key)
