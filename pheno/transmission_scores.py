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
_OTARI_TABLE = "variant_effects_comprehensive.tsv"
_OTARI_GENES = "interpretability_analysis.tsv"
_SKIPPED_OTARI = {"variant_id", "transcript_id", "max_effect", "mean_effect"}
_LABEL_COLUMNS = ("chrom", "pos", "ref", "alt")
_FIXED_FIELDS = ["CHROM", "gene", "gene_name", "gene_type", "n_variants", "n_families", "n_scored"]
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
) -> list[dict[str, object]]:
    """One wide row per gene. ``p`` is the two-sided normal tail of ``z``."""
    if not records:
        return []
    records.sort(key=lambda item: (item[0], item[2], item[3]))
    intervals = _intervals_for(chrom, gene_index)
    if not intervals:
        return []
    assigned = _genes_at([pos for pos, _text, _ref, _alt, _dosages in records], intervals)
    counts: dict[str, int] = defaultdict(int)
    scored: dict[str, int] = defaultdict(int)
    carriers: dict[str, set[str]] = defaultdict(set)
    weights: dict[str, dict[str, np.ndarray]] = {}
    first_pos: dict[str, int] = {}
    width = len(source.columns)
    for (pos, pos_text, ref, alt, dosages), gene_ids in zip(records, assigned):
        for gene_id in gene_ids:
            counts[gene_id] += 1
            if gene_id not in first_pos or pos < first_pos[gene_id]:
                first_pos[gene_id] = pos
            carriers[gene_id].update(dosages)
            values = score_vector(loaded, source, chrom, pos_text, ref, alt, gene_id)
            if values is None:
                continue
            scored[gene_id] += 1
            bucket = weights.setdefault(gene_id, {})
            for person, dosage in dosages.items():
                current = bucket.get(person)
                if current is None:
                    current = np.zeros(width, dtype=np.float64)
                    bucket[person] = current
                current += dosage * values
    rows: list[dict[str, object]] = []
    for gene_id in sorted(counts, key=lambda gene: (first_pos[gene], gene)):
        summary = _summarize(carriers[gene_id], weights.get(gene_id, {}), people, members, width)
        if summary is None:
            continue
        n_families, informative, z_values = summary
        row: dict[str, object] = {
            "CHROM": chrom,
            "gene": gene_id,
            "gene_name": gene_names.get(gene_id, ""),
            "gene_type": gene_types.get(gene_id, ""),
            "n_variants": counts[gene_id],
            "n_families": n_families,
            "n_scored": scored[gene_id],
        }
        for index, column in enumerate(source.columns):
            count = int(informative[index])
            row[f"{column}_n_informative"] = count
            if count == 0:
                row[f"{column}_z"] = ""
                row[f"{column}_p"] = format(1.0, ".12g")
            else:
                z_value = float(z_values[index])
                row[f"{column}_z"] = format(z_value, ".12g")
                row[f"{column}_p"] = format(float(2.0 * norm.sf(abs(z_value))), ".12g")
        rows.append(row)
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
    collapsed: dict[tuple[str, str], np.ndarray] = {}
    for key, rows in grouped.items():
        collapsed[key] = _collapse_transcripts(rows)
    return collapsed


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
