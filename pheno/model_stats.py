"""Column statistics for model score tables."""

from __future__ import annotations

import csv
import math
import re
import shutil
from fractions import Fraction
from pathlib import Path

import numpy as np

_CHROM_DIR = re.compile(r"^chr(\d+|X|Y)$")
_SHARD_DIR = re.compile(r"^(\d+)_(\d+)$")
_OTARI_NAME = "variant_effects_comprehensive.tsv"
_ID_COLUMNS = ("variant_id", "transcript_id")
_BINS = 50
_FLUSH = 100_000
_QUANTILES = (
    ("q005", Fraction(5, 1000)),
    ("q01", Fraction(1, 100)),
    ("q02", Fraction(2, 100)),
    ("q025", Fraction(25, 1000)),
    ("q05", Fraction(5, 100)),
    ("q50", Fraction(1, 2)),
    ("q95", Fraction(95, 100)),
    ("q975", Fraction(975, 1000)),
    ("q98", Fraction(98, 100)),
    ("q99", Fraction(99, 100)),
    ("q995", Fraction(995, 1000)),
)
_SHOWN = {
    "positive": frozenset({"q50", "q95", "q98", "q99"}),
    "negative": frozenset({"q01", "q02", "q05", "q50"}),
    "both": frozenset({"q005", "q01", "q025", "q50", "q975", "q99", "q995"}),
}
_TRIMS = (
    ("1", {"positive": "q99", "negative": "q01", "both": ("q005", "q995")}),
    ("2", {"positive": "q98", "negative": "q02", "both": ("q01", "q99")}),
    ("5", {"positive": "q95", "negative": "q05", "both": ("q025", "q975")}),
)
_HIST_FIELDS = ["column", "bin_left", "bin_right", "n"]


def compute_model_stats(
    input_dir: Path,
    output_dir: Path,
    *,
    otari: bool = False,
    seqweaver: bool = False,
    column_names: Path | None = None,
) -> int:
    """Write per-chromosome column statistics for model score tables.

    ``chrN/columns.by_row.tsv`` uses each input row as one observation.
    For Otari, ``chrN/columns.by_variant.tsv`` first takes the plain mean of
    the stored values for each variant. A Seqweaver matrix has one row per
    variant, so those two tables match. Histograms are written beside them.
    A chromosome with no data rows is omitted.
    """
    if otari == seqweaver:
        raise ValueError("set exactly one of --otari and --seqweaver")
    if seqweaver and column_names is None:
        raise ValueError("--column-names is required with --seqweaver")
    if column_names is not None and not seqweaver:
        raise ValueError("--column-names requires --seqweaver")
    if not input_dir.is_dir():
        raise ValueError(f"input directory not found: {input_dir}")
    names: list[str] = []
    if otari:
        chromosomes = _discover(input_dir)
        if not chromosomes:
            raise ValueError(f"no {_OTARI_NAME} files under {input_dir}")
    else:
        names = _load_column_names(column_names)
        chromosomes = _discover_seqweaver(input_dir)
        if not chromosomes:
            raise ValueError(f"no seqweaver annot.*diffs.h5 files under {input_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    written = 0
    try:
        for chrom, files in chromosomes:
            partial = output_dir / ".partial" / chrom
            shutil.rmtree(partial, ignore_errors=True)
            partial.mkdir(parents=True)
            try:
                loaded = _spill(files, partial) if otari else _spill_seqweaver(files, partial, names)
                if loaded[1] == 0:
                    continue
                try:
                    _emit(output_dir / chrom, loaded)
                except Exception:
                    shutil.rmtree(output_dir / chrom, ignore_errors=True)
                    raise
                written += 1
                print(f"finished processing {chrom}", flush=True)
            finally:
                shutil.rmtree(partial, ignore_errors=True)
    finally:
        shutil.rmtree(output_dir / ".partial", ignore_errors=True)
    return written


def _fields(count_name: str) -> list[str]:
    names = ["column", "n", count_name, "sign", "min", "max"]
    names.extend(name for name, _probability in _QUANTILES)
    names.extend(["mean", "std", "mean_1", "std_1", "mean_2", "std_2", "mean_5", "std_5"])
    return names


def _discover(input_dir: Path) -> list[tuple[str, list[Path]]]:
    found: list[tuple[str, list[Path]]] = []
    for chrom_dir in _chrom_dirs(input_dir):
        files = []
        for shard in chrom_dir.iterdir():
            if not shard.is_dir():
                continue
            path = shard / _OTARI_NAME
            if path.is_file():
                files.append(path)
        if files:
            found.append((chrom_dir.name, sorted(files, key=_shard_key)))
    return found


def _discover_seqweaver(input_dir: Path) -> list[tuple[str, list[tuple[Path, Path]]]]:
    found: list[tuple[str, list[tuple[Path, Path]]]] = []
    for chrom_dir in _chrom_dirs(input_dir):
        pairs: list[tuple[Path, Path]] = []
        for path in chrom_dir.rglob("annot.*diffs.h5"):
            if path.is_file() and path.parent.name == "seqweaver_human":
                pairs.append((path, _row_labels_path(path)))
        if pairs:
            found.append((chrom_dir.name, sorted(pairs, key=lambda item: item[0].as_posix())))
    return found


def _row_labels_path(diffs: Path) -> Path:
    name = diffs.name
    if not name.endswith("_diffs.h5"):
        raise ValueError(f"{diffs} is not an annot.*diffs.h5 file")
    labels = diffs.with_name(name[: -len("_diffs.h5")] + "_row_labels.txt")
    if not labels.is_file():
        raise ValueError(f"{diffs} is missing {labels.name}")
    return labels


def _load_column_names(path: Path | None) -> list[str]:
    if path is None or not path.is_file():
        raise ValueError(f"column names file not found: {path}")
    names = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not names:
        raise ValueError(f"{path} has no column names")
    if len(names) != len(set(names)):
        raise ValueError(f"{path} has duplicate column names")
    return names


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


def _shard_key(path: Path) -> tuple[int, int, int, str]:
    match = _SHARD_DIR.fullmatch(path.parent.name)
    if match is None:
        return (1, 0, 0, path.parent.name)
    return (0, int(match.group(1)), int(match.group(2)), path.parent.name)


def _spill(
    files: list[Path],
    partial: Path,
) -> tuple[list[str], int, int, list[Path], list[Path]]:
    columns: list[str] | None = None
    value_handles = []
    index_handles = []
    pending_values: list[list[float]] = []
    pending_index: list[list[int]] = []
    variants: dict[str, int] = {}
    n_rows = 0
    try:
        for path in files:
            with path.open(encoding="utf-8", newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader, None)
                if not header:
                    raise ValueError(f"{path} is empty")
                header = [name.strip() for name in header]
                if len(header) != len(set(header)):
                    raise ValueError(f"{path} has duplicate columns")
                missing = [name for name in _ID_COLUMNS if name not in header]
                if missing:
                    raise ValueError(f"{path} is missing {', '.join(missing)}")
                score_columns = [name for name in header if name not in _ID_COLUMNS]
                if columns is None:
                    columns = score_columns
                    for index, _name in enumerate(columns):
                        value_handles.append((partial / f"{index}.values").open("wb"))
                        index_handles.append((partial / f"{index}.index").open("wb"))
                        pending_values.append([])
                        pending_index.append([])
                elif score_columns != columns:
                    raise ValueError(f"{path} columns differ from {files[0]}")
                variant_at = header.index("variant_id")
                score_at = [header.index(name) for name in columns]
                for line_number, row in enumerate(reader, start=2):
                    if not row or all(not cell.strip() for cell in row):
                        continue
                    if len(row) != len(header):
                        raise ValueError(f"{path}:{line_number} expected {len(header)} columns")
                    try:
                        variant_id = _standard_variant_id(row[variant_at])
                    except ValueError as exc:
                        raise ValueError(f"{path}:{line_number} {exc}") from exc
                    variant_index = variants.setdefault(variant_id, len(variants))
                    n_rows += 1
                    for column_index, source in enumerate(score_at):
                        text = row[source].strip()
                        if not text:
                            continue
                        try:
                            value = float(text)
                        except ValueError as exc:
                            raise ValueError(
                                f"{path}:{line_number} column {columns[column_index]} value {text!r} is not a finite number"
                            ) from exc
                        if not math.isfinite(value):
                            raise ValueError(
                                f"{path}:{line_number} column {columns[column_index]} value {text!r} is not a finite number"
                            )
                        pending_values[column_index].append(value)
                        pending_index[column_index].append(variant_index)
                        if len(pending_values[column_index]) >= _FLUSH:
                            _flush_column(
                                value_handles[column_index],
                                index_handles[column_index],
                                pending_values[column_index],
                                pending_index[column_index],
                            )
        for index, handle in enumerate(value_handles):
            _flush_column(handle, index_handles[index], pending_values[index], pending_index[index])
    finally:
        for handle in value_handles:
            handle.close()
        for handle in index_handles:
            handle.close()
    if columns is None:
        raise ValueError(f"{files[0]} is empty")
    value_paths = [partial / f"{index}.values" for index in range(len(columns))]
    index_paths = [partial / f"{index}.index" for index in range(len(columns))]
    return columns, n_rows, len(variants), value_paths, index_paths


def _spill_seqweaver(
    files: list[tuple[Path, Path]],
    partial: Path,
    columns: list[str],
) -> tuple[list[str], int, int, list[Path], list[Path]]:
    try:
        import h5py
    except ImportError as exc:
        raise ValueError("h5py is required for --seqweaver") from exc
    n_rows = 0
    widths: list[tuple[Path, int]] = []
    for diffs, labels in files:
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
            width = int(data.shape[0])
        _check_row_labels(labels, width)
        widths.append((diffs, width))
        n_rows += width
    value_handles = [(partial / f"{index}.values").open("wb") for index in range(len(columns))]
    index_handles = [(partial / f"{index}.index").open("wb") for index in range(len(columns))]
    try:
        row_base = 0
        for diffs, width in widths:
            with h5py.File(diffs, "r") as handle:
                data = handle["data"]
                for start in range(0, width, _FLUSH):
                    block = np.asarray(data[start : start + _FLUSH], dtype=np.float64)
                    for column_index in range(len(columns)):
                        values = block[:, column_index]
                        finite = np.isfinite(values)
                        if not finite.any():
                            continue
                        chosen = np.flatnonzero(finite)
                        values[chosen].astype(np.float64, copy=False).tofile(value_handles[column_index])
                        (row_base + start + chosen).astype(np.int32).tofile(index_handles[column_index])
            row_base += width
    finally:
        for handle in value_handles:
            handle.close()
        for handle in index_handles:
            handle.close()
    value_paths = [partial / f"{index}.values" for index in range(len(columns))]
    index_paths = [partial / f"{index}.index" for index in range(len(columns))]
    return columns, n_rows, n_rows, value_paths, index_paths


def _check_row_labels(path: Path, n_rows: int) -> None:
    lines = sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line.strip())
    if lines != n_rows and lines != n_rows + 1:
        raise ValueError(f"{path} has {lines} label rows for {n_rows} variants")


def _flush_column(values_handle, index_handle, values: list[float], indices: list[int]) -> None:
    if not values:
        return
    np.asarray(values, dtype=np.float64).tofile(values_handle)
    np.asarray(indices, dtype=np.int32).tofile(index_handle)
    values.clear()
    indices.clear()


def _emit(
    chrom_dir: Path,
    loaded: tuple[list[str], int, int, list[Path], list[Path]],
) -> None:
    columns, n_rows, n_variants, value_paths, index_paths = loaded
    chrom_dir.mkdir(parents=True, exist_ok=True)
    row_stats = []
    variant_stats = []
    row_hist = []
    variant_hist = []
    for name, value_path, index_path in zip(columns, value_paths, index_paths):
        values = np.fromfile(value_path, dtype=np.float64)
        indices = np.fromfile(index_path, dtype=np.int32)
        if values.size != indices.size:
            raise ValueError(f"column {name} spill is truncated")
        row_summary, row_bins = _summarize(values)
        row_stats.append(_stat_record(name, n_rows, row_summary))
        row_hist.extend(_hist_records(name, row_bins))
        means = _variant_means(values, indices, n_variants)
        variant_summary, variant_bins = _summarize(means)
        variant_stats.append(_stat_record(name, n_variants, variant_summary))
        variant_hist.extend(_hist_records(name, variant_bins))
    _write_table(chrom_dir / "columns.by_row.tsv", _fields("n_rows"), row_stats)
    _write_table(chrom_dir / "columns.by_variant.tsv", _fields("n_variants"), variant_stats)
    _write_table(chrom_dir / "columns.by_row.hist.tsv", _HIST_FIELDS, row_hist)
    _write_table(chrom_dir / "columns.by_variant.hist.tsv", _HIST_FIELDS, variant_hist)


def _variant_means(values: np.ndarray, indices: np.ndarray, n_variants: int) -> np.ndarray:
    if values.size == 0:
        return values
    sums = np.bincount(indices, weights=values, minlength=n_variants)
    counts = np.bincount(indices, minlength=n_variants)
    kept = counts > 0
    return sums[kept] / counts[kept]


def _summarize(values: np.ndarray) -> tuple[dict[str, object], list[tuple[float | None, float | None, int]]]:
    summary: dict[str, object] = {"n": int(values.size)}
    if values.size == 0:
        summary["sign"] = ""
        return summary, []
    ordered = np.sort(values)
    sign = _sign(float(ordered[0]), float(ordered[-1]))
    quantiles = {name: float(ordered[_quantile_index(ordered.size, probability)]) for name, probability in _QUANTILES}
    summary["sign"] = sign
    summary["min"] = float(ordered[0])
    summary["max"] = float(ordered[-1])
    summary["quantiles"] = quantiles
    summary["mean"] = float(ordered.mean())
    summary["std"] = _std(ordered)
    for label, bounds in _TRIMS:
        kept = _trimmed(ordered, sign, quantiles, bounds)
        summary[f"mean_{label}"] = float(kept.mean()) if kept.size else None
        summary[f"std_{label}"] = _std(kept)
    return summary, _histogram(ordered, sign, quantiles)


def _sign(low: float, high: float) -> str:
    if low < 0 < high:
        return "both"
    if low < 0:
        return "negative"
    return "positive"


def _quantile_index(n: int, probability: Fraction) -> int:
    """Smallest index whose empirical rank is at least ``probability``."""
    ceil = (probability.numerator * n + probability.denominator - 1) // probability.denominator
    return min(max(ceil - 1, 0), n - 1)


def _trimmed(
    ordered: np.ndarray,
    sign: str,
    quantiles: dict[str, float],
    bounds: dict[str, str | tuple[str, str]],
) -> np.ndarray:
    rule = bounds[sign]
    if sign == "positive":
        assert isinstance(rule, str)
        return ordered[ordered <= quantiles[rule]]
    if sign == "negative":
        assert isinstance(rule, str)
        return ordered[ordered >= quantiles[rule]]
    assert isinstance(rule, tuple)
    low, high = rule
    return ordered[(ordered >= quantiles[low]) & (ordered <= quantiles[high])]


def _std(values: np.ndarray) -> float | None:
    if values.size < 2:
        return None
    return float(np.std(values, ddof=1))


def _histogram(
    ordered: np.ndarray,
    sign: str,
    quantiles: dict[str, float],
) -> list[tuple[float | None, float | None, int]]:
    if sign == "positive":
        left, right = 0.0, quantiles["q99"]
        open_low, open_high = False, True
    elif sign == "negative":
        left, right = quantiles["q01"], 0.0
        open_low, open_high = True, False
    else:
        left, right = quantiles["q01"], quantiles["q99"]
        open_low, open_high = True, True
    rows: list[tuple[float | None, float | None, int]] = []
    if left == right:
        inside = int(np.sum(ordered == left))
        bulk = [(left, right, inside)]
    else:
        counts, edges = np.histogram(ordered, bins=_BINS, range=(left, right))
        bulk = [
            (float(edges[index]), float(edges[index + 1]), int(counts[index]))
            for index in range(counts.size)
        ]
    if open_low:
        rows.append((float("-inf"), left, int(np.sum(ordered < left))))
    rows.extend(bulk)
    if open_high:
        rows.append((right, float("inf"), int(np.sum(ordered > right))))
    return rows


def _stat_record(column: str, count: int, summary: dict[str, object]) -> dict[str, object]:
    record: dict[str, object] = {
        "column": column,
        "n": summary["n"],
        "count": count,
        "sign": summary.get("sign", ""),
        "min": summary.get("min"),
        "max": summary.get("max"),
        "mean": summary.get("mean"),
        "std": summary.get("std"),
    }
    shown = _SHOWN.get(str(summary.get("sign", "")), frozenset())
    quantiles = summary.get("quantiles", {})
    assert isinstance(quantiles, dict)
    for name, _probability in _QUANTILES:
        record[name] = quantiles.get(name) if name in shown else None
    for label, _bounds in _TRIMS:
        record[f"mean_{label}"] = summary.get(f"mean_{label}")
        record[f"std_{label}"] = summary.get(f"std_{label}")
    return record


def _hist_records(
    column: str,
    bins: list[tuple[float | None, float | None, int]],
) -> list[dict[str, object]]:
    return [
        {"column": column, "bin_left": left, "bin_right": right, "n": count}
        for left, right, count in bins
    ]


def _write_table(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    count_name = fieldnames[2] if fieldnames[2] in {"n_rows", "n_variants"} else None
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            output = {}
            for name in fieldnames:
                if name == count_name:
                    output[name] = row["count"]
                elif name in {"column", "n", "sign"}:
                    output[name] = "" if row.get(name) is None else row[name]
                elif name == "bin_left" or name == "bin_right":
                    output[name] = _format_edge(row[name])
                else:
                    output[name] = _format_number(row.get(name))
            writer.writerow(output)


def _format_number(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
        return str(int(value))
    return format(float(value), ".16g")


def _format_edge(value: object) -> str:
    if value is None:
        return ""
    number = float(value)
    if math.isinf(number):
        return "inf" if number > 0 else "-inf"
    return format(number, ".16g")


def _standard_variant_id(raw: str) -> str:
    """Map ``21_10414931_C_T_hg38`` to ``chr21_10414931_C_T``."""
    text = raw.strip()
    if text.endswith("_hg38"):
        text = text[: -len("_hg38")]
    contig, separator, rest = text.partition("_")
    if not separator or not rest or not contig:
        raise ValueError(f"variant id {raw!r} is not contig_pos_ref_alt")
    if not contig.startswith("chr"):
        contig = "chr" + contig
    return f"{contig}_{rest}"
