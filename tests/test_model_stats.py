import csv
import statistics
from pathlib import Path

import pytest

from pheno.cli_model_stats import main
from pheno.model_stats import compute_model_stats


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _otari(root: Path, chrom: str, shard: str, text: str) -> None:
    _write(root / chrom / shard / "variant_effects_comprehensive.tsv", text)


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _by_column(path: Path) -> dict[str, dict[str, str]]:
    return {row["column"]: row for row in _rows(path)}


def test_plain_mean_across_transcripts(tmp_path: Path) -> None:
    source = tmp_path / "otari"
    _otari(
        source,
        "chr21",
        "1_100",
        "\n".join(
            [
                "variant_id\ttranscript_id\tmax_effect\tBrain\tLiver",
                "21_10_A_G_hg38\tt1\t9\t2\t-2",
                "21_10_A_G_hg38\tt2\t9\t0\t-4",
                "21_20_C_T_hg38\tt1\t1\t1\t1",
                "",
            ]
        ),
    )
    output = tmp_path / "stats"
    assert compute_model_stats(source, output, otari=True) == 1
    by_row = _by_column(output / "chr21" / "columns.by_row.tsv")
    by_variant = _by_column(output / "chr21" / "columns.by_variant.tsv")
    assert list(by_row) == ["max_effect", "Brain", "Liver"]
    assert by_row["Brain"]["n"] == "3"
    assert by_row["Brain"]["n_rows"] == "3"
    assert by_row["Brain"]["sign"] == "positive"
    assert float(by_row["Brain"]["mean"]) == pytest.approx(1)
    assert by_row["Brain"]["q01"] == ""
    assert by_row["max_effect"]["n_rows"] == "3"
    assert float(by_row["max_effect"]["mean"]) == pytest.approx(19 / 3)
    assert by_row["Liver"]["sign"] == "both"
    assert float(by_row["Liver"]["mean"]) == pytest.approx(-5 / 3)
    assert by_variant["Brain"]["n"] == "2"
    assert by_variant["Brain"]["n_variants"] == "2"
    assert float(by_variant["Brain"]["mean"]) == pytest.approx(1)
    assert float(by_variant["max_effect"]["mean"]) == pytest.approx(5)
    assert float(by_variant["Liver"]["mean"]) == pytest.approx(-1)
    assert by_variant["Liver"]["sign"] == "both"
    assert not (output / ".partial").exists()


def test_positive_trim_and_histogram(tmp_path: Path) -> None:
    lines = ["variant_id\ttranscript_id\tscore"]
    lines.extend(f"21_{value}_A_C_hg38\tt1\t{value}" for value in range(1, 101))
    source = tmp_path / "otari"
    _otari(source, "chr21", "1_100", "\n".join(lines) + "\n")
    output = tmp_path / "stats"
    compute_model_stats(source, output, otari=True)
    row = _by_column(output / "chr21" / "columns.by_row.tsv")["score"]
    variant = _by_column(output / "chr21" / "columns.by_variant.tsv")["score"]
    assert row["n_rows"] == variant["n_variants"] == "100"
    assert row["sign"] == "positive"
    assert float(row["min"]) == 1
    assert float(row["max"]) == 100
    assert float(row["q50"]) == 50
    assert float(row["q95"]) == 95
    assert float(row["q98"]) == 98
    assert float(row["q99"]) == 99
    assert row["q01"] == ""
    values = list(range(1, 101))
    assert float(row["mean"]) == pytest.approx(statistics.fmean(values))
    assert float(row["std"]) == pytest.approx(statistics.stdev(values))
    assert float(row["mean_1"]) == pytest.approx(statistics.fmean(range(1, 100)))
    assert float(row["std_1"]) == pytest.approx(statistics.stdev(range(1, 100)))
    assert float(row["mean_2"]) == pytest.approx(statistics.fmean(range(1, 99)))
    assert float(row["mean_5"]) == pytest.approx(statistics.fmean(range(1, 96)))
    assert row["mean"] == variant["mean"]
    hist = _rows(output / "chr21" / "columns.by_row.hist.tsv")
    assert len(hist) == 51
    assert hist[-1]["bin_right"] == "inf"
    assert int(hist[-1]["n"]) == 1
    assert sum(int(item["n"]) for item in hist) == 100
    assert all(item["bin_left"] != "-inf" for item in hist)


def test_two_sided_trim_uses_both_tails(tmp_path: Path) -> None:
    values = list(range(-100, 100))
    lines = ["variant_id\ttranscript_id\tscore"]
    lines.extend(f"21_{index}_A_C_hg38\tt1\t{value}" for index, value in enumerate(values))
    source = tmp_path / "otari"
    _otari(source, "chr2", "10_20", "\n".join(lines) + "\n")
    output = tmp_path / "stats"
    compute_model_stats(source, output, otari=True)
    row = _by_column(output / "chr2" / "columns.by_row.tsv")["score"]
    assert row["sign"] == "both"
    assert row["q95"] == ""
    assert float(row["q005"]) == -100
    assert float(row["q01"]) == -99
    assert float(row["q50"]) == -1
    assert float(row["q99"]) == 97
    assert float(row["q995"]) == 98
    kept = [value for value in values if -99 <= value <= 97]
    assert float(row["mean_2"]) == pytest.approx(statistics.fmean(kept))
    assert float(row["std_2"]) == pytest.approx(statistics.stdev(kept))
    hist = _rows(output / "chr2" / "columns.by_row.hist.tsv")
    assert hist[0]["bin_left"] == "-inf"
    assert hist[-1]["bin_right"] == "inf"
    assert len(hist) == 52
    assert sum(int(item["n"]) for item in hist) == 200


def test_negative_column_trims_the_low_tail(tmp_path: Path) -> None:
    values = list(range(-100, 0))
    lines = ["variant_id\ttranscript_id\tscore"]
    lines.extend(f"21_{index}_A_C_hg38\tt1\t{value}" for index, value in enumerate(values))
    source = tmp_path / "otari"
    _otari(source, "chr21", "1_2", "\n".join(lines) + "\n")
    output = tmp_path / "stats"
    compute_model_stats(source, output, otari=True)
    row = _by_column(output / "chr21" / "columns.by_row.tsv")["score"]
    assert row["sign"] == "negative"
    assert row["q99"] == ""
    assert float(row["q01"]) == -100
    assert float(row["q02"]) == -99
    assert float(row["q05"]) == -96
    assert float(row["q50"]) == -51
    kept = [value for value in values if value >= -96]
    assert float(row["mean_5"]) == pytest.approx(statistics.fmean(kept))
    hist = _rows(output / "chr21" / "columns.by_row.hist.tsv")
    assert hist[0]["bin_left"] == "-inf"
    assert hist[-1]["bin_right"] != "inf"
    assert len(hist) == 51
    assert sum(int(item["n"]) for item in hist) == 100


def test_mapped_variant_ids_are_averaged_across_shards(tmp_path: Path) -> None:
    source = tmp_path / "otari"
    _otari(
        source,
        "chrX",
        "2_20",
        "variant_id\ttranscript_id\tscore\nchrX_5_A_C\tt2\t4\n",
    )
    _otari(
        source,
        "chrX",
        "1_10",
        "variant_id\ttranscript_id\tscore\nX_5_A_C_hg38\tt1\t2\n",
    )
    output = tmp_path / "stats"
    compute_model_stats(source, output, otari=True)
    row = _by_column(output / "chrX" / "columns.by_row.tsv")["score"]
    variant = _by_column(output / "chrX" / "columns.by_variant.tsv")["score"]
    assert row["n"] == "2"
    assert row["n_rows"] == "2"
    assert variant["n"] == "1"
    assert variant["n_variants"] == "1"
    assert float(variant["mean"]) == pytest.approx(3)
    assert variant["std"] == ""


def test_blank_cell_keeps_the_row_count(tmp_path: Path) -> None:
    source = tmp_path / "otari"
    _otari(
        source,
        "chr21",
        "1_5",
        "\n".join(
            [
                "variant_id\ttranscript_id\tBrain\tLiver",
                "21_1_A_C_hg38\tt1\t1\t2",
                "21_2_A_C_hg38\tt1\t\t4",
                "",
            ]
        ),
    )
    output = tmp_path / "stats"
    compute_model_stats(source, output, otari=True)
    by_row = _by_column(output / "chr21" / "columns.by_row.tsv")
    assert by_row["Brain"]["n"] == "1"
    assert by_row["Brain"]["n_rows"] == "2"
    assert by_row["Brain"]["std"] == ""
    assert by_row["Liver"]["n"] == "2"
    assert by_row["Liver"]["n_rows"] == "2"
    hist_columns = {row["column"] for row in _rows(output / "chr21" / "columns.by_row.hist.tsv")}
    assert hist_columns == {"Brain", "Liver"}


def test_empty_chromosome_is_omitted(tmp_path: Path) -> None:
    source = tmp_path / "otari"
    _otari(
        source,
        "chr21",
        "1_5",
        "variant_id\ttranscript_id\tscore\n21_1_A_C_hg38\tt1\t1\n",
    )
    _otari(source, "chr22", "1_5", "variant_id\ttranscript_id\tscore\n")
    _otari(source, "chrM", "1_5", "variant_id\ttranscript_id\tscore\nM_1_A_C_hg38\tt1\t5\n")
    output = tmp_path / "stats"
    assert compute_model_stats(source, output, otari=True) == 1
    assert (output / "chr21" / "columns.by_row.tsv").is_file()
    assert not (output / "chr22").exists()
    assert not (output / "chrM").exists()


def test_bad_value_writes_nothing(tmp_path: Path) -> None:
    source = tmp_path / "otari"
    _otari(
        source,
        "chr21",
        "1_5",
        "variant_id\ttranscript_id\tscore\n21_1_A_C_hg38\tt1\tNA\n",
    )
    output = tmp_path / "stats"
    with pytest.raises(ValueError, match="not a finite number"):
        compute_model_stats(source, output, otari=True)
    assert not (output / "chr21" / "columns.by_row.tsv").exists()
    assert not (output / ".partial").exists()


def test_differing_columns_are_rejected(tmp_path: Path) -> None:
    source = tmp_path / "otari"
    _otari(source, "chr21", "1_5", "variant_id\ttranscript_id\tBrain\n21_1_A_C_hg38\tt1\t1\n")
    _otari(
        source,
        "chr21",
        "6_9",
        "variant_id\ttranscript_id\tBrain\tLiver\n21_2_A_C_hg38\tt1\t1\t2\n",
    )
    output = tmp_path / "stats"
    with pytest.raises(ValueError, match="columns differ"):
        compute_model_stats(source, output, otari=True)
    assert not (output / "chr21" / "columns.by_row.tsv").exists()


def test_cli_requires_otari(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    code = main(["--input-dir", str(tmp_path), "--output-dir", str(tmp_path / "out")])
    assert code == 1
    assert "error: set exactly one of --otari and --seqweaver" in capsys.readouterr().err


def test_cli_writes_statistics(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    source = tmp_path / "otari"
    _otari(source, "chr21", "1_5", "variant_id\ttranscript_id\tscore\n21_1_A_C_hg38\tt1\t1\n")
    output = tmp_path / "stats"
    code = main(["--input-dir", str(source), "--output-dir", str(output), "--otari"])
    assert code == 0
    assert "Wrote column statistics for 1 chromosomes" in capsys.readouterr().out
    assert (output / "chr21" / "columns.by_variant.hist.tsv").is_file()


def _seqweaver(path: Path, matrix: list[list[float]], labels: list[str]) -> None:
    import h5py

    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as handle:
        handle.create_dataset("data", data=matrix)
    path.with_name(path.name.replace("_diffs.h5", "_row_labels.txt")).write_text(
        "".join(f"{label}\n" for label in labels),
        encoding="utf-8",
    )


def test_seqweaver_concatenates_slices_and_matches_row_and_variant(tmp_path: Path) -> None:
    source = tmp_path / "seqweaver"
    root = source / "chr21" / "variants_00000"
    _seqweaver(
        root / "variants_00000_slice_024" / "seqweaver_human" / "annot.variants_00000_slice_024_diffs.h5",
        [[1.0, -2.0], [3.0, float("nan")]],
        ["v1", "v2"],
    )
    _seqweaver(
        root / "variants_00000_slice_025" / "seqweaver_human" / "annot.variants_00000_slice_025_diffs.h5",
        [[5.0, 6.0]],
        ["v3"],
    )
    names = tmp_path / "seqweaver.colnames"
    names.write_text("score_a\nscore_b\n", encoding="utf-8")
    output = tmp_path / "stats"
    assert compute_model_stats(source, output, seqweaver=True, column_names=names) == 1
    by_row = _by_column(output / "chr21" / "columns.by_row.tsv")
    by_variant = _by_column(output / "chr21" / "columns.by_variant.tsv")
    assert list(by_row) == ["score_a", "score_b"]
    assert by_row["score_a"]["n"] == "3"
    assert by_row["score_a"]["n_rows"] == "3"
    assert by_row["score_a"]["sign"] == "positive"
    assert float(by_row["score_a"]["mean"]) == pytest.approx(3)
    assert by_row["score_b"]["n"] == "2"
    assert by_row["score_b"]["sign"] == "both"
    assert float(by_row["score_b"]["mean"]) == pytest.approx(2)
    assert by_variant["score_a"]["n_variants"] == "3"
    assert by_variant["score_a"]["mean"] == by_row["score_a"]["mean"]
    assert by_variant["score_b"]["n"] == "2"


def test_seqweaver_requires_column_names(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    code = main(["--input-dir", str(tmp_path), "--output-dir", str(tmp_path / "out"), "--seqweaver"])
    assert code == 1
    assert "error: --column-names is required with --seqweaver" in capsys.readouterr().err
