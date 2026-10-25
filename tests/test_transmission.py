import csv
import math
from pathlib import Path

import pytest

from pheno.cli_transmission import main
from pheno.gene_burden import load_analysis_groups
from pheno.mapping import load_column_map
from pheno.transmission import compute_transmission

FAMILY = """\
person	family	mother	father	sex	asd
A	F1	mom1	dad1	male	2
B	F1	mom1	dad1	female	1
C	F2	mom2	dad2	male	2
D	F2	mom2	dad2	female	1
E	F3	mom3	dad3	male	1
U	F1	mom1	dad1	male	1
"""
HEADER = "#CHROM\tPOS\tID\tREF\tALT\tTRIO_CALLS\n"


def _write(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HEADER + body, encoding="utf-8")


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _run(tmp_path: Path, body: str) -> list[dict[str, str]]:
    source = tmp_path / "results"
    _write(source / "chr21" / "inherited_00000.tsv", body)
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    output = tmp_path / "transmission.tsv"
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    compute_transmission(source, "inherited", people, output)
    return _rows(output)


def test_pair_and_asymmetric_family(tmp_path: Path):
    rows = _run(
        tmp_path,
        "chr21\t10\t.\tG\tA\tA=0/0|0/1|0/1|99;B=0/1|0/0|0/0|99;"
        "C=0/0|0/1|0/1|99;D=0/1|0/0|0/0|99\n"
        "chr21\t20\t.\tT\tC\tA=0/1|0/1|0/0|99;B=0/0|0/1|0/1|99;U=0/0|0/1|0/1|99\n"
        "chr21\t30\t.\tC\tG\tA=0/0|0/1|0/1|99;B=0/1|0/0|0/1|99\n",
    )
    by_pos = {row["POS"]: row for row in rows}
    assert set(by_pos) == {"10", "20", "30"}

    pair = by_pos["10"]
    assert pair["n_families"] == "2"
    assert pair["n_informative"] == "2"
    assert float(pair["z"]) == pytest.approx(math.sqrt(2))
    assert float(pair["p_exact"]) == pytest.approx(0.5)

    uneven = by_pos["20"]
    assert uneven["n_families"] == "1"
    assert uneven["n_informative"] == "1"
    assert float(uneven["z"]) == pytest.approx(-math.sqrt(2))
    assert float(uneven["p_exact"]) == pytest.approx(1 / 3)

    tied = by_pos["30"]
    assert tied["n_families"] == "1"
    assert tied["n_informative"] == "0"
    assert tied["z"] == ""
    assert float(tied["p_exact"]) == pytest.approx(1)


def test_missing_sibling_is_not_reference(tmp_path: Path):
    rows = _run(tmp_path, "chr21\t10\t.\tG\tA\tA=0/0|0/1|0/1|99\n")
    assert rows == []


def test_multiallelic_variant_is_skipped(tmp_path: Path):
    rows = _run(
        tmp_path,
        "chr21\t10\t.\tG\tA,C\tA=0/0|0/1|0/1|99;B=0/1|0/0|0/0|99\n"
        "chr21\t20\t.\tG\tA\tA=0/0|0/1|1/2|99;B=0/1|0/0|0/0|99\n"
        "chr21\t30\t.\tG\tA\tA=0/0|1/1|0/1|99;B=1/1|0/0|0/0|99\n",
    )
    assert [row["POS"] for row in rows] == ["30"]
    assert float(rows[0]["z"]) == pytest.approx(1)
    assert float(rows[0]["p_exact"]) == pytest.approx(1)


def test_three_field_genotype_uses_child(tmp_path: Path):
    rows = _run(
        tmp_path,
        "chr21\t10\t.\tG\tA\tA=0/1|0/1|99;B=0/0|0/0|99\n",
    )
    assert rows[0]["n_informative"] == "1"
    assert float(rows[0]["z"]) == pytest.approx(1)


def test_patients_column_is_rejected(tmp_path: Path, capsys):
    source = tmp_path / "results" / "chr21"
    source.mkdir(parents=True)
    (source / "inherited_00000.tsv").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\nchr21\t10\t.\tG\tA\tA;B\n",
        encoding="utf-8",
    )
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    output = tmp_path / "transmission.tsv"
    code = main(
        [
            "--family-file",
            str(family),
            "--input-dir",
            str(tmp_path / "results"),
            "--file-pattern",
            "inherited",
            "--output",
            str(output),
        ]
    )
    assert code == 1
    assert "requires a TRIO_CALLS column" in capsys.readouterr().err
    assert not output.exists()


def test_output_dir_writes_one_sorted_file_per_chromosome(tmp_path: Path):
    source = tmp_path / "results"
    pair = "A=0/0|0/1|0/1|99;B=0/1|0/0|0/0|99"
    _write(
        source / "chr21" / "inherited_200_300.tsv",
        f"chr21\t250\t.\tG\tA\t{pair}\nchr21\t30\t.\tC\tT\t{pair}\n",
    )
    _write(
        source / "chr22" / "inherited_10_20.tsv",
        "chr22\t15\t.\tA\tG\tC=0/0|0/1|0/1|99;D=0/1|0/0|0/0|99\n",
    )
    _write(source / "chrX" / "inherited_1_2.tsv", "chrX\t5\t.\tA\tG\tA=0/0|0/1|0/1|99\n")
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    output_dir = tmp_path / "by-chrom"
    written = compute_transmission(source, "inherited", people, output_dir=output_dir)
    assert written == 3
    assert [row["POS"] for row in _rows(output_dir / "chr21" / "inherited.tsv")] == ["30", "250"]
    assert [row["POS"] for row in _rows(output_dir / "chr22" / "inherited.tsv")] == ["15"]
    assert not (output_dir / "chrX").exists()


def test_output_and_output_dir_are_exclusive(tmp_path: Path, capsys):
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    shared = [
        "--family-file",
        str(family),
        "--input-dir",
        str(tmp_path),
        "--file-pattern",
        "inherited",
    ]
    neither = main(shared)
    assert neither == 1
    assert "exactly one of --output and --output-dir" in capsys.readouterr().err
    both = main([*shared, "--output", "out.tsv", "--output-dir", "out"])
    assert both == 1
    assert "exactly one of --output and --output-dir" in capsys.readouterr().err
