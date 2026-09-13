from pathlib import Path

import pandas as pd

from pheno.counts import count_directory, count_variants, discover_files
from pheno.mapping import load_column_map
from pheno.pedigree import load_cohort

FAMILY = """\
person	family	mother	father	sex	asd
A	F1	mom1	dad1	Male	2
B	F1	mom1	dad1	Female	1
C	F2	mom2	dad2	male	True
D	F2	mom2	dad2	1	False
E	F3	mom3	dad3	female	2
F	F3	mom3	dad3	Male	2
G	F3	mom3	dad3	Female	1
"""


def _cohort(tmp_path: Path):
    fam = tmp_path / "family.tsv"
    fam.write_text(FAMILY)
    return load_cohort(fam, load_column_map(None))


def test_discover_files(tmp_path: Path):
    (tmp_path / "inherited.tsv").write_text("x")
    (tmp_path / "inherited_1000_2000.tsv").write_text("x")
    (tmp_path / "inherited_0001.tsv").write_text("x")
    (tmp_path / "other.tsv").write_text("x")
    names = [p.name for p in discover_files(tmp_path, "inherited")]
    assert names == ["inherited.tsv", "inherited_0001.tsv", "inherited_1000_2000.tsv"]


def test_carrier_and_mcnemar_counts(tmp_path: Path):
    cohort = _cohort(tmp_path)
    variants = pd.DataFrame(
        {
            "#CHROM": ["chr21"] * 5,
            "POS": ["1", "2", "3", "4", "5"],
            "ID": ["."] * 5,
            "REF": ["C"] * 5,
            "ALT": ["T"] * 5,
            "PATIENTS": [
                "A",
                "A;B",
                "B",
                "E",
                "UNKNOWN;A",
            ],
        }
    )
    out = count_variants(variants, cohort).set_index("POS")

    # A only: 1 ASD male, 0 ctrl. F1 only-ASD.
    assert list(out.loc["1", ["a_all", "b_all", "c_all", "d_all"]]) == [1, 0, 3, 3]
    assert list(out.loc["1", ["a_male", "b_male", "a_female", "b_female"]]) == [1, 0, 0, 0]
    assert out.loc["1", "n_only_asd"] == 1
    assert out.loc["1", "n_only_ctrl"] == 0

    # A and B: both sides of F1 carry → uninformative
    assert list(out.loc["2", ["a_all", "b_all"]]) == [1, 1]
    assert out.loc["2", "n_only_asd"] == 0
    assert out.loc["2", "n_only_ctrl"] == 0

    # B only: F1 only-ctrl
    assert list(out.loc["3", ["a_all", "b_all"]]) == [0, 1]
    assert out.loc["3", "n_only_asd"] == 0
    assert out.loc["3", "n_only_ctrl"] == 1

    # E only: F3 has two ASD and one ctrl; ASD side carries, ctrl does not.
    # No pair picking: one family contribution.
    assert list(out.loc["4", ["a_all", "b_all", "a_female"]]) == [1, 0, 1]
    assert out.loc["4", "n_only_asd"] == 1
    assert out.loc["4", "n_only_ctrl"] == 0

    # UNKNOWN skipped; A counted
    assert out.loc["5", "a_all"] == 1
    assert out.loc["5", "n_only_asd"] == 1


def test_empty_patients_keeps_row(tmp_path: Path):
    cohort = _cohort(tmp_path)
    variants = pd.DataFrame(
        {
            "CHROM": ["chr21"],
            "POS": ["9"],
            "ID": ["."],
            "REF": ["A"],
            "ALT": ["G"],
            "PATIENTS": [""],
        }
    )
    out = count_variants(variants, cohort)
    assert len(out) == 1
    assert out.loc[0, "a_all"] == 0
    assert out.loc[0, "c_all"] == 4
    assert out.loc[0, "n_only_asd"] == 0


def test_count_directory(tmp_path: Path):
    cohort = _cohort(tmp_path)
    header = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"
    (tmp_path / "inherited.tsv").write_text(header + "chr21\t10\t.\tC\tT\tA\n")
    (tmp_path / "inherited_1000_2000.tsv").write_text(header + "chr21\t20\t.\tC\tA\tB;C\n")
    out = count_directory(tmp_path, "inherited", cohort)
    assert list(out["POS"]) == ["10", "20"]
    assert list(out["a_all"]) == [1, 1]
    assert list(out["b_all"]) == [0, 1]
