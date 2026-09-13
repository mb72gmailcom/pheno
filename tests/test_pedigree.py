from pathlib import Path

import pytest

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
H	F4	0	dad4	Male	2
I	F5	mom5	-	Female	1
J	F6	mom6	dad6	Male	
mom1	F0	-	-	Female	1
"""


def _write_family(tmp_path: Path, text: str = FAMILY) -> Path:
    path = tmp_path / "family.tsv"
    path.write_text(text)
    return path


def test_cohort_membership(tmp_path: Path):
    cohort = load_cohort(_write_family(tmp_path), load_column_map(None))
    assert set(cohort.person_to_idx) == {"A", "B", "C", "D", "E", "F", "G"}
    assert cohort.n_asd == 4  # A, C, E, F
    assert cohort.n_ctrl == 3  # B, D, G
    assert cohort.n_male_asd == 3  # A, C, F
    assert cohort.n_female_asd == 1  # E
    assert cohort.n_male_ctrl == 1  # D (sex 1)
    assert cohort.n_female_ctrl == 2  # B, G
    assert cohort.n_informative_families == 3  # F1, F2, F3
    assert cohort.n_children == 7


def test_missing_asd_and_bad_parents_excluded(tmp_path: Path):
    cohort = load_cohort(_write_family(tmp_path), load_column_map(None))
    assert "H" not in cohort.person_to_idx
    assert "I" not in cohort.person_to_idx
    assert "J" not in cohort.person_to_idx
    assert "mom1" not in cohort.person_to_idx


def test_sequenced_filter(tmp_path: Path):
    text = """\
person	family	mother	father	mother_sequenced	father_sequenced	sex	asd
A	F1	mom1	dad1	True	True	Male	2
B	F1	mom1	dad1	False	True	Female	1
"""
    path = tmp_path / "family.tsv"
    path.write_text(text)
    cmap_path = tmp_path / "map.json"
    cmap_path.write_text(
        """
        {"columns": {"mother_sequenced": "mother_sequenced", "father_sequenced": "father_sequenced"}}
        """
    )
    cohort = load_cohort(path, load_column_map(cmap_path))
    assert set(cohort.person_to_idx) == {"A"}


def test_empty_cohort_raises(tmp_path: Path):
    path = tmp_path / "family.tsv"
    path.write_text("person\tfamily\tmother\tfather\tsex\tasd\nX\tF\t0\t-\tMale\t2\n")
    with pytest.raises(ValueError, match="no children"):
        load_cohort(path, load_column_map(None))
