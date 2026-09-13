from pathlib import Path

from pheno.mapping import load_column_map


def test_default_map():
    cmap = load_column_map(None)
    assert cmap.person == "person"
    assert cmap.mother_sequenced is None
    assert cmap.asd_status("2") == "asd"
    assert cmap.asd_status("True") == "asd"
    assert cmap.asd_status("1") == "ctrl"
    assert cmap.asd_status("False") == "ctrl"
    assert cmap.asd_status("") is None
    assert cmap.sex_status("Male") == "male"
    assert cmap.sex_status("1") == "male"
    assert cmap.sex_status("Female") == "female"
    assert cmap.sex_status("2") == "female"
    assert cmap.is_bad_id("0")
    assert cmap.is_bad_id("-")
    assert cmap.is_bad_id("false")
    assert cmap.is_bad_id("")
    assert not cmap.is_bad_id("mom1")


def test_custom_map(tmp_path: Path):
    path = tmp_path / "map.json"
    path.write_text(
        """
        {
          "columns": {
            "person": "IID",
            "family": "FID",
            "mother": "MothID",
            "father": "FathID",
            "mother_sequenced": "-",
            "father_sequenced": "dad_seq",
            "sex": "SEX",
            "asd": "AFF"
          }
        }
        """
    )
    cmap = load_column_map(path)
    assert cmap.person == "IID"
    assert cmap.mother == "MothID"
    assert cmap.mother_sequenced is None
    assert cmap.father_sequenced == "dad_seq"
    assert cmap.is_sequenced("TRUE")
