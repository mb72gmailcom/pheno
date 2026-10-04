import json
from pathlib import Path

import pytest

from pheno.burden import compute_burden, load_people
from pheno.cli_burden import main
from pheno.mapping import load_column_map

FAMILY = """\
person	family	mother	father	sex	asd
A	F1	mom1	dad1	male	2
B	F1	mom1	dad1	female	1
C	F2	mom2	dad2	male	2
"""

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"
OTARI_HEADER = "variant_id\ttranscript_id\tmax_effect\tmean_effect\tBrain\n"


def _write(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")


def _cohort(tmp_path: Path) -> tuple[Path, Path, Path]:
    source = tmp_path / "results"
    otari = tmp_path / "otari"
    _write(
        source / "chr21" / "inherited_25000001_27500000.tsv",
        HEADER
        + "chr21\t25019786\t.\tG\tT\tA;Z\n"
        + "chr21\t25019787\t.\tC\tG\tA\n"
        + "chr21\t25019800\t.\tC\tT\tA;B\n",
    )
    _write(
        source / "chr21" / "denovo_25000001_27500000.tsv",
        HEADER + "chr21\t25019786\t.\tG\tT\tA\n",
    )
    _write(
        otari / "inherited" / "chr21" / "25000001_27500000" / "variant_effects_comprehensive.tsv",
        OTARI_HEADER
        + "21_25019786_G_T_hg38\tENST00000441009\t0.2\t0.1\t-0.2\n"
        + "21_25019786_G_T_hg38\tENST00000999999\t0.8\t0.2\t-0.9\n"
        + "21_25019787_C_G_hg38\tENST00000441009\t0.4\t0.4\t0.1\n",
    )
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    return source, otari, family


def test_burden_collapses_transcripts_and_counts_unscored(tmp_path: Path):
    source, otari, family = _cohort(tmp_path)
    people = load_people(family, load_column_map(None))
    payload = compute_burden(
        source,
        otari,
        "inherited",
        people,
        ["max_effect", "Brain"],
        threshold=0.5,
        use_abs=True,
        transcripts="max",
    )
    assert payload["transcripts"] == "max"
    assert payload["abs"] is True
    assert set(payload["patients"]) == {"A", "B"}
    person_a = payload["patients"]["A"]
    assert person_a["status"] == "asd"
    assert person_a["family_id"] == "F1"
    assert person_a["n_variants"] == 3
    assert person_a["n_unscored"] == 1
    assert person_a["max_effect"]["n_scored"] == 2
    assert person_a["max_effect"]["n_damaging"] == 1
    assert person_a["max_effect"]["fraction_damaging"] == pytest.approx(0.5)
    assert person_a["max_effect"]["sum_effect"] == pytest.approx(1.2)
    assert person_a["max_effect"]["sum_damaging"] == pytest.approx(0.8)
    assert person_a["max_effect"]["max_score"] == pytest.approx(0.8)
    assert person_a["Brain"]["n_damaging"] == 1
    assert person_a["Brain"]["sum_effect"] == pytest.approx(1.0)
    assert person_a["Brain"]["sum_damaging"] == pytest.approx(0.9)
    assert person_a["Brain"]["max_score"] == pytest.approx(0.9)
    person_b = payload["patients"]["B"]
    assert person_b["status"] == "unaffected"
    assert person_b["n_variants"] == 1
    assert person_b["n_unscored"] == 1
    assert person_b["max_effect"]["n_scored"] == 0
    assert person_b["max_effect"]["fraction_damaging"] is None
    assert person_b["max_effect"]["max_score"] is None
    assert person_b["max_effect"]["sum_effect"] == 0.0
    assert person_b["max_effect"]["sum_damaging"] == 0.0


def test_burden_mean_transcripts(tmp_path: Path):
    source, otari, family = _cohort(tmp_path)
    people = load_people(family, load_column_map(None))
    payload = compute_burden(
        source,
        otari,
        "inherited",
        people,
        ["Brain"],
        threshold=0.5,
        use_abs=True,
        transcripts="mean",
    )
    brain = payload["patients"]["A"]["Brain"]
    assert brain["sum_effect"] == pytest.approx(0.65)
    assert brain["n_damaging"] == 1
    assert brain["sum_damaging"] == pytest.approx(0.55)
    assert brain["max_score"] == pytest.approx(0.55)


def test_burden_cli_no_abs(tmp_path: Path):
    source, otari, family = _cohort(tmp_path)
    output_dir = tmp_path / "burden"
    rc = main(
        [
            "--family-file",
            str(family),
            "--input-dir",
            str(source),
            "--otari-dir",
            str(otari),
            "--output-dir",
            str(output_dir),
            "--file-pattern",
            "inherited",
            "--otari-columns",
            "Brain",
            "--threshold",
            "0.5",
            "--no-abs",
            "--max-transcripts",
        ]
    )
    assert rc == 0
    payload = json.loads((output_dir / "burden.json").read_text(encoding="utf-8"))
    assert payload["abs"] is False
    brain = payload["patients"]["A"]["Brain"]
    assert brain["max_score"] == pytest.approx(0.1)
    assert brain["sum_effect"] == pytest.approx(-0.1)
    assert brain["n_damaging"] == 0
    chromosome = json.loads((output_dir / "chr21" / "burden.json").read_text(encoding="utf-8"))
    assert chromosome["patients"]["A"]["Brain"]["max_score"] == pytest.approx(0.1)


def test_burden_sums_chromosomes(tmp_path: Path):
    source, otari, family = _cohort(tmp_path)
    _write(
        source / "chr22" / "inherited_1_100.tsv",
        HEADER + "chr22\t10\t.\tA\tG\tA;B\n" + "chr22\t11\t.\tT\tC\tC\n",
    )
    _write(
        otari / "inherited" / "chr22" / "1_100" / "variant_effects_comprehensive.tsv",
        OTARI_HEADER
        + "22_10_A_G_hg38\tENST00000441009\t0.2\t0.2\t0.2\n"
        + "22_11_T_C_hg38\tENST00000441009\t0.9\t0.9\t0.9\n",
    )
    output_dir = tmp_path / "out"
    people = load_people(family, load_column_map(None))
    payload = compute_burden(
        source,
        otari,
        "inherited",
        people,
        ["max_effect"],
        threshold=0.5,
        use_abs=True,
        transcripts="max",
        output_dir=output_dir,
    )

    chr21 = json.loads((output_dir / "chr21" / "burden.json").read_text(encoding="utf-8"))
    chr22 = json.loads((output_dir / "chr22" / "burden.json").read_text(encoding="utf-8"))
    assert set(chr21["patients"]) == {"A", "B"}
    assert chr21["patients"]["A"]["n_variants"] == 3
    assert chr21["patients"]["A"]["max_effect"]["fraction_damaging"] == pytest.approx(0.5)
    assert set(chr22["patients"]) == {"A", "B", "C"}
    assert chr22["patients"]["A"]["n_variants"] == 1
    assert chr22["patients"]["A"]["max_effect"]["n_damaging"] == 0
    assert chr22["patients"]["A"]["max_effect"]["fraction_damaging"] == pytest.approx(0.0)
    assert chr22["patients"]["C"]["max_effect"]["max_score"] == pytest.approx(0.9)

    person_a = payload["patients"]["A"]
    assert person_a["n_variants"] == 4
    assert person_a["n_unscored"] == 1
    effect = person_a["max_effect"]
    assert effect["n_scored"] == 3
    assert effect["n_damaging"] == 1
    assert effect["fraction_damaging"] == pytest.approx(1 / 3)
    assert effect["sum_effect"] == pytest.approx(1.4)
    assert effect["sum_damaging"] == pytest.approx(0.8)
    assert effect["max_score"] == pytest.approx(0.8)

    person_b = payload["patients"]["B"]["max_effect"]
    assert payload["patients"]["B"]["n_variants"] == 2
    assert payload["patients"]["B"]["n_unscored"] == 1
    assert person_b["n_scored"] == 1
    assert person_b["n_damaging"] == 0
    assert person_b["fraction_damaging"] == pytest.approx(0.0)
    assert person_b["max_score"] == pytest.approx(0.2)
    assert set(payload["patients"]) == {"A", "B", "C"}

    on_disk = json.loads((output_dir / "burden.json").read_text(encoding="utf-8"))
    assert on_disk["patients"]["A"]["max_effect"]["fraction_damaging"] == pytest.approx(1 / 3)


def test_burden_otari_prefix_can_differ(tmp_path: Path):
    source, otari, family = _cohort(tmp_path)
    inherited = (source / "chr21" / "inherited_25000001_27500000.tsv").read_text(encoding="utf-8")
    _write(source / "chr21" / "inherited_asd_25000001_27500000.tsv", inherited)
    people = load_people(family, load_column_map(None))
    payload = compute_burden(
        source,
        otari,
        "inherited_asd",
        people,
        ["max_effect"],
        threshold=0.5,
        use_abs=True,
        transcripts="max",
        otari_prefix="inherited",
    )
    effect = payload["patients"]["A"]["max_effect"]
    assert effect["n_scored"] == 2
    assert effect["max_score"] == pytest.approx(0.8)
