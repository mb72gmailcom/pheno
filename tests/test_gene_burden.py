import csv
import gzip
import json
from pathlib import Path

import pytest

from pheno.cli_gene_burden import main
from pheno.gene_burden import compute_gene_burden, load_analysis_groups
from pheno.mapping import load_column_map

FAMILY = """\
person	family	mother	father	sex	asd
A	F1	mom1	dad1	male	2
B	F1	mom1	dad1	female	1
C	F2	mom2	dad2	male	2
D	F2	mom2	dad2	female	1
E	F3	mom3	dad3	male	1
F	F4	mom4	dad4	female	2
P	F1	0	0	male	1
"""

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"
OTARI_HEADER = "variant_id\ttranscript_id\tmax_effect\n"
GENE_HEADER = "variant_id\tgene_id\ttranscript_id\tmost_affected_node\ttop_features\n"
ANNOTATION = """\
feature	chr	start	end	strand	name
gene	chr21	1000	2000	+	GENE1
gene	chr21	100000	101000	+	GENE2
gene	chr21	3000	5000	+	GENE3
"""


def _load_genes(path: Path) -> dict[str, dict[str, dict[str, object]]]:
    genes: dict[str, dict[str, dict[str, object]]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        score_columns = [name[: -len("_n_scored")] for name in fieldnames if name.endswith("_n_scored")]
        for row in reader:
            record: dict[str, object] = {
                "status": row["status"],
                "family_id": row["family_id"],
                "n_variants": int(row["n_variants"]),
                "n_unscored": int(row["n_unscored"]),
            }
            for column in score_columns:
                fraction = row[f"{column}_fraction_damaging"]
                var_mean = row[f"{column}_var_mean_effect"]
                max_score = row[f"{column}_max_score"]
                record[column] = {
                    "n_scored": int(row[f"{column}_n_scored"]),
                    "n_damaging": int(row[f"{column}_n_damaging"]),
                    "fraction_damaging": None if fraction == "" else float(fraction),
                    "sum_effect": float(row[f"{column}_sum_effect"]),
                    "var_mean_effect": None if var_mean == "" else float(var_mean),
                    "sum_damaging": float(row[f"{column}_sum_damaging"]),
                    "max_score": None if max_score == "" else float(max_score),
                }
            genes.setdefault(row["gene"], {})[row["person"]] = record
    return {"genes": genes}


def _write(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")


def _cohort(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    source = tmp_path / "results"
    otari = tmp_path / "otari"
    _write(
        source / "chr21" / "inherited_1000_5000.tsv",
        HEADER
        + "chr21\t1500\t.\tA\tT\tA\n"
        + "chr21\t1800\t.\tA\tT\tA\n"
        + "chr21\t1900\t.\tA\tT\tB\n"
        + "chr21\t2000\t.\tA\tT\tD\n"
        + "chr21\t4001\t.\tA\tT\tA\n"
        + "chr21\t3500\t.\tC\tG\tA\n"
        + "chr21\t100050\t.\tC\tG\tC;P\n",
    )
    shard = otari / "inherited" / "chr21" / "1000_5000"
    _write(
        shard / "variant_effects_comprehensive.tsv",
        OTARI_HEADER
        + "21_1500_A_T_hg38\tENST1\t0.2\n"
        + "21_1500_A_T_hg38\tENST2\t0.8\n"
        + "21_1800_A_T_hg38\tENST1\t0.2\n"
        + "21_1900_A_T_hg38\tENST1\t0.1\n"
        + "21_4001_A_T_hg38\tENST1\t5.0\n"
        + "21_3500_C_G_hg38\tENST1\t0.4\n"
        + "21_3500_C_G_hg38\tENST3\t0.7\n"
        + "21_100050_C_G_hg38\tENST9\t0.9\n",
    )
    _write(
        shard / "interpretability_analysis.tsv",
        GENE_HEADER
        + "21_1500_A_T_hg38\tGENE1\tENST1\t0\t\n"
        + "21_1500_A_T_hg38\tGENE1\tENST2\t0\t\n"
        + "21_1800_A_T_hg38\tGENE1\tENST1\t0\t\n"
        + "21_1900_A_T_hg38\tGENE1\tENST1\t0\t\n"
        + "21_4001_A_T_hg38\tGENE1\tENST1\t0\t\n"
        + "21_3500_C_G_hg38\tGENE1\tENST1\t0\t\n"
        + "21_3500_C_G_hg38\tGENE3\tENST3\t0\t\n"
        + "21_100050_C_G_hg38\tGENE2\tENST9\t0\t\n",
    )
    annotation = tmp_path / "genes.tsv"
    annotation.write_text(ANNOTATION, encoding="utf-8")
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    return source, otari, annotation, family


def test_gene_burden_cohort_and_per_gene_scores(tmp_path: Path):
    source, otari, annotation, family = _cohort(tmp_path)
    people, n_asd, n_siblings = load_analysis_groups(family, load_column_map(None))
    assert n_asd == 3
    assert n_siblings == 2
    assert set(people) == {"A", "B", "C", "D", "F"}
    output = tmp_path / "out"
    summary = compute_gene_burden(
        source,
        otari,
        annotation,
        "inherited",
        people,
        ["max_effect"],
        threshold=0.5,
        use_abs=True,
        transcripts="max",
        n_asd=n_asd,
        n_unaffected_sibling=n_siblings,
        output_dir=output,
    )
    genes = _load_genes(output / "chr21" / "genes.tsv.gz")
    assert not (output / "chr21" / ".partial").exists()

    gene1 = genes["genes"]["GENE1"]
    assert set(gene1) == {"A", "B", "D"}
    person_a = gene1["A"]
    assert person_a["n_variants"] == 3
    assert person_a["n_unscored"] == 0
    effect = person_a["max_effect"]
    assert effect["n_scored"] == 3
    assert effect["n_damaging"] == 1
    assert effect["fraction_damaging"] == pytest.approx(1 / 3)
    assert effect["sum_effect"] == pytest.approx(1.4)
    assert effect["var_mean_effect"] == pytest.approx(1.4 / 3)
    assert effect["sum_damaging"] == pytest.approx(0.8)
    assert effect["max_score"] == pytest.approx(0.8)
    assert gene1["B"]["max_effect"]["sum_effect"] == pytest.approx(0.1)
    assert gene1["B"]["max_effect"]["var_mean_effect"] == pytest.approx(0.1)
    assert gene1["D"]["n_variants"] == 1
    assert gene1["D"]["n_unscored"] == 1
    assert gene1["D"]["max_effect"]["n_scored"] == 0
    assert gene1["D"]["max_effect"]["var_mean_effect"] is None
    assert "P" not in genes["genes"]["GENE2"]
    assert genes["genes"]["GENE2"]["C"]["max_effect"]["sum_effect"] == pytest.approx(0.9)
    assert set(genes["genes"]["GENE3"]) == {"A", "B", "D"}
    assert genes["genes"]["GENE3"]["A"]["max_effect"]["sum_effect"] == pytest.approx(0.7)
    assert genes["genes"]["GENE3"]["A"]["max_effect"]["n_scored"] == 1
    assert genes["genes"]["GENE3"]["B"]["n_unscored"] == 1
    assert genes["genes"]["GENE3"]["B"]["max_effect"]["n_scored"] == 0

    assert summary["n_asd"] == 3
    assert summary["n_unaffected_sibling"] == 2
    asd = summary["genes"]["GENE1"]["asd"]
    assert asd["n_people"] == 3
    assert asd["n_carriers"] == 1
    assert asd["max_effect"]["mean_sum_effect"] == pytest.approx(1.4)
    assert asd["max_effect"]["mean_var_mean_effect"] == pytest.approx(1.4 / 3)
    assert asd["max_effect"]["n_scored_carriers"] == 1
    assert asd["max_effect"]["n_scored"] == 3
    assert asd["max_effect"]["n_damaging"] == 1
    assert asd["max_effect"]["fraction_damaging"] == pytest.approx(1 / 3)
    siblings = summary["genes"]["GENE1"]["unaffected_sibling"]
    assert siblings["n_people"] == 2
    assert siblings["n_carriers"] == 2
    assert siblings["max_effect"]["mean_sum_effect"] == pytest.approx(0.05)
    assert siblings["max_effect"]["mean_var_mean_effect"] == pytest.approx(0.1)
    assert siblings["max_effect"]["n_scored_carriers"] == 1
    assert siblings["max_effect"]["fraction_damaging"] == pytest.approx(0.0)
    gene2_asd = summary["genes"]["GENE2"]["asd"]
    assert gene2_asd["n_carriers"] == 1
    assert gene2_asd["max_effect"]["mean_sum_effect"] == pytest.approx(0.9)
    gene2_siblings = summary["genes"]["GENE2"]["unaffected_sibling"]
    assert gene2_siblings["n_carriers"] == 0
    assert gene2_siblings["n_people"] == 2
    assert gene2_siblings["max_effect"]["mean_sum_effect"] is None
    assert gene2_siblings["max_effect"]["mean_var_mean_effect"] is None
    assert gene2_siblings["max_effect"]["n_scored_carriers"] == 0
    assert gene2_siblings["max_effect"]["fraction_damaging"] is None


def test_gene_burden_cli_writes_summary(tmp_path: Path):
    source, otari, annotation, family = _cohort(tmp_path)
    output_dir = tmp_path / "out"
    rc = main(
        [
            "--family-file",
            str(family),
            "--input-dir",
            str(source),
            "--otari-dir",
            str(otari),
            "--annotation",
            str(annotation),
            "--output-dir",
            str(output_dir),
            "--file-pattern",
            "inherited",
            "--threshold",
            "0.5",
        ]
    )
    assert rc == 0
    summary = json.loads((output_dir / "gene_summary.json").read_text(encoding="utf-8"))
    chromosome = _load_genes(output_dir / "chr21" / "genes.tsv.gz")
    chromosome_summary = json.loads((output_dir / "chr21" / "gene_summary.json").read_text(encoding="utf-8"))
    assert not (output_dir / "genes.tsv.gz").exists()
    assert chromosome["genes"]["GENE2"]["C"]["status"] == "asd"
    assert chromosome["genes"]["GENE1"]["A"]["n_variants"] == 3
    assert chromosome_summary["genes"]["GENE1"]["asd"]["n_people"] == 3
    assert summary["n_asd"] == 3
    assert summary["genes"]["GENE1"]["asd"]["n_people"] == 3
    assert not (output_dir / "chr21" / ".partial").exists()


def test_gene_burden_merges_person_across_shards(tmp_path: Path):
    source, otari, annotation, family = _cohort(tmp_path)
    single = source / "chr21" / "inherited_1000_5000.tsv"
    body = single.read_text(encoding="utf-8")
    lines = [line for line in body.splitlines() if line and not line.startswith("#")]
    early = [line for line in lines if int(line.split("\t")[1]) <= 1800]
    late = [line for line in lines if int(line.split("\t")[1]) > 1800]
    single.unlink()
    _write(source / "chr21" / "inherited_1000_1800.tsv", HEADER + "".join(f"{line}\n" for line in early))
    _write(source / "chr21" / "inherited_1801_200000.tsv", HEADER + "".join(f"{line}\n" for line in late))
    shard = otari / "inherited" / "chr21" / "1000_5000"
    effects = (shard / "variant_effects_comprehensive.tsv").read_text(encoding="utf-8")
    genes_map = (shard / "interpretability_analysis.tsv").read_text(encoding="utf-8")
    early_ids = {"21_1500_A_T_hg38", "21_1800_A_T_hg38"}

    def _split(text: str, keep_early: bool) -> str:
        rows = text.splitlines()
        kept = [rows[0]]
        for row in rows[1:]:
            variant_id = row.split("\t", 1)[0]
            if (variant_id in early_ids) == keep_early:
                kept.append(row)
        return "\n".join(kept) + "\n"

    _write(otari / "inherited" / "chr21" / "1000_1800" / "variant_effects_comprehensive.tsv", _split(effects, True))
    _write(otari / "inherited" / "chr21" / "1000_1800" / "interpretability_analysis.tsv", _split(genes_map, True))
    _write(otari / "inherited" / "chr21" / "1801_200000" / "variant_effects_comprehensive.tsv", _split(effects, False))
    _write(otari / "inherited" / "chr21" / "1801_200000" / "interpretability_analysis.tsv", _split(genes_map, False))
    shutil_shard = shard
    for path in shutil_shard.iterdir():
        path.unlink()
    shutil_shard.rmdir()

    people, n_asd, n_siblings = load_analysis_groups(family, load_column_map(None))
    output = tmp_path / "out"
    compute_gene_burden(
        source,
        otari,
        annotation,
        "inherited",
        people,
        ["max_effect"],
        threshold=0.5,
        use_abs=True,
        transcripts="max",
        n_asd=n_asd,
        n_unaffected_sibling=n_siblings,
        output_dir=output,
    )
    genes = _load_genes(output / "chr21" / "genes.tsv.gz")
    person_a = genes["genes"]["GENE1"]["A"]
    assert person_a["n_variants"] == 3
    assert person_a["max_effect"]["n_scored"] == 3
    assert person_a["max_effect"]["sum_effect"] == pytest.approx(1.4)
    assert person_a["max_effect"]["var_mean_effect"] == pytest.approx(1.4 / 3)
    assert person_a["max_effect"]["max_score"] == pytest.approx(0.8)
    summary = json.loads((output / "chr21" / "gene_summary.json").read_text(encoding="utf-8"))
    assert summary["genes"]["GENE1"]["asd"]["n_carriers"] == 1
    assert summary["genes"]["GENE1"]["asd"]["max_effect"]["mean_var_mean_effect"] == pytest.approx(1.4 / 3)
    assert not (output / "chr21" / ".partial").exists()
