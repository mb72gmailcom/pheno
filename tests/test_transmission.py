import csv
import math
from pathlib import Path

import pytest
from scipy.stats import norm

from pheno.cli_transmission import main
from pheno.gene_burden import load_analysis_groups
from pheno.mapping import load_column_map
from pheno.transmission import compute_transmission
from pheno.transmission_scores import gene_shard_covers

FAMILY = """\
person	family	mother	father	sex	asd
A	F1	mom1	dad1	male	2
B	F1	mom1	dad1	female	1
C	F2	mom2	dad2	male	2
D	F2	mom2	dad2	female	1
E	F3	mom3	dad3	male	1
U	F1	mom1	dad1	male	1
P	F5	mom5	dad5	male	2
Q	F5	mom5	dad5	female	1
S	F6	mom6	dad6	male	2
"""
HEADER = "#CHROM\tPOS\tID\tREF\tALT\tTRIO_CALLS\n"


def _write(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HEADER + body, encoding="utf-8")


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _run(tmp_path: Path, body: str, fcap: int | None = None) -> list[dict[str, str]]:
    source = tmp_path / "results"
    _write(source / "chr21" / "inherited_00000.tsv", body)
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    output = tmp_path / "transmission.tsv"
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    compute_transmission(source, "inherited", people, output, fcap=fcap)
    return _rows(output)


def test_pair_and_asymmetric_family(tmp_path: Path):
    rows = _run(
        tmp_path,
        "chr21\t10\t.\tG\tA\tP=0/0|0/1|0/1|99;Q=0/1|0/0|0/0|99;"
        "C=0/0|0/1|0/1|99;D=0/1|0/0|0/0|99\n"
        "chr21\t20\t.\tT\tC\tA=0/1|0/1|0/0|99;B=0/0|0/1|0/1|99;U=0/0|0/1|0/1|99\n"
        "chr21\t30\t.\tC\tG\tA=0/0|0/1|0/1|99;B=0/1|0/0|0/1|99;U=0/0|0/1|0/1|99\n",
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


def test_absent_sibling_is_reference(tmp_path: Path):
    rows = _run(
        tmp_path,
        "chr21\t10\t.\tG\tA\tC=0/0|0/1|0/1|99\n"
        "chr21\t20\t.\tG\tA\tD=0/1|0/0|0/1|99\n",
    )
    by_pos = {row["POS"]: row for row in rows}
    assert by_pos["10"]["n_families"] == "1"
    assert by_pos["10"]["n_informative"] == "1"
    assert float(by_pos["10"]["z"]) == pytest.approx(1)
    assert float(by_pos["10"]["p_exact"]) == pytest.approx(1)
    assert float(by_pos["20"]["z"]) == pytest.approx(-1)


def test_multiallelic_variant_is_skipped(tmp_path: Path):
    rows = _run(
        tmp_path,
        "chr21\t10\t.\tG\tA,C\tA=0/0|0/1|0/1|99;B=0/1|0/0|0/0|99\n"
        "chr21\t20\t.\tG\tA\tA=0/0|0/1|1/2|99;B=0/1|0/0|0/0|99\n"
        "chr21\t30\t.\tG\tA\tC=0/0|1/1|0/1|99;D=1/1|0/0|0/0|99\n",
    )
    assert [row["POS"] for row in rows] == ["30"]
    assert float(rows[0]["z"]) == pytest.approx(1)
    assert float(rows[0]["p_exact"]) == pytest.approx(1)


def test_three_field_genotype_uses_child(tmp_path: Path):
    rows = _run(
        tmp_path,
        "chr21\t10\t.\tG\tA\tC=0/1|0/1|99;D=0/0|0/0|99\n",
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
    _write(source / "chrX" / "inherited_1_2.tsv", "chrX\t5\t.\tA\tG\tS=0/0|0/1|0/1|99\n")
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    output_dir = tmp_path / "by-chrom"
    written, gene_count = compute_transmission(source, "inherited", people, output_dir=output_dir)
    assert written == 3
    assert gene_count == 0
    assert [row["POS"] for row in _rows(output_dir / "chr21" / "inherited.tsv")] == ["30", "250"]
    assert [row["POS"] for row in _rows(output_dir / "chr22" / "inherited.tsv")] == ["15"]
    assert not (output_dir / "chrX").exists()


def test_gene_sums_three_variants_and_shuffles_once(tmp_path: Path):
    source = tmp_path / "results"
    _write(
        source / "chr21" / "inherited_1000_2000.tsv",
        "chr21\t900\t.\tG\tA\tA=0/0|0/1|1/1|99\n"
        "chr21\t1100\t.\tG\tA\tA=0/0|0/1|0/1|99;C=0/0|0/1|0/1|99\n"
        "chr21\t1500\t.\tT\tC\tA=0/1|0/1|1/1|99;B=0/0|0/1|0/1|99;"
        "D=0/0|0/1|0/1|99;E=0/0|0/1|0/1|99\n"
        "chr21\t1800\t.\tC\tG\tB=0/0|0/1|0/1|99;C=0/0|0/1|0/1|99\n",
    )
    family = tmp_path / "family.tsv"
    family.write_text(
        "person\tfamily\tmother\tfather\tsex\tasd\n"
        "A\tF1\tmom1\tdad1\tmale\t2\n"
        "B\tF1\tmom1\tdad1\tfemale\t1\n"
        "C\tF2\tmom2\tdad2\tmale\t2\n"
        "D\tF2\tmom2\tdad2\tfemale\t1\n"
        "E\tF2\tmom2\tdad2\tmale\t1\n",
        encoding="utf-8",
    )
    annotation = tmp_path / "genes.tsv"
    annotation.write_text(
        "feature\tchr\tstart\tend\tstrand\tname\tattributes\n"
        'gene\tchr21\t1000\t2000\t+\tGENE1\tgene_id "ENSG00000000001.2"; '
        'gene_type "protein_coding"; gene_name "AAA";\n',
        encoding="utf-8",
    )
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    output_dir = tmp_path / "by-chrom"
    _variants, gene_count = compute_transmission(
        source,
        "inherited",
        people,
        output_dir=output_dir,
        genes=annotation,
        window=0,
    )
    assert gene_count == 1
    rows = _rows(output_dir / "chr21" / "inherited.genes.tsv")
    assert len(rows) == 1
    gene = rows[0]
    assert gene["gene_id"] == "ENSG00000000001"
    assert gene["gene_name"] == "AAA"
    assert gene["gene_type"] == "protein_coding"
    assert gene["n_variants"] == "3"
    assert gene["n_families"] == "2"
    assert gene["n_informative"] == "2"
    assert float(gene["z"]) == pytest.approx(2 / math.sqrt(1.5))
    assert float(gene["p_exact"]) == pytest.approx(1 / 6)


def test_fcap_uses_the_normal_tail_above_the_cutoff(tmp_path: Path):
    body = (
        "chr21\t10\t.\tG\tA\tP=0/0|0/1|0/1|99;Q=0/1|0/0|0/0|99;"
        "C=0/0|0/1|0/1|99;D=0/1|0/0|0/0|99\n"
        "chr21\t30\t.\tC\tG\tA=0/0|0/1|0/1|99;B=0/1|0/0|0/1|99;U=0/0|0/1|0/1|99\n"
    )
    capped = {row["POS"]: row for row in _run(tmp_path, body, fcap=1)}
    z = math.sqrt(2)
    assert capped["10"]["n_informative"] == "2"
    assert float(capped["10"]["z"]) == pytest.approx(z)
    assert float(capped["10"]["p_exact"]) == pytest.approx(2 * norm.sf(z))
    assert float(capped["30"]["p_exact"]) == pytest.approx(1)

    exact = {row["POS"]: row for row in _run(tmp_path, body, fcap=2)}
    assert float(exact["10"]["p_exact"]) == pytest.approx(0.5)


def test_fcap_applies_to_gene_rows(tmp_path: Path):
    source = tmp_path / "results"
    _write(
        source / "chr21" / "inherited_1000_2000.tsv",
        "chr21\t1100\t.\tG\tA\tA=0/0|0/1|0/1|99;C=0/0|0/1|0/1|99\n"
        "chr21\t1500\t.\tT\tC\tA=0/1|0/1|1/1|99;B=0/0|0/1|0/1|99;"
        "D=0/0|0/1|0/1|99;E=0/0|0/1|0/1|99\n"
        "chr21\t1800\t.\tC\tG\tB=0/0|0/1|0/1|99;C=0/0|0/1|0/1|99\n",
    )
    family = tmp_path / "family.tsv"
    family.write_text(
        "person\tfamily\tmother\tfather\tsex\tasd\n"
        "A\tF1\tmom1\tdad1\tmale\t2\n"
        "B\tF1\tmom1\tdad1\tfemale\t1\n"
        "C\tF2\tmom2\tdad2\tmale\t2\n"
        "D\tF2\tmom2\tdad2\tfemale\t1\n"
        "E\tF2\tmom2\tdad2\tmale\t1\n",
        encoding="utf-8",
    )
    annotation = tmp_path / "genes.tsv"
    annotation.write_text(
        "feature\tchr\tstart\tend\tstrand\tname\tattributes\n"
        'gene\tchr21\t1000\t2000\t+\tGENE1\tgene_id "ENSG00000000001.2"; '
        'gene_type "protein_coding"; gene_name "AAA";\n',
        encoding="utf-8",
    )
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    output_dir = tmp_path / "by-chrom"
    compute_transmission(
        source,
        "inherited",
        people,
        output_dir=output_dir,
        genes=annotation,
        window=0,
        fcap=1,
    )
    gene = _rows(output_dir / "chr21" / "inherited.genes.tsv")[0]
    z = 2 / math.sqrt(1.5)
    assert gene["n_informative"] == "2"
    assert float(gene["z"]) == pytest.approx(z)
    assert float(gene["p_exact"]) == pytest.approx(2 * norm.sf(z))


def test_negative_fcap_is_rejected(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    code = main(
        [
            "--family-file",
            str(family),
            "--input-dir",
            str(tmp_path),
            "--file-pattern",
            "inherited",
            "--output",
            str(tmp_path / "out.tsv"),
            "--fcap",
            "-1",
        ]
    )
    assert code == 1
    assert "error: --fcap must be at least 0" in capsys.readouterr().err


def _annotation(path: Path) -> None:
    path.write_text(
        "feature\tchr\tstart\tend\tstrand\tname\tattributes\n"
        'gene\tchr21\t1\t5000\t+\tGENE1\tgene_id "ENSG00000000001.2"; '
        'gene_type "protein_coding"; gene_name "AAA";\n',
        encoding="utf-8",
    )


def _pair_family(path: Path) -> None:
    path.write_text(
        "person\tfamily\tmother\tfather\tsex\tasd\n"
        "A\tF1\tmom\tdad\tmale\t2\n"
        "B\tF1\tmom\tdad\tfemale\t1\n",
        encoding="utf-8",
    )


def test_otari_weights_keep_the_strongest_signed_transcript(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    source = tmp_path / "results"
    _write(
        source / "chr21" / "inherited_1_5000.tsv",
        "chr21\t100\t.\tG\tA\tA=0/0|0/1|0/1|99\n"
        "chr21\t200\t.\tC\tT\tB=0/0|0/1|0/1|99\n"
        "chr21\t300\t.\tA\tG\tA=0/0|0/1|0/1|99\n",
    )
    family = tmp_path / "family.tsv"
    _pair_family(family)
    annotation = tmp_path / "genes.tsv"
    _annotation(annotation)
    shard = tmp_path / "otari" / "chr21" / "1_5000"
    shard.mkdir(parents=True)
    (shard / "variant_effects_comprehensive.tsv").write_text(
        "variant_id\ttranscript_id\tBrain\tLiver\tHeart\tCancel\n"
        "21_100_G_A_hg38\tt1\t-0.75\t0.10\t-1.5\t2\n"
        "21_100_G_A_hg38\tt2\t0.25\t0.30\t1.5\t0\n"
        "21_100_G_A_hg38\tt3\t9\t9\t9\t9\n"
        "21_200_C_T_hg38\tt1\t0\t0\t0\t2\n",
        encoding="utf-8",
    )
    (shard / "interpretability_analysis.tsv").write_text(
        "variant_id\tgene_id\ttranscript_id\n"
        "21_100_G_A_hg38\tENSG00000000001.5\tt1\n"
        "21_100_G_A_hg38\tENSG00000000001.5\tt2\n"
        "21_100_G_A_hg38\tENSG00000000002.1\tt3\n"
        "21_200_C_T_hg38\tENSG00000000001.5\tt1\n",
        encoding="utf-8",
    )
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    output_dir = tmp_path / "by-chrom"
    compute_transmission(
        source,
        "inherited",
        people,
        output_dir=output_dir,
        genes=annotation,
        window=0,
        otari_dir=tmp_path / "otari",
        debug_gene="ENSG00000000001.5",
        chrom="21",
    )
    assert not (output_dir / "chr21" / "inherited.genes.tsv").exists()
    rows = _rows(output_dir / "chr21" / "inherited.gene_scores.tsv")
    assert len(rows) == 1
    gene = rows[0]
    assert gene["gene_id"] == "ENSG00000000001"
    assert gene["n_variants"] == "3"
    assert gene["n_informative"] == "1"
    assert float(gene["z"]) == pytest.approx(1)
    assert float(gene["p_exact"]) == pytest.approx(1)
    assert gene["n_scored"] == "2"
    assert gene["n_families"] == "1"
    normal_p = 2 * norm.sf(1)
    for column, expected in (("Brain", -1), ("Liver", 1), ("Heart", 1)):
        assert gene[f"{column}_n_informative"] == "1"
        assert float(gene[f"{column}_z"]) == pytest.approx(expected)
        assert float(gene[f"{column}_p"]) == pytest.approx(normal_p)
    assert gene["Cancel_n_informative"] == "0"
    assert gene["Cancel_z"] == ""
    assert float(gene["Cancel_p"]) == pytest.approx(1)
    variants = {row["POS"]: row for row in _rows(output_dir / "chr21" / "inherited.tsv")}
    assert float(variants["100"]["p_exact"]) == pytest.approx(1)
    printed = capsys.readouterr().out
    assert "debug ENSG00000000001: variant chr21 100 G A shard 1_5000; dosages A=1; score found" in printed
    assert "debug ENSG00000000001: variant chr21 300 A G shard 1_5000; dosages A=1; no score" in printed
    assert "writing n_variants=3 n_families=1 n_scored=2 n_informative=1 z=" in printed


def test_debug_gene_reports_a_missing_gene_and_a_dropped_row(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    source = tmp_path / "results"
    _write(source / "chr21" / "inherited_1_20.tsv", "chr21\t10\t.\tG\tA\tA=0/0|0/1|0/1|99\n")
    family = tmp_path / "family.tsv"
    family.write_text(
        "person\tfamily\tmother\tfather\tsex\tasd\nA\tF1\tmom\tdad\tmale\t2\n",
        encoding="utf-8",
    )
    annotation = tmp_path / "genes.tsv"
    _annotation(annotation)
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    output_dir = tmp_path / "by-chrom"
    compute_transmission(
        source,
        "inherited",
        people,
        output_dir=output_dir,
        genes=annotation,
        window=0,
        debug_gene="ENSG00000000999",
        chrom="chr21",
    )
    assert "debug ENSG00000000999: chr21: not in the annotation" in capsys.readouterr().out
    compute_transmission(
        source,
        "inherited",
        people,
        output_dir=output_dir,
        genes=annotation,
        window=0,
        debug_gene="ENSG00000000001",
        chrom="21",
    )
    printed = capsys.readouterr().out
    assert "dropping; no carrier family has both an ASD child and an unaffected sibling" in printed
    assert "A asd family F1 asd=1 unaffected=0" in printed
    assert not (output_dir / "chr21" / "inherited.genes.tsv").exists()


def test_debug_gene_requires_a_chromosome(tmp_path: Path):
    source = tmp_path / "results"
    _write(source / "chr21" / "inherited_1_20.tsv", "chr21\t10\t.\tG\tA\tA=0/0|0/1|0/1|99\n")
    _write(source / "chr22" / "inherited_1_20.tsv", "chr22\t10\t.\tG\tA\tA=0/0|0/1|0/1|99\n")
    family = tmp_path / "family.tsv"
    _pair_family(family)
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    with pytest.raises(ValueError, match="--debug-gene requires --chr"):
        compute_transmission(
            source,
            "inherited",
            people,
            output_dir=tmp_path / "by-chrom",
            debug_gene="ENSG00000000001",
        )
    output_dir = tmp_path / "one-chrom"
    compute_transmission(source, "inherited", people, output_dir=output_dir, chrom="21")
    assert (output_dir / "chr21" / "inherited.tsv").is_file()
    assert not (output_dir / "chr22").exists()
    code = main(
        [
            "--family-file",
            str(family),
            "--input-dir",
            str(source),
            "--file-pattern",
            "inherited",
            "--output-dir",
            str(tmp_path / "cli"),
            "--debug-gene",
            "ENSG00000000001",
        ]
    )
    assert code == 1


def _score_shard(root: Path, chrom: str, shard: str, rows: str) -> None:
    directory = root / chrom / shard
    directory.mkdir(parents=True)
    (directory / "variant_effects_comprehensive.tsv").write_text(
        "variant_id\ttranscript_id\tBrain\n" + rows,
        encoding="utf-8",
    )
    interpreted = "".join(
        f"{line.split(chr(9))[0]}\tENSG00000000002\t{line.split(chr(9))[1]}\n"
        if "21_150_" not in line
        else f"{line.split(chr(9))[0]}\tENSG00000000001\t{line.split(chr(9))[1]}\n"
        for line in rows.splitlines()
        if line
    )
    (directory / "interpretability_analysis.tsv").write_text(
        "variant_id\tgene_id\ttranscript_id\n" + interpreted,
        encoding="utf-8",
    )


def test_otari_keeps_a_gene_open_until_the_next_shard(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    source = tmp_path / "results"
    _write(source / "chr21" / "inherited_1_1000.tsv", "chr21\t150\t.\tG\tA\tA=0/0|0/1|0/1|99\nchr21\t900\t.\tC\tT\tA=0/0|0/1|0/1|99\n")
    _write(source / "chr21" / "inherited_1001_2000.tsv", "chr21\t1200\t.\tA\tG\tA=0/0|0/1|0/1|99\n")
    family = tmp_path / "family.tsv"
    _pair_family(family)
    annotation = tmp_path / "genes.tsv"
    annotation.write_text(
        "feature\tchr\tstart\tend\tstrand\tname\tattributes\n"
        'gene\tchr21\t100\t200\t+\tEARLY\tgene_id "ENSG00000000001.1"; '
        'gene_type "protein_coding"; gene_name "EARLY";\n'
        'gene\tchr21\t500\t1500\t+\tSPAN\tgene_id "ENSG00000000002.1"; '
        'gene_type "protein_coding"; gene_name "SPAN";\n',
        encoding="utf-8",
    )
    otari = tmp_path / "otari"
    _score_shard(otari, "chr21", "1_1000", "21_150_G_A_hg38\tt1\t1\n21_900_C_T_hg38\tt1\t1\n")
    _score_shard(otari, "chr21", "1001_2000", "21_1200_A_G_hg38\tt1\t1\n")
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    output_dir = tmp_path / "by-chrom"
    compute_transmission(
        source,
        "inherited",
        people,
        output_dir=output_dir,
        genes=annotation,
        window=0,
        otari_dir=otari,
    )
    rows = _rows(output_dir / "chr21" / "inherited.gene_scores.tsv")
    assert [row["gene_name"] for row in rows] == ["EARLY", "SPAN"]
    assert rows[1]["n_variants"] == "2"
    assert rows[1]["n_scored"] == "2"
    assert float(rows[1]["Brain_z"]) == pytest.approx(1)
    printed = capsys.readouterr().out
    assert "chr21 shard 1_1000: variants 2, assigned 2, gene rows 1, still open 1" in printed
    assert "chr21 shard 1001_2000: variants 1, assigned 1, gene rows 1, still open 0" in printed


def test_gene_shard_cover_counts_the_outer_span():
    shards = [(1, 1000), (1001, 2000), (2001, 3000)]
    short = gene_shard_covers([(100, 200, "EARLY"), (500, 1500, "SPAN")], shards)
    by_gene = {gene_id: held for gene_id, _lo, _hi, held in short}
    assert [item[0] for item in by_gene["EARLY"]] == [1]
    assert [item[0] for item in by_gene["SPAN"]] == [1, 1001]
    long = gene_shard_covers([(100, 2600, "LONG")], shards)
    assert [item[0] for item in long[0][3]] == [1, 1001, 2001]


def test_a_gene_spanning_three_shards_is_printed(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    source = tmp_path / "results"
    _write(source / "chr21" / "inherited_1_1000.tsv", "chr21\t500\t.\tC\tT\tA=0/0|0/1|0/1|99\n")
    _write(source / "chr21" / "inherited_1001_2000.tsv", "chr21\t1500\t.\tA\tG\tA=0/0|0/1|0/1|99\n")
    _write(source / "chr21" / "inherited_2001_3000.tsv", "chr21\t2500\t.\tG\tC\tA=0/0|0/1|0/1|99\n")
    family = tmp_path / "family.tsv"
    _pair_family(family)
    annotation = tmp_path / "genes.tsv"
    annotation.write_text(
        "feature\tchr\tstart\tend\tstrand\tname\tattributes\n"
        'gene\tchr21\t100\t200\t+\tEARLY\tgene_id "ENSG00000000001.1"; '
        'gene_type "protein_coding"; gene_name "EARLY";\n'
        'gene\tchr21\t100\t2600\t+\tLONG\tgene_id "ENSG00000000002.1"; '
        'gene_type "protein_coding"; gene_name "LONG";\n',
        encoding="utf-8",
    )
    otari = tmp_path / "otari"
    _score_shard(otari, "chr21", "1_1000", "21_500_C_T_hg38\tt1\t1\n")
    _score_shard(otari, "chr21", "1001_2000", "21_1500_A_G_hg38\tt1\t1\n")
    _score_shard(otari, "chr21", "2001_3000", "21_2500_G_C_hg38\tt1\t1\n")
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    output_dir = tmp_path / "by-chrom"
    compute_transmission(
        source,
        "inherited",
        people,
        output_dir=output_dir,
        genes=annotation,
        window=0,
        otari_dir=otari,
    )
    printed = capsys.readouterr().out
    assert "ENSG00000000002 LONG (100-2600) spans 3 shards and stays open across 1_1000, 1001_2000, 2001_3000" in printed
    assert "chr21: 1 gene spans 3 or more shards" in printed
    assert "chr21 shard 1_1000: variants 1, assigned 1, gene rows 0, still open 1" in printed
    assert "chr21 shard 2001_3000: variants 1, assigned 1, gene rows 1, still open 0" in printed
    assert "EARLY" not in printed
    rows = _rows(output_dir / "chr21" / "inherited.gene_scores.tsv")
    assert [row["gene_name"] for row in rows] == ["LONG"]
    assert rows[0]["n_variants"] == "3"


def test_otari_shard_splits_must_match(tmp_path: Path):
    source = tmp_path / "results" / "chr21"
    source.mkdir(parents=True)
    (source / "inherited_1_1000.tsv").write_text(HEADER, encoding="utf-8")
    (source / "inherited_1001_2000.tsv").write_text(HEADER, encoding="utf-8")
    otari = tmp_path / "otari"
    _score_shard(otari, "chr21", "1_1000", "")
    family = tmp_path / "family.tsv"
    _pair_family(family)
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    output_dir = tmp_path / "by-chrom"
    with pytest.raises(ValueError, match="chr21 shard splits differ"):
        compute_transmission(
            source.parent,
            "inherited",
            people,
            output_dir=output_dir,
            genes=tmp_path / "genes.tsv",
            otari_dir=otari,
        )
    assert not output_dir.exists()


def test_seqweaver_uses_label_coordinates_and_the_first_duplicate(tmp_path: Path):
    source = tmp_path / "results"
    _write(source / "chr21" / "inherited_1_20.tsv", "chr21\t10\t.\tG\tA\tA=0/0|0/1|0/1|99\n")
    family = tmp_path / "family.tsv"
    _pair_family(family)
    annotation = tmp_path / "genes.tsv"
    _annotation(annotation)
    root = tmp_path / "seqweaver" / "chr21"

    def write_h5(directory: str, name: str, matrix: object, labels: str) -> None:
        import h5py

        path = root / directory / "seqweaver_human" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(path, "w") as handle:
            handle.create_dataset("data", data=matrix)
        path.with_name(name.replace("_diffs.h5", "_row_labels.txt")).write_text(labels, encoding="utf-8")

    label = (
        "chrom pos name ref alt strand ref_match contains_unk\n"
        "chr21 10 . G A + True False\n"
    )
    write_h5("variants_1_20/a", "annot.variants_00000_diffs.h5", [[1.0]], label)
    write_h5("variants_1_20/b", "annot.variants_00001_diffs.h5", [[-1.0]], label)
    import numpy as np

    write_h5(
        "variants_1_20/c",
        "annot.variants_00002_diffs.h5",
        np.zeros((0, 1)),
        "chrom\tpos\tname\tref\talt\tstrand\tref_match\tcontains_unk\n",
    )
    names = tmp_path / "seqweaver.colnames"
    names.write_text("ago_adult_brain.BA4.hg19\n", encoding="utf-8")
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    output = tmp_path / "transmission.tsv"
    compute_transmission(
        source,
        "inherited",
        people,
        output,
        genes=annotation,
        window=0,
        seqweaver_dir=tmp_path / "seqweaver",
        column_names=names,
    )
    assert not gene_table_path_exists(output)
    rows = _rows(output.with_name("transmission.gene_scores.tsv"))
    gene = rows[0]
    column = "ago_adult_brain.BA4.hg19"
    assert gene["n_variants"] == "1"
    assert gene["n_scored"] == "1"
    assert gene[f"{column}_n_informative"] == "1"
    assert float(gene[f"{column}_z"]) == pytest.approx(1)
    assert float(gene[f"{column}_p"]) == pytest.approx(2 * norm.sf(1))


def gene_table_path_exists(output: Path) -> bool:
    suffix = output.suffix or ".tsv"
    return output.with_name(f"{output.stem}.genes{suffix}").exists()


def test_seqweaver_header_only_file_must_match_the_matrix(tmp_path: Path):
    import h5py

    source = tmp_path / "results"
    _write(source / "chr21" / "inherited_1_20.tsv", "chr21\t10\t.\tG\tA\tA=0/0|0/1|0/1|99\n")
    family = tmp_path / "family.tsv"
    _pair_family(family)
    annotation = tmp_path / "genes.tsv"
    _annotation(annotation)
    diffs = (
        tmp_path
        / "seqweaver"
        / "chr21"
        / "variants_1_20"
        / "seqweaver_human"
        / "annot.variants_00000_diffs.h5"
    )
    diffs.parent.mkdir(parents=True)
    with h5py.File(diffs, "w") as handle:
        handle.create_dataset("data", data=[[1.0]])
    diffs.with_name("annot.variants_00000_row_labels.txt").write_text(
        "chrom\tpos\tname\tref\talt\tstrand\tref_match\tcontains_unk\n",
        encoding="utf-8",
    )
    names = tmp_path / "seqweaver.colnames"
    names.write_text("track\n", encoding="utf-8")
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    with pytest.raises(ValueError, match="0 label rows for 1 variants"):
        compute_transmission(
            source,
            "inherited",
            people,
            tmp_path / "transmission.tsv",
            genes=annotation,
            window=0,
            seqweaver_dir=tmp_path / "seqweaver",
            column_names=names,
        )


def test_seqweaver_keeps_a_gene_open_until_the_next_shard(tmp_path: Path):
    import h5py

    source = tmp_path / "results"
    _write(source / "chr21" / "inherited_1_1000.tsv", "chr21\t900\t.\tC\tT\tA=0/0|0/1|0/1|99\n")
    _write(source / "chr21" / "inherited_1001_2000.tsv", "chr21\t1200\t.\tA\tG\tA=0/0|0/1|0/1|99\n")
    family = tmp_path / "family.tsv"
    _pair_family(family)
    annotation = tmp_path / "genes.tsv"
    annotation.write_text(
        "feature\tchr\tstart\tend\tstrand\tname\tattributes\n"
        'gene\tchr21\t500\t1500\t+\tSPAN\tgene_id "ENSG00000000002.1"; '
        'gene_type "protein_coding"; gene_name "SPAN";\n',
        encoding="utf-8",
    )

    def write_h5(shard: str, pos: str, ref: str, alt: str) -> None:
        path = (
            tmp_path
            / "seqweaver"
            / "chr21"
            / shard
            / "seqweaver_human"
            / f"annot.{shard}_diffs.h5"
        )
        path.parent.mkdir(parents=True)
        with h5py.File(path, "w") as handle:
            handle.create_dataset("data", data=[[1.0]])
        path.with_name(f"annot.{shard}_row_labels.txt").write_text(
            "chrom pos name ref alt strand ref_match contains_unk\n"
            f"chr21 {pos} . {ref} {alt} + True False\n",
            encoding="utf-8",
        )

    write_h5("variants_1_1000", "900", "C", "T")
    write_h5("variants_1001_2000", "1200", "A", "G")
    names = tmp_path / "seqweaver.colnames"
    names.write_text("track\n", encoding="utf-8")
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    output_dir = tmp_path / "by-chrom"
    compute_transmission(
        source,
        "inherited",
        people,
        output_dir=output_dir,
        genes=annotation,
        window=0,
        seqweaver_dir=tmp_path / "seqweaver",
        column_names=names,
    )
    rows = _rows(output_dir / "chr21" / "inherited.gene_scores.tsv")
    assert len(rows) == 1
    assert rows[0]["n_variants"] == "2"
    assert rows[0]["n_scored"] == "2"
    assert float(rows[0]["track_z"]) == pytest.approx(1)


def test_seqweaver_shard_splits_must_match(tmp_path: Path):
    source = tmp_path / "results" / "chr21"
    source.mkdir(parents=True)
    (source / "inherited_1_1000.tsv").write_text(HEADER, encoding="utf-8")
    (source / "inherited_1001_2000.tsv").write_text(HEADER, encoding="utf-8")
    matrix = (
        tmp_path
        / "seqweaver"
        / "chr21"
        / "variants_1_1000"
        / "seqweaver_human"
        / "annot.variants_1_1000_diffs.h5"
    )
    matrix.parent.mkdir(parents=True)
    matrix.write_bytes(b"")
    names = tmp_path / "seqweaver.colnames"
    names.write_text("track\n", encoding="utf-8")
    family = tmp_path / "family.tsv"
    _pair_family(family)
    people, _n_asd, _n_siblings = load_analysis_groups(family, load_column_map(None))
    output_dir = tmp_path / "by-chrom"
    with pytest.raises(ValueError, match="seqweaver has"):
        compute_transmission(
            source.parent,
            "inherited",
            people,
            output_dir=output_dir,
            genes=tmp_path / "genes.tsv",
            seqweaver_dir=tmp_path / "seqweaver",
            column_names=names,
        )
    assert not output_dir.exists()


def test_model_score_flags_are_rejected(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    shared = [
        "--family-file",
        str(family),
        "--input-dir",
        str(tmp_path),
        "--file-pattern",
        "inherited",
        "--output",
        str(tmp_path / "out.tsv"),
    ]
    code = main([*shared, "--otari-dir", str(tmp_path), "--seqweaver", str(tmp_path)])
    assert code == 1
    assert "set only one of --otari-dir and --seqweaver" in capsys.readouterr().err

    code = main([*shared, "--seqweaver", str(tmp_path)])
    assert code == 1
    assert "--column-names is required with --seqweaver" in capsys.readouterr().err

    code = main([*shared, "--column-names", str(tmp_path / "names")])
    assert code == 1
    assert "--column-names requires --seqweaver" in capsys.readouterr().err

    code = main([*shared, "--otari-dir", str(tmp_path / "otari")])
    assert code == 1
    assert "--genes is required with --otari-dir or --seqweaver" in capsys.readouterr().err


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
