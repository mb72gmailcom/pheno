import csv
from pathlib import Path

import pytest

from pheno.cli_tss import main
from pheno.tss import allele_interval, write_tss_distances

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"
GTF = """\
chr21\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG1"; gene_type "protein_coding"; gene_name "PLUS";
chr21\tHAVANA\ttranscript\t1000\t1800\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; transcript_type "protein_coding"; gene_name "PLUS";
chr21\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1B"; transcript_type "retained_intron"; gene_type "protein_coding"; gene_name "PLUS";
chr21\tHAVANA\ttranscript\t500\t1000\t.\t-\t.\tgene_id "ENSG2"; transcript_id "ENST2"; transcript_type "protein_coding"; gene_name "MINUS";
chr21\tHAVANA\ttranscript\t5000\t6000\t.\t+\t.\tgene_id "ENSG3"; transcript_id "ENST3"; transcript_type "protein_coding"; gene_name "FAR";
chr22\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "ENSG4"; transcript_id "ENST4"; transcript_type "protein_coding"; gene_name "OTHER";
chrX\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "ENSGX"; transcript_id "ENSTX"; transcript_type "protein_coding"; gene_name "CHRX";
"""


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_run(root: Path) -> tuple[Path, Path]:
    gtf = root / "genes.gtf"
    gtf.write_text(GTF, encoding="utf-8")
    source = root / "run"
    (source / "chr21").mkdir(parents=True)
    (source / "chr22").mkdir()
    (source / "chrX").mkdir()
    (source / "chr21" / "inherited_100_200.tsv").write_text(
        HEADER
        + "21\t990\t.\tC\tT\tA\n"
        + "21\t1000\t.\tA\tG\tB\n"
        + "21\t2000\t.\tC\tT\tC\n"
        + "21\t2001\t.\tC\tT\tD\n"
        + "21\t999\t.\tAT\tA\tE\n"
        + "21\t1000\t.\tAT\tA\tF\n"
        + "21\t999\t.\tA\tAT\tG\n"
        + "21\t999\t.\tAA\tGG\tH\n"
        + "21\t1010\t.\tC\tT\tI\n"
        + "21\t1000\t.\tT\tAT\tJ\n",
        encoding="utf-8",
    )
    (source / "chr21" / "inherited_females_00000.tsv").write_text(
        HEADER + "21\t1005\t.\tG\tA\tK\n",
        encoding="utf-8",
    )
    (source / "chr21" / "denovo_00000.tsv").write_text(
        HEADER + "21\t1000\t.\tC\tT\tNO\n",
        encoding="utf-8",
    )
    (source / "chr21" / "inherited_per_variant.json").write_text("{", encoding="utf-8")
    (source / "chr22" / "inherited_00000.tsv").write_text(
        HEADER + "22\t1000\t.\tA\tG\tL\n",
        encoding="utf-8",
    )
    (source / "chrX" / "inherited_males_nonPar_00001.tsv").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tTRIO_CALLS\nX\t1000\t.\tA\tG\tchild=0/1|0/0|0/1|30\n",
        encoding="utf-8",
    )
    return source, gtf


def test_allele_interval_anchor_and_alignment():
    assert allele_interval(100, "A", "G") == (100, 100, True)
    assert allele_interval(100, "AA", "GG") == (100, 101, True)
    assert allele_interval(100, "AT", "A") == (101, 101, True)
    assert allele_interval(100, "TA", "A") == (100, 100, True)
    assert allele_interval(100, "A", "AT") == (101, 101, False)
    assert allele_interval(100, "T", "AT") == (100, 100, False)
    assert allele_interval(100, "ATG", "GC") == (100, 102, True)


def test_signed_distance_and_tss_coverage(tmp_path: Path):
    source, gtf = _write_run(tmp_path)
    output = tmp_path / "out"
    written = write_tss_distances(source, output, "inherited", gtf, window=1000, chrom="21")
    rows = _rows(output / "chr21" / "tss.tsv")
    assert written == len(rows)
    by_key = {(row["pos"], row["ref"], row["alt"], row["transcript_id"]): row for row in rows}

    upstream = by_key[("990", "C", "T", "ENST1")]
    assert upstream["distance"] == "-10"
    assert upstream["covers_tss"] == "0"
    assert upstream["strand"] == "+"
    assert upstream["patients"] == "A"

    at_tss = by_key[("1000", "A", "G", "ENST1")]
    assert at_tss["distance"] == "0"
    assert at_tss["covers_tss"] == "1"

    edge = by_key[("2000", "C", "T", "ENST1")]
    assert edge["distance"] == "1000"
    assert ("2001", "C", "T", "ENST1") not in by_key

    deletion = by_key[("999", "AT", "A", "ENST1")]
    assert deletion["distance"] == "0"
    assert deletion["covers_tss"] == "1"
    downstream_deletion = by_key[("1000", "AT", "A", "ENST1")]
    assert downstream_deletion["distance"] == "1"
    assert downstream_deletion["covers_tss"] == "0"

    insertion = by_key[("999", "A", "AT", "ENST1")]
    assert insertion["distance"] == "0"
    assert insertion["covers_tss"] == "0"
    right_insertion = by_key[("1000", "T", "AT", "ENST1")]
    assert right_insertion["distance"] == "0"
    assert right_insertion["covers_tss"] == "0"

    mnv = by_key[("999", "AA", "GG", "ENST1")]
    assert mnv["covers_tss"] == "1"

    minus = by_key[("1010", "C", "T", "ENST2")]
    assert minus["distance"] == "-10"
    assert minus["strand"] == "-"
    assert minus["tss"] == "1000"
    assert "ENST1B" not in {row["transcript_id"] for row in rows}
    assert not (output / "chr22").exists()


def test_indels_flag_drops_equal_length_alleles(tmp_path: Path):
    source, gtf = _write_run(tmp_path)
    output = tmp_path / "out"
    write_tss_distances(source, output, "inherited", gtf, window=1000, chrom="21", indels_only=True)
    alleles = {(row["ref"], row["alt"]) for row in _rows(output / "chr21" / "tss.tsv")}
    assert alleles == {("AT", "A"), ("A", "AT"), ("T", "AT")}


def test_prefix_chr_and_sex_bucket(tmp_path: Path):
    source, gtf = _write_run(tmp_path)
    output = tmp_path / "out"
    write_tss_distances(source, output, "inherited", gtf)
    assert (output / "chr21" / "tss.tsv").is_file()
    assert (output / "chr22" / "tss.tsv").is_file()
    sex = _rows(output / "chrX" / "tss.tsv")
    assert sex[0]["patients"] == "child=0/1|0/0|0/1|30"
    assert sex[0]["chrom"] == "chrX"
    assert "NO" not in {row["patients"] for row in _rows(output / "chr21" / "tss.tsv")}

    limited = tmp_path / "one"
    write_tss_distances(source, limited, "denovo", gtf, chrom="chr21")
    denovo = _rows(limited / "chr21" / "tss.tsv")
    assert {row["patients"] for row in denovo} == {"NO"}
    assert not (limited / "chr22").exists()


def test_cli_rejects_a_chromosome_list(tmp_path: Path):
    source, gtf = _write_run(tmp_path)
    rc = main(
        [
            "--input-dir",
            str(source),
            "--output-dir",
            str(tmp_path / "out"),
            "--prefix",
            "inherited",
            "--gtf",
            str(gtf),
            "--chr",
            "21,22",
        ]
    )
    assert rc == 1


def test_missing_prefix_is_an_error(tmp_path: Path):
    source, gtf = _write_run(tmp_path)
    with pytest.raises(ValueError, match="no missing\\*\\.tsv"):
        write_tss_distances(source, tmp_path / "out", "missing", gtf, chrom="21")
