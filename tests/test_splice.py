import csv
from pathlib import Path

import pytest

from pheno.cli_splice import main
from pheno.splice import load_splice_sites, write_splice_distances

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"
GTF = """\
chr21\tHAVANA\ttranscript\t100\t500\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; transcript_type "protein_coding"; gene_name "PLUS";
chr21\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; transcript_type "protein_coding"; gene_name "PLUS";
chr21\tHAVANA\texon\t300\t400\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; transcript_type "protein_coding"; gene_name "PLUS";
chr21\tHAVANA\texon\t450\t500\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; transcript_type "protein_coding"; gene_name "PLUS";
chr21\tHAVANA\ttranscript\t100\t400\t.\t-\t.\tgene_id "ENSG2"; transcript_id "ENST2"; transcript_type "protein_coding"; gene_name "MINUS";
chr21\tHAVANA\texon\t100\t150\t.\t-\t.\tgene_id "ENSG2"; transcript_id "ENST2"; transcript_type "protein_coding"; gene_name "MINUS";
chr21\tHAVANA\texon\t250\t400\t.\t-\t.\tgene_id "ENSG2"; transcript_id "ENST2"; transcript_type "protein_coding"; gene_name "MINUS";
chr21\tHAVANA\ttranscript\t1000\t1100\t.\t+\t.\tgene_id "ENSG3"; transcript_id "ENST3"; transcript_type "protein_coding"; gene_name "ONE";
chr21\tHAVANA\texon\t1000\t1100\t.\t+\t.\tgene_id "ENSG3"; transcript_id "ENST3"; transcript_type "protein_coding"; gene_name "ONE";
chr21\tHAVANA\ttranscript\t100\t400\t.\t+\t.\tgene_id "ENSG4"; transcript_id "ENST4"; transcript_type "retained_intron"; gene_type "protein_coding"; gene_name "SKIP";
chr21\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id "ENSG4"; transcript_id "ENST4"; transcript_type "retained_intron"; gene_name "SKIP";
chr21\tHAVANA\texon\t300\t400\t.\t+\t.\tgene_id "ENSG4"; transcript_id "ENST4"; transcript_type "retained_intron"; gene_name "SKIP";
chr22\tHAVANA\ttranscript\t100\t400\t.\t+\t.\tgene_id "ENSG5"; transcript_id "ENST5"; transcript_type "protein_coding"; gene_name "OTHER";
chr22\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id "ENSG5"; transcript_id "ENST5"; transcript_type "protein_coding"; gene_name "OTHER";
chr22\tHAVANA\texon\t300\t400\t.\t+\t.\tgene_id "ENSG5"; transcript_id "ENST5"; transcript_type "protein_coding"; gene_name "OTHER";
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
    (source / "chr21" / "inherited_100_200.tsv").write_text(
        HEADER
        + "21\t201\t.\tC\tT\tA\n"
        + "21\t200\t.\tC\tT\tB\n"
        + "21\t299\t.\tA\tG\tC\n"
        + "21\t401\t.\tG\tA\tD\n"
        + "21\t249\t.\tC\tT\tE\n"
        + "21\t250\t.\tC\tT\tF\n"
        + "21\t200\t.\tAT\tA\tG\n"
        + "21\t200\t.\tA\tAT\tH\n"
        + "21\t1000\t.\tC\tT\tI\n",
        encoding="utf-8",
    )
    (source / "chr22" / "inherited_00000.tsv").write_text(
        HEADER + "22\t201\t.\tA\tG\tJ\n",
        encoding="utf-8",
    )
    return source, gtf


def test_junctions_follow_transcription_direction(tmp_path: Path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GTF, encoding="utf-8")
    sites = load_splice_sites(gtf)
    plus = {(site.site_type, site.site) for site in sites["chr21"] if site.transcript_id == "ENST1"}
    assert plus == {("donor", 201), ("acceptor", 299), ("donor", 401), ("acceptor", 449)}
    minus = {(site.site_type, site.site) for site in sites["chr21"] if site.transcript_id == "ENST2"}
    assert minus == {("donor", 249), ("acceptor", 151)}
    assert all(site.transcript_id != "ENST3" for site in sites["chr21"])
    assert all(site.transcript_id != "ENST4" for site in sites["chr21"])


def test_distance_and_coverage(tmp_path: Path):
    source, gtf = _write_run(tmp_path)
    output = tmp_path / "out"
    written = write_splice_distances(source, output, "inherited", gtf, window=30, chrom="21")
    rows = _rows(output / "chr21" / "splice.tsv")
    assert written == len(rows)
    by_key = {(row["pos"], row["ref"], row["alt"], row["transcript_id"], row["site_type"]): row for row in rows}

    donor = by_key[("201", "C", "T", "ENST1", "donor")]
    assert donor["site"] == "201"
    assert donor["distance"] == "0"
    assert donor["covers_site"] == "1"
    assert donor["strand"] == "+"

    upstream = by_key[("200", "C", "T", "ENST1", "donor")]
    assert upstream["distance"] == "-1"
    assert upstream["covers_site"] == "0"

    acceptor = by_key[("299", "A", "G", "ENST1", "acceptor")]
    assert acceptor["distance"] == "0"
    assert acceptor["covers_site"] == "1"

    second_donor = by_key[("401", "G", "A", "ENST1", "donor")]
    assert second_donor["site"] == "401"

    minus_donor = by_key[("249", "C", "T", "ENST2", "donor")]
    assert minus_donor["distance"] == "0"
    assert minus_donor["covers_site"] == "1"
    minus_upstream = by_key[("250", "C", "T", "ENST2", "donor")]
    assert minus_upstream["distance"] == "-1"

    deletion = by_key[("200", "AT", "A", "ENST1", "donor")]
    assert deletion["distance"] == "0"
    assert deletion["covers_site"] == "1"
    insertion = by_key[("200", "A", "AT", "ENST1", "donor")]
    assert insertion["distance"] == "0"
    assert insertion["covers_site"] == "0"

    assert ("1000", "C", "T", "ENST3", "donor") not in by_key
    assert not (output / "chr22").exists()


def test_indels_and_chromosome_limit(tmp_path: Path):
    source, gtf = _write_run(tmp_path)
    output = tmp_path / "out"
    write_splice_distances(source, output, "inherited", gtf, window=30, chrom="21", indels_only=True)
    alleles = {(row["ref"], row["alt"]) for row in _rows(output / "chr21" / "splice.tsv")}
    assert alleles == {("AT", "A"), ("A", "AT")}

    both = tmp_path / "both"
    write_splice_distances(source, both, "inherited", gtf, window=30)
    assert (both / "chr22" / "splice.tsv").is_file()


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
        write_splice_distances(source, tmp_path / "out", "missing", gtf, chrom="21")
