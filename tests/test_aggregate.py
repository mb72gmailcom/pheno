from pathlib import Path

import pandas as pd

from pheno.aggregate import aggregate_genes
from pheno.cli_gene import main
from pheno.mapping import load_column_map
from pheno.pedigree import load_cohort

FAMILY = """\
person	family	mother	father	sex	asd
A	F1	mom1	dad1	Male	2
B	F1	mom1	dad1	Female	1
C	F2	mom2	dad2	male	True
D	F2	mom2	dad2	Female	False
"""

GTF = """\
chr21\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG1"; gene_type "protein_coding"; gene_name "PLUS";
chr21\tHAVANA\tgene\t3000\t4000\t.\t+\t.\tgene_id "ENSG2"; gene_type "protein_coding"; gene_name "OTHER";
"""


def _cohort(tmp_path: Path):
    fam = tmp_path / "family.tsv"
    fam.write_text(FAMILY)
    return load_cohort(fam, load_column_map(None))


def test_union_of_patients_and_overlap(tmp_path: Path):
    cohort = _cohort(tmp_path)
    gtf = tmp_path / "chr21.gtf"
    gtf.write_text(GTF)
    idir = tmp_path / "chr21"
    idir.mkdir()
    header = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"
    (idir / "inherited.tsv").write_text(header + "chr21\t1500\t.\tC\tT\tA\n")
    (idir / "inherited_1000_2000.tsv").write_text(header + "chr21\t1600\t.\tC\tA\tA;B\n")
    (idir / "inherited_9000_9100.tsv").write_text(header + "chr21\t9000\t.\tG\tT\tC\n")

    out = aggregate_genes(idir, "inherited", gtf, cohort, upstream=0, downstream=0).set_index(
        "gene_id"
    )
    assert list(out.index) == ["ENSG1"]
    assert out.loc["ENSG1", "n_variants"] == 2
    assert out.loc["ENSG1", "a_all"] == 1
    assert out.loc["ENSG1", "b_all"] == 1
    assert out.loc["ENSG1", "n_only_asd"] == 0
    assert out.loc["ENSG1", "n_only_ctrl"] == 0


def test_cli_gene(tmp_path: Path):
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY)
    gtf = tmp_path / "chr21.gtf"
    gtf.write_text(GTF)
    idir = tmp_path / "chr21"
    idir.mkdir()
    (idir / "inherited.tsv").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\nchr21\t1500\t.\tC\tT\tA\n"
    )
    output = tmp_path / "genes.tsv"
    rc = main(
        [
            "--family-file",
            str(family),
            "--idir",
            str(idir),
            "--file-pattern",
            "inherited",
            "--gtf",
            str(gtf),
            "--upstream",
            "0",
            "--downstream",
            "0",
            "--output",
            str(output),
        ]
    )
    assert rc == 0
    out = pd.read_csv(output, sep="\t")
    assert list(out["gene_id"]) == ["ENSG1"]
    assert out.loc[0, "a_all"] == 1
    assert out.loc[0, "n_only_asd"] == 1
    assert "mcnemar_log10p" in out.columns
    assert "gene_name" in out.columns
