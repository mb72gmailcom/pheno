from pheno.gtf import expand_flanks, load_genes, normalize_chrom, overlap_variants
import pandas as pd


def test_normalize_chrom():
    assert normalize_chrom("21") == "chr21"
    assert normalize_chrom("chr21") == "chr21"
    assert normalize_chrom("CHR21") == "chr21"


def test_expand_flanks_strand():
    assert expand_flanks(1000, 2000, "+", 100, 10) == (900, 2010)
    assert expand_flanks(1000, 2000, "-", 100, 10) == (990, 2100)
    assert expand_flanks(50, 80, "+", 100, 10) == (1, 90)


def test_load_genes_and_overlap(tmp_path):
    gtf = tmp_path / "chr21.gtf"
    gtf.write_text(
        """\
##gff-version 2
chr21\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG1"; gene_type "protein_coding"; gene_name "PLUS";
chr21\tHAVANA\tgene\t1000\t2000\t.\t-\t.\tgene_id "ENSG2"; gene_type "protein_coding"; gene_name "MINUS";
chr21\tHAVANA\tgene\t5000\t6000\t.\t+\t.\tgene_id "ENSG3"; gene_type "lncRNA"; gene_name "LNCR";
chr21\tHAVANA\texon\t1000\t1100\t.\t+\t.\tgene_id "ENSG1";
"""
    )
    genes = load_genes(gtf, upstream=100, downstream=10)
    assert set(genes["gene_id"]) == {"ENSG1", "ENSG2"}
    plus = genes.set_index("gene_id").loc["ENSG1"]
    assert plus["region_start"] == 900
    assert plus["region_end"] == 2010

    variants = pd.DataFrame(
        {
            "vid": [0, 1, 2],
            "CHROM": ["21", "chr21", "chr21"],
            "POS": ["950", "2050", "9000"],
        }
    )
    hits = overlap_variants(variants, genes)
    by_var = hits.groupby("vid")["gene_id"].apply(set)
    assert by_var[0] == {"ENSG1"}
    assert by_var[1] == {"ENSG2"}
    assert 2 not in by_var.index
