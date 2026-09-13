from pathlib import Path

import pandas as pd

from pheno.cli import main

FAMILY = """\
person	family	mother	father	sex	asd
A	F1	mom1	dad1	Male	2
B	F1	mom1	dad1	Female	1
C	F2	mom2	dad2	male	True
D	F2	mom2	dad2	Female	False
"""

VARIANTS = """\
#CHROM	POS	ID	REF	ALT	PATIENTS
chr21	7927554	.	C	T	A
chr21	7927704	.	C	A	A;B
"""


def test_cli_end_to_end(tmp_path: Path):
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY)
    idir = tmp_path / "chr21"
    idir.mkdir()
    (idir / "inherited.tsv").write_text(VARIANTS)
    (idir / "inherited_1000_2000.tsv").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\nchr21\t7928720\t.\tC\tT\tC\n"
    )
    output = tmp_path / "chr21.tsv"
    rc = main(
        [
            "--family-file",
            str(family),
            "--idir",
            str(idir),
            "--file-pattern",
            "inherited",
            "--output",
            str(output),
        ]
    )
    assert rc == 0
    out = pd.read_csv(output, sep="\t")
    assert len(out) == 3
    assert list(out["POS"].astype(str)) == ["7927554", "7927704", "7928720"]
    assert out.loc[0, "a_all"] == 1
    assert out.loc[0, "n_only_asd"] == 1
    assert out.loc[1, "n_only_asd"] == 0
    assert "mcnemar_log10p" in out.columns
    assert "fisher_p_female" in out.columns
