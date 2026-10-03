import json
from pathlib import Path

from pheno.cli_patient_counts import main
from pheno.mapping import load_column_map
from pheno.patient_counts import count_patient_variants
from pheno.split import load_asd_status

FAMILY = """\
person	family	mother	father	sex	asd
A	F1	mom1	dad1	male	2
B	F1	mom1	dad1	female	1
C	F2	mom2	dad2	male	2
U	F3	mom3	dad3	female	9
"""

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"


def _write(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HEADER + body, encoding="utf-8")


def test_counts_include_zeros_and_sum_chromosomes(tmp_path: Path):
    source = tmp_path / "results"
    _write(
        source / "chr21" / "inherited_10000001_12500000.tsv",
        "chr21\t10\t.\tA\tG\tA\n"
        "chr21\t20\t.\tC\tT\tA;B\n"
        "chr21\t30\t.\tG\tA\tZ\n",
    )
    _write(source / "chr22" / "inherited.tsv", "chr22\t40\t.\tT\tC\tA\n")
    _write(source / "chr22" / "denovo.tsv", "chr22\t50\t.\tA\tC\tB\n")
    _write(source / "notes" / "inherited.tsv", "chr21\t60\t.\tA\tT\tA\n")
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")

    counts = count_patient_variants(
        source, "inherited", load_asd_status(family, load_column_map(None))
    )

    assert counts == {
        "asd": {"A": 3, "C": 0},
        "unaffected": {"B": 1},
    }


def test_counts_cli_writes_json(tmp_path: Path):
    source = tmp_path / "results"
    _write(source / "chrX" / "denovo_00000.tsv", "chrX\t10\t.\tA\tG\tC;B;C\n")
    family = tmp_path / "family.tsv"
    family.write_text("spid\taffected\nC\ttrue\nB\tfalse\nD\t2\n", encoding="utf-8")
    column_map = tmp_path / "columns.json"
    column_map.write_text(
        json.dumps({"columns": {"person": "spid", "asd": "affected"}}),
        encoding="utf-8",
    )
    output = tmp_path / "counts.json"
    rc = main(
        [
            "--family-file",
            str(family),
            "--input-dir",
            str(source),
            "--output",
            str(output),
            "--file-pattern",
            "denovo",
            "--column-map",
            str(column_map),
        ]
    )
    assert rc == 0
    assert json.loads(output.read_text(encoding="utf-8")) == {
        "asd": {"C": 1, "D": 0},
        "unaffected": {"B": 1},
    }
