import json
from pathlib import Path

from pheno.cli_split import main
from pheno.mapping import load_column_map
from pheno.split import TSV_HEADER, load_asd_status, split_directory

FAMILY = """\
person\tfamily\tmother\tfather\tsex\tasd
A\tF1\tmom1\tdad1\tmale\t2
B\tF1\tmom1\tdad1\tfemale\t1
C\tF2\tmom2\tdad2\tmale\t2
U\tF3\tmom3\tdad3\tfemale\t9
"""

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"


def _write(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HEADER + body, encoding="utf-8")


def test_split_assigns_each_variant_to_one_class(tmp_path: Path):
    source = tmp_path / "results"
    _write(
        source / "chr21" / "inherited_10000001_12500000.tsv",
        "chr21\t10\t.\tA\tG\tA\n"
        "chr21\t20\t.\tC\tT\tB\n"
        "chr21\t30\t.\tG\tA\tA;B\n"
        "chr21\t40\t.\tT\tC\tA;U\n"
        "chr21\t50\t.\tA\tC\tU\n"
        "chr21\t60\t.\tC\tG\tZ\n",
    )
    _write(source / "chr21" / "denovo_1_2.tsv", "chr21\t70\t.\tA\tT\tA\n")
    _write(source / "notes" / "inherited_1_2.tsv", "chr21\t80\t.\tA\tT\tA\n")
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    output = tmp_path / "out"

    split_directory(source, output, "inherited", load_asd_status(family, load_column_map(None)))

    asd = output / "chr21" / "inherited_asd_10000001_12500000.tsv"
    unaffected = output / "chr21" / "inherited_unaffected_10000001_12500000.tsv"
    both = output / "chr21" / "inherited_both_10000001_12500000.tsv"
    assert asd.read_text(encoding="utf-8") == (
        TSV_HEADER + "chr21\t10\t.\tA\tG\tA\n" + "chr21\t40\t.\tT\tC\tA\n"
    )
    assert unaffected.read_text(encoding="utf-8") == (
        TSV_HEADER + "chr21\t20\t.\tC\tT\tB\n"
    )
    assert both.read_text(encoding="utf-8") == (
        TSV_HEADER + "chr21\t30\t.\tG\tA\tA;B\n"
    )
    assert not (output / "chr21" / "denovo_asd_1_2.tsv").exists()
    assert not (output / "notes").exists()


def test_split_skips_empty_class_files(tmp_path: Path):
    source = tmp_path / "results"
    _write(source / "chr22" / "inherited.tsv", "chr22\t10\t.\tA\tG\tA\n")
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    output = tmp_path / "out"

    split_directory(source, output, "inherited", load_asd_status(family, load_column_map(None)))

    assert (output / "chr22" / "inherited_asd.tsv").is_file()
    assert not (output / "chr22" / "inherited_unaffected.tsv").exists()
    assert not (output / "chr22" / "inherited_both.tsv").exists()


def test_split_uses_column_map(tmp_path: Path):
    source = tmp_path / "results"
    _write(source / "chrX" / "denovo_00000.tsv", "chrX\t10\t.\tA\tG\tC;B\n")
    family = tmp_path / "family.tsv"
    family.write_text(
        "spid\taffected\nC\ttrue\nB\tfalse\n",
        encoding="utf-8",
    )
    column_map = tmp_path / "columns.json"
    column_map.write_text(
        json.dumps({"columns": {"person": "spid", "asd": "affected"}}),
        encoding="utf-8",
    )
    output = tmp_path / "out"
    rc = main(
        [
            "--family-file",
            str(family),
            "--input-dir",
            str(source),
            "--output-dir",
            str(output),
            "--file-pattern",
            "denovo",
            "--column-map",
            str(column_map),
        ]
    )
    assert rc == 0
    text = (output / "chrX" / "denovo_both_00000.tsv").read_text(encoding="utf-8")
    assert text == TSV_HEADER + "chrX\t10\t.\tA\tG\tC;B\n"


def test_split_rejects_full_format(tmp_path: Path):
    source = tmp_path / "results" / "chr21"
    source.mkdir(parents=True)
    (source / "inherited.tsv").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tTRIO_CALLS\nchr21\t10\t.\tA\tG\tA=0/1|0/0|0/1|30\n",
        encoding="utf-8",
    )
    family = tmp_path / "family.tsv"
    family.write_text(FAMILY, encoding="utf-8")
    rc = main(
        [
            "--family-file",
            str(family),
            "--input-dir",
            str(tmp_path / "results"),
            "--output-dir",
            str(tmp_path / "out"),
            "--file-pattern",
            "inherited",
        ]
    )
    assert rc == 1
