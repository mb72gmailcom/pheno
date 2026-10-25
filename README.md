# pheno

Count inherited-variant carriers among sequenced trio children and test association with ASD. Input is one chromosome directory of variant TSVs; output is one TSV for that chromosome.

The primary test is **McNemar exact** (within-family, discordant sibships). **Fisher exact** 2×2 tables for all / male / female children are extra columns.

## Install

Python 3.10+. From the repo root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Run

`--idir` must already be split by chromosome. File names do not include the chromosome; they match `{pattern}.tsv` and `{pattern}_*.tsv` (for example `inherited.tsv`, `inherited_1000_2000.tsv`).

```bash
pheno-asd \
  --family-file family.tsv \
  --idir /path/to/chr21 \
  --file-pattern inherited \
  --column-map examples/column_map.json \
  --output chr21.tsv
```

`--column-map` is optional if the family file already uses the default column names (`person`, `family`, `mother`, `father`, `sex`, `asd`).

Gene-level burden (separate command): union of carriers for every rare variant whose `POS` falls in a GTF gene plus strand-aware flanks, then the same Fisher / McNemar tests. One output row per gene that has at least one overlapping variant. Default flanks are 5 kb upstream of the TSS and 1 kb downstream of the TES; default biotype is `protein_coding`.

```bash
pheno-asd-gene \
  --family-file family.tsv \
  --idir /path/to/chr21 \
  --file-pattern inherited \
  --gtf ~/resources/annotations_by_chrm/chr21/chr21.gencode.v47.basic.annotation.gtf \
  --column-map examples/column_map.json \
  --upstream 5000 \
  --downstream 1000 \
  --output chr21.genes.tsv
```

Use `--all-gene-types` to keep every GTF `gene` feature, or `--gene-types protein_coding,lncRNA` to choose types.

Split each variant into one of three short-format TSVs. A row is `asd` when every remaining carrier is affected, `unaffected` when every remaining carrier is unaffected, and `both` when both groups are present. People missing from the family file, or with an unrecognized ASD value, are dropped before that decision. A class with no variants is not written.

```bash
pheno-asd-split \
  --family-file family.tsv \
  --input-dir /path/to/results \
  --file-pattern inherited \
  --column-map examples/column_map.json \
  --output-dir /path/to/split
```

`--input-dir` is the parent of the `chrN` directories. For `chr21/inherited_10000001_12500000.tsv` the outputs are `inherited_asd_10000001_12500000.tsv`, `inherited_unaffected_10000001_12500000.tsv`, and `inherited_both_10000001_12500000.tsv` under `output-dir/chr21/`. The `PATIENTS` column keeps only the carriers who belong in that file.

Count variants per patient. Only people who appear in the selected TSVs are included. Their ASD status comes from the family file. Ids missing from the family file, or with an unrecognized ASD value, are not counted.

```bash
pheno-asd-counts \
  --family-file family.tsv \
  --input-dir /path/to/results \
  --file-pattern inherited \
  --column-map examples/column_map.json \
  --output counts.json
```

`--input-dir` is the parent of the `chrN` directories and must contain short-format TSVs with a `PATIENTS` column. The JSON object has two maps, `asd` and `unaffected`, from person id to the number of matching rows across all chromosomes.

Within-family dosage contrast for variants that appear in a family with both an ASD child and an unaffected sibling. Dosage is 0 for `0/0`, 1 for `0/1` or `1/0`, and 2 for `1/1`. For each such family, `d` is the mean ASD dosage minus the mean unaffected dosage. `T` is the sum of those differences. A family is informative when the siblings do not all have the same dosage. The z-score is `T` divided by the square root of the sum of the within-family null variances. The exact p-value is the two-sided probability of an absolute sum at least as large as `|T|` when the ASD labels are reassigned inside each family. A sibling who belongs to a family touched by the variant and is absent from the row is dosage 0. Multiallelic variants are skipped.

```bash
pheno-asd-transmission \
  --family-file family.tsv \
  --input-dir /path/to/results \
  --file-pattern inherited \
  --column-map examples/column_map.json \
  --output transmission.tsv
```

`--output` writes every chromosome into that one file. `--output-dir` writes one file per chromosome, `chrN/inherited.tsv`, ordered by position. Set exactly one of the two. The columns are `CHROM`, `POS`, `ID`, `REF`, `ALT`, `n_families`, `n_informative`, `z`, and `p_exact`. `n_families` counts sibships with both groups genotyped on that row. `z` is blank when `n_informative` is 0, and `p_exact` is then 1. A chromosome with no qualifying variant is omitted from `--output-dir`.

Per-person Otari burden. For each carried variant, transcript scores are made absolute unless `--no-abs` is set, then collapsed across transcripts. A scored variant is damaging when that collapsed score is greater than `--threshold`.

```bash
pheno-asd-burden \
  --family-file family.tsv \
  --input-dir /path/to/results \
  --otari-dir /path/to/otari/inherited \
  --file-pattern inherited_asd \
  --column-map examples/column_map.json \
  --otari-columns max_effect,Brain \
  --threshold 0.5 \
  --output-dir /path/to/burden
```

Variant TSVs must be named `{file-pattern}_{start}_{end}.tsv`. The Otari file for that shard is `{otari-dir}/chrN/{start}_{end}/variant_effects_comprehensive.tsv`. A row `chr21 25019786 G T` matches Otari id `21_25019786_G_T_hg38`. Each chromosome is written to `{output-dir}/chrN/burden.json` as it finishes. `{output-dir}/burden.json` is the sum of those files. Counts and sums are added, `max_score` is the maximum across chromosomes, and `fraction_damaging` is recomputed from the summed `n_damaging` and `n_scored`. A person is left out of a chromosome file when they carry no variant on that chromosome. Each person has `n_variants`, `n_unscored`, and, for every selected column, `n_scored`, `n_damaging`, `fraction_damaging`, `sum_effect`, `sum_damaging`, and `max_score`.

Per-gene Otari burden for ASD children and their unaffected siblings. A variant is assigned to a gene when it falls in the gene body or within `--window` base pairs of either end. The default window is 2000. Transcripts are collapsed within that gene.

```bash
pheno-asd-burden-gene \
  --family-file family.tsv \
  --input-dir /path/to/results \
  --otari-dir /path/to/otari/inherited \
  --annotation resources/gencode.v47.basic.annotation.clean.gtf.gz \
  --file-pattern inherited \
  --column-map examples/column_map.json \
  --otari-columns max_effect \
  --threshold 0.5 \
  --output-dir /path/to/gene-burden
```

An unaffected sibling is an unaffected child who shares a family id with an ASD child. Parents are left out. While a chromosome is processed, each shard is written to a temporary file and removed after `{output-dir}/chrN/genes.tsv.gz` and `{output-dir}/chrN/gene_summary.json` are complete. `genes.tsv.gz` has one row per person and gene. The `gene` column is the gene id with the version removed. `gene_name` and `gene_type` are copied from the annotation `attributes` field, and either is left blank when that attribute is absent. `{output-dir}/gene_summary.json` joins the chromosome summaries. `n_carriers` counts people with at least one variant in the gene. `mean_sum_effect` is the mean of their `sum_effect` values. `var_mean_effect` is `sum_effect / n_scored` for that person, and it is empty when `n_scored` is 0. `mean_var_mean_effect` is the mean of those values among carriers who have a score. `n_scored_carriers` is that count. The gene for each Otari transcript is read from `interpretability_analysis.tsv` in the same Otari shard.

The same command reads a variant score table when `--score-dir` is set instead of `--otari-dir`. `--score-prefix` is the filename prefix, and `--score-columns` names the numeric columns. A carrier row matches a score row on chromosome, position, ref, and alt. A blank cell is unscored. The stored number is used as it is: there is no absolute value and no transcript collapse. For `ag_PHRED`, pass that column alone.

```bash
pheno-asd-burden-gene \
  --family-file family.tsv \
  --input-dir /path/to/results \
  --score-dir /path/to/scores \
  --annotation resources/gencode.v47.basic.annotation.clean.gtf.gz \
  --file-pattern inherited \
  --score-prefix variants_enrich \
  --score-columns ag_PHRED \
  --threshold 1 \
  --output-dir /path/to/gene-burden
```

## Family file

TSV with one row per person. A **child** is anyone with both parents defined (IDs not `0`, `false`, `-`, or empty) and a known ASD status. Missing ASD excludes the person. Missing sex keeps them in the `*_all` columns only.

Default value codes:

| Field | Affected / male | Unaffected / female |
| --- | --- | --- |
| ASD | `true`, `2` | `false`, `1` |
| Sex | `male`, `1` | `female`, `2` |

If `mother_sequenced` / `father_sequenced` are mapped, both must be true (`true`, `1`, `yes`). If those columns are absent, set them to `null` or `"-"` in the JSON map.

## Variant files

```
#CHROM  POS  ID  REF  ALT  PATIENTS
chr21   7927554  .  C  T  SP0201052
chr21   7927704  .  C  A  SP0155907;SP0294915
```

`PATIENTS` is a `;`-separated list of person ids. Counts are **carrier presence**, not allele count. Ids not in the family analysis set are skipped.

## Output columns

`pheno-asd` writes one row per variant, keyed by `chr_pos_ref_alt`. `pheno-asd-gene` writes one row per gene (`gene_id`, `gene_name`, `n_variants`, locus and flank coordinates) and the same count / p-value columns. A child is a gene carrier if they carry **any** overlapping variant (union of `PATIENTS`, each person once).

For each stratum (`all`, `male`, `female`):

| Columns | Meaning |
| --- | --- |
| `a_*` | ASD carriers (in `PATIENTS`) |
| `b_*` | Unaffected carriers |
| `c_*` | ASD non-carriers = N(ASD) − `a` |
| `d_*` | Unaffected non-carriers = N(unaffect) − `b` |
| `or_*` | Haldane-Anscombe odds ratio |
| `fisher_p_*`, `fisher_log10p_*` | Two-sided Fisher exact |

`N(ASD)` and `N(unaffect)` are the cohort sizes from the family file (children with both parents and known ASD), not the number of people who appear in the variant files.

**McNemar** (all children, primary):

| Columns | Meaning |
| --- | --- |
| `n_only_asd` | Families with ≥1 ASD child and ≥1 unaffect child where the ASD side carries and the unaffect side does not |
| `n_only_ctrl` | Same families, only the unaffect side carries |
| `mcnemar_p`, `mcnemar_log10p` | Exact two-sided binomial test of those two counts (null p = 0.5) |

Families where both sides carry, or neither carries, are uninformative. Each family contributes at most 1. There is no pair picking when a family has extra siblings.

## Tests

```bash
pytest
```

## License

MIT
