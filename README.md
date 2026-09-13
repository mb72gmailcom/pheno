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

Each row is one variant, keyed by `chr_pos_ref_alt`.

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
