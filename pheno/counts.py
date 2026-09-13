"""Vectorized carrier counts from per-chromosome variant TSVs."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from pheno.pedigree import Cohort

logger = logging.getLogger(__name__)

VARIANT_COLUMNS = ("CHROM", "POS", "REF", "ALT", "PATIENTS")
COUNT_COLUMNS = (
    "a_all",
    "b_all",
    "c_all",
    "d_all",
    "a_male",
    "b_male",
    "c_male",
    "d_male",
    "a_female",
    "b_female",
    "c_female",
    "d_female",
    "n_only_asd",
    "n_only_ctrl",
)


def discover_files(idir: str | Path, file_pattern: str) -> list[Path]:
    directory = Path(idir)
    if not directory.is_dir():
        raise FileNotFoundError(f"input directory does not exist: {directory}")
    files = sorted(
        set(directory.glob(f"{file_pattern}.tsv")) | set(directory.glob(f"{file_pattern}_*.tsv"))
    )
    if not files:
        raise FileNotFoundError(
            f"no files matching {file_pattern}.tsv or {file_pattern}_*.tsv in {directory}"
        )
    return files


def _normalize_variant_columns(frame: pd.DataFrame) -> pd.DataFrame:
    renamed = frame.rename(columns=lambda col: str(col).lstrip("#").strip())
    missing = [name for name in VARIANT_COLUMNS if name not in renamed.columns]
    if missing:
        raise ValueError(f"variant file is missing columns: {missing}")
    return renamed


def counts_from_patient_lists(
    base: pd.DataFrame, patients: pd.Series, cohort: Cohort
) -> pd.DataFrame:
    """Count unique child carriers per row of `base`; rows align with `patients`."""
    empty = _empty_counts(base, cohort)
    if len(base) == 0:
        return empty

    exploded = base.assign(person=patients.fillna("").astype(str).str.split(";")).explode(
        "person", ignore_index=True
    )
    exploded["person"] = exploded["person"].str.strip()
    exploded = exploded[exploded["person"].ne("")]
    exploded["idx"] = exploded["person"].map(cohort.person_to_idx)
    n_unknown = int(exploded["idx"].isna().sum())
    if n_unknown:
        logger.warning("skipped %s patient ids not present in the analysis cohort", n_unknown)
    exploded = exploded.dropna(subset=["idx"])
    if exploded.empty:
        return empty

    exploded["idx"] = exploded["idx"].astype(np.int64)
    exploded = exploded.drop_duplicates(["vid", "idx"], keep="first")
    idx = exploded["idx"].to_numpy()
    exploded = exploded.assign(
        is_asd=cohort.is_asd[idx],
        is_ctrl=cohort.is_ctrl[idx],
        is_male_asd=cohort.is_male_asd[idx],
        is_male_ctrl=cohort.is_male_ctrl[idx],
        is_female_asd=cohort.is_female_asd[idx],
        is_female_ctrl=cohort.is_female_ctrl[idx],
        family_id=cohort.family_id[idx],
        informative=cohort.informative[idx],
    )

    sums = exploded.groupby("vid", sort=False)[
        ["is_asd", "is_ctrl", "is_male_asd", "is_male_ctrl", "is_female_asd", "is_female_ctrl"]
    ].sum()

    out = base.set_index("vid").join(sums).fillna(0)
    out["a_all"] = out["is_asd"].astype(np.int64)
    out["b_all"] = out["is_ctrl"].astype(np.int64)
    out["a_male"] = out["is_male_asd"].astype(np.int64)
    out["b_male"] = out["is_male_ctrl"].astype(np.int64)
    out["a_female"] = out["is_female_asd"].astype(np.int64)
    out["b_female"] = out["is_female_ctrl"].astype(np.int64)
    out["c_all"] = cohort.n_asd - out["a_all"]
    out["d_all"] = cohort.n_ctrl - out["b_all"]
    out["c_male"] = cohort.n_male_asd - out["a_male"]
    out["d_male"] = cohort.n_male_ctrl - out["b_male"]
    out["c_female"] = cohort.n_female_asd - out["a_female"]
    out["d_female"] = cohort.n_female_ctrl - out["b_female"]

    info = exploded.loc[
        exploded["informative"] & (exploded["family_id"] >= 0),
        ["vid", "family_id", "is_asd", "is_ctrl"],
    ]
    asd_pairs = info.loc[info["is_asd"], ["vid", "family_id"]].drop_duplicates()
    ctrl_pairs = info.loc[info["is_ctrl"], ["vid", "family_id"]].drop_duplicates()
    merged = asd_pairs.merge(ctrl_pairs, on=["vid", "family_id"], how="outer", indicator=True)
    only_asd = merged.loc[merged["_merge"] == "left_only"].groupby("vid").size().rename("n_only_asd")
    only_ctrl = (
        merged.loc[merged["_merge"] == "right_only"].groupby("vid").size().rename("n_only_ctrl")
    )
    out = out.join(only_asd).join(only_ctrl)
    out["n_only_asd"] = out["n_only_asd"].fillna(0).astype(np.int64)
    out["n_only_ctrl"] = out["n_only_ctrl"].fillna(0).astype(np.int64)
    drop_helper = [
        "is_asd",
        "is_ctrl",
        "is_male_asd",
        "is_male_ctrl",
        "is_female_asd",
        "is_female_ctrl",
    ]
    return out.drop(columns=[c for c in drop_helper if c in out.columns]).reset_index(drop=True)


def count_variants(frame: pd.DataFrame, cohort: Cohort) -> pd.DataFrame:
    """Count unique child carriers per variant; one output row per input row."""
    frame = _normalize_variant_columns(frame)
    n_var = len(frame)
    base = pd.DataFrame(
        {
            "CHROM": frame["CHROM"].astype(str).to_numpy(),
            "POS": frame["POS"].astype(str).to_numpy(),
            "REF": frame["REF"].astype(str).to_numpy(),
            "ALT": frame["ALT"].astype(str).to_numpy(),
            "vid": np.arange(n_var, dtype=np.int64),
        }
    )
    base["key"] = base["CHROM"] + "_" + base["POS"] + "_" + base["REF"] + "_" + base["ALT"]
    out = counts_from_patient_lists(base, frame["PATIENTS"], cohort)
    return out[["CHROM", "POS", "REF", "ALT", "key", *COUNT_COLUMNS]]


def _empty_counts(base: pd.DataFrame, cohort: Cohort) -> pd.DataFrame:
    out = base.drop(columns=["vid"]).copy() if "vid" in base.columns else base.copy()
    zeros = {
        "a_all": 0,
        "b_all": 0,
        "c_all": cohort.n_asd,
        "d_all": cohort.n_ctrl,
        "a_male": 0,
        "b_male": 0,
        "c_male": cohort.n_male_asd,
        "d_male": cohort.n_male_ctrl,
        "a_female": 0,
        "b_female": 0,
        "c_female": cohort.n_female_asd,
        "d_female": cohort.n_female_ctrl,
        "n_only_asd": 0,
        "n_only_ctrl": 0,
    }
    for name, value in zeros.items():
        out[name] = np.int64(value)
    return out


def count_directory(idir: str | Path, file_pattern: str, cohort: Cohort) -> pd.DataFrame:
    frames = []
    files = discover_files(idir, file_pattern)
    logger.info("reading %s variant files from %s", len(files), idir)
    for path in files:
        logger.info("  %s", path.name)
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        frames.append(count_variants(frame, cohort))
    return pd.concat(frames, ignore_index=True)
