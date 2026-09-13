"""Build analysis sets from a family TSV."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from pheno.mapping import ColumnMap

logger = logging.getLogger(__name__)


@dataclass
class Cohort:
    """Children with both parents defined and a known ASD status."""

    person_to_idx: dict[str, int]
    is_asd: np.ndarray
    is_ctrl: np.ndarray
    is_male_asd: np.ndarray
    is_male_ctrl: np.ndarray
    is_female_asd: np.ndarray
    is_female_ctrl: np.ndarray
    family_id: np.ndarray
    informative: np.ndarray
    n_asd: int
    n_ctrl: int
    n_male_asd: int
    n_male_ctrl: int
    n_female_asd: int
    n_female_ctrl: int
    n_informative_families: int
    n_children: int


def _require_columns(frame: pd.DataFrame, names: list[str], source: str) -> None:
    missing = [name for name in names if name not in frame.columns]
    if missing:
        raise ValueError(f"{source} is missing columns: {missing}")


def load_cohort(family_file: str | Path, column_map: ColumnMap) -> Cohort:
    frame = pd.read_csv(family_file, sep="\t", dtype=str, keep_default_na=False)
    needed = [
        column_map.person,
        column_map.family,
        column_map.mother,
        column_map.father,
        column_map.sex,
        column_map.asd,
    ]
    if column_map.mother_sequenced:
        needed.append(column_map.mother_sequenced)
    if column_map.father_sequenced:
        needed.append(column_map.father_sequenced)
    _require_columns(frame, needed, "family file")

    person = frame[column_map.person].str.strip()
    family = frame[column_map.family].str.strip()
    mother = frame[column_map.mother].str.strip()
    father = frame[column_map.father].str.strip()
    sex = frame[column_map.sex]
    asd = frame[column_map.asd]

    has_parents = ~(mother.map(column_map.is_bad_id) | father.map(column_map.is_bad_id))
    if column_map.mother_sequenced:
        has_parents &= frame[column_map.mother_sequenced].map(column_map.is_sequenced)
    if column_map.father_sequenced:
        has_parents &= frame[column_map.father_sequenced].map(column_map.is_sequenced)

    asd_status = asd.map(column_map.asd_status)
    sex_status = sex.map(column_map.sex_status)
    keep = has_parents & asd_status.notna()

    kids = pd.DataFrame(
        {
            "person": person[keep].to_numpy(),
            "family": family[keep].to_numpy(),
            "asd": asd_status[keep].to_numpy(),
            "sex": sex_status[keep].to_numpy(),
        }
    )
    n_before = len(kids)
    kids = kids.drop_duplicates(subset="person", keep="first")
    if len(kids) < n_before:
        logger.warning("dropped %s duplicate person ids in family file", n_before - len(kids))

    if kids.empty:
        raise ValueError("no children with both parents defined and known ASD status")

    person_to_idx = {pid: i for i, pid in enumerate(kids["person"].tolist())}
    is_asd = (kids["asd"].to_numpy() == "asd")
    is_ctrl = (kids["asd"].to_numpy() == "ctrl")
    is_male = (kids["sex"].to_numpy() == "male")
    is_female = (kids["sex"].to_numpy() == "female")

    family_codes, _uniques = pd.factorize(kids["family"].replace("", np.nan), use_na_sentinel=True)
    family_id = family_codes.astype(np.int64)

    informative = np.zeros(len(kids), dtype=bool)
    valid_fam = family_id >= 0
    if valid_fam.any():
        asd_by_fam = pd.Series(is_asd[valid_fam]).groupby(family_id[valid_fam]).any()
        ctrl_by_fam = pd.Series(is_ctrl[valid_fam]).groupby(family_id[valid_fam]).any()
        both = asd_by_fam & ctrl_by_fam
        informative_ids = set(both[both].index.tolist())
        informative = np.array([fid in informative_ids for fid in family_id], dtype=bool)
        n_informative = len(informative_ids)
    else:
        n_informative = 0

    cohort = Cohort(
        person_to_idx=person_to_idx,
        is_asd=is_asd,
        is_ctrl=is_ctrl,
        is_male_asd=is_asd & is_male,
        is_male_ctrl=is_ctrl & is_male,
        is_female_asd=is_asd & is_female,
        is_female_ctrl=is_ctrl & is_female,
        family_id=family_id,
        informative=informative,
        n_asd=int(is_asd.sum()),
        n_ctrl=int(is_ctrl.sum()),
        n_male_asd=int((is_asd & is_male).sum()),
        n_male_ctrl=int((is_ctrl & is_male).sum()),
        n_female_asd=int((is_asd & is_female).sum()),
        n_female_ctrl=int((is_ctrl & is_female).sum()),
        n_informative_families=n_informative,
        n_children=len(kids),
    )
    logger.info(
        "cohort: %s children (%s ASD, %s unaffect); %s informative families",
        cohort.n_children,
        cohort.n_asd,
        cohort.n_ctrl,
        cohort.n_informative_families,
    )
    return cohort
