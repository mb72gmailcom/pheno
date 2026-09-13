"""Fisher exact (extra columns) and McNemar exact (primary) p-values."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import binom, fisher_exact

STRATA = ("all", "male", "female")


def haldane_odds_ratio(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> np.ndarray:
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def neglog10(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    finite = np.isfinite(p) & (p > 0)
    out[finite] = -np.log10(p[finite])
    out[np.isfinite(p) & (p == 0)] = np.inf
    return out


def fisher_p_vectorized(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> np.ndarray:
    """Two-sided Fisher p-values, cached over unique 2x2 tables."""
    tables = np.column_stack([a, b, c, d]).astype(np.int64, copy=False)
    uniq, inverse = np.unique(tables, axis=0, return_inverse=True)
    cached = np.empty(len(uniq), dtype=float)
    for i, (aa, bb, cc, dd) in enumerate(uniq):
        if min(aa, bb, cc, dd) < 0:
            cached[i] = np.nan
            continue
        cached[i] = fisher_exact([[int(aa), int(bb)], [int(cc), int(dd)]], alternative="two-sided")[
            1
        ]
    return cached[inverse]


def mcnemar_p_vectorized(n_only_asd: np.ndarray, n_only_ctrl: np.ndarray) -> np.ndarray:
    """Exact two-sided McNemar / binomial test, p = 0.5 under the null."""
    n_only_asd = np.asarray(n_only_asd, dtype=np.int64)
    n_only_ctrl = np.asarray(n_only_ctrl, dtype=np.int64)
    n = n_only_asd + n_only_ctrl
    p = np.ones(n.shape, dtype=float)
    mask = n > 0
    k = np.minimum(n_only_asd[mask], n_only_ctrl[mask])
    p[mask] = np.minimum(1.0, 2.0 * binom.cdf(k, n[mask], 0.5))
    return p


def add_stats(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for stratum in STRATA:
        a = out[f"a_{stratum}"].to_numpy()
        b = out[f"b_{stratum}"].to_numpy()
        c = out[f"c_{stratum}"].to_numpy()
        d = out[f"d_{stratum}"].to_numpy()
        out[f"or_{stratum}"] = haldane_odds_ratio(a, b, c, d)
        out[f"fisher_p_{stratum}"] = fisher_p_vectorized(a, b, c, d)
        out[f"fisher_log10p_{stratum}"] = neglog10(out[f"fisher_p_{stratum}"].to_numpy())

    out["mcnemar_p"] = mcnemar_p_vectorized(
        out["n_only_asd"].to_numpy(), out["n_only_ctrl"].to_numpy()
    )
    out["mcnemar_log10p"] = neglog10(out["mcnemar_p"].to_numpy())
    return out[
        [
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "key",
            "a_all",
            "b_all",
            "c_all",
            "d_all",
            "or_all",
            "fisher_p_all",
            "fisher_log10p_all",
            "a_male",
            "b_male",
            "c_male",
            "d_male",
            "or_male",
            "fisher_p_male",
            "fisher_log10p_male",
            "a_female",
            "b_female",
            "c_female",
            "d_female",
            "or_female",
            "fisher_p_female",
            "fisher_log10p_female",
            "n_only_asd",
            "n_only_ctrl",
            "mcnemar_p",
            "mcnemar_log10p",
        ]
    ]
