import numpy as np
from scipy.stats import binomtest, fisher_exact

from pheno.stats import add_stats, fisher_p_vectorized, mcnemar_p_vectorized, neglog10


def test_mcnemar_matches_binomtest():
    n_asd = np.array([0, 1, 5, 3, 0])
    n_ctrl = np.array([0, 0, 0, 1, 4])
    got = mcnemar_p_vectorized(n_asd, n_ctrl)
    assert got[0] == 1.0
    for i in range(1, len(n_asd)):
        n = int(n_asd[i] + n_ctrl[i])
        expected = binomtest(int(n_asd[i]), n, 0.5, alternative="two-sided").pvalue
        assert abs(got[i] - expected) < 1e-12


def test_fisher_matches_scipy():
    a, b, c, d = np.array([1]), np.array([0]), np.array([3]), np.array([3])
    got = fisher_p_vectorized(a, b, c, d)[0]
    expected = fisher_exact([[1, 0], [3, 3]], alternative="two-sided")[1]
    assert abs(got - expected) < 1e-12


def test_neglog10():
    vals = neglog10(np.array([1.0, 0.01, 0.0]))
    assert vals[0] == 0.0
    assert abs(vals[1] - 2.0) < 1e-12
    assert np.isinf(vals[2])


def test_add_stats_column_order():
    import pandas as pd

    frame = pd.DataFrame(
        {
            "CHROM": ["chr21"],
            "POS": ["1"],
            "REF": ["C"],
            "ALT": ["T"],
            "key": ["chr21_1_C_T"],
            "a_all": [1],
            "b_all": [0],
            "c_all": [3],
            "d_all": [3],
            "a_male": [1],
            "b_male": [0],
            "c_male": [2],
            "d_male": [1],
            "a_female": [0],
            "b_female": [0],
            "c_female": [1],
            "d_female": [2],
            "n_only_asd": [1],
            "n_only_ctrl": [0],
        }
    )
    out = add_stats(frame)
    assert out.columns[0] == "CHROM"
    assert "mcnemar_p" in out.columns
    assert "fisher_p_all" in out.columns
    assert out.loc[0, "mcnemar_p"] == 1.0
