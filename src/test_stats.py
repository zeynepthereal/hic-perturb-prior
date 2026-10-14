import numpy as np
import pandas as pd

from .stats import paired_table, signflip_test


def _cases(effect, n_pert=300, n_split=5, seed=0):
    rng = np.random.RandomState(seed)
    rows = []
    base = rng.randn(n_pert) * 0.5                                   # perturbation-level shared component (induces dependence)
    for s in range(n_split):
        te = rng.choice(n_pert, 80, replace=False)
        for i in te:
            rows.append((i, f"seed{s}", effect + base[i] + 0.3 * rng.randn()))
    return pd.DataFrame(rows, columns=["pert_i", "split", "d"])


def test_signflip_null_and_effect():
    null_p = [signflip_test(_cases(0.0, seed=k), n_perm=5000, n_boot=500)["p_signflip"] for k in range(20)]
    assert np.mean(np.array(null_p) < 0.05) <= 0.15                  # roughly calibrated under H0 (20 reps, loose bound)
    r = signflip_test(_cases(0.25), n_perm=5000, n_boot=500)
    assert r["p_signflip"] < 0.01 and r["D"] > 0 and r["ci95"][0] > 0


def test_aggregation_formula_and_paired_table():
    d = pd.DataFrame({"pert_i": [1, 2, 1, 3], "split": ["a", "a", "b", "b"], "d": [1.0, 3.0, 5.0, 7.0]})
    r = signflip_test(d, n_perm=2000, n_boot=100)
    assert abs(r["D"] - np.mean([np.mean([1.0, 3.0]), np.mean([5.0, 7.0])])) < 1e-12
    pp = pd.DataFrame({"pert_i": [1, 1], "split": ["a", "a"], "model": ["ctl", "hic"], "eval_set": "common", "target": "excluded",
                       "subset": "top1000_expr", "l2": [4.0, 3.5], "pearson": [0.1, 0.2]})
    assert paired_table(pp, "ctl", "hic")["d"].iloc[0] == 0.5                         # L2 lower = better -> positive
    assert paired_table(pp, "ctl", "hic", metric="pearson")["d"].iloc[0] == 0.1       # correlation higher = better -> positive
