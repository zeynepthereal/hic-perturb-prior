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


def test_tost_equivalence_and_decision_rule():
    from .stats import classify, tost_signflip
    big = _cases(0.0, n_pert=600, seed=3)
    big["d"] = big["d"] * 0.2                                         # precise, centred at 0
    tt = tost_signflip(big, delta=0.1, n_perm=5000)
    assert tt["equivalent"] and tt["p_tost"] < 0.05
    noisy = _cases(0.0, n_pert=90, n_split=2, seed=4)
    assert not tost_signflip(noisy, delta=0.02, n_perm=5000)["equivalent"]
    eff = _cases(0.4, seed=5)
    assert not tost_signflip(eff, delta=0.1, n_perm=5000)["equivalent"]
    assert classify(0.001, 0.3, 0.9) == "SIGNIFICANT: Hi-C better" and classify(0.001, -0.3, 0.9) == "SIGNIFICANT: Hi-C worse"
    assert classify(0.5, 0.0, 0.01) == "EQUIVALENT" and classify(0.5, 0.0, 0.4) == "INCONCLUSIVE"
    assert classify(0.01, 0.01, 0.01).startswith("SIGNIFICANT but within")


def test_control_averaged_over_models_and_primary_analysis():
    from .stats import primary_analysis
    rng = np.random.RandomState(0)
    rows = []
    for s in range(3):
        for i in rng.choice(100, 40, replace=False):
            for m, off in (("hic", 0.0), ("c0", 0.05), ("c1", 0.07), ("c2", 0.06)):
                rows.append((i, f"seed{s}", m, "common", "excluded", "top1000_expr", 4.0 + off + 0.01 * rng.randn(), 0.3))
    pp = pd.DataFrame(rows, columns=["pert_i", "split", "model", "eval_set", "target", "subset", "l2", "pearson"])
    d = paired_table(pp, ["c0", "c1", "c2"], "hic")
    assert abs(d["d"].mean() - 0.06) < 0.01                                          # control = average of the three control models
    r = primary_analysis(pp, "hic", ["c0", "c1", "c2"], n_perm=2000)
    assert r["decision"] == "SIGNIFICANT: Hi-C better" and abs(r["control_mean"] - 4.06) < 0.02
