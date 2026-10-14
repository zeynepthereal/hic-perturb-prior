"""Paired comparison of two models over (perturbation, split) test cases -- the analysis code for the pre-registered primary test.
Written and tested BEFORE any Hi-C model is run (see PREREGISTRATION.md).

Per (perturbation i, split s) paired difference  d_is = metric(reference) - metric(new)  for L2 (positive = `new` better);
for correlation metrics the sign is reversed so that positive always means `new` is better.
Aggregate effect (equal weight per split):  D = (1/S) sum_s mean_{i in test_s} d_is = sum_i a_i,   a_i = (1/S) sum_{s: i in test_s} d_is / n_s .
Test: two-sided sign-flip permutation test. All d_is of one perturbation get the SAME random sign (a perturbation appears in several
splits, so cases of one perturbation are not independent); under H0 (symmetric around 0) D -> sum_i eps_i a_i.
CI: percentile cluster bootstrap over perturbations."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps


def paired_table(per_pert: pd.DataFrame, ref: str, new: str, *, eval_set="common", target="excluded", subset="top1000_expr", metric="l2") -> pd.DataFrame:
    """per_pert = per_perturbation_metrics.csv of a run. Returns columns pert_i, split, d (positive = `new` better)."""
    q = per_pert[(per_pert.eval_set == eval_set) & (per_pert.target == target) & (per_pert["subset"] == subset)]
    a = q[q.model == ref].set_index(["pert_i", "split"])[metric]
    b = q[q.model == new].set_index(["pert_i", "split"])[metric]
    j = pd.concat([a.rename("ref"), b.rename("new")], axis=1, join="inner").dropna()
    sign = 1.0 if metric == "l2" else -1.0
    out = (sign * (j["ref"] - j["new"])).rename("d").reset_index()
    return out


def signflip_test(d: pd.DataFrame, n_perm: int = 100_000, seed: int = 0, n_boot: int = 10_000) -> dict:
    d = d.copy()
    n_s = d.groupby("split").pert_i.transform("count")
    S = d.split.nunique()
    d["w"] = d["d"] / n_s / S
    a = d.groupby("pert_i").w.sum()                       # a_i
    D = float(a.sum())
    rng = np.random.RandomState(seed)
    av = a.values
    cnt, chunk = 0, 5000
    for _ in range(0, n_perm, chunk):
        signs = rng.choice([-1.0, 1.0], size=(chunk, len(av)))
        cnt += int((np.abs(signs @ av) >= abs(D) - 1e-15).sum())
    p = (cnt + 1) / (n_perm + 1)
    boots = np.array([av[rng.randint(len(av), size=len(av))].sum() for _ in range(n_boot)])
    per_split = d.groupby("split")["d"].agg(["mean", "count"])
    tt = {s: sps.ttest_1samp(g["d"], 0.0).pvalue for s, g in d.groupby("split")}
    wx = {s: (sps.wilcoxon(g["d"]).pvalue if (g["d"] != 0).any() else 1.0) for s, g in d.groupby("split")}
    return dict(D=D, p_signflip=p, ci95=(float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))), n_unique_perts=len(av),
                n_cases=len(d), n_splits=S, per_split_mean=per_split["mean"].to_dict(), per_split_n=per_split["count"].to_dict(),
                per_split_paired_t_p=tt, per_split_wilcoxon_p=wx, frac_perts_new_better=float((av > 0).mean()))
