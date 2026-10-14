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


def _per_case(q: pd.DataFrame, models, metric: str) -> pd.Series:
    """metric of one model, or the AVERAGE OF THE PER-CASE METRIC VALUES of several models (e.g. the 10 random-trans control seeds),
    per (pert_i, split). This is the mean of ten L2 values, each computed from that seed's OWN prediction; it is deliberately NOT the L2
    of an averaged (ensembled) prediction, which would be lower by Jensen's inequality and would hand the control an ensemble advantage.
    Cases for which not every listed model has a value are dropped."""
    models = [models] if isinstance(models, str) else list(models)
    cols = [q[q.model == m].set_index(["pert_i", "split"])[metric].rename(m) for m in models]
    return pd.concat(cols, axis=1, join="inner").dropna().mean(axis=1)


def paired_table(per_pert: pd.DataFrame, ref, new, *, eval_set="common", target="excluded", subset="top1000_expr", metric="l2") -> pd.DataFrame:
    """per_pert = per_perturbation_metrics.csv of a run. `ref`/`new` = model name or list of names (averaged per case).
    Returns columns pert_i, split, d (positive = `new` better) and ref_value (the reference/control metric, for the SESOI)."""
    q = per_pert[(per_pert.eval_set == eval_set) & (per_pert.target == target) & (per_pert["subset"] == subset)]
    a, b = _per_case(q, ref, metric), _per_case(q, new, metric)
    j = pd.concat([a.rename("ref"), b.rename("new")], axis=1, join="inner").dropna()
    sign = 1.0 if metric == "l2" else -1.0
    out = pd.DataFrame({"d": sign * (j["ref"] - j["new"]), "ref_value": j["ref"]}).reset_index()
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


def _one_sided_p(d: pd.DataFrame, shift: float, upper: bool, n_perm: int, seed: int) -> float:
    """Sign-flip p-value for H0: true mean = `shift`... implemented by shifting the cases so the boundary null is centred at 0."""
    dd = d.copy(); dd["d"] = dd["d"] - shift
    n_s = dd.groupby("split").pert_i.transform("count"); S = dd.split.nunique()
    av = (dd["d"] / n_s / S).groupby(dd.pert_i).sum().values
    D = av.sum(); rng = np.random.RandomState(seed); cnt = 0
    for _ in range(0, n_perm, 5000):
        sd = rng.choice([-1.0, 1.0], size=(5000, len(av))) @ av
        cnt += int((sd <= D + 1e-15).sum()) if upper else int((sd >= D - 1e-15).sum())
    return (cnt + 1) / (n_perm + 1)


def tost_signflip(d: pd.DataFrame, delta: float, n_perm: int = 100_000, seed: int = 1, alpha: float = 0.05) -> dict:
    """Two one-sided tests for equivalence of the aggregate effect D to 0 within (-delta, +delta), with the same perturbation-level
    sign-flip machinery. Lower test: H0 D <= -delta; upper test: H0 D >= +delta. TOST p = max of the two; equivalent if < alpha
    (equivalently the 1-2*alpha bootstrap/permutation interval lies inside the bounds)."""
    p_lo = _one_sided_p(d, -delta, upper=False, n_perm=n_perm, seed=seed)
    p_hi = _one_sided_p(d, +delta, upper=True, n_perm=n_perm, seed=seed + 1)
    return dict(delta=delta, p_lower=p_lo, p_upper=p_hi, p_tost=max(p_lo, p_hi), equivalent=bool(max(p_lo, p_hi) < alpha))


def sesoi(d: pd.DataFrame, frac: float = 0.005) -> float:
    """Smallest effect size of interest = `frac` of the control's mean metric over the primary test cases (pre-registered: 0.5 %)."""
    return float(frac * d["ref_value"].mean())


def classify(p_primary: float, D: float, p_tost: float, alpha: float = 0.05) -> str:
    """Pre-registered decision rule."""
    sig, eq = p_primary < alpha, p_tost < alpha
    if sig and eq:
        return "SIGNIFICANT but within the equivalence bounds (negligible effect)"
    if sig:
        return "SIGNIFICANT: Hi-C better" if D > 0 else "SIGNIFICANT: Hi-C worse"
    return "EQUIVALENT" if eq else "INCONCLUSIVE"


def primary_analysis(per_pert: pd.DataFrame, new: str, control_models, *, n_perm: int = 100_000, seed: int = 0, **kw) -> dict:
    """The complete pre-registered primary analysis: paired sign-flip test + TOST at the SESOI + decision."""
    d = paired_table(per_pert, control_models, new, **kw)
    res = signflip_test(d, n_perm=n_perm, seed=seed)
    delta = sesoi(d)
    tost = tost_signflip(d, delta, n_perm=n_perm, seed=seed + 100)
    ranks = seed_rank_summary(per_pert, new, control_models, **kw)
    return {**res, **tost, "control_mean": float(d["ref_value"].mean()), "decision": classify(res["p_signflip"], res["D"], tost["p_tost"]),
            "seed_ranks_descriptive": ranks}


def seed_rank_summary(per_pert: pd.DataFrame, new: str, control_models, *, eval_set="common", target="excluded", subset="top1000_expr",
                      metric: str = "l2") -> dict:
    """DESCRIPTIVE companion to the primary test (no p-value): rank of the Hi-C arm's mean test metric among the control seeds.
    Per split: mean over that split's test cases of the metric for the Hi-C arm and for each control seed; rank 1 = best (lowest L2 /
    highest correlation) among the 1 + n_seeds values, and the number of control seeds beaten strictly. Pooled: mean of the split means
    (equal weight per split, the same aggregation as the primary effect D). Cases need a value for the Hi-C arm and every control seed."""
    controls = [control_models] if isinstance(control_models, str) else list(control_models)
    q = per_pert[(per_pert.eval_set == eval_set) & (per_pert.target == target) & (per_pert["subset"] == subset)]
    cols = {m: q[q.model == m].set_index(["pert_i", "split"])[metric] for m in [new] + controls}
    T = pd.concat(cols, axis=1, join="inner").dropna()
    lower_better = metric == "l2"
    per_split = T.groupby(level="split").mean()
    pooled = per_split.mean()

    def rank(v):
        c = v[controls]
        beaten = int((c > v[new]).sum()) if lower_better else int((c < v[new]).sum())
        return dict(hic=float(v[new]), control_mean=float(c.mean()), control_min=float(c.min()), control_max=float(c.max()),
                    rank=1 + (len(controls) - beaten - int((c == v[new]).sum())), n_controls_beaten=beaten, n_controls=len(controls))
    return dict(per_split={s: rank(per_split.loc[s]) for s in per_split.index}, pooled=rank(pooled), n_cases=int(len(T)))
