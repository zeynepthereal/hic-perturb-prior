"""python -m scripts.run_s4_vs_mean   -> results/<ts>_prereg_s4_vs_mean/{s4_results.json,s4_results.md}
PRE-REGISTERED S4 (and the S4/S9 PCA reference): each arm vs the MEAN BASELINE on the primary test cases (common set, GEARS splits 1-5, L2 on the top-1000
expressed genes, target gene excluded), with the same paired sign-flip test and the same D definition as the primary but with the reference = mean baseline:
    D = (1/5) sum_s mean_{i in test_s} [ L2_mean(i,s) - L2_arm(i,s) ]      (positive = the arm is better than the mean baseline).
Arms: b500_hic_trans; b500_lin_randtrans (mean of the 10 per-seed case-level L2 values, as in the primary); the 1,000 random-Gaussian arms (summary);
at K=10 (paper hp, lambda 0.1) and K=50 (fixed, lambda 0.1). Reference: the PCA embedding vs the mean baseline on its own predictable subset.
All tests are `src.stats.signflip_test` (frozen). The 1,000 random arms use 10,000 permutations / 1,000 bootstrap draws per arm (the named arms use the frozen 100,000 / 10,000)."""
import glob
import json

import numpy as np
import pandas as pd

from src.data import PerturbationDataset, split_from_conditions
from src.evaluate import fixed_name, lm_name
from src.hic.embeddings import GraphSpectralEmbedding
from src.nulls import NullEvaluator
from src.stats import paired_table, signflip_test
from src.utils import make_run_dir

PRIMARY = sorted(glob.glob("results/*_prereg_primary"))[-1]
PCA_RUN = sorted(glob.glob("results/*_sanity1_pca_only"))[-1]
KS = (10, 50)
name = lambda emb, K: lm_name(emb, False) if K == 10 else fixed_name(emb, 50, 0.1)


def wmean(d):                                   # equal weight per split, as in D
    return d.groupby("split").mean().mean()


def summarize(d, n_perm=100_000, n_boot=10_000):
    r = signflip_test(d, n_perm=n_perm, n_boot=n_boot)
    ref_w = wmean(d.set_index(["pert_i", "split"])["ref_value"])
    return {"mean_L2_arm": float(ref_w - r["D"]), "mean_L2_mean_baseline": float(ref_w), "D_vs_mean_baseline": r["D"], "p_two_sided": r["p_signflip"],
            "ci95": r["ci95"], "n_cases": r["n_cases"], "n_unique_perts": r["n_unique_perts"], "frac_perts_arm_better": r["frac_perts_new_better"],
            "per_split_D": r["per_split_mean"]}


pp = pd.read_csv(f"{PRIMARY}/per_perturbation_metrics.csv")
ds = PerturbationDataset.load("data/processed/k562_gears.npz")
splits = {f"gears_seed{s}": split_from_conditions(ds, json.load(open(f"configs/splits/gears_k562_seed{s}.json"))) for s in range(1, 6)}
out = {"primary_run": PRIMARY, "pca_run": PCA_RUN, "K": {}}
for K in KS:
    res = {}
    res["b500_hic_trans"] = summarize(paired_table(pp, "mean", name("b500_hic_trans", K)))
    res["b500_lin_randtrans (10-seed mean of case-level L2)"] = summarize(paired_table(pp, "mean", [name(f"b500_lin_randtrans_s{k}", K) for k in range(10)]))
    base = pp[(pp.eval_set == "common") & (pp.target == "excluded") & (pp["subset"] == "top1000_expr") & (pp.model == "mean")].set_index(["pert_i", "split"])["l2"]
    ev = NullEvaluator(ds, splits, K=K, lam=0.1)
    rows = []
    for j in range(1000):
        arm = ev.case_l2(GraphSpectralEmbedding("b500_hic_trans", "random", j).perturbation_matrix(ds, K))
        j_ = pd.concat([base.rename("ref_value"), arm.rename("arm")], axis=1, join="inner").dropna()
        d = pd.DataFrame({"d": j_["ref_value"] - j_["arm"], "ref_value": j_["ref_value"]}).reset_index()
        r = signflip_test(d, n_perm=10_000, n_boot=1_000)
        rows.append({"seed": j, "mean_L2": float(wmean(j_["arm"])), "D": r["D"], "p": r["p_signflip"], "ci_lo": r["ci95"][0], "ci_hi": r["ci95"][1]})
    R = pd.DataFrame(rows)
    q = lambda c: {"min": float(R[c].min()), "q05": float(R[c].quantile(.05)), "median": float(R[c].median()), "q95": float(R[c].quantile(.95)), "max": float(R[c].max())}
    res["random_gaussian (1,000 arms; summary)"] = {"n_arms": 1000, "mean_L2": q("mean_L2"), "D_vs_mean_baseline": q("D"), "p_two_sided": q("p"),
                                                    "n_D_positive": int((R.D > 0).sum()), "n_p_lt_0.05_arm_better": int(((R.p < .05) & (R.D > 0)).sum()),
                                                    "n_p_lt_0.05_arm_worse": int(((R.p < .05) & (R.D < 0)).sum()), "median_ci95": [float(R.ci_lo.median()), float(R.ci_hi.median())],
                                                    "mean_L2_mean_baseline": float(wmean(base)), "permutations_per_arm": 10_000}
    out["K"][str(K)] = res
    R.to_csv(f"/tmp/s4_random_K{K}.csv", index=False) if False else None
# S4/S9 reference: PCA embedding vs the mean baseline on its own predictable subset (the `common` set of the PCA-only run)
pc = pd.read_csv(f"{PCA_RUN}/per_perturbation_metrics.csv")
out["pca_reference"] = {str(K): summarize(paired_table(pc, "mean", name("train_pca_expr", K))) for K in KS}
run_dir = make_run_dir({"analysis": "PRE-REGISTERED S4 / S9: arms vs the mean baseline", "primary_run": PRIMARY, "pca_run": PCA_RUN}, "prereg_s4_vs_mean")
(run_dir / "s4_results.json").write_text(json.dumps(out, indent=1, default=float))
print(run_dir)
