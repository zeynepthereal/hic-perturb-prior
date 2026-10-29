"""Study 2 variance summary: GEARS (3 splits x 3 training seeds = 9 runs) vs the deterministic baselines (3 splits), median and spread, all three metrics side by side.
python -m scripts.study2_aggregate [--splits seed1 seed2 seed3] [--seeds 1 2 3]     (scores any finished run that has no summary yet; run in env hic-perturb-prior)
GEARS rows: median over the 9 runs, [min, max], IQR, and the range of the 3 per-split medians (between-split) vs the mean within-split range (between-seed).
Baseline rows (mean, no-change, linear PCA; fitted on each split's training list, deterministic): median over the 3 splits and [min, max].
Shuffled controls = the same predictions permuted across the evaluated perturbations (mean of 20 permutations) -> the floor of each metric. NO rescaling of any prediction is applied anywhere."""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(); ap.add_argument("--splits", nargs="+", default=["seed1", "seed2", "seed3"]); ap.add_argument("--seeds", nargs="+", default=["1", "2", "3"]); a = ap.parse_args()
rows, amps, missing = [], [], []
for sp in a.splits:
    for sd in a.seeds:
        label = f"gears_{sp}_s{sd}"; d = Path(f"results/study2/gears_k562_{label}")
        if not (d / "predictions_test.npz").exists():
            missing.append(label); continue
        if not (d / "gears_eval_amplitude.csv").exists():
            subprocess.run([sys.executable, "-m", "scripts.study2_gears_eval", "--label", label, "--split", sp], check=True, stdout=subprocess.DEVNULL)
        rows.append(pd.read_csv(d / "gears_eval_summary.csv").assign(split=sp, seed=sd)); amps.append(pd.read_csv(d / "gears_eval_amplitude.csv").assign(split=sp, seed=sd))
if not rows:
    sys.exit("no finished runs")
S = pd.concat(rows); A = pd.concat(amps)
S = S.merge(A, on=["set", "target", "model", "subset", "split", "seed"], how="left")
S = S[S["subset"] == "top1000_expr"]
metrics = ["L2", "Pearson", "Pearson_centred", "mean_delta_ratio", "pert_norm_ratio", "collapse_pred", "collapse_obs"]


def fmt(x):
    x = x.dropna()
    return "n.a." if len(x) == 0 else f"{x.median():.4f} [{x.min():.4f}, {x.max():.4f}]"


out = [f"# Study 2: variance summary (finished runs: {S[['split', 'seed']].drop_duplicates().shape[0]} of {len(a.splits) * len(a.seeds)}; missing: {missing or 'none'})\n",
       "Top-1000 control-expressed genes; cells are `median [min, max]` over runs (GEARS: up to 9 runs; deterministic baselines: one value per split, up to 3). Pearson_centred = Pearson after subtracting the mean perturbation effect across the evaluated test perturbations (NaN for constant predictors). Amplitude ratios: >1 overshoot, <1 undershoot.\n"]
for (sset, tgt), g in S.groupby(["set", "target"]):
    out.append(f"\n## {sset} - target {tgt}\n\n| model | n runs | L2 | Pearson | Pearson (centred) | amp. ratio, mean delta | amp. ratio, per-pert norm | collapse(pred) | collapse(obs) |\n|---|---|---|---|---|---|---|---|---|")
    for m, gm in g.groupby("model"):
        gm = gm.drop_duplicates(["split", "seed"]) if not m.startswith("GEARS") else gm            # deterministic models: identical across seeds -> one row per split
        out.append(f"| {m} | {len(gm)} | " + " | ".join(fmt(gm[c]) for c in metrics) + " |")
G = S[(S.model == "GEARS") & S.set.str.startswith("A") & (S.target == "included")]
out.append("\n## GEARS spread decomposition (set A, target included)\n\n| metric | IQR over runs | between-split range (of per-split medians) | mean within-split range (between seeds) |\n|---|---|---|---|")
for c in metrics:
    ps = G.groupby("split")[c]
    out.append(f"| {c} | {G[c].quantile(.75) - G[c].quantile(.25):.4f} | {ps.median().max() - ps.median().min():.4f} | {(ps.max() - ps.min()).mean():.4f} |")
if "mean_delta_ratio" in G:
    out.append(f"\nCollapse index = ||mean_p X_p|| / mean_p ||X_p|| (1 = all perturbations along one shared direction); GEARS collapse(pred) > collapse(obs) in {int((G.collapse_pred > G.collapse_obs).sum())} of {len(G)} runs.")
    out.append(f"\nGEARS mean-delta amplitude ratio > 1 in {int((G.mean_delta_ratio > 1).sum())} of {len(G)} runs; per-perturbation norm ratio < 1 in {int((G.pert_norm_ratio < 1).sum())} of {len(G)}.")
Path("results/study2").mkdir(exist_ok=True)
open("results/study2/variance_summary.md", "w").write("\n".join(out) + "\n"); S.to_csv("results/study2/variance_all_runs.csv", index=False)
print("\n".join(out))
