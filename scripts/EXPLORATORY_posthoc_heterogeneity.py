"""EXPLORATORY, POST HOC, NOT PRE-REGISTERED.  python -m scripts.EXPLORATORY_posthoc_heterogeneity
Descriptive look at the heterogeneity of the primary effect (b500_hic_trans vs b500_lin_randtrans, K=10): distribution of the per-perturbation contributions a_i
(D = sum_i a_i, a_i = (1/5) sum_s d_is / n_s, d = control L2 - Hi-C L2), concentration of D in the top perturbations, and the largest contributors in split 1.
NO tests and NO p-values."""
import glob
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.baselines import EMBEDDINGS
from src.data import PerturbationDataset, split_from_conditions
from src.evaluate import lm_name
from src.hic.embeddings import gene_tables
from src.stats import paired_table
from src.utils import make_run_dir

PRIMARY = sorted(glob.glob("results/*_prereg_primary"))[-1]
pp = pd.read_csv(f"{PRIMARY}/per_perturbation_metrics.csv")
ds = PerturbationDataset.load("data/processed/k562_gears.npz")
nodes, sym, _ = gene_tables(); node = nodes.set_index("gene_id")
d = paired_table(pp, [lm_name(f"b500_lin_randtrans_s{k}", False) for k in range(10)], lm_name("b500_hic_trans", False))   # d = control - Hi-C
hic = paired_table(pp, "mean", lm_name("b500_hic_trans", False)).set_index(["pert_i", "split"])["ref_value"]               # (only for the Hi-C L2 column)
S = d.split.nunique(); n_s = d.groupby("split").pert_i.transform("count"); d["w"] = d["d"] / n_s / S
a = d.groupby("pert_i").w.sum().sort_values(ascending=False); N = len(a); D = float(a.sum())
run_dir = make_run_dir({"analysis": "EXPLORATORY POST HOC (not pre-registered): heterogeneity of the primary effect", "primary_run": PRIMARY}, "EXPLORATORY_posthoc_primary_heterogeneity")
(run_dir / "README.md").write_text("# EXPLORATORY, POST HOC, NOT PRE-REGISTERED\n\nDescriptive heterogeneity of the primary effect, produced after the pre-registered results were known. No tests, no p-values. "
                                   "See DEVIATIONS.md (D7). Source: `scripts/EXPLORATORY_posthoc_heterogeneity.py`.\n")
out = {"D": D, "n_unique_perts": N, "sum_a_equals_D": bool(abs(a.sum() - D) < 1e-12)}
qs = [0, .01, .05, .10, .25, .5, .75, .90, .95, .99, 1.0]
out["a_i_quantiles"] = {f"{q:g}": float(a.quantile(q)) for q in qs}
out["N_times_a_i_quantiles (mean = D)"] = {f"{q:g}": float(N * a.quantile(q)) for q in qs}
out["frac_a_i_positive"] = float((a > 0).mean())
conc = {}
for pct in (5, 10):
    k = math.ceil(N * pct / 100)
    top, bot = float(a.iloc[:k].sum()), float(a.iloc[-k:].sum())
    ab = a.abs().sort_values(ascending=False)
    conc[f"top_{pct}pct (k={k})"] = {"sum_of_largest_a_i": top, "share_of_D": top / D, "sum_of_smallest_a_i": bot, "share_of_D_smallest": bot / D,
                                      "D_without_largest": D - top, "share_of_sum_abs_a_i_top_by_abs": float(ab.iloc[:k].sum() / ab.sum())}
out["concentration"] = conc
a.rename("a_i").reset_index().assign(target=lambda x: ds.perts[x.pert_i.values], N_times_a_i=lambda x: N * x.a_i).to_csv(run_dir / "a_i_per_perturbation.csv", index=False)
fig, ax = plt.subplots(figsize=(7.5, 3.8)); fig.patch.set_facecolor("#fcfcfb"); ax.set_facecolor("#fcfcfb")
v = a.values * 1e4; hi = np.quantile(v, .995); lo = np.quantile(v, .005)
ax.hist(np.clip(v, lo, hi), bins=60, color="#2a78d6", edgecolor="#fcfcfb", linewidth=0.5)
ax.axvline(0, color="#52514e", lw=1); ax.axvline(D / N * 1e4, color="#eb6834", lw=2)
ax.set_xlabel("a_i × 10⁴   (end bins clip values beyond the 0.5 / 99.5 % quantiles)"); ax.set_ylabel("perturbations")
ax.set_title(f"Per-perturbation contributions a_i (N = {N}); orange = mean of a_i = D/N", loc="left", fontsize=10, color="#52514e")
ax.grid(color="#e4e3df", lw=0.6); ax.set_axisbelow(True)
for sp_ in ("top", "right"): ax.spines[sp_].set_visible(False)
fig.tight_layout(); fig.savefig(run_dir / "a_i_histogram.png", dpi=150); plt.close(fig)
hist_counts, edges = np.histogram(np.clip(v, lo, hi), bins=60); out["histogram_counts_a_i_x1e4"] = {"edges": edges.tolist(), "counts": hist_counts.tolist()}

# ---- split 1: top contributors
okx = EMBEDDINGS["b500x_hic_trans"].predictable(ds)          # False = target bin lies in an S10 rearranged/amplified region (removed before balancing)
splits = {f"gears_seed{s}": split_from_conditions(ds, json.load(open(f"configs/splits/gears_k562_seed{s}.json"))) for s in range(1, 6)}
rank_all = {p: r + 1 for r, p in enumerate(a.index)}
per_split = {}
for lab in sorted(d.split.unique()):
    s = d[d.split == lab].sort_values("d", ascending=False); Ds = float(s["d"].mean()); ns = len(s)
    k5 = math.ceil(ns * .05)
    per_split[lab] = {"D_s": Ds, "n": ns, "top10_share_of_D_s": float(s["d"].iloc[:10].sum() / ns / Ds), "top5pct_k": k5, "top5pct_share_of_D_s": float(s["d"].iloc[:k5].sum() / ns / Ds),
                      "D_s_without_top10": float(s["d"].iloc[10:].sum() / ns), "n_cases_d_gt_0": int((s["d"] > 0).sum())}
out["per_split"] = per_split
s1 = d[d.split == "gears_seed1"].sort_values("d", ascending=False).reset_index(drop=True); n1 = len(s1); D1 = float(s1["d"].mean())
rows = []; cum = 0.0
for r in range(10):
    p = int(s1.pert_i[r]); gid = sym.get(ds.perts[p]); g = node.loc[gid]
    others = d[(d.pert_i == p) & (d.split != "gears_seed1")]
    cum += s1["d"][r] / n1
    rows.append({"rank_in_split1": r + 1, "pert_i": p, "target_gene": ds.perts[p], "chromosome": g.chr, "tss": int(g.tss), "target_in_5000_readout_genes": bool(ds.target_idx[p] >= 0),
                 "in_S10_rearranged_region": bool(not okx[p]), "d_split1 (control L2 - Hi-C L2)": float(s1["d"][r]), "contribution_to_D_split1 (d/n)": float(s1["d"][r] / n1),
                 "share_of_D_split1": float(s1["d"][r] / n1 / D1), "cumulative_share": float(cum / D1), "rank_of_a_i_overall": rank_all[p],
                 "n_other_splits_in_test": int(len(others)), "mean_d_in_other_splits": float(others["d"].mean()) if len(others) else None})
T = pd.DataFrame(rows); T.to_csv(run_dir / "split1_top10_contributors.csv", index=False)
out["split1_top10"] = rows
out["split1_D"] = D1
out["split1_top10_in_S10_region"] = int(T.in_S10_rearranged_region.sum())
out["split1_chromosomes_top10"] = T.chromosome.value_counts().to_dict()
(run_dir / "heterogeneity.json").write_text(json.dumps(out, indent=1, default=float))
print(run_dir)
