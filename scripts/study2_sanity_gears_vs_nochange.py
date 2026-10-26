"""Study 2 sanity check: is 'GEARS L2 > no-change L2' (set A, target included, top-1000) real, or a mean-centering / scaling mismatch between the GEARS prediction and the control used for L2?
python -m scripts.study2_sanity_gears_vs_nochange --label gears_seed1        (reads only the saved predictions and the data; no training)"""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np

from src.data import PerturbationDataset, split_from_conditions
from src.evaluate import gene_masks

ap = argparse.ArgumentParser(); ap.add_argument("--label", default="gears_seed1"); a = ap.parse_args()
d = Path(f"results/study2/gears_k562_{a.label}")
ds = PerturbationDataset.load("data/processed/k562_gears.npz")
split = split_from_conditions(ds, json.load(open(f"configs/splits/gears_k562_{a.label.replace('gears_', '')}.json")))
z = np.load(d / "predictions_test.npz")
out = {}
# 1) control expression: recompute from the GEARS h5ad exactly as GEARS does (mean of the ctrl cells of adata.X) and compare with the control used for the deltas
h = ad.read_h5ad("data/gears/replogle_k562_essential/perturb_processed.h5ad", backed="r")
ctrl_rows = np.where((h.obs.condition == "ctrl").values)[0]
ctrl_gears = np.asarray(h.X[ctrl_rows].mean(0)).ravel()
out["ctrl_mean_max_abs_diff_vs_dataset (GEARS-style ctrl mean vs the control used for L2)"] = float(np.abs(ctrl_gears - ds.ctrl_mean).max())
out["gene_order_identical"] = bool(list(z["gene_names"]) == list(ds.genes))
out["corr(mean predicted absolute expression over test perts, ctrl mean)"] = float(np.corrcoef(z["pred_expression"].mean(0), ctrl_gears)[0, 1])
# 2) set A (GEARS & PCA predictable), top-1000, target included
pos = {p: i for i, p in enumerate(ds.perts)}
idx = np.array([pos[c.replace("+ctrl", "")] for c in z["conditions"]])
P = z["pred_expression"].astype(np.float64) - ds.ctrl_mean                        # predicted delta
sel = ds.target_idx[idx] >= 0
O = ds.delta[idx]; P, O, idx = P[sel], O[sel], idx[sel]
M = gene_masks(O, ds.ctrl_mean, ds.target_idx[idx], False, 20, 1000)["top1000_expr"]
Pm, Om = P * M, O * M
l2 = lambda X_: np.sqrt(((X_ - Om) ** 2).sum(1))
no_change = np.sqrt((Om ** 2).sum(1)); g = l2(Pm)
out["n_set_A"] = int(len(idx))
out["mean L2: GEARS / no_change / mean_baseline(train)"] = [float(g.mean()), float(no_change.mean()), float(l2(np.tile(ds.delta[split['train']].mean(0), (len(idx), 1)) * M).mean())]
out["fraction of perturbations where GEARS L2 > no-change L2"] = float((g > no_change).mean())
nP, nO = np.sqrt((Pm ** 2).sum(1)), np.sqrt((Om ** 2).sum(1))
cos = (Pm * Om).sum(1) / (nP * nO)
out["mean ||pred delta||, mean ||obs delta||, ratio"] = [float(nP.mean()), float(nO.mean()), float(nP.mean() / nO.mean())]
out["mean cosine(pred delta, obs delta) (uncentred)"] = float(cos.mean())
alpha = float((Pm * Om).sum() / (Pm * Pm).sum())                                  # global amplitude that minimises the summed squared error
out["optimal global scale alpha* (pred*alpha*)"] = alpha
out["mean L2 with pred*alpha*"] = float(l2(alpha * Pm).mean())
out["no-change beats GEARS iff ||p||^2 > 2<p,o>: fraction of perturbations satisfying it"] = float(((Pm ** 2).sum(1) > 2 * (Pm * Om).sum(1)).mean())
# 3) mean-centering: mean predicted delta across test perturbations vs the mean observed delta and the training mean delta
pbar, obar, tbar = Pm.mean(0), Om.mean(0), ds.delta[split["train"]].mean(0) * M[0]
out["||mean pred delta||, ||mean obs delta||, ||train mean delta||"] = [float(np.linalg.norm(pbar)), float(np.linalg.norm(obar)), float(np.linalg.norm(tbar))]
out["cosine(mean pred, mean obs), cosine(mean pred, train mean)"] = [float(pbar @ obar / np.linalg.norm(pbar) / np.linalg.norm(obar)), float(pbar @ tbar / np.linalg.norm(pbar) / np.linalg.norm(tbar))]
out["mean L2 if the mean predicted delta is replaced by the observed mean"] = float(l2(Pm - pbar + obar).mean())
out["mean L2 of the centred parts only: (pred - mean pred) vs (obs - mean obs)"] = float(np.sqrt((((Pm - pbar) - (Om - obar)) ** 2).sum(1)).mean())
out["mean L2 of the centred observed deltas alone (no-change after removing the mean response)"] = float(np.sqrt(((Om - obar) ** 2).sum(1)).mean())
json.dump(out, open(d / "sanity_gears_vs_nochange.json", "w"), indent=1)
for k, v in out.items():
    print(f"{k}: {v}")
