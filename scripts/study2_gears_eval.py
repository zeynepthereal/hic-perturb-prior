"""STUDY 2 / Phase 1: score GEARS predictions on the Study 1 metric. python -m scripts.study2_gears_eval --label gears_seed1
Metric (as Study 1 / the paper): per-perturbation L2 between predicted and observed delta (perturbation - control pseudobulk) and Pearson delta, on the 1,000 genes most expressed in
control cells; 'included' = target gene kept (the paper's convention), 'excluded' = target removed first (Study 1 primary convention). Baselines are fitted on the training list of the split
exactly as in Study 1 (mean baseline = mean training delta; linear = PCA, K=10, lambda=0.1) so the comparison is on identical test cases.
Sets: A = test perturbations GEARS can predict (GO graph) AND the PCA linear model can predict (target is a read-out gene) -> the paper's `common` convention;
      B = all test perturbations GEARS can predict. NO Hi-C of any kind enters this script."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from src.baselines import EMBEDDINGS
from src.data import PerturbationDataset, split_from_conditions
from src.evaluate import fit_predict_linear, gene_masks, per_pert_metrics
from src.study2_metrics import amplitude_ratios, centred_pearson

N_SHUFFLE = 20                                 # random permutations averaged for the shuffled control
PAPER = {"L2": 3.991, "Pearson": 0.387}          # Ahlmann-Eltze et al. 2025, GEARS, Replogle K562 (mean over their 2 splits), top-1000 expressed genes, target included


def score(ds, split, gears_delta: dict, n_top_expr=1000, n_top_de=20):
    """gears_delta: {ds perturbation index -> predicted delta (n_gene)} for the test perturbations GEARS can predict. Returns (per-perturbation metrics for sets A and B, amplitude-ratio table); shuffled-perturbation controls are extra model rows."""
    tr, te = split["train"], split["test"]
    g_ok = np.array([i in gears_delta for i in te])
    pca_ok = EMBEDDINGS["train_pca_expr"].predictable(ds)[te]
    sets = {"A: GEARS & PCA predictable": g_ok & pca_ok, "B: GEARS predictable": g_ok}
    mean_pred = np.tile(ds.delta[tr].mean(0), (len(te), 1))
    lin_pred = fit_predict_linear(ds, tr, te, 10, 0.1, "train_pca_expr")
    gears_pred = np.stack([gears_delta.get(i, np.full(ds.delta.shape[1], np.nan)) for i in te])
    models = {"GEARS": gears_pred, "mean": mean_pred, "no_change": np.zeros_like(mean_pred), "linear[PCA] (K=10, lam=0.1)": lin_pred}
    out, amp = [], []
    rng = np.random.RandomState(0)
    for sname, sel in sets.items():
        for tmode, excl in (("included", False), ("excluded", True)):
            masks = gene_masks(ds.delta[te][sel], ds.ctrl_mean, ds.target_idx[te][sel], excl, n_top_de, n_top_expr)
            for mname, P in models.items():
                if mname.startswith("linear") and not sname.startswith("A"):
                    continue                                               # the PCA model cannot predict targets outside the read-out set
                d = per_pert_metrics(P[sel], ds.delta[te][sel], masks); d["pert_i"] = te[sel][d["pert_i"].values]
                d["pearson_centred"] = np.nan                                  # third metric: Pearson after subtracting the mean perturbation effect across the evaluated test perturbations
                for sub, m_ in masks.items():
                    d.loc[d["subset"] == sub, "pearson_centred"] = centred_pearson(P[sel], ds.delta[te][sel], m_)
                out.append(d.assign(model=mname, set=sname, target=tmode))
                for sub, m_ in masks.items():
                    if mname != "no_change":
                        amp.append({"set": sname, "target": tmode, "model": mname, "subset": sub, **amplitude_ratios(P[sel], ds.delta[te][sel], m_)})
                if mname == "GEARS" or mname.startswith("linear"):                 # shuffled-perturbation control: the same predictions permuted across the evaluated perturbations (floor)
                    reps = []
                    for r in range(N_SHUFFLE):
                        Ps = P[sel][rng.permutation(int(sel.sum()))]
                        dd = per_pert_metrics(Ps, ds.delta[te][sel], masks); dd["pert_i"] = te[sel][dd["pert_i"].values]; dd["pearson_centred"] = np.nan
                        for sub, m_ in masks.items():
                            dd.loc[dd["subset"] == sub, "pearson_centred"] = centred_pearson(Ps, ds.delta[te][sel], m_)
                        reps.append(dd)
                    dd = pd.concat(reps).groupby(["pert_i", "subset"], as_index=False)[["l2", "pearson", "pearson_centred"]].mean()      # mean over the permutations, per perturbation
                    out.append(dd.assign(model=f"{mname} [shuffled across perts]", set=sname, target=tmode))
    return pd.concat(out), pd.DataFrame(amp)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--label", default="gears_seed1"); ap.add_argument("--dir", default=None); ap.add_argument("--split", default=None, help="split name, e.g. seed1 (default: derived from --label)"); a = ap.parse_args()
    d = Path(a.dir or f"results/study2/gears_k562_{a.label}")
    ds = PerturbationDataset.load("data/processed/k562_gears.npz")
    split = split_from_conditions(ds, json.load(open(f"configs/splits/gears_k562_{a.split or a.label.replace('gears_', '')}.json")))
    z = np.load(d / "predictions_test.npz", allow_pickle=False)
    assert list(z["gene_names"]) == list(ds.genes), "gene order of the GEARS predictions differs from the dataset"
    pos = {p: i for i, p in enumerate(ds.perts)}
    gd = {pos[c.replace("+ctrl", "")]: z["pred_expression"][k].astype(np.float64) - ds.ctrl_mean for k, c in enumerate(z["conditions"])}
    M, AMP = score(ds, split, gd); AMP.to_csv(d / "gears_eval_amplitude.csv", index=False); M.to_csv(d / "gears_eval_per_perturbation.csv", index=False)
    S = M.groupby(["set", "target", "model", "subset"]).agg(L2=("l2", "mean"), Pearson=("pearson", "mean"), Pearson_centred=("pearson_centred", "mean"), n=("pert_i", "nunique")).reset_index()
    S.to_csv(d / "gears_eval_summary.csv", index=False)
    print(S[S["subset"] == "top1000_expr"].round(4).to_string(index=False))
    like = S[(S["set"].str.startswith("A")) & (S.target == "included") & (S["subset"] == "top1000_expr") & (S.model == "GEARS")].iloc[0]
    json.dump({"label": a.label, "gears_like_for_like (set A, target included, top-1000)": {"L2": float(like.L2), "Pearson": float(like.Pearson), "Pearson_centred": float(like.Pearson_centred), "n": int(like.n)}, "paper_gears": PAPER}, open(d / "gears_eval_headline.json", "w"), indent=1)


if __name__ == "__main__":
    main()
