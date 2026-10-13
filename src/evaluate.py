"""Metrics, model selection on the validation split, and the multi-seed experiment loop."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .baselines import EMBEDDINGS, LinearModel, MeanBaseline, NoChange
from .data import PerturbationDataset, make_split, split_from_conditions
from .utils import set_seed

PAPER_LM = "linear (paper K=10, lam=0.1)"
TUNED_LM = "linear (val-tuned)"


def gene_masks(obs_delta, ctrl_mean, target_idx, exclude_target, n_de=20, n_expr=1000) -> dict:
    """Boolean (n_pert, n_gene) masks per gene subset.
    all:          every read-out gene
    top20_de:     per-perturbation largest |observed delta| (oracle, as the paper's Extended Data Fig. 8b)
    top1000_expr: the paper's headline subset, most expressed genes in control
    With exclude_target, each perturbation's own target gene is removed *before* the top-n selection."""
    n_pert, n_gene = obs_delta.shape
    allowed = np.ones((n_pert, n_gene), dtype=bool)
    if exclude_target:
        r = np.where(target_idx >= 0)[0]
        allowed[r, target_idx[r]] = False

    def top(score, n):
        m = np.zeros((n_pert, n_gene), dtype=bool)
        idx = np.argsort(-np.where(allowed, score, -np.inf), axis=1)[:, :min(n, n_gene)]
        np.put_along_axis(m, idx, True, axis=1)
        return m & allowed

    return {"all": allowed, "top20_de": top(np.abs(obs_delta), n_de),
            "top1000_expr": top(np.tile(ctrl_mean, (n_pert, 1)), n_expr)}


def _pearson(a, b):
    a, b = a - a.mean(), b - b.mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return np.nan if den < 1e-12 else float((a * b).sum() / den)  # NaN for constant predictions (no-change)


def per_pert_metrics(pred, obs, masks: dict) -> pd.DataFrame:
    """L2 and Pearson (predicted vs observed delta) per perturbation per gene subset."""
    rows = []
    for name, m in masks.items():
        for i in range(obs.shape[0]):
            p, o = pred[i, m[i]], obs[i, m[i]]
            rows.append((i, name, float(np.sqrt(((p - o) ** 2).sum())), _pearson(p, o)))
    return pd.DataFrame(rows, columns=["pert_i", "subset", "l2", "pearson"])


def fit_predict_linear(ds, train_idx, test_idx, K, lam, embedding):
    Y = ds.delta[train_idx]
    G, P_all = EMBEDDINGS[embedding](Y, train_idx, ds, K)
    return LinearModel(lam).fit(Y, G, P_all[train_idx]).predict(P_all[test_idx])


def select_hparams(ds, split, Ks, lams, embedding, exclude_target):
    """Pick (K, lambda) minimising mean L2 over all read-out genes on the validation perturbations."""
    obs = ds.delta[split["val"]]
    allowed = gene_masks(obs, ds.ctrl_mean, ds.target_idx[split["val"]], exclude_target)["all"]
    best, grid = None, []
    for K in Ks:
        for lam in lams:
            pred = fit_predict_linear(ds, split["train"], split["val"], K, lam, embedding)
            err = float(np.sqrt((((pred - obs) ** 2) * allowed).sum(axis=1)).mean())
            grid.append({"K": K, "lambda": lam, "val_l2": err})
            if best is None or err < best[0]:
                best = (err, K, lam)
    return best[1], best[2], pd.DataFrame(grid)


def run_split(ds: PerturbationDataset, cfg: dict, split: dict, label, seed: int):
    set_seed(seed)
    e = cfg["eval"]
    emb, excl = e["embedding"], e["exclude_target_primary"]
    tr, va, te = split["train"], split["val"], split["test"]
    K, lam, grid = select_hparams(ds, split, e["K_grid"], e["lambda_grid"], emb, excl)
    trva = np.concatenate([tr, va])   # models are refit on train+val before scoring on test
    obs = ds.delta[te]
    preds = {
        "no_change": NoChange().fit(ds.delta[trva]).predict(len(te)),
        "mean": MeanBaseline().fit(ds.delta[trva]).predict(len(te)),
        PAPER_LM: fit_predict_linear(ds, trva, te, e["paper_K"], e["paper_lambda"], emb),
        TUNED_LM: fit_predict_linear(ds, trva, te, K, lam, emb),
    }
    out = []
    for target_mode, excl_flag in (("excluded", True), ("included", False)):
        masks = gene_masks(obs, ds.ctrl_mean, ds.target_idx[te], excl_flag, e["n_top_de"], e["n_top_expr"])
        for name, p in preds.items():
            df = per_pert_metrics(p, obs, masks)
            df["model"], df["split"], df["target"] = name, label, target_mode
            out.append(df)
    info = {"split": label, "K": K, "lambda": lam, "n_train": len(tr), "n_val": len(va), "n_test": len(te)}
    return pd.concat(out), info, grid


def get_splits(ds, cfg):
    """Yield (label, seed, split) for the random-seed runs and any GEARS/paper split files in the config."""
    e = cfg["eval"]
    for s in e["seeds"]:
        yield f"seed{s}", s, make_split(len(ds.perts), s, tuple(e["split_fracs"]))
    for label, path in (e.get("split_files") or {}).items():
        with open(path) as f:
            yield label, 0, split_from_conditions(ds, json.load(f))


def summarize(per_pert: pd.DataFrame, group_col="split") -> pd.DataFrame:
    """Mean over test perturbations within a split, then mean +- sd over splits."""
    per_split = per_pert.groupby(["target", "model", "subset", group_col])[["l2", "pearson"]].mean().reset_index()
    g = per_split.groupby(["target", "model", "subset"])[["l2", "pearson"]].agg(["mean", "std"])
    g.columns = [f"{a}_{b}" for a, b in g.columns]
    return g.reset_index()


def format_table(summary: pd.DataFrame, target: str) -> str:
    def f(m, s): return "n/a" if np.isnan(m) else f"{m:.3f} ± {0 if np.isnan(s) else s:.3f}"
    wide = {}
    for _, r in summary[summary.target == target].iterrows():
        wide.setdefault(r["model"], {})[f"{r['subset']} L2"] = f(r["l2_mean"], r["l2_std"])
        wide[r["model"]][f"{r['subset']} Pearson"] = f(r["pearson_mean"], r["pearson_std"])
    return pd.DataFrame(wide).T.to_markdown()
