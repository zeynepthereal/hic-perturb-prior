"""Metrics, model selection on the validation split, and the multi-seed experiment loop."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .baselines import EMBEDDINGS, LinearModel, MeanBaseline, NoChange
from .data import PerturbationDataset, make_split, split_from_conditions
from .utils import REPO_ROOT, set_seed



def fixed_name(emb: str, K: int, lam: float) -> str:
    return f"linear[{emb}] (fixed K={K}, lam={lam})"


def lm_name(emb: str, tuned: bool) -> str:
    return f"linear[{emb}] ({'val-tuned' if tuned else 'paper hp K=10, lam=0.1'})"


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
    """G is fit on ALL training perturbations (as in the paper's code); W and b only on those the embedding can embed.
    Returns (len(test_idx), n_gene) with NaN rows for test perturbations the embedding cannot embed."""
    emb = EMBEDDINGS[embedding]
    G, P_all = emb.fit(ds.delta[train_idx], train_idx, ds, K)
    ok = emb.predictable(ds)
    tr = train_idx[ok[train_idx]]
    m = LinearModel(lam).fit(ds.delta[tr], G, P_all[tr])
    out = np.full((len(test_idx), ds.delta.shape[1]), np.nan)
    t_ok = ok[test_idx]
    out[t_ok] = m.predict(P_all[test_idx[t_ok]])
    return out


def select_hparams(ds, split, Ks, lams, embedding, exclude_target):
    """Pick (K, lambda) minimising mean L2 (all read-out genes) on the validation perturbations the embedding can predict."""
    va = split["val"][EMBEDDINGS[embedding].predictable(ds)[split["val"]]]
    obs = ds.delta[va]
    allowed = gene_masks(obs, ds.ctrl_mean, ds.target_idx[va], exclude_target)["all"]
    best, grid = None, []
    for K in Ks:
        for lam in lams:
            pred = fit_predict_linear(ds, split["train"], va, K, lam, embedding)
            err = float(np.sqrt((((pred - obs) ** 2) * allowed).sum(axis=1)).mean())
            grid.append({"embedding": embedding, "K": K, "lambda": lam, "val_l2": err, "n_val": len(va)})
            if best is None or err < best[0]:
                best = (err, K, lam)
    return best[1], best[2], pd.DataFrame(grid)


def run_split(ds: PerturbationDataset, cfg: dict, split: dict, label, seed: int):
    """Fit all models on one train/val/test split and score them on the eval sets:
      common            test perturbations predictable by ALL compared embeddings  (primary comparison)
      extended:<emb>    all test perturbations predictable by <emb> (superset of common), <emb> vs mean baseline only
      <name>            (optional, eval.extra_eval_sets) common ∩ a named perturbation list, e.g. the perturbations kept by GEARS
    """
    set_seed(seed)
    e = cfg["eval"]
    embs, excl = e["embeddings"], e["exclude_target_primary"]
    tr, va, te = split["train"], split["val"], split["test"]
    fit_idx = np.concatenate([tr, va]) if e.get("refit_trainval", True) else tr
    obs = ds.delta[te]
    preds = {"no_change": NoChange().fit(ds.delta[fit_idx]).predict(len(te)),
             "mean": MeanBaseline().fit(ds.delta[fit_idx]).predict(len(te))}
    ok, info, grids, emb_models = {}, {"split": label, "n_train": len(tr), "n_val": len(va), "n_test": len(te)}, [], {}
    for emb in embs:
        K, lam, grid = select_hparams(ds, split, e["K_grid"], e["lambda_grid"], emb, excl)
        preds[lm_name(emb, False)] = fit_predict_linear(ds, fit_idx, te, e["paper_K"], e["paper_lambda"], emb)
        preds[lm_name(emb, True)] = fit_predict_linear(ds, fit_idx, te, K, lam, emb)
        emb_models[emb] = [lm_name(emb, False), lm_name(emb, True)]
        for K2, lam2 in e.get("extra_fixed_hp", []):                  # e.g. [[20, 0.1], [50, 0.1]]: fixed (not tuned) secondary settings
            preds[fixed_name(emb, K2, lam2)] = fit_predict_linear(ds, fit_idx, te, K2, lam2, emb)
            emb_models[emb].append(fixed_name(emb, K2, lam2))
        ok[emb] = EMBEDDINGS[emb].predictable(ds)[te]
        info[f"{emb}_K"], info[f"{emb}_lambda"] = K, lam
        grids.append(grid)
    common = np.logical_and.reduce([ok[m] for m in embs])
    sets = {"common": (common, list(preds))}
    for emb in embs:
        if (ok[emb] != common).any():
            sets[f"extended:{emb}"] = (ok[emb], ["mean"] + emb_models[emb])
    for set_name, path in (e.get("extra_eval_sets") or {}).items():   # restrict `common` to a named perturbation list
        names = set(pd.read_csv(REPO_ROOT / path)["perturbation"])
        sets[set_name] = (common & np.isin(ds.perts[te], list(names)), list(preds))
    out = []
    for set_name, (sel, models) in sets.items():
        info[f"n_eval_{set_name}"] = int(sel.sum())
        if sel.sum() == 0:
            continue
        for target_mode, excl_flag in (("excluded", True), ("included", False)):
            masks = gene_masks(obs[sel], ds.ctrl_mean, ds.target_idx[te][sel], excl_flag, e["n_top_de"], e["n_top_expr"])
            for name in models:
                df = per_pert_metrics(preds[name][sel], obs[sel], masks)
                df["pert_i"] = te[sel][df["pert_i"].values]
                df["model"], df["split"], df["target"], df["eval_set"] = name, label, target_mode, set_name
                out.append(df)
    return pd.concat(out), info, pd.concat(grids)


def get_splits(ds, cfg):
    """Yield (label, seed, split) for random-seed runs plus any GEARS-style split files ({train,val,test: [conditions]})."""
    e = cfg["eval"]
    for s in e.get("seeds") or []:
        yield f"seed{s}", s, make_split(len(ds.perts), s, tuple(e["split_fracs"]))
    for label, path in (e.get("split_files") or {}).items():
        with open(REPO_ROOT / path) as f:
            yield label, int(label.split("seed")[-1]) if "seed" in label else 0, split_from_conditions(ds, json.load(f))


KEYS = ["eval_set", "target", "model", "subset"]


def summarize(per_pert: pd.DataFrame) -> pd.DataFrame:
    """Mean over test perturbations within a split, then mean +- sd over splits; n_perts = mean #perts per split."""
    per_split = per_pert.groupby(KEYS + ["split"]).agg(l2=("l2", "mean"), pearson=("pearson", "mean"),
                                                       n=("pert_i", "nunique")).reset_index()
    g = per_split.groupby(KEYS).agg(l2_mean=("l2", "mean"), l2_std=("l2", "std"), pearson_mean=("pearson", "mean"),
                                    pearson_std=("pearson", "std"), n_perts=("n", "mean"), n_splits=("split", "nunique"))
    return g.reset_index()


def format_table(summary: pd.DataFrame, eval_set: str, target: str) -> str:
    def f(m, s): return "n/a" if np.isnan(m) else f"{m:.3f} ± {0 if np.isnan(s) else s:.3f}"
    sub = summary[(summary.eval_set == eval_set) & (summary.target == target)]
    wide = {}
    for _, r in sub.iterrows():
        wide.setdefault(r["model"], {})[f"{r['subset']} L2"] = f(r["l2_mean"], r["l2_std"])
        wide[r["model"]][f"{r['subset']} Pearson"] = f(r["pearson_mean"], r["pearson_std"])
    n, k = sub.n_perts.mean(), sub.n_splits.max()
    title = f"**eval_set = {eval_set}**, target gene {target}; {n:.0f} test perturbations per split (mean over {k} splits)"
    return title + "\n\n" + pd.DataFrame(wide).T.to_markdown()
