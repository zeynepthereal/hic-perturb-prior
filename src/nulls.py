"""Lean evaluator and driver for the pre-registered null distributions (S2 / H1 support condition).

Why a second evaluator: the standard path (`evaluate.run_split`) fits ~8 models per embedding and split (~3 s); 2 x 1,000 null seeds x 5 splits
would take hours. Everything that does not depend on the embedding's P (the training PCA G, G'(Y-b), observed deltas, gene masks) is cached;
per seed only the P-dependent algebra and the L2 remain (~4-5 ms per seed and split). The arithmetic is identical to
`evaluate.fit_predict_linear` + `evaluate.per_pert_metrics` (L2, top-N expressed genes, target gene excluded), which a unit test verifies
against the standard path on synthetic data. Only the primary fixed-hyper-parameter model (K, lambda) is evaluated.

CLI (NOT to be run before the pre-registration is frozen):  python -m src.nulls --config configs/prereg/s2_nulls.yaml"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from .baselines import pca_gene_embedding
from .evaluate import gene_masks


class NullEvaluator:
    def __init__(self, ds, splits: dict, K: int, lam: float = 0.1, n_top_expr: int = 1000, n_top_de: int = 20, exclude_target: bool = True):
        self.ds, self.splits, self.K, self.lam = ds, splits, K, lam
        self.n_top_expr, self.n_top_de, self.excl = n_top_expr, n_top_de, exclude_target
        self.G = {lab: pca_gene_embedding(ds.delta[s["train"]] + ds.ctrl_mean, K) for lab, s in splits.items()}   # fit on training perturbations only
        self._prep_cache: dict = {}

    def _prep(self, label, ok):
        key = (label, ok.tobytes())
        if key not in self._prep_cache:
            s, ds, K = self.splits[label], self.ds, self.K
            tr = s["train"][ok[s["train"]]]                              # W and b only from training perturbations the embedding can embed
            te = s["test"][ok[s["test"]]]
            Yt, G = ds.delta[tr].T, self.G[label]
            b = Yt.mean(axis=1, keepdims=True)
            left = np.linalg.solve(G.T @ G + self.lam * np.eye(K), G.T @ (Yt - b))
            obs = ds.delta[te]
            M = gene_masks(obs, ds.ctrl_mean, ds.target_idx[te], self.excl, self.n_top_de, self.n_top_expr)["top1000_expr"]
            self._prep_cache[key] = (tr, te, G, b, left, obs, M)
        return self._prep_cache[key]

    def case_l2(self, P_all: np.ndarray) -> pd.Series:
        """Case-level L2 indexed by (pert_i, split) for all test perturbations whose P row is available (not NaN)."""
        ok = ~np.isnan(P_all).any(1)
        out = {}
        for label in self.splits:
            tr, te, G, b, left, obs, M = self._prep(label, ok)
            P_tr = P_all[tr]
            W = np.linalg.solve((P_tr.T @ P_tr + self.lam * np.eye(self.K)).T, (left @ P_tr).T).T
            pred = (G @ W @ P_all[te].T + b).T
            l2 = np.sqrt((((pred - obs) ** 2) * M).sum(1))
            out.update({(int(i), label): float(v) for i, v in zip(te, l2)})
        s = pd.Series(out, name="l2"); s.index.names = ["pert_i", "split"]
        return s


def _mean_case(series_list) -> pd.Series:
    """per-case mean of several case-level series (the control = mean of per-seed L2 values, never an averaged prediction)."""
    return pd.concat(series_list, axis=1, join="inner").mean(axis=1)


def run_null_analysis(ds, splits, K, lam, hic_P, control_Ps, perm_Ps, rand_Ps, rewire_Ps, n_top_expr=1000, n_top_de=20) -> dict:
    """hic_P / control_Ps ... : P_all matrices (or iterables yielding them) for the Hi-C arm, the control seeds and each null draw."""
    from .stats import h1_support, null_effect
    ev = NullEvaluator(ds, splits, K, lam, n_top_expr, n_top_de)
    ctrl = _mean_case([ev.case_l2(P) for P in control_Ps])
    hic = ev.case_l2(hic_P)
    D_obs = null_effect(ctrl, hic)
    Ds = {name: np.array([null_effect(ctrl, ev.case_l2(P)) for P in Ps]) for name, Ps in (("permutation", perm_Ps), ("random_gaussian", rand_Ps), ("rewire", rewire_Ps))}
    return {"K": K, "D_obs": D_obs, "n_cases": int(len(ctrl)), **h1_support(D_obs, Ds["permutation"], Ds["random_gaussian"], Ds["rewire"]), "null_D": {k: v.tolist() for k, v in Ds.items()}}


def main():
    from .baselines import EMBEDDINGS
    from .data import PerturbationDataset, split_from_conditions
    from .hic.embeddings import GraphSpectralEmbedding
    from .utils import REPO_ROOT, load_config, make_run_dir
    ap = argparse.ArgumentParser(); ap.add_argument("--config", required=True); a = ap.parse_args()
    cfg = load_config(a.config); n = cfg["nulls"]
    ds = PerturbationDataset.load(REPO_ROOT / cfg["data"]["processed_path"])
    splits = {lab: split_from_conditions(ds, json.load(open(REPO_ROOT / p))) for lab, p in cfg["eval"]["split_files"].items()}
    run_dir = make_run_dir(cfg, cfg["name"]); base = cfg["hic"]; res = {}
    for K in cfg["K"]:
        gen = lambda kind, seeds: (GraphSpectralEmbedding(base, kind, s).perturbation_matrix(ds, K) for s in seeds)
        res[K] = run_null_analysis(ds, splits, K, cfg["lambda"], EMBEDDINGS[base].perturbation_matrix(ds, K),
                                   (EMBEDDINGS[c].perturbation_matrix(ds, K) for c in cfg["control"]),
                                   gen("perm", range(n["permutation"]["n"])), gen("random", range(n["random_gaussian"]["n"])), gen("rewire", range(n["rewire"]["n"])))
        print(f"K={K}: D_obs={res[K]['D_obs']:.4f}  H1 null-support: {res[K]['support']}  p={res[K]['p']}")
    (run_dir / "nulls.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
