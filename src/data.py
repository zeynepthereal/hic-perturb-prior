"""Data loading, preprocessing and perturbation-level splits.

Output of preprocessing is a `PerturbationDataset`: one row per perturbation,
one column per read-out gene, values = pseudobulk(perturbation) - pseudobulk(controls)
on log1p-normalised expression.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class PerturbationDataset:
    delta: np.ndarray        # (n_pert, n_genes) pseudobulk(pert) - pseudobulk(control)
    genes: np.ndarray        # (n_genes,) read-out gene symbols (the selected HVGs)
    perts: np.ndarray        # (n_pert,) perturbation labels
    target_idx: np.ndarray   # (n_pert,) column index of the perturbed gene in `genes`; -1 if not a read-out gene
    n_cells: np.ndarray      # (n_pert,) cells per perturbation
    ctrl_mean: np.ndarray    # (n_genes,) control pseudobulk (log1p-normalised)
    meta: dict               # preprocessing choices, source file, checksum

    def subset(self, mask_or_idx) -> "PerturbationDataset":
        i = np.asarray(mask_or_idx)
        return PerturbationDataset(self.delta[i], self.genes, self.perts[i], self.target_idx[i],
                                   self.n_cells[i], self.ctrl_mean, self.meta)

    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, delta=self.delta, genes=self.genes.astype(str), perts=self.perts.astype(str),
                            target_idx=self.target_idx, n_cells=self.n_cells, ctrl_mean=self.ctrl_mean,
                            meta=json.dumps(self.meta))

    @classmethod
    def load(cls, path: str | Path) -> "PerturbationDataset":
        z = np.load(path, allow_pickle=False)
        return cls(z["delta"], z["genes"], z["perts"], z["target_idx"], z["n_cells"], z["ctrl_mean"],
                   json.loads(str(z["meta"])))


def md5(path: str | Path, chunk: int = 1 << 24) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def seurat_hvg_from_moments(mean, var, n_top_genes, n_bins=20):
    """Replica of scanpy's flavor='seurat' HVG selection (sc.pp.highly_variable_genes, n_top_genes set), but fed
    with streamed per-gene mean/variance of the *normalised, non-logged* values (= expm1 of the log1p data, as scanpy
    does). Verified against scanpy in src/test_synthetic.py. Uses all cells (no perturbation labels)."""
    mean = np.where(mean == 0, 1e-12, mean)
    disp = var / mean
    disp[disp == 0] = np.nan
    df = pd.DataFrame({"m": np.log1p(mean), "d": np.log(disp)})
    df["bin"] = pd.cut(df["m"], bins=n_bins)
    g = df.groupby("bin", observed=False)["d"]
    mu, sd = g.mean(), g.std(ddof=1)
    one = sd.isnull()
    sd[one] = mu[one]
    mu[one] = 0
    dn = ((df["d"].values - mu[df["bin"].values].values) / sd[df["bin"].values].values)
    cut = np.sort(dn[~np.isnan(dn)])[::-1][n_top_genes - 1]
    return dn >= cut


def preprocess(h5ad_path: str | Path, cfg: dict):
    """Streaming (the file is a dense 310k x 8.5k matrix, too big to hold comfortably):
    per cell normalize_total(target_sum) -> log1p; HVGs = top-n by scanpy-seurat dispersion over ALL cells;
    pseudobulk mean of log1p-normalised expression per perturbation; delta = pseudobulk(pert) - pseudobulk(control).
    Returns (PerturbationDataset restricted to the HVGs, QC DataFrame). No perturbation is filtered by quality."""
    import anndata as ad

    d = cfg["data"]
    a = ad.read_h5ad(h5ad_path, backed="r")
    labels = a.obs[d["pert_col"]].astype(str).values
    ctrl = labels == d["control_label"]
    if ctrl.sum() == 0:
        raise ValueError(f"no cells with {d['pert_col']}=={d['control_label']!r}")
    names, inv = np.unique(labels, return_inverse=True)       # groups incl. control
    ctrl_g = int(np.where(names == d["control_label"])[0][0])
    n_g, n_genes = len(names), a.n_vars
    S_log = np.zeros((n_g, n_genes)); S_norm = np.zeros((n_g, n_genes))
    gsum = np.zeros(n_genes); gsq = np.zeros(n_genes)
    for i in range(0, a.n_obs, d.get("chunk", 20000)):
        c = np.asarray(a.X[i:i + d.get("chunk", 20000)], dtype=np.float64)
        if i == 0 and not np.allclose(c[:2000], np.round(c[:2000])):
            raise ValueError("X is not integer-valued counts; inspect the file")
        norm = c / c.sum(1, keepdims=True) * d["target_sum"]
        lg = np.log1p(norm)
        np.add.at(S_log, inv[i:i + len(c)], lg); np.add.at(S_norm, inv[i:i + len(c)], norm)
        gsum += norm.sum(0); gsq += (norm ** 2).sum(0)
    n = np.bincount(inv, minlength=n_g).astype(float)
    N = a.n_obs
    g_mean = gsum / N
    g_var = (gsq / N - g_mean ** 2) * N / (N - 1)
    hvg = seurat_hvg_from_moments(g_mean, g_var, d["n_hvg"])

    all_genes = a.var_names.values.astype(str)
    pb_log = S_log / n[:, None]                      # (group, all genes)
    ctrl_log = pb_log[ctrl_g]
    keep = (np.arange(n_g) != ctrl_g) & (n >= d["min_cells"])
    perts = names[keep]
    pos_all = {g: i for i, g in enumerate(all_genes)}
    hv_idx = np.where(hvg)[0]
    pos_hv = {all_genes[j]: k for k, j in enumerate(hv_idx)}
    target_idx = np.array([pos_hv.get(p, -1) for p in perts])
    delta = (pb_log - ctrl_log)[keep][:, hv_idx]

    # QC (no filtering): knockdown of the target gene in the full gene matrix, pseudocount-0.01 log2 FC of mean
    # normalised expression, plus the log1p-scale difference.
    m_norm, c_norm = S_norm / n[:, None], (S_norm / n[:, None])[ctrl_g]
    rows = []
    for gi in np.where(keep)[0]:
        t = names[gi]; j = pos_all.get(t, -1)
        rows.append({"perturbation": t, "n_cells": int(n[gi]), "target_in_matrix": j >= 0,
                     "target_in_readout": t in pos_hv,
                     "target_ctrl_mean_norm": c_norm[j] if j >= 0 else np.nan,
                     "target_log2fc": np.log2((m_norm[gi, j] + 0.01) / (c_norm[j] + 0.01)) if j >= 0 else np.nan,
                     "target_log1p_diff": (pb_log[gi, j] - ctrl_log[j]) if j >= 0 else np.nan})
    qc = pd.DataFrame(rows)
    meta = {
        "source_file": str(h5ad_path), "source_md5": md5(h5ad_path), "n_hvg": d["n_hvg"],
        "hvg_method": "scanpy-seurat replica on all cells (no perturbation labels used)",
        "target_sum": d["target_sum"], "min_cells": d["min_cells"],
        "n_control_cells": int(n[ctrl_g]), "n_perts_total": int(keep.sum()),
        "n_perts_target_in_readout": int((target_idx >= 0).sum()),
    }
    ds = PerturbationDataset(delta, all_genes[hv_idx], perts, target_idx, n[keep].astype(int),
                             ctrl_log[hv_idx], meta)
    return ds, qc


def make_split(n: int, seed: int, fracs=(0.7, 0.1, 0.2)) -> dict[str, np.ndarray]:
    """Random split of *perturbations* (never cells) into train/val/test index arrays."""
    perm = np.random.RandomState(seed).permutation(n)
    n_tr, n_va = int(round(fracs[0] * n)), int(round(fracs[1] * n))
    return {"train": np.sort(perm[:n_tr]), "val": np.sort(perm[n_tr:n_tr + n_va]), "test": np.sort(perm[n_tr + n_va:])}


def predictable_subset(ds: PerturbationDataset, require_target_in_readout: bool = True) -> PerturbationDataset:
    """Linear model with P = rows of G can only handle perturbations whose target is a read-out gene
    (same restriction as the paper). Apply to ALL methods so they are evaluated on the same perturbations."""
    return ds.subset(np.where(ds.target_idx >= 0)[0]) if require_target_in_readout else ds


def split_from_conditions(ds: PerturbationDataset, set2conditions: dict) -> dict[str, np.ndarray]:
    """Map GEARS-style {'train': [...], 'val': [...], 'test': [...]} (conditions like 'GENE+ctrl') onto ds rows.
    Conditions not present in ds (e.g. dropped by our filters) are skipped and counted in ds.meta by the caller."""
    pos = {p: i for i, p in enumerate(ds.perts)}
    out = {}
    for k, conds in set2conditions.items():
        names = [c.replace("+ctrl", "").replace("ctrl+", "") for c in conds if c != "ctrl"]
        out[k] = np.array(sorted(pos[n] for n in names if n in pos), dtype=int)
    return out
