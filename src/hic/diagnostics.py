"""Per-variant diagnostic metrics shared by the diagnostics scripts: structure, R² on chromosome/position, A/B correlation."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..baselines import EMBEDDINGS
from .embeddings import ART, load_spectral
from .position import position_basis, r2_columns

OUT = None  # set by callers: results/hic_v2


def structure(v: str) -> dict:
    z = load_spectral(v); W = sp.load_npz(ART / f"{v}.npz"); d = np.diff(W.tocsr().indptr); d = d[d > 0]
    return dict(variant=v, nodes_embedded=int(z["embedded"].sum()), edges=int(W.nnz // 2), components=int(np.sum(z["comp_sizes"] > 0)),
                deg_median=float(np.median(d)), deg_mean=float(d.mean()), deg_p99=float(np.percentile(d, 99)), deg_max=int(d.max()),
                eigval_1=float(z["eigvals"][0]), eigval_10=float(z["eigvals"][9]), eigval_50=float(z["eigvals"][49]))


def r2_rows(v: str, out_dir) -> list[dict]:
    e = EMBEDDINGS[v]; z = load_spectral(v.split(":")[0])
    E, ok = e._node_matrix(); ci, pos = e._node_coords(z)
    X = {k: position_basis(ci, pos, kind=k) for k in ("chromosome", "spline", "cubic")}
    rows = []
    for K in (10, 50):
        Y = E[ok][:, :K]
        rows.append(dict(variant=v, dims=K, R2_chromosome=float(r2_columns(Y, X["chromosome"][ok]).mean()), R2_chrom_cubic=float(r2_columns(Y, X["cubic"][ok]).mean()),
                         R2_chrom_spline=float(r2_columns(Y, X["spline"][ok]).mean()), frac_nodes_nonzero_in_a_dim=float((np.abs(Y) > 1e-8).any(1).mean())))
    return rows


def ab_rows(v: str, out_dir) -> list[dict]:
    e = EMBEDDINGS[v]; z = load_spectral(v.split(":")[0])
    E, ok = e._node_matrix(); ci, pos = e._node_coords(z)
    t = pd.read_csv(f"{out_dir}/ab_{int(z['res']) // 1000}k.csv")
    Xs = position_basis(ci, pos, kind="spline")
    rows = []
    for name, a in (("AB_PC1", t.ab_pc1.values), ("AB_best_of_3", t.ab_best_of_3.values)):
        m = ok & np.isfinite(a)
        for K in (10, 50):
            Y = E[m][:, :K]
            R2 = float(r2_columns(a[m][:, None], np.hstack([Y, np.ones((m.sum(), 1))]))[0])
            ar = a[m] - Xs[m] @ np.linalg.lstsq(Xs[m], a[m], rcond=None)[0]
            Yr = Y - Xs[m] @ np.linalg.lstsq(Xs[m], Y, rcond=None)[0]
            R2b = float(r2_columns(ar[:, None], np.hstack([Yr, np.ones((m.sum(), 1))]))[0])
            rows.append(dict(variant=v, ab=name, dims=K, multiple_R=float(np.sqrt(max(R2, 0))), R2_beyond_position=R2b))
    return rows
