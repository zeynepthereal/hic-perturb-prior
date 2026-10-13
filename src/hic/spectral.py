"""Spectral embedding of a weighted symmetric gene graph.

Normalized Laplacian L = I - D^-1/2 W D^-1/2. Nodes with degree 0 are removed (they get no embedding and are NOT predictable).
For the remaining graph with c connected components, L has exactly c zero eigenvalues whose eigenvectors are (scaled) component
indicators; these `c` trivial vectors are dropped explicitly (and checked to have eigenvalue ~0) and the first K_max non-trivial
eigenvectors are returned. Columns are scaled to unit RMS over embedded nodes so that the ridge penalty of the linear model
acts on a scale comparable between embeddings (raw unit-norm eigenvectors would have entries ~1/sqrt(n))."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.sparse.csgraph import connected_components


def spectral_embedding(W: sp.spmatrix, k_max: int = 50, seed: int = 0, tol: float = 1e-9):
    """Returns dict(emb (n, k_max) with NaN rows for isolated nodes, eigvals (k_max,) of L, n_components, comp_sizes,
    trivial_eigvals (c,), embedded (bool n), component (int n, -1 for isolated))."""
    W = sp.csr_matrix(W, dtype=np.float64)
    n = W.shape[0]
    deg = np.asarray(W.sum(1)).ravel()
    keep = deg > 0
    idx = np.where(keep)[0]
    Wk = W[idx][:, idx]
    dk = deg[idx]
    ncomp, lab = connected_components(Wk, directed=False)
    S = sp.diags(dk ** -0.5) @ Wk @ sp.diags(dk ** -0.5)
    S = (S + S.T) / 2
    n_eig = ncomp + k_max + 5
    rng = np.random.RandomState(seed)
    vals, vecs = sla.eigsh(S, k=min(n_eig, Wk.shape[0] - 2), which="LA", tol=tol, v0=rng.rand(Wk.shape[0]) - 0.5, ncv=min(2 * n_eig + 20, Wk.shape[0]))
    order = np.argsort(-vals)
    vals, vecs = vals[order], vecs[:, order]
    lam = 1.0 - vals                                   # eigenvalues of L, ascending
    assert np.all(lam[:ncomp] < 1e-6), f"expected {ncomp} zero eigenvalues of L, got {lam[:ncomp+2]}"
    lam_nt, vec_nt = lam[ncomp: ncomp + k_max], vecs[:, ncomp: ncomp + k_max]
    rms = np.sqrt((vec_nt ** 2).mean(0, keepdims=True))
    emb = np.full((n, k_max), np.nan)
    emb[idx] = vec_nt / rms
    component = np.full(n, -1)
    component[idx] = lab
    sizes = np.bincount(lab)
    return dict(emb=emb, eigvals=lam_nt, n_components=int(ncomp), comp_sizes=np.sort(sizes)[::-1], trivial_eigvals=lam[:ncomp],
                embedded=keep, component=component, n_isolated=int((~keep).sum()), all_eigvals=lam)
