"""Spectral embedding of a weighted symmetric gene graph.

Normalized Laplacian L = I - D^-1/2 W D^-1/2. Nodes with degree 0 are removed (they get no embedding and are NOT predictable).
The remaining graph has c connected components, so L has exactly c zero eigenvalues (component indicators scaled by sqrt(degree)).
A single Lanczos run cannot resolve such a degenerate eigenspace (verified: ARPACK returned 3 of 39), so every component is
decomposed on its own (dense for small, shift-invert eigsh for large), its trivial zero mode is dropped (and checked), and the
non-trivial eigenpairs of all components are merged by eigenvalue; the first K_max overall are the embedding.
Columns are scaled to unit RMS over embedded nodes so that the ridge penalty acts on a comparable scale across embeddings
(unit-norm eigenvectors would have entries ~1/sqrt(n))."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.sparse.csgraph import connected_components

DENSE_MAX = 1500


def _component_eigs(Wc: sp.csr_matrix, k: int):
    """Smallest k+1 eigenpairs of L for one connected component; returns (lam[1:], vec[:,1:]) with the trivial mode removed."""
    n = Wc.shape[0]
    d = np.asarray(Wc.sum(1)).ravel()
    Dm = sp.diags(d ** -0.5)
    L = sp.identity(n, format="csr") - Dm @ Wc @ Dm
    L = (L + L.T) / 2
    if n <= DENSE_MAX:
        lam, vec = np.linalg.eigh(L.toarray())
        lam, vec = lam[: k + 1], vec[:, : k + 1]
    else:
        lam, vec = sla.eigsh(L.tocsc(), k=k + 1, sigma=-0.01, which="LM", tol=1e-10,
                             v0=np.random.RandomState(0).rand(n) - 0.5)
        o = np.argsort(lam); lam, vec = lam[o], vec[:, o]
    assert abs(lam[0]) < 1e-6, f"component of size {n}: smallest eigenvalue {lam[0]} is not ~0"
    return lam[1:], vec[:, 1:]


def spectral_embedding(W: sp.spmatrix, k_max: int = 50):
    """Returns dict(emb (n, k_max) NaN rows for isolated nodes, eigvals (k_max,), n_components, comp_sizes, trivial_eigvals,
    embedded (bool n), component (int n, -1 isolated), n_isolated, all_eigvals)."""
    W = sp.csr_matrix(W, dtype=np.float64)
    n = W.shape[0]
    deg = np.asarray(W.sum(1)).ravel()
    keep = deg > 0
    idx = np.where(keep)[0]
    Wk = W[idx][:, idx]
    ncomp, lab = connected_components(Wk, directed=False)
    cands = []                                              # (eigenvalue, component, column within that component's vectors)
    store = {}
    for c in range(ncomp):
        members = np.where(lab == c)[0]
        if len(members) < 2:
            continue
        lam, vec = _component_eigs(Wk[members][:, members], min(k_max, len(members) - 2) if len(members) > 2 else 1)
        store[c] = (members, vec)
        cands += [(lam[j], c, j) for j in range(len(lam))]
    cands.sort(key=lambda t: t[0])
    chosen = cands[:k_max]
    emb_k = np.zeros((len(idx), len(chosen)))
    for col, (_, c, j) in enumerate(chosen):
        members, vec = store[c]
        emb_k[members, col] = vec[:, j]
    rms = np.sqrt((emb_k ** 2).mean(0, keepdims=True))
    emb = np.full((n, k_max), np.nan)
    emb[idx, : emb_k.shape[1]] = emb_k / rms
    component = np.full(n, -1); component[idx] = lab
    lam_sel = np.array([t[0] for t in chosen])
    return dict(emb=emb, eigvals=lam_sel, n_components=int(ncomp), comp_sizes=np.sort(np.bincount(lab))[::-1],
                trivial_eigvals=np.zeros(ncomp), embedded=keep, component=component, n_isolated=int((~keep).sum()),
                all_eigvals=np.concatenate([np.zeros(ncomp), lam_sel]),
                source_component=np.array([t[1] for t in chosen]))
