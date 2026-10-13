"""Synthetic checks for the Hi-C graph code (pytest src/test_hic.py)."""
import numpy as np
import scipy.sparse as sp

from .hic.embeddings import double_edge_swap
from .hic.graphs import ice_coo, symmetrize_union, topk_rows
from .hic.spectral import spectral_embedding


def test_ice_equalizes_marginals():
    rng = np.random.RandomState(0)
    n = 60
    true = rng.rand(n, n); true = true + true.T
    bias = rng.uniform(0.3, 3, n)
    M = rng.poisson(200 * true * np.outer(bias, bias) / 4).astype(float)
    M = np.triu(M) ; M = M + np.triu(M, 1).T
    iu = np.triu_indices(n)
    b, valid = ice_coo(iu[0], iu[1], np.triu(M)[iu], n)
    bal = M / np.outer(b, b)
    r = bal.sum(1)
    assert valid.all() and r.std() / r.mean() < 1e-3


def test_spectral_handles_components_and_isolated_nodes():
    # two rings of 30 nodes (2 components) + 3 isolated nodes
    def ring(n, off):
        a = np.arange(n) + off; b = (np.arange(n) + 1) % n + off
        return a, b
    a1, b1 = ring(30, 0); a2, b2 = ring(30, 30)
    a, b = np.concatenate([a1, a2]), np.concatenate([b1, b2])
    W = sp.csr_matrix((np.ones(len(a)), (a, b)), shape=(63, 63)); W = W + W.T
    r = spectral_embedding(W, k_max=6)
    assert r["n_components"] == 2 and r["n_isolated"] == 3
    assert np.all(r["trivial_eigvals"] < 1e-6) and r["eigvals"][0] > 1e-6          # trivial vectors dropped
    assert np.isnan(r["emb"][60:]).all() and not np.isnan(r["emb"][:60]).any()   # isolated nodes unembedded
    assert np.allclose(np.sqrt((r["emb"][:60] ** 2).mean(0)), 1)


def test_double_edge_swap_preserves_degrees_and_weights():
    rng = np.random.RandomState(1)
    n = 80
    iu = np.triu_indices(n, 1)
    m = rng.rand(len(iu[0])) < 0.1
    W = sp.csr_matrix((rng.rand(m.sum()) + 0.5, (iu[0][m], iu[1][m])), shape=(n, n)); W = W + W.T
    R = double_edge_swap(W, seed=3)
    assert (np.asarray((W > 0).sum(1)).ravel() == np.asarray((R > 0).sum(1)).ravel()).all()
    assert R.diagonal().sum() == 0
    assert np.allclose(np.sort(sp.triu(W, 1).data), np.sort(sp.triu(R, 1).data))
    assert (R != W).nnz > 0


def test_topk_and_union():
    W = sp.csr_matrix(np.array([[0, 3, 2, 1], [3, 0, 0, 0], [2, 0, 0, 5], [1, 0, 5, 0.]]))
    D = topk_rows(W, 1)
    S = symmetrize_union(D)
    assert S[0, 1] == 3 and S[2, 3] == 5 and S[0, 2] == 0 and (S != S.T).nnz == 0
