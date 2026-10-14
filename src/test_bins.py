"""Synthetic checks for the bin-level graph code (pytest src/test_bins.py)."""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

from .hic.bins import _msf_mask, _topk_mask, _to_graph
from .hic.position import position_basis, r2_columns, residualize


def test_msf_is_spanning_and_chain_for_distance_weights():
    pos = np.sort(np.random.RandomState(0).rand(40)) * 1e7
    W = np.exp(-np.abs(pos[:, None] - pos[None, :]) / 5e6); np.fill_diagonal(W, 0)
    M = _msf_mask(W); M = M | M.T
    assert M.sum() // 2 == 39
    ncomp, _ = connected_components(sp.csr_matrix(M))
    assert ncomp == 1
    # for a monotone-in-distance weight the maximum spanning tree is the path of consecutive bins
    r, c = np.where(np.triu(M, 1))
    assert set(c - r) == {1}


def test_topk_mask_and_graph_normalisation():
    W = np.random.RandomState(1).rand(30, 30); W = W + W.T; np.fill_diagonal(W, 0)
    m = _topk_mask(W, 4)
    assert (m.sum(1) == 4).all()
    r, c = np.where(np.triu(m | m.T, 1))
    A = _to_graph([r], [c], [W[r, c]], 30, 0.999)
    assert abs(A.data.mean() - 1) < 0.05 and (A != A.T).nnz == 0


def test_residualize_removes_position_signal():
    rng = np.random.RandomState(2)
    chr_idx = np.repeat(np.arange(23), 40)
    pos = np.tile(np.arange(40) * 5e5, 23)
    X = position_basis(chr_idx, pos)
    E = np.stack([np.sin(pos / 6e6) + chr_idx * 0.1, rng.randn(len(pos))], 1)
    ok = np.ones(len(pos), bool)
    R = residualize(E, ok, chr_idx, pos)
    assert r2_columns(E, X)[0] > 0.95            # first column is positional
    assert r2_columns(R, X).max() < 1e-6          # residual has no chromosome/position content
    assert np.allclose(np.sqrt((R ** 2).mean(0)), 1)
