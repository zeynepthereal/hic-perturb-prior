"""Synthetic checks of the graph-neighbour model: it must recover a planted graph and must not beat the mean baseline on a random graph."""
import numpy as np
import scipy.sparse as sp

from .data import PerturbationDataset
from .graphnn import GraphNeighbourModel, GraphSource, fit_predict_graph_neighbour


def _modules(n_mod=5, per=40, n_gene=60, noise=0.3, seed=0):
    rng = np.random.RandomState(seed)
    n = n_mod * per
    mod = np.repeat(np.arange(n_mod), per)
    sig = rng.randn(n_mod, n_gene) * 2.0                                              # module-specific response signature
    delta = sig[mod] + noise * rng.randn(n, n_gene)
    return n, mod, sig, delta, rng


def _planted_graph(mod, rng, k=6, p_cross=0.02):
    n = len(mod); rows, cols = [], []
    for i in range(n):
        same = np.where(mod == mod[i])[0]; same = same[same != i]
        for j in rng.choice(same, k, replace=False): rows.append(i); cols.append(j)
        if rng.rand() < p_cross: rows.append(i); cols.append(rng.randint(n))
    W = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n)); return (W + W.T).tocsr()


def _random_graph(n, rng, k=6):
    rows = np.repeat(np.arange(n), k); cols = rng.randint(n, size=n * k)
    W = sp.csr_matrix((np.ones(n * k), (rows, cols)), shape=(n, n)); W.setdiag(0); W.eliminate_zeros(); return (W + W.T).tocsr()


def _l2(pred, obs): return float(np.sqrt(((pred - obs) ** 2).sum(1)).mean())


def test_recovers_planted_graph_and_beats_mean_baseline():
    n, mod, sig, delta, rng = _modules()
    perm = rng.permutation(n); tr, te = perm[: int(.7 * n)], perm[int(.7 * n):]
    m = GraphNeighbourModel(_planted_graph(mod, rng)).fit(delta[tr], tr)             # node = perturbation index
    pred = m.predict(te)
    base = np.tile(delta[tr].mean(0), (len(te), 1))
    assert _l2(pred, delta[te]) < 0.45 * _l2(base, delta[te])                        # large gain over the mean baseline
    # recovers the planted structure: the prediction is close to the module signature of the held-out perturbation
    assert np.mean([np.corrcoef(pred[i], sig[mod[te[i]]])[0, 1] for i in range(len(te))]) > 0.95


def test_random_graph_does_not_beat_mean_baseline():
    n, mod, sig, delta, rng = _modules(seed=1)
    losses = []
    for s in range(5):
        r = np.random.RandomState(10 + s); perm = r.permutation(n); tr, te = perm[: int(.7 * n)], perm[int(.7 * n):]
        pred = GraphNeighbourModel(_random_graph(n, r)).fit(delta[tr], tr).predict(te)
        base = np.tile(delta[tr].mean(0), (len(te), 1))
        losses.append((_l2(pred, delta[te]), _l2(base, delta[te])))
    g, b = np.mean([x[0] for x in losses]), np.mean([x[1] for x in losses])
    assert g >= b                                                                    # random neighbours average fewer perturbations -> never better than the global mean
    assert all(x[0] >= x[1] - 1e-9 for x in losses)


def test_weights_fallback_and_nan_rows():
    # path 0-1-2-3-4 with weights; hops=1 => node 2 averages nodes 1 and 3 with P-normalised weights
    W = sp.lil_matrix((5, 5)); W[0, 1] = W[1, 0] = 1.0; W[1, 2] = W[2, 1] = 2.0; W[2, 3] = W[3, 2] = 6.0; W[3, 4] = W[4, 3] = 1.0
    delta = np.arange(5)[:, None] * np.array([[1.0, 10.0]])                           # delta_i = i * (1, 10)
    m = GraphNeighbourModel(W.tocsr(), hops=1, same_node_weight=0.0).fit(delta[[0, 1, 3, 4]], np.array([0, 1, 3, 4]))
    p = m.predict(np.array([2]))[0]
    assert np.allclose(p, (2 / 8) * delta[1] + (6 / 8) * delta[3])                    # edge-weighted: 2:6
    # two hops: manual S
    m2 = GraphNeighbourModel(W.tocsr(), hops=2, gamma=0.5, same_node_weight=0.0).fit(delta[[0, 1, 3, 4]], np.array([0, 1, 3, 4]))
    Wd = W.toarray(); P = Wd / Wd.sum(1, keepdims=True); S = P + 0.5 * P @ P
    wts = S[2, [0, 1, 3, 4]]
    assert np.allclose(m2.predict(np.array([2]))[0], (wts @ delta[[0, 1, 3, 4]]) / wts.sum())
    # isolated test node -> exactly the mean baseline; unmapped target -> NaN row
    W6 = sp.block_diag([W.tocsr(), sp.csr_matrix((1, 1))]).tocsr()
    m3 = GraphNeighbourModel(W6, hops=2).fit(delta[[0, 1, 3, 4]], np.array([0, 1, 3, 4]))
    out = m3.predict(np.array([5, -1]))
    assert np.allclose(out[0], delta[[0, 1, 3, 4]].mean(0)) and np.isnan(out[1]).all() and bool(m3.last_fallback_[0])
    # same-node training perturbation (two genes in one bin) gets the same_node_weight
    m4 = GraphNeighbourModel(W.tocsr(), hops=1, same_node_weight=100.0).fit(delta[[1, 3]], np.array([2, 1]))
    assert np.linalg.norm(m4.predict(np.array([2]))[0] - delta[1]) < 0.1 * np.linalg.norm(delta[1] - delta[3])


def test_dataset_interface_mirrors_linear_path():
    n, mod, sig, delta, rng = _modules(n_mod=3, per=30, n_gene=40, seed=2)
    ds = PerturbationDataset(delta, np.array([f"g{i}" for i in range(40)]), np.array([f"p{i}" for i in range(n)]), np.zeros(n, int), np.full(n, 50), np.zeros(40), {})
    node = np.arange(n); node[::11] = -1                                              # some targets have no node
    src = GraphSource(_planted_graph(mod, rng), node)
    perm = rng.permutation(n); tr, te = perm[:60], perm[60:]
    pred = fit_predict_graph_neighbour(ds, tr, te, src)
    assert pred.shape == (len(te), 40) and np.isnan(pred[node[te] < 0]).all() and not np.isnan(pred[node[te] >= 0]).any()
    assert (src.predictable()[te] == (node[te] >= 0)).all()
