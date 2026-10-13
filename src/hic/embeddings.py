"""Perturbation embeddings (P) derived from gene graphs, as `Embedding` classes plugging into the eval_set mechanism.

All classes keep G = PCA of the training expression pseudobulks (identical to `train_pca_expr`), so only P changes.
`predictable(ds)`: perturbation target maps to a graph node that has an embedding (degree > 0). Unlike the PCA embedding, the
target does NOT need to be a read-out gene -> larger `extended:<emb>` eval sets.

Names understood by `make_embedding`:
  hic_cis | hic_cis_trans | linear_distance                      graph spectral embeddings
  <variant>:perm:<seed>      gene_label_permutation: same graph/spectrum, which gene sits on which node is permuted
  hic_cis_trans:rewire:<seed> degree_preserving_rewire of the graph, re-embedded (computed lazily, cached)
  <variant>:random:<seed>    random Gaussian embedding with matched dimensions and the same embedded node set
"""
from __future__ import annotations

import hashlib
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..baselines import Embedding, pca_gene_embedding
from .spectral import spectral_embedding

ART = Path(__file__).resolve().parents[2] / "data" / "hic" / "graphs"
TABLES = Path(__file__).resolve().parents[2] / "results" / "hic"


@lru_cache(maxsize=None)
def gene_tables():
    nodes = pd.read_csv(TABLES / "gene_table.csv")
    sym = pd.read_csv(TABLES / "symbol_map.csv").set_index("symbol").gene_id.to_dict()
    node_of_id = {g: i for i, g in enumerate(nodes.gene_id)}
    return nodes, sym, node_of_id


@lru_cache(maxsize=None)
def load_spectral(name: str) -> dict:
    z = np.load(ART / f"{name}_spectral.npz", allow_pickle=False)
    return {k: z[k] for k in z.files}


def double_edge_swap(W: sp.csr_matrix, n_swap_factor: int = 10, seed: int = 0) -> sp.csr_matrix:
    """Degree-preserving rewiring of an undirected weighted graph: pick two edges (a,b,w1), (c,d,w2) and swap partners to
    (a,d,w1), (c,b,w2) if that creates no self-loop or duplicate edge. Node degrees (counts) and the multiset of weights are kept."""
    rng = np.random.RandomState(seed)
    U = sp.triu(W, 1).tocoo()
    a, b, w = U.row.copy(), U.col.copy(), U.data.copy()
    m = len(a)
    edges = set(zip(a.tolist(), b.tolist()))
    key = lambda x, y: (x, y) if x < y else (y, x)
    done = 0
    for _ in range(n_swap_factor * m):
        e1, e2 = rng.randint(m, size=2)
        if e1 == e2:
            continue
        p, q, r, s = a[e1], b[e1], a[e2], b[e2]
        if rng.rand() < 0.5:
            r, s = s, r
        if len({p, q, r, s}) < 4:
            continue
        n1, n2 = key(p, s), key(r, q)
        if n1 in edges or n2 in edges:
            continue
        edges.discard(key(p, q)); edges.discard(key(r, s)); edges.add(n1); edges.add(n2)
        a[e1], b[e1] = n1; a[e2], b[e2] = n2
        done += 1
    R = sp.csr_matrix((w, (a, b)), shape=W.shape)
    return (R + R.T).tocsr()


class GraphSpectralEmbedding(Embedding):
    def __init__(self, variant: str, null: str | None = None, seed: int = 0):
        self.variant, self.null, self.seed = variant, null, seed
        self.name = variant + (f":{null}:{seed}" if null else "")

    # -- node-level embedding (n_nodes, 50) and its row availability
    def _node_matrix(self):
        nodes, _, _ = gene_tables()
        if self.null == "rewire":
            nm = f"{self.variant}__rewire{self.seed}"
            f = ART / f"{nm}_spectral.npz"
            if not f.exists():                      # lazily built the first time it is needed ("run later")
                W = sp.load_npz(ART / f"{self.variant}.npz")
                R = double_edge_swap(W, seed=self.seed)
                sp.save_npz(ART / f"{nm}.npz", R)
                np.savez_compressed(f, **{k: v for k, v in spectral_embedding(R).items()})
            sp_ = load_spectral(nm)
        else:
            sp_ = load_spectral(self.variant)
        E, ok = sp_["emb"].copy(), sp_["embedded"].astype(bool)
        if self.null == "random":
            rng = np.random.RandomState(self.seed)
            E = np.where(ok[:, None], rng.randn(*E.shape), np.nan)
        return E, ok

    def _node_index(self, ds):
        """dataset perturbation -> node whose embedding row is used (-1 if none)."""
        _, ok = self._node_matrix()
        own = _own_index(ds)
        avail = (own >= 0) & ok[np.maximum(own, 0)]
        idx = np.where(avail, own, -1)
        if self.null == "perm":      # gene_label_permutation: same graph, permute which gene sits on which (embedded) node
            emb_nodes = np.where(ok)[0]
            perm = np.random.RandomState(self.seed).permutation(len(emb_nodes))
            pos = np.full(len(ok), 0); pos[emb_nodes] = np.arange(len(emb_nodes))
            idx = np.where(avail, emb_nodes[perm[pos[np.maximum(own, 0)]]], -1)
        return idx

    def predictable(self, ds):
        return self._node_index(ds) >= 0

    def fit(self, Y_train, train_idx, ds, K):
        G = pca_gene_embedding(Y_train + ds.ctrl_mean, K)      # identical G to `train_pca_expr`
        E, _ = self._node_matrix()
        idx = self._node_index(ds)
        P = np.full((len(ds.perts), K), np.nan)
        P[idx >= 0] = E[idx[idx >= 0], :K]
        return G, P


def _own_index(ds):
    _, sym, node_of_id = gene_tables()
    return np.array([node_of_id.get(sym.get(p), -1) for p in ds.perts])


def make_embedding(name: str) -> Embedding:
    parts = name.split(":")
    base = parts[0]
    if base in ("hic_cis", "hic_cis_trans", "linear_distance"):
        return GraphSpectralEmbedding(base, parts[1] if len(parts) > 1 else None, int(parts[2]) if len(parts) > 2 else 0)
    raise KeyError(name)
