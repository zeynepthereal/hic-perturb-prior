"""Perturbation embeddings (P) derived from Hi-C graphs, as `Embedding` classes plugging into the eval_set mechanism.

All classes keep G = PCA of the training expression pseudobulks (identical to `train_pca_expr`), so only P changes.
`predictable(ds)`: the perturbation target maps to a node that has an embedding. Unlike the PCA embedding the target does NOT need to
be a read-out gene -> larger `extended:<emb>` eval sets.

Graph names
  b<res>_<cis>_<trans>   bin-level whole-genome graphs: cis in {hic, lin}, trans in {trans (Hi-C), cis (= none), randtrans_s<seed> (random trans edges
                         matched to the Hi-C trans graph in count and degree sequence)}; e.g. b500_hic_trans (primary Hi-C),
                         b500_lin_randtrans_s0..s9 (primary controls, averaged over seeds in the analysis). b500x_* = rearranged regions excluded.
                         A gene gets the embedding of its TSS bin; genes in one bin share it.
  hic_local25k           gene-level 25 kb local-cis graph (band <= 10 Mb), the earlier 'hic_cis'
Suffixes
  :percomp               (cis-only graphs) top eigenvectors chosen PER CONNECTED COMPONENT, number allocated proportionally to component size
  :resid                 embedding residualized on chromosome + per-chromosome spline of position (secondary analysis), columns re-scaled
  :perm:<seed>           gene_label_permutation: same graph, which gene sits on which (embedded) node is permuted
  :rewire:<seed>         degree_preserving_rewire of the graph, re-embedded (built lazily, cached)
  :random:<seed>         random Gaussian embedding with matched dimensions and the same embedded node set
Example: 'b500_hic_trans:rewire:3', 'b500_hic_trans:resid', 'b500_lin_cis:perm:1'.
"""
from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..baselines import Embedding, pca_gene_embedding
from .bins import Genome
from .genes import HIC_CHROMS
from .position import residualize
from .spectral import percomp_matrix, spectral_embedding

ROOT = Path(__file__).resolve().parents[2]
ART, TABLES = ROOT / "data" / "hic" / "graphs", ROOT / "results" / "hic"
BASE_RE = re.compile(r"^(hic_local25k|b\d+x?_(hic|lin)_(trans|cis|randtrans_s\d+))$")


@lru_cache(maxsize=None)
def gene_tables():
    nodes = pd.read_csv(TABLES / "gene_table.csv")
    sym = pd.read_csv(TABLES / "symbol_map.csv").set_index("symbol").gene_id.to_dict()
    return nodes, sym, {g: i for i, g in enumerate(nodes.gene_id)}


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
    R = sp.csr_matrix((w, (a, b)), shape=W.shape)
    return (R + R.T).tocsr()


class GraphSpectralEmbedding(Embedding):
    def __init__(self, base: str, null: str | None = None, seed: int = 0, resid: bool = False, percomp: bool = False):
        assert BASE_RE.match(base), base
        assert not (percomp and (null or resid)), "percomp is a stand-alone variant"
        self.base, self.null, self.seed, self.resid, self.percomp = base, null, seed, resid, percomp
        self.name = base + (f":{null}:{seed}" if null else "") + (":resid" if resid else "") + (":percomp" if percomp else "")

    # ---- node level
    def _artifact(self) -> str:
        if self.null != "rewire":
            return self.base
        nm = f"{self.base}__rewire{self.seed}"
        if not (ART / f"{nm}_spectral.npz").exists():              # lazily built the first time it is needed
            W = sp.load_npz(ART / f"{self.base}.npz")
            R = double_edge_swap(W, seed=self.seed)
            sp.save_npz(ART / f"{nm}.npz", R)
            meta = load_spectral(self.base)
            r = spectral_embedding(R, k_max=50)
            np.savez_compressed(ART / f"{nm}_spectral.npz", emb=r["emb"], eigvals=r["eigvals"], all_eigvals=r["all_eigvals"],
                                embedded=r["embedded"], component=r["component"], comp_sizes=r["comp_sizes"],
                                **{k: meta[k] for k in ("node_chr", "node_pos", "valid", "res", "level") if k in meta})
        return nm

    def _node_matrix(self, K: int = 50):
        if self.percomp:
            pc = np.load(ART / f"{self.base}_percomp.npz")
            return percomp_matrix(pc["modes"], pc["component"], pc["sizes"], pc["avail"], K), pc["embedded"].astype(bool)
        z = load_spectral(self._artifact())
        E, ok = z["emb"].copy(), z["embedded"].astype(bool)
        if self.null == "random":
            E = np.where(ok[:, None], np.random.RandomState(self.seed).randn(*E.shape), np.nan)
        if self.resid:
            chr_idx, pos = self._node_coords(z)
            E = residualize(E, ok, chr_idx, pos)
        return E, ok

    def _node_coords(self, z):
        if "node_chr" in z:
            return z["node_chr"], z["node_pos"]
        nodes, _, _ = gene_tables()                                # gene-level graph: nodes are gene-table rows
        return nodes.chr.map({c: i for i, c in enumerate(HIC_CHROMS)}).values, nodes.tss.values.astype(float)

    def _own_index(self, ds):
        """dataset perturbation -> index of the node carrying its target gene (-1 if none)."""
        nodes, sym, node_of_id = gene_tables()
        gi = np.array([node_of_id.get(sym.get(p), -1) for p in ds.perts])
        z = load_spectral(self.base)
        if "res" not in z:                                         # gene-level graph
            return gi
        gb = Genome(int(z["res"])).gene_bins(nodes.chr.values, nodes.tss.values)
        return np.where(gi >= 0, gb[np.maximum(gi, 0)], -1)

    def _node_index(self, ds):
        _, ok = self._node_matrix()
        own = self._own_index(ds)
        avail = (own >= 0) & ok[np.maximum(own, 0)]
        idx = np.where(avail, own, -1)
        if self.null == "perm":      # gene_label_permutation: same graph, permute which gene sits on which (embedded) node
            emb_nodes = np.where(ok)[0]
            perm = np.random.RandomState(self.seed).permutation(len(emb_nodes))
            pos = np.zeros(len(ok), dtype=int); pos[emb_nodes] = np.arange(len(emb_nodes))
            idx = np.where(avail, emb_nodes[perm[pos[np.maximum(own, 0)]]], -1)
        return idx

    def predictable(self, ds):
        return self._node_index(ds) >= 0

    def perturbation_matrix(self, ds, K: int) -> np.ndarray:
        """P rows for every dataset perturbation (NaN rows = not predictable)."""
        E, _ = self._node_matrix(K)
        idx = self._node_index(ds)
        P = np.full((len(ds.perts), K), np.nan)
        P[idx >= 0] = E[idx[idx >= 0], :K]
        return P

    def fit(self, Y_train, train_idx, ds, K):
        G = pca_gene_embedding(Y_train + ds.ctrl_mean, K)          # identical G to `train_pca_expr`
        return G, self.perturbation_matrix(ds, K)


def make_embedding(name: str) -> Embedding:
    parts = name.split(":")
    resid, percomp = parts[-1] == "resid", parts[-1] == "percomp"
    if resid or percomp:
        parts = parts[:-1]
    base = parts[0]
    null = parts[1] if len(parts) > 1 else None
    seed = int(parts[2]) if len(parts) > 2 else 0
    if null not in (None, "perm", "rewire", "random"):
        raise KeyError(name)
    return GraphSpectralEmbedding(base, null, seed, resid, percomp)
