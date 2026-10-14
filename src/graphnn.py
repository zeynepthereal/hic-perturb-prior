"""Graph-neighbour model (Study 2, family 2).

Predict the delta (perturbation minus control) of a test perturbation as the weighted average of the deltas of TRAINING perturbations whose target genes are neighbours
of the test target on a gene/bin graph; fall back to the mean baseline (mean training delta) when no training perturbation is reachable.

Neighbour weights (k-hop, edge-weighted):  S = w0 * I  +  sum_{h=1..hops} gamma^(h-1) * P^h,   P = D^-1 W  (row-stochastic transition matrix of the weighted graph)
  S[node_p, node_q] is the weight of training perturbation q for test perturbation p; I is the "same node" term (two genes in the same bin: w0 = same_node_weight).
  prediction_p = sum_q S[node_p, node_q] * delta_q / sum_q S[node_p, node_q]      if the denominator is > 0, else the mean training delta.
No parameters are fitted: (hops, gamma, same_node_weight) are fixed by the user. A perturbation whose target has no node (or an isolated node) is not graph-predictable:
  no node -> NaN row (not predictable, as in the eval_set mechanism); isolated node -> falls back to the mean baseline unless a training perturbation shares its node.
Interface: `fit_predict_graph_neighbour(ds, train_idx, test_idx, source)` mirrors `evaluate.fit_predict_linear` (NaN rows for unpredictable perturbations)."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp


class GraphNeighbourModel:
    def __init__(self, W: sp.spmatrix, hops: int = 2, gamma: float = 0.5, same_node_weight: float = 1.0):
        W = sp.csr_matrix(W, dtype=np.float64)
        n = W.shape[0]
        deg = np.asarray(W.sum(1)).ravel()
        P = sp.diags(np.where(deg > 0, 1.0 / np.maximum(deg, 1e-300), 0.0)) @ W            # isolated nodes: zero rows
        S = same_node_weight * sp.identity(n, format="csr")
        Ph = sp.identity(n, format="csr")
        for h in range(1, hops + 1):
            Ph = Ph @ P
            S = S + (gamma ** (h - 1)) * Ph
        self.S, self.n_nodes, self.hops, self.gamma, self.w0 = S.tocsr(), n, hops, gamma, same_node_weight

    def fit(self, delta_train: np.ndarray, node_train: np.ndarray):
        keep = node_train >= 0
        self.D, self.node_tr = np.asarray(delta_train)[keep], np.asarray(node_train)[keep]
        self.mean = np.asarray(delta_train).mean(axis=0)                                  # mean baseline over ALL training perturbations
        return self

    def weights(self, node_test: np.ndarray) -> np.ndarray:
        """(n_test, n_train) neighbour weights; zero rows for test perturbations without a node."""
        node_test = np.asarray(node_test)
        out = np.zeros((len(node_test), len(self.node_tr)))
        ok = node_test >= 0
        if ok.any():
            out[ok] = self.S[node_test[ok]][:, self.node_tr].toarray()
        return out

    def predict(self, node_test: np.ndarray) -> np.ndarray:
        """(n_test, n_gene) predicted deltas; NaN rows where the target has no node."""
        node_test = np.asarray(node_test)
        Wt = self.weights(node_test)
        mass = Wt.sum(axis=1)
        pred = np.tile(self.mean, (len(node_test), 1))
        has = mass > 0
        if has.any():
            pred[has] = (Wt[has] @ self.D) / mass[has, None]
        pred[node_test < 0] = np.nan
        self.last_mass_, self.last_fallback_ = mass, (~has) & (node_test >= 0)
        return pred


class GraphSource:
    """A named information source for the graph-neighbour model: a graph over nodes + the mapping perturbation -> node (-1 = none)."""
    def __init__(self, W: sp.spmatrix, node_of_pert: np.ndarray):
        self.W, self.node_of_pert = sp.csr_matrix(W), np.asarray(node_of_pert)

    def predictable(self, ds=None) -> np.ndarray:
        return self.node_of_pert >= 0

    @classmethod
    def from_hic_graph(cls, name: str, ds):
        """Hi-C / control graphs built in Study 1 (e.g. 'b500_hic_trans', 'b500_lin_randtrans_s0', 'hic_local25k'); target -> TSS bin (or gene) node;
        targets whose node has no edges at all are mapped to -1."""
        import scipy.sparse as sp_
        from .hic.embeddings import ART, GraphSpectralEmbedding, load_spectral
        emb = GraphSpectralEmbedding(name)
        own = emb._own_index(ds)
        W = sp_.load_npz(ART / f"{name}.npz")
        embedded = load_spectral(name)["embedded"].astype(bool)
        return cls(W, np.where((own >= 0) & embedded[np.maximum(own, 0)], own, -1))


def fit_predict_graph_neighbour(ds, train_idx, test_idx, source: GraphSource, hops: int = 2, gamma: float = 0.5, same_node_weight: float = 1.0) -> np.ndarray:
    """Mirror of evaluate.fit_predict_linear: (len(test_idx), n_gene) with NaN rows for perturbations the source cannot place on the graph."""
    m = GraphNeighbourModel(source.W, hops, gamma, same_node_weight).fit(ds.delta[train_idx], source.node_of_pert[train_idx])
    return m.predict(source.node_of_pert[test_idx])
