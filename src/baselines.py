"""Baselines. Convention everywhere: Y is (n_pert, n_gene) deltas (perturbation minus control).

An `Embedding` (see below) maps (Y_train, train_idx, ds, K) -> (G, P_all):
    G      (n_gene, K)  gene embedding
    P_all  (n_pert, L)  perturbation embedding for ALL perturbations (train and unseen)
Swap `P_all` (GO spectral, Hi-C spectral, shuffled Hi-C, random) by registering a new function in EMBEDDINGS;
LinearModel and evaluation code stay unchanged.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


class NoChange:
    def fit(self, Y_train): self.n_genes = Y_train.shape[1]; return self
    def predict(self, n): return np.zeros((n, self.n_genes))


class MeanBaseline:
    """Predict the mean training delta for every unseen perturbation (the paper's `b`)."""
    def fit(self, Y_train): self.b = Y_train.mean(axis=0); return self
    def predict(self, n): return np.tile(self.b, (n, 1))


class LinearModel:
    """Y_hat = G W P~^T + b  (paper's notation, transposed to pert x gene).

    W = (G'G + lam I)^-1 G' (Y - b) P (P'P + lam I)^-1, with Y (gene x pert) and b = row means of Y_train.
    """
    def __init__(self, lam: float = 0.1): self.lam = lam

    def fit(self, Y_train, G, P_train):
        Yt = Y_train.T                                   # (gene, pert)
        self.G, self.b = G, Yt.mean(axis=1, keepdims=True)
        K, L = G.shape[1], P_train.shape[1]
        left = np.linalg.solve(G.T @ G + self.lam * np.eye(K), G.T @ (Yt - self.b))      # (K, pert)
        self.W = np.linalg.solve((P_train.T @ P_train + self.lam * np.eye(L)).T, (left @ P_train).T).T  # (K, L)
        return self

    def predict(self, P_test):
        return (self.G @ self.W @ P_test.T + self.b).T   # (pert_test, gene)


def pca_gene_embedding(Y_train, K):
    """Top-K PCs (scores) of the gene x pert matrix, columns (perturbations) centred, as R's prcomp()."""
    Yt = Y_train.T
    U, S, _ = np.linalg.svd(Yt - Yt.mean(axis=0, keepdims=True), full_matrices=False)
    return U[:, :K] * S[:K]


@dataclass
class Embedding:
    """A (G, P) provider. `fit` -> (G (n_gene, K), P_all (n_pert, L)) with NaN rows for perturbations it cannot embed;
    `predictable(ds)` -> bool (n_pert,) saying which perturbations it can embed (independent of the split).
    Register new ones (GO spectral, Hi-C spectral, shuffled Hi-C, random) in EMBEDDINGS."""
    fit: Callable
    predictable: Callable


def _train_pca_rows(Y_train, ds, K, on_expression):
    # Paper (run_linear_pretrained_model.R): prcomp on the training *expression* pseudobulks, not on deltas.
    X = Y_train + ds.ctrl_mean if on_expression else Y_train
    G = pca_gene_embedding(X, K)          # fit on training perturbations ONLY (caller passes Y_train)
    P = np.full((len(ds.perts), K), np.nan)
    ok = ds.target_idx >= 0
    P[ok] = G[ds.target_idx[ok]]          # P = rows of G at the perturbed genes
    return G, P


def _pca_ok(ds): return ds.target_idx >= 0


EMBEDDINGS = {
    "train_pca_expr": Embedding(lambda Y, tr, ds, K: _train_pca_rows(Y, ds, K, True), _pca_ok),
    "train_pca_delta": Embedding(lambda Y, tr, ds, K: _train_pca_rows(Y, ds, K, False), _pca_ok),
}
