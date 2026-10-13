"""Baselines. Convention everywhere: Y is (n_pert, n_gene) deltas (perturbation minus control).

An *embedding function* maps (Y_train, train_idx, ds, K) -> (G, P_all):
    G      (n_gene, K)  gene embedding
    P_all  (n_pert, L)  perturbation embedding for ALL perturbations (train and unseen)
Swap `P_all` (GO spectral, Hi-C spectral, shuffled Hi-C, random) by registering a new function in EMBEDDINGS;
LinearModel and evaluation code stay unchanged.
"""
from __future__ import annotations

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


def _train_pca_rows(Y_train, ds, K, on_expression):
    # Paper (run_linear_pretrained_model.R): prcomp on the training *expression* pseudobulks, not on deltas.
    X = Y_train + ds.ctrl_mean if on_expression else Y_train
    G = pca_gene_embedding(X, K)          # fit on training perturbations ONLY (caller passes Y_train)
    return G, G[ds.target_idx]            # P = rows of G at the perturbed genes (all perts, incl. unseen)


def train_pca_expr_rows(Y_train, train_idx, ds, K): return _train_pca_rows(Y_train, ds, K, True)
def train_pca_delta_rows(Y_train, train_idx, ds, K): return _train_pca_rows(Y_train, ds, K, False)


EMBEDDINGS = {"train_pca_expr": train_pca_expr_rows, "train_pca_delta": train_pca_delta_rows}
