"""Chromosome + position design matrices, used (a) to measure how much of an embedding is explained by genomic position and
(b) to residualize embeddings on it (secondary analysis). Unsupervised: uses only node coordinates, never expression."""
from __future__ import annotations

import numpy as np
from sklearn.preprocessing import SplineTransformer

SPLINE_DF = 8        # cubic B-spline basis functions per chromosome (n_knots=6, degree=3, with bias -> 8 columns)


def position_basis(chr_idx: np.ndarray, pos: np.ndarray, n_chr: int = 23, kind: str = "spline") -> np.ndarray:
    """Rows = nodes. kind='chromosome': chromosome indicators; 'spline': chromosome-specific cubic B-spline of normalised position
    (the basis sums to 1, so it also contains the chromosome indicator); 'cubic': chromosome indicators x {1, x, x^2, x^3}."""
    chr_idx = np.asarray(chr_idx)
    pos = np.asarray(pos, dtype=float)
    oh = np.eye(n_chr)[chr_idx]
    if kind == "chromosome":
        return oh
    x = np.zeros(len(pos))
    for c in range(n_chr):
        m = chr_idx == c
        if m.any():
            lo, hi = pos[m].min(), pos[m].max()
            x[m] = (pos[m] - lo) / max(hi - lo, 1.0)
    if kind == "cubic":
        return np.hstack([oh] + [oh * ((x[:, None] - 0.5) ** p) for p in (1, 2, 3)])
    B = SplineTransformer(n_knots=SPLINE_DF - 2, degree=3, knots="uniform", extrapolation="constant", include_bias=True).fit_transform(x[:, None])
    return np.hstack([oh[:, [c]] * B for c in range(n_chr)])


def r2_columns(Y: np.ndarray, X: np.ndarray) -> np.ndarray:
    B, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return 1 - ((Y - X @ B) ** 2).sum(0) / ((Y - Y.mean(0)) ** 2).sum(0)


def residualize(E: np.ndarray, ok: np.ndarray, chr_idx, pos, rescale: bool = True) -> np.ndarray:
    """Regress each embedding column on chromosome + per-chromosome spline of position (fit on embedded nodes) and keep the
    residual; columns are re-scaled to unit RMS (so the ridge penalty acts on the same scale as for the raw embedding)."""
    X = position_basis(chr_idx, pos)
    Eo = E[ok]
    B, *_ = np.linalg.lstsq(X[ok], Eo, rcond=None)
    R = Eo - X[ok] @ B
    if rescale:
        R = R / np.maximum(np.sqrt((R ** 2).mean(0, keepdims=True)), 1e-12)
    out = np.full_like(E, np.nan)
    out[ok] = R
    return out
