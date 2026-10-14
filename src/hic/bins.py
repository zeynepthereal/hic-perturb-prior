"""Bin-level whole-genome Hi-C graphs (nodes = genomic bins, e.g. 500 kb / 250 kb, chr1-22,X).

Normalisation  : genome-wide ICE on the raw map (all cis bins up to the FULL chromosome length + all trans), cooler-style low-coverage filter.
Observed/exp.  : cis  = balanced / mean balanced value of the same diagonal (distance-stratified, per chromosome; diagonals with fewer
                 than `min_pairs` valid pairs use the farthest well-sampled diagonal's expectation);
                 trans = balanced / global mean balanced trans value.
Sparsification : the SAME rule for every graph (edge classes are treated separately):
                 cis   class: kNN edges (k_cis strongest partners per bin, union of both directions, each bin keeping at most deg_cap
                              incident edges = 'capped_union', which limits hub bins) U the maximum spanning forest per chromosome
                              (guarantees exactly one cis component per chromosome with data),
                 trans class: the same kNN rule with k_trans strongest trans partners (no spanning forest).
                 Union, symmetric (max). Each class's weights are clipped at the clip_q percentile and rescaled to mean 1 so that the
                 relative cis/trans weighting does not depend on the weight source.
Cis weight source: 'hic' = cis O/E;  'lin' = exp(-d/d0) of the 1D distance between bin centres (a chain along each chromosome).
Trans          : 'hic' = trans O/E partners (identical edges in hic_trans and lin_trans);  'none'.
Factorial cells: {hic,lin} x {trans,none}  ->  hic_trans (primary Hi-C), hic_cis, lin_trans, lin_cis (primary control).
All four cells share the same node set (bins that pass the Hi-C coverage filter), so predictable perturbation sets are identical."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import minimum_spanning_tree

from .genes import CHROM_LEN, HIC_CHROMS
from .graphs import valid_bins

BIN_CFG = dict(k_cis=10, k_trans=10, clip_q=0.999, min_pairs=20, ice_iters=200, ice_tol=1e-5, d0_lin=5_000_000, knn_mode="capped_union", deg_cap=20)
NCH = len(HIC_CHROMS)


class Genome:
    def __init__(self, res: int):
        self.res = res
        self.nb = np.array([-(-CHROM_LEN[c] // res) for c in HIC_CHROMS])
        self.off = np.concatenate([[0], np.cumsum(self.nb)])
        self.n = int(self.off[-1])
        self.chr_idx = np.repeat(np.arange(NCH), self.nb)
        self.local = np.concatenate([np.arange(k) for k in self.nb])
        self.pos = (self.local + 0.5) * res

    def sl(self, c: int) -> slice:
        return slice(int(self.off[c]), int(self.off[c + 1]))

    def gene_bins(self, chr_names, tss) -> np.ndarray:
        """global bin index of each gene's TSS bin (-1 for chromosomes without Hi-C)."""
        ci = np.array([HIC_CHROMS.index(c) if c in HIC_CHROMS else -1 for c in chr_names])
        return np.where(ci >= 0, self.off[np.maximum(ci, 0)] + (np.asarray(tss) - 1) // self.res, -1)


def load_raw(G: Genome, pairs_dir: str | Path) -> np.ndarray:
    """Assemble the genome-wide raw (float32) matrix from the per-chromosome-pair blocks written by scripts/fetch_hic.py."""
    R = np.zeros((G.n, G.n), dtype=np.float32)
    for ia, a in enumerate(HIC_CHROMS):
        for ib in range(ia, NCH):
            b = HIC_CHROMS[ib]
            M = np.load(Path(pairs_dir) / f"{a}__{b}.npy")
            if ia == ib:
                R[G.sl(ia), G.sl(ia)] = np.triu(M) + np.triu(M, 1).T
            else:
                R[G.sl(ia), G.sl(ib)] = M
                R[G.sl(ib), G.sl(ia)] = M.T
    return R


def ice(R: np.ndarray, iters: int, tol: float):
    """Genome-wide ICE. Modifies R in place (invalid rows/cols zeroed). Returns (x = 1/bias, valid)."""
    valid = valid_bins(R.sum(1, dtype=np.float64))
    R[~valid, :] = 0
    R[:, ~valid] = 0
    b = np.ones(len(valid))
    for _ in range(iters):
        x = 1.0 / b
        r = (R @ x.astype(np.float32)).astype(np.float64) * x
        s = r / r[valid].mean()
        s[~valid | (s == 0)] = 1.0
        b *= s
        if np.abs(s[valid] - 1).max() < tol:
            break
    return 1.0 / b, valid


def oe_matrix(R: np.ndarray, x: np.ndarray, valid: np.ndarray, G: Genome, min_pairs: int):
    """Returns (OE float32 (n, n) with 0 for invalid bins and the diagonal, m_trans, per-chromosome expected arrays)."""
    x32 = x.astype(np.float32)
    OE = np.zeros_like(R)
    s, c = 0.0, 0
    for a in range(NCH):
        for b in range(a + 1, NCH):
            blk = R[G.sl(a), G.sl(b)] * x32[G.sl(a), None] * x32[None, G.sl(b)]
            vm = valid[G.sl(a), None] & valid[None, G.sl(b)]
            s += float(blk[vm].sum(dtype=np.float64)); c += int(vm.sum())
    m_trans = s / c
    for a in range(NCH):
        for b in range(a + 1, NCH):
            blk = R[G.sl(a), G.sl(b)] * x32[G.sl(a), None] * x32[None, G.sl(b)] / np.float32(m_trans)
            vm = valid[G.sl(a), None] & valid[None, G.sl(b)]
            blk = np.where(vm, blk, 0)
            OE[G.sl(a), G.sl(b)] = blk
            OE[G.sl(b), G.sl(a)] = blk.T
    expected = []
    for a in range(NCH):
        sa = G.sl(a); n = sa.stop - sa.start
        bal = R[sa, sa] * x32[sa, None] * x32[None, sa]
        v = valid[sa]
        ex, npairs = np.full(n, np.nan), np.zeros(n, dtype=int)
        for d in range(n):
            vm = v[: n - d] & v[d:]
            npairs[d] = vm.sum()
            if npairs[d] > 0:
                ex[d] = np.diagonal(bal, d)[vm].mean(dtype=np.float64)
        good = np.where(npairs >= min_pairs)[0]
        last = good.max() if len(good) else 0
        ex[last + 1:] = ex[last]                                          # far diagonals: plateau of the last well-sampled one
        ex = pd.Series(ex).ffill().bfill().values
        i = np.arange(n)
        D = np.abs(i[:, None] - i[None, :])
        o = bal / ex[D].astype(np.float32)
        o[~(v[:, None] & v[None, :])] = 0
        np.fill_diagonal(o, 0)
        OE[sa, sa] = o
        expected.append(ex)
    return OE, m_trans, expected


# ------------------------------------------------------------------ compartments
def ab_compartments(OE: np.ndarray, valid: np.ndarray, G: Genome, density: np.ndarray):
    """A/B compartment track per bin: PC1 of the Pearson correlation matrix of the cis O/E matrix, per chromosome, oriented so that
    it correlates positively with gene density (A = gene-rich), z-scored per chromosome. Also reports PC2/PC3 correlations with density
    and the best-of-3 PC (the usual remedy when PC1 captures the chromosome arms instead of compartments)."""
    ab = np.full(G.n, np.nan)
    ab_best = np.full(G.n, np.nan)
    rows = []
    for c in range(NCH):
        sa = G.sl(c)
        v = valid[sa]
        if v.sum() < 20:
            continue
        M = OE[sa, sa][np.ix_(v, v)].astype(np.float64)
        np.fill_diagonal(M, 1.0)                  # O/E of a bin with itself is its expectation
        C = np.corrcoef(M)
        C = np.nan_to_num(C)
        w, U = np.linalg.eigh(C)
        U = U[:, ::-1][:, :3]
        dens = density[sa][v]
        cors = [np.corrcoef(U[:, k], dens)[0, 1] for k in range(3)]
        def orient(k):
            u = U[:, k] * (1 if cors[k] >= 0 else -1)
            return (u - u.mean()) / u.std()
        t = np.full(len(v), np.nan); t[v] = orient(0); ab[sa] = t
        kb = int(np.argmax(np.abs(cors)))
        t = np.full(len(v), np.nan); t[v] = orient(kb); ab_best[sa] = t
        rows.append(dict(chrom=HIC_CHROMS[c], n_valid_bins=int(v.sum()), corr_pc1_density=cors[0], corr_pc2_density=cors[1],
                         corr_pc3_density=cors[2], best_pc=kb + 1))
    return ab, ab_best, pd.DataFrame(rows)


# ------------------------------------------------------------------ graph cells
def _topk_mask(W: np.ndarray, k: int) -> np.ndarray:
    """Directed top-k mask per row of a dense weight matrix (entries <= 0 are not edges)."""
    n = W.shape[0]
    k = min(k, n - 1)
    m = np.zeros(W.shape, dtype=bool)
    if k <= 0:
        return m
    idx = np.argpartition(-W, k - 1, axis=1)[:, :k]
    np.put_along_axis(m, idx, True, axis=1)
    return m & (W > 0)


def _msf_mask(W: np.ndarray) -> np.ndarray:
    """Maximum spanning forest of a dense non-negative weight matrix (0 = no edge) as a boolean mask."""
    m = np.zeros(W.shape, dtype=bool)
    if not (W > 0).any():
        return m
    cost = np.where(W > 0, W.max() + 1.0 - W, 0.0)
    T = minimum_spanning_tree(sp.csr_matrix(cost)).tocoo()
    m[T.row, T.col] = True
    return m


def _select(W: np.ndarray, k: int, mode: str, cap: int) -> np.ndarray:
    """Symmetric boolean edge mask of the kNN rule. 'union': kNN in either direction; 'mutual': in both; 'capped_union': union, then
    every node keeps at most `cap` incident edges (its strongest by weight) and an edge survives only if BOTH endpoints keep it
    (limits hub bins without making the graph sparser for ordinary bins)."""
    K = _topk_mask(W, k)
    if mode == "mutual":
        return K & K.T
    U = K | K.T
    if mode == "union":
        return U
    deg = U.sum(1)
    keep = np.ones_like(U)
    for i in np.where(deg > cap)[0]:
        cand = np.where(U[i])[0]
        drop = cand[np.argsort(-W[i, cand], kind="stable")[cap:]]
        keep[i, drop] = False
    return U & keep & keep.T


def _to_graph(rows, cols, vals, n, q) -> sp.csr_matrix:
    """symmetric graph from undirected edge lists; weights clipped at the q-th percentile and rescaled to mean 1."""
    rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
    vals = np.minimum(vals, np.quantile(vals, q))
    vals = vals / vals.mean()
    A = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    return A.maximum(A.T).tocsr()


def build_cells(OE: np.ndarray, valid: np.ndarray, G: Genome, cfg=BIN_CFG) -> dict:
    """Returns {'hic_trans', 'hic_cis', 'lin_trans', 'lin_cis'} symmetric csr graphs over all G.n bins (invalid bins isolated)."""
    n = G.n
    out = {}
    for src in ("hic", "lin"):
        R_, C_, V_ = [], [], []
        for c in range(NCH):
            sa = G.sl(c)
            v = np.where(valid[sa])[0]
            if len(v) < 3:
                continue
            pos = G.pos[sa][v]
            if src == "hic":
                W = OE[sa, sa][np.ix_(v, v)].astype(np.float64)
            else:
                W = np.exp(-np.abs(pos[:, None] - pos[None, :]) / cfg["d0_lin"])
                np.fill_diagonal(W, 0)
            M = _select(W, cfg["k_cis"], cfg["knn_mode"], cfg["deg_cap"])
            F = _msf_mask(W); M = M | F | F.T          # spanning forest is exempt from the cap (guarantees one cis component)
            r, c2 = np.where(np.triu(M, 1))
            R_.append(v[r] + sa.start); C_.append(v[c2] + sa.start); V_.append(W[r, c2])
        out[src] = _to_graph(R_, C_, V_, n, cfg["clip_q"])
    # trans class: top-k_trans partners among the other chromosomes (identical for hic_trans and lin_trans)
    Wt = OE.copy()
    for c in range(NCH):
        Wt[G.sl(c), G.sl(c)] = 0
    M = _select(Wt, cfg["k_trans"], cfg["knn_mode"], cfg["deg_cap"])
    r, c2 = np.where(np.triu(M, 1))
    T = _to_graph([r], [c2], [Wt[r, c2].astype(np.float64)], n, cfg["clip_q"])
    del Wt
    return {"hic_trans": (out["hic"] + T).tocsr(), "hic_cis": out["hic"], "lin_trans": (out["lin"] + T).tocsr(), "lin_cis": out["lin"],
            "_trans_only": T}
