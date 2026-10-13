"""Gene x gene graphs from Hi-C (all parameters in GRAPH_CFG, documented in README).

Pipeline per variant
  raw counts (remote, io.py) -> ICE balancing (own implementation, KR vectors are absent in the file) -> observed/expected
  -> gene-level weights via TSS bins -> kNN sparsification (same rule for every variant) -> symmetric weighted graph.

Sparsification rule (identical for all variants): every node keeps its `k_cis` strongest cis edges (cis edges only exist within
`max_edge_dist` and need >= `min_raw` raw reads); variant hic_cis_trans additionally lets every node keep its `k_trans`
strongest trans partners (distinct 500 kb partner bins, one gene drawn at random per bin because trans weights are only
defined at bin resolution). The graph is the union of the kept edges (symmetric; weight = O/E, clipped at the q-th percentile).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .genes import HIC_CHROMS

GRAPH_CFG = dict(cis_res=25_000, trans_res=500_000, max_fetch_dist=10_000_000, max_edge_dist=5_000_000, min_raw=3,
                 k_cis=10, k_trans=3, clip_q=0.999, ice_iters=200, ice_tol=1e-5, seed=0)


# ------------------------------------------------------------------ balancing
def valid_bins(cov: np.ndarray, mad_max: float = 5.0) -> np.ndarray:
    """cooler-style low-coverage filter: drop empty bins and bins with log-coverage < median - mad_max*MAD."""
    nz = cov > 0
    if nz.sum() == 0:
        return nz
    lc = np.log(cov[nz])
    med = np.median(lc)
    mad = np.median(np.abs(lc - med))
    return nz & (np.log(np.where(nz, cov, 1)) > med - mad_max * 1.4826 * mad)


def ice_coo(i, j, c, n, iters=200, tol=1e-5):
    """Iterative correction on a symmetric matrix given as upper-triangle COO (i<=j). Returns (bias, valid)."""
    off = i != j
    ii, jj, cc = np.concatenate([i, j[off]]), np.concatenate([j, i[off]]), np.concatenate([c, c[off]]).astype(np.float64)
    cov = np.bincount(ii, weights=cc, minlength=n)
    valid = valid_bins(cov)
    keep = valid[ii] & valid[jj]
    ii, jj, cc = ii[keep], jj[keep], cc[keep]
    b = np.ones(n)
    for _ in range(iters):
        r = np.bincount(ii, weights=cc / (b[ii] * b[jj]), minlength=n)
        s = r / r[valid].mean()
        s[~valid | (s == 0)] = 1.0
        b *= s
        if np.abs(s[valid] - 1).max() < tol:
            break
    b = b / np.sqrt(np.mean(b[valid] ** 2)) if valid.any() else b
    return b, valid


# ------------------------------------------------------------------ cis observed/expected
def cis_oe_bins(npz_path: str | Path, cfg=GRAPH_CFG):
    """O/E for one chromosome. Returns dict(i, j, oe, raw, n_bins, valid) with i<=j, distance<=max_edge_dist, raw>=min_raw."""
    z = np.load(npz_path)
    i, j, c, n = z["i"].astype(np.int64), z["j"].astype(np.int64), z["c"].astype(np.float64), int(z["n_bins"])
    b, valid = ice_coo(i, j, c, n, cfg["ice_iters"], cfg["ice_tol"])
    ok = valid[i] & valid[j]
    i, j, c = i[ok], j[ok], c[ok]
    bal = c / (b[i] * b[j])
    d = j - i
    maxd = int(d.max()) if len(d) else 0
    # expected(d): mean balanced value over all valid bin pairs at distance d (zeros included)
    s = np.bincount(d, weights=bal, minlength=maxd + 1)
    nv = np.array([np.sum(valid[: n - k] & valid[k:]) for k in range(maxd + 1)], dtype=np.float64)
    expd = np.where(nv > 0, s / np.maximum(nv, 1), np.nan)
    oe = bal / expd[d]
    keep = (c >= cfg["min_raw"]) & (d * cfg["cis_res"] <= cfg["max_edge_dist"]) & np.isfinite(oe)
    return dict(i=i[keep], j=j[keep], oe=oe[keep], raw=c[keep], n_bins=n, valid=valid, expected=expd)


def topk_rows(W: sp.csr_matrix, k_per_row) -> sp.csr_matrix:
    """Keep the k_per_row[r] largest entries of every row (directed selection)."""
    W = W.tocsr()
    rows, cols, vals = [], [], []
    for r in range(W.shape[0]):
        a, b = W.indptr[r], W.indptr[r + 1]
        if b == a:
            continue
        k = int(k_per_row if np.isscalar(k_per_row) else k_per_row[r])
        if k <= 0:
            continue
        v = W.data[a:b]
        idx = np.argsort(-v, kind="stable")[:k]
        rows.append(np.full(len(idx), r)); cols.append(W.indices[a:b][idx]); vals.append(v[idx])
    if not rows:
        return sp.csr_matrix(W.shape)
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=W.shape)


def symmetrize_union(D: sp.csr_matrix) -> sp.csr_matrix:
    """Union of directed kNN edges; weight = max of the two directions (equal for O/E weights, kept for generality)."""
    D = D.tocsr()
    S = D.maximum(D.T)
    S.setdiag(0)
    S.eliminate_zeros()
    return S.tocsr()


def gene_cis_candidates(nodes: pd.DataFrame, cis_dir: str | Path, cfg=GRAPH_CFG) -> sp.csr_matrix:
    """Gene x gene cis O/E candidate weights (all pairs within max_edge_dist with >= min_raw reads). Shape (n_nodes, n_nodes)."""
    n = len(nodes)
    blocks = []
    for chrom in HIC_CHROMS:
        idx = np.where(nodes.chr.values == chrom)[0]
        if len(idx) == 0:
            continue
        o = cis_oe_bins(Path(cis_dir) / f"{chrom}.npz", cfg)
        nb = o["n_bins"]
        Wb = sp.coo_matrix((o["oe"], (o["i"], o["j"])), shape=(nb, nb)).tocsr()
        Wb = Wb + sp.triu(Wb, 1).T                      # symmetric (diagonal kept once)
        A = sp.csr_matrix((np.ones(len(idx)), (np.arange(len(idx)), nodes.bin25k.values[idx])), shape=(len(idx), nb))
        Wg = (A @ Wb @ A.T).tocoo()
        m = Wg.row != Wg.col
        blocks.append(sp.coo_matrix((Wg.data[m], (idx[Wg.row[m]], idx[Wg.col[m]])), shape=(n, n)))
    return sum(blocks[1:], blocks[0]).tocsr()


# ------------------------------------------------------------------ trans (500 kb, genome-wide)
def load_genome_matrix(pairs_dir: str | Path, chrom_bins: dict):
    """Assemble the genome-wide raw 500 kb matrix from the per-pair blocks. Returns (R, offsets)."""
    offs, o = {}, 0
    for c in HIC_CHROMS:
        offs[c] = o; o += chrom_bins[c]
    R = np.zeros((o, o), dtype=np.float64)
    for ia, a in enumerate(HIC_CHROMS):
        for b in HIC_CHROMS[ia:]:
            M = np.load(Path(pairs_dir) / f"{a}__{b}.npy")
            R[offs[a]: offs[a] + M.shape[0], offs[b]: offs[b] + M.shape[1]] = M
            if a != b:
                R[offs[b]: offs[b] + M.shape[1], offs[a]: offs[a] + M.shape[0]] = M.T
            else:
                R[offs[a]: offs[a] + M.shape[0], offs[a]: offs[a] + M.shape[0]] = np.triu(M) + np.triu(M, 1).T
    return R, offs


def trans_oe(R: np.ndarray, offs: dict, chrom_bins: dict, iters=200, tol=1e-5):
    """Genome-wide ICE on the dense 500 kb matrix, then trans O/E = balanced / mean balanced of that chromosome pair."""
    n = R.shape[0]
    cov = R.sum(1)
    valid = valid_bins(cov)
    Rv = R * (valid[:, None] & valid[None, :])
    b = np.ones(n)
    for _ in range(iters):
        r = (Rv / np.outer(b, b)).sum(1)
        s = r / r[valid].mean()
        s[~valid | (s == 0)] = 1.0
        b *= s
        if np.abs(s[valid] - 1).max() < tol:
            break
    bal = Rv / np.outer(b, b)
    OE = np.zeros_like(bal)
    for a in HIC_CHROMS:
        sa = slice(offs[a], offs[a] + chrom_bins[a])
        for c in HIC_CHROMS:
            if a == c:
                continue
            sc_ = slice(offs[c], offs[c] + chrom_bins[c])
            blk = bal[sa, sc_]
            vm = np.outer(valid[sa], valid[sc_])
            e = blk[vm].mean() if vm.any() else np.nan
            if e and np.isfinite(e) and e > 0:
                OE[sa, sc_] = blk / e
    return OE, valid


def gene_trans_edges(nodes: pd.DataFrame, OE: np.ndarray, offs: dict, k_trans: int, seed: int) -> sp.csr_matrix:
    """Directed trans kNN: every gene keeps its k_trans best partner bins (other chromosomes, O/E > 0); in each partner bin one
    gene is drawn at random (seeded). Weight = bin-pair trans O/E."""
    rng = np.random.RandomState(seed)
    n = len(nodes)
    gbin = np.array([offs[c] for c in nodes.chr.values]) + nodes.bin500k.values
    members = {}
    for g, b in enumerate(gbin):
        members.setdefault(b, []).append(g)
    has_gene = np.zeros(OE.shape[0], dtype=bool)
    has_gene[list(members)] = True
    rows, cols, vals = [], [], []
    for b in np.unique(gbin):
        w = np.where(has_gene, OE[b], 0.0)
        top = np.argsort(-w)[:k_trans]
        top = top[w[top] > 0]
        for g in members[b]:
            for t in top:
                h = members[t][rng.randint(len(members[t]))]
                rows.append(g); cols.append(h); vals.append(w[t])
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


# ------------------------------------------------------------------ variants
def build_variants(nodes: pd.DataFrame, cis_dir, pairs_dir, chrom_bins: dict, cfg=GRAPH_CFG) -> dict:
    """Returns {'hic_cis', 'hic_cis_trans', 'linear_distance'} -> symmetric csr graphs (rows/cols follow `nodes`) + extras."""
    cand = gene_cis_candidates(nodes, cis_dir, cfg)
    k = cfg["k_cis"]
    D_cis = topk_rows(cand, k)
    cap = np.quantile(D_cis.data, cfg["clip_q"])
    D_cis.data = np.minimum(D_cis.data, cap)
    G_cis = symmetrize_union(D_cis)

    R, offs = load_genome_matrix(pairs_dir, chrom_bins)
    OE, _ = trans_oe(R, offs, chrom_bins, cfg["ice_iters"], cfg["ice_tol"])
    D_tr = gene_trans_edges(nodes, OE, offs, cfg["k_trans"], cfg["seed"])
    D_tr.data = np.minimum(D_tr.data, np.quantile(D_tr.data, cfg["clip_q"]))
    G_ct = symmetrize_union(D_cis + D_tr)

    # linear-distance control: same per-node out-degree as the hic_cis directed selection, nearest genes by TSS distance on the
    # same chromosome and inside the same candidate window; weight exp(-d/d0), d0 = median TSS distance of the retained hic_cis edges.
    kdeg = np.diff(D_cis.tocsr().indptr)
    tss = nodes.tss.values.astype(np.float64)
    chr_ = nodes.chr.values
    cc = D_cis.tocoo()
    d0 = float(np.median(np.abs(tss[cc.row] - tss[cc.col])))
    rows, cols, vals = [], [], []
    for c in HIC_CHROMS:
        idx = np.where(chr_ == c)[0]                       # nodes are sorted by (chr, tss)
        t = tss[idx]
        for a, g in enumerate(idx):
            kk = int(kdeg[g])
            if kk == 0:
                continue
            lo, hi = max(0, a - kk - 1), min(len(idx), a + kk + 2)
            cand_i = [x for x in range(lo, hi) if x != a]
            dist = np.abs(t[cand_i] - t[a])
            ok = dist <= cfg["max_edge_dist"]
            order = np.argsort(dist[ok], kind="stable")[:kk]
            for o in order:
                h = idx[np.array(cand_i)[ok][o]]
                rows.append(g); cols.append(h); vals.append(float(np.exp(-dist[ok][o] / d0)))
    G_lin = symmetrize_union(sp.csr_matrix((vals, (rows, cols)), shape=D_cis.shape))
    return {"hic_cis": G_cis, "hic_cis_trans": G_ct, "linear_distance": G_lin,
            "_extras": {"d0": d0, "cis_clip": float(cap), "trans_OE": None}}
