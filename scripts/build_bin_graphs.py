"""python -m scripts.build_bin_graphs <res_bp> [--randtrans N] [--only-randtrans] [--exclude BED --tag x]
   e.g.  500000 | 250000 | 500000 --randtrans 10 --only-randtrans | 500000 --exclude configs/k562_rearranged_regions_hg38.bed --tag x
Needs data/hic/pairs_<res>kb (scripts/fetch_hic.py trans <res>) and results/hic/gene_table.csv.
Writes data/hic/graphs/b<res>k<tag>_<cell>{.npz,_spectral.npz}; (untagged) results/hic_v2/ab_<res>k.csv, ab_<res>k_chrom.csv; bin<res>k<tag>_stats.json
Cells: {hic,lin}_{trans,cis}; with --randtrans N also {hic,lin}_randtrans_s<seed> (matched random trans edges, seeds 0..N-1)."""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from src.hic.bins import BIN_CFG, Genome, ab_compartments, build_cells, ice, load_raw, oe_matrix, rewire_trans
from src.hic.genes import HIC_CHROMS
from src.hic.spectral import spectral_embedding

ap = argparse.ArgumentParser()
ap.add_argument("res", type=int); ap.add_argument("--randtrans", type=int, default=0); ap.add_argument("--only-randtrans", action="store_true")
ap.add_argument("--exclude", default=None); ap.add_argument("--tag", default="")
a = ap.parse_args()
res = a.res; tag = f"b{res // 1000}{a.tag}"
ART = Path("data/hic/graphs"); ART.mkdir(parents=True, exist_ok=True)
OUT = Path("results/hic_v2"); OUT.mkdir(parents=True, exist_ok=True)
t0 = time.time(); log = lambda m: print(f"[{time.time()-t0:5.0f}s] {m}", flush=True)

G = Genome(res)
R = load_raw(G, f"data/hic/pairs_{res // 1000}kb"); log(f"raw matrix {R.shape}")
n_excl = 0
if a.exclude:                                   # bins overlapping the exclusion BED are removed BEFORE balancing
    mask = np.zeros(G.n, dtype=bool)
    for line in open(a.exclude):
        if line.startswith("#") or not line.strip():
            continue
        c, s, e = line.split("\t")[:3]
        if c in HIC_CHROMS:
            ci = HIC_CHROMS.index(c)
            lo, hi = G.off[ci] + int(s) // res, G.off[ci] + (int(e) - 1) // res
            mask[lo: min(hi, G.off[ci + 1] - 1) + 1] = True
    R[mask, :] = 0; R[:, mask] = 0; n_excl = int(mask.sum()); log(f"excluded {n_excl} bins ({100*n_excl/G.n:.1f}%) from {a.exclude}")
x, valid = ice(R, BIN_CFG["ice_iters"], BIN_CFG["ice_tol"]); log(f"ICE done; valid bins {valid.sum()}/{G.n}")
OE, m_trans, expected = oe_matrix(R, x, valid, G, BIN_CFG["min_pairs"]); del R; log(f"O/E done (global trans mean {m_trans:.3e})")

nodes = pd.read_csv("results/hic/gene_table.csv")
if not a.tag and not a.only_randtrans:
    gb = G.gene_bins(nodes.chr.values, nodes.tss.values)
    density = np.bincount(gb[gb >= 0], minlength=G.n).astype(float)
    ab, ab_best, ab_tab = ab_compartments(OE, valid, G, density)
    pd.DataFrame(dict(chr=[HIC_CHROMS[i] for i in G.chr_idx], start=(G.local * res), valid=valid, gene_density=density, ab_pc1=ab,
                      ab_best_of_3=ab_best)).to_csv(OUT / f"ab_{res // 1000}k.csv", index=False)
    ab_tab.to_csv(OUT / f"ab_{res // 1000}k_chrom.csv", index=False)
    log(f"A/B: median corr(PC1, gene density) = {ab_tab.corr_pc1_density.median():.2f}; PC1 best of 3 on {(ab_tab.best_pc == 1).sum()}/{len(ab_tab)} chromosomes")

cells = build_cells(OE, valid, G); T = cells.pop("_trans_only"); log("cells built")
stats = {"res": res, "tag": a.tag, "cfg": BIN_CFG, "n_bins": G.n, "valid_bins": int(valid.sum()), "excluded_bins": n_excl, "m_trans": float(m_trans),
         "trans_only_edges": int(T.nnz // 2)}


def save(name, W):
    sp.save_npz(ART / f"{tag}_{name}.npz", W.astype(np.float32))
    r = spectral_embedding(W, k_max=50)
    np.savez_compressed(ART / f"{tag}_{name}_spectral.npz", emb=r["emb"], eigvals=r["eigvals"], all_eigvals=r["all_eigvals"],
                        embedded=r["embedded"], component=r["component"], comp_sizes=r["comp_sizes"],
                        node_chr=G.chr_idx, node_pos=G.pos, valid=valid, res=res, level="bin")
    deg = np.diff(W.tocsr().indptr)
    return dict(nodes_embedded=int(r["embedded"].sum()), isolated=r["n_isolated"], edges=int(W.nnz // 2), components=r["n_components"],
                largest_component=int(r["comp_sizes"][0]), median_degree=float(np.median(deg[deg > 0])), mean_degree=float(deg[deg > 0].mean()),
                max_degree=int(deg.max()), mean_weight=float(W.data.mean()))


if not a.only_randtrans:
    for name, W in cells.items():
        stats[name] = save(name, W); log(f"{tag}_{name}: {stats[name]}")
if a.randtrans:
    tdeg = np.diff(T.indptr)
    stats["randtrans"] = {}
    for seed in range(a.randtrans):
        Tr = rewire_trans(T, G.chr_idx, seed)
        assert (np.diff(Tr.indptr) == tdeg).all() and np.allclose(np.sort(sp.triu(Tr, 1).data), np.sort(sp.triu(T, 1).data))
        frac = 1 - ((T > 0).multiply(Tr > 0)).nnz / T.nnz
        for cis in ("hic", "lin"):
            nm = f"{cis}_randtrans_s{seed}"
            stats["randtrans"][nm] = {**save(nm, (cells[f"{cis}_cis"] + Tr).tocsr()), "trans_edges_changed_frac": float(frac)}
        log(f"seed {seed}: lin_randtrans {stats['randtrans'][f'lin_randtrans_s{seed}']}")
suffix = "_randtrans" if a.only_randtrans else ""
Path(OUT / f"bin{res // 1000}{a.tag}{suffix}_stats.json").write_text(json.dumps(stats, indent=2, default=float))
