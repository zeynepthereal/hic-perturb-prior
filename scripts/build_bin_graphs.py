"""python -m scripts.build_bin_graphs <res_bp>   e.g. 500000 | 250000
Needs data/hic/pairs_<res>kb (scripts/fetch_hic.py trans <res>) and results/hic/gene_table.csv.
Writes data/hic/graphs/b<res>k_<cell>{.npz,_spectral.npz}, results/hic_v2/ab_<res>k.csv, ab_<res>k_chrom.csv, bin<res>k_stats.json"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from src.hic.bins import BIN_CFG, Genome, ab_compartments, build_cells, ice, load_raw, oe_matrix
from src.hic.genes import HIC_CHROMS
from src.hic.spectral import spectral_embedding

res = int(sys.argv[1]); tag = f"b{res // 1000}"
ART = Path("data/hic/graphs"); ART.mkdir(parents=True, exist_ok=True)
OUT = Path("results/hic_v2"); OUT.mkdir(parents=True, exist_ok=True)
t0 = time.time(); log = lambda m: print(f"[{time.time()-t0:5.0f}s] {m}", flush=True)

G = Genome(res)
R = load_raw(G, f"data/hic/pairs_{res // 1000}kb"); log(f"raw matrix {R.shape}, total reads {R.sum(dtype=np.float64)/2:.3e}")
x, valid = ice(R, BIN_CFG["ice_iters"], BIN_CFG["ice_tol"]); log(f"ICE done; valid bins {valid.sum()}/{G.n}")
OE, m_trans, expected = oe_matrix(R, x, valid, G, BIN_CFG["min_pairs"]); del R; log(f"O/E done (global trans mean {m_trans:.3e})")

nodes = pd.read_csv("results/hic/gene_table.csv")
gb = G.gene_bins(nodes.chr.values, nodes.tss.values)
density = np.bincount(gb[gb >= 0], minlength=G.n).astype(float)
ab, ab_best, ab_tab = ab_compartments(OE, valid, G, density)
pd.DataFrame(dict(chr=[HIC_CHROMS[i] for i in G.chr_idx], start=(G.local * res), valid=valid, gene_density=density, ab_pc1=ab,
                  ab_best_of_3=ab_best)).to_csv(OUT / f"ab_{res // 1000}k.csv", index=False)
ab_tab.to_csv(OUT / f"ab_{res // 1000}k_chrom.csv", index=False)
log(f"A/B: median corr(PC1, gene density) = {ab_tab.corr_pc1_density.median():.2f}; chromosomes where PC1 is best of 3: {(ab_tab.best_pc == 1).sum()}/{len(ab_tab)}")

cells = build_cells(OE, valid, G); T = cells.pop("_trans_only"); log("cells built")
stats = {"res": res, "cfg": BIN_CFG, "n_bins": G.n, "valid_bins": int(valid.sum()), "m_trans": float(m_trans),
         "trans_only_edges": int(T.nnz // 2)}
for name, W in cells.items():
    sp.save_npz(ART / f"{tag}_{name}.npz", W.astype(np.float32))
    r = spectral_embedding(W, k_max=50)
    np.savez_compressed(ART / f"{tag}_{name}_spectral.npz", emb=r["emb"], eigvals=r["eigvals"], all_eigvals=r["all_eigvals"],
                        embedded=r["embedded"], component=r["component"], comp_sizes=r["comp_sizes"],
                        node_chr=G.chr_idx, node_pos=G.pos, valid=valid, res=res, level="bin")
    deg = np.diff(W.tocsr().indptr)
    stats[name] = dict(nodes_embedded=int(r["embedded"].sum()), isolated=r["n_isolated"], edges=int(W.nnz // 2), components=r["n_components"],
                       largest_component=int(r["comp_sizes"][0]), median_degree=float(np.median(deg[deg > 0])), mean_degree=float(deg[deg > 0].mean()),
                       max_degree=int(deg.max()), mean_weight=float(W.data.mean()))
    log(f"{tag}_{name}: {stats[name]}")
Path(OUT / f"bin{res // 1000}k_stats.json").write_text(json.dumps(stats, indent=2, default=float))
