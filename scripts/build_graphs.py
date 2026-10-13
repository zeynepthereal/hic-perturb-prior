"""python -m scripts.build_graphs  -> data/hic/graphs/{hic_cis,hic_cis_trans,linear_distance}{.npz,_spectral.npz}, results/hic/graph_stats.json
Requires scripts/build_gene_table.py and scripts/fetch_hic.py (cis, trans) to have run."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr

from src.hic.genes import CHROM_LEN, HIC_CHROMS
from src.hic.graphs import GRAPH_CFG, build_variants, gene_cis_candidates
from src.hic.spectral import spectral_embedding

ART = Path("data/hic/graphs"); ART.mkdir(parents=True, exist_ok=True)
nodes = pd.read_csv("results/hic/gene_table.csv")
chrom_bins = {c: -(-CHROM_LEN[c] // GRAPH_CFG["trans_res"]) for c in HIC_CHROMS}
V = build_variants(nodes, "data/hic/cis_25kb", "data/hic/pairs_500kb", chrom_bins)
extras = V.pop("_extras")
stats = {"graph_cfg": GRAPH_CFG, "n_nodes_total": len(nodes), **extras}
for name, W in V.items():
    sp.save_npz(ART / f"{name}.npz", W.astype(np.float32))
    r = spectral_embedding(W, k_max=50)
    np.savez_compressed(ART / f"{name}_spectral.npz", emb=r["emb"], eigvals=r["eigvals"], all_eigvals=r["all_eigvals"],
                        trivial_eigvals=r["trivial_eigvals"], embedded=r["embedded"], component=r["component"],
                        comp_sizes=r["comp_sizes"])
    deg = np.diff(W.tocsr().indptr)
    stats[name] = {"nodes_embedded": int(r["embedded"].sum()), "isolated_nodes": r["n_isolated"], "edges": int(W.nnz // 2),
                   "n_components": r["n_components"], "largest_component": int(r["comp_sizes"][0]),
                   "median_degree": float(np.median(deg[deg > 0])), "mean_degree": float(deg[deg > 0].mean()),
                   "max_degree": int(deg.max()), "median_weight": float(np.median(W.data)),
                   "n_zero_eigenvalues_dropped": r["n_components"]}
    print(name, stats[name], flush=True)

# hic_cis vs linear_distance after O/E: (i) overlap/correlation of the sparsified graphs, (ii) O/E vs 1D distance on ALL candidate pairs
A, B = sp.triu(V["hic_cis"], 1).tocoo(), sp.triu(V["linear_distance"], 1).tocoo()
ea, eb = set(zip(A.row.tolist(), A.col.tolist())), set(zip(B.row.tolist(), B.col.tolist()))
shared = ea & eb
tss = nodes.tss.values.astype(float)
cand = gene_cis_candidates(nodes, "data/hic/cis_25kb").tocoo()
m = cand.row < cand.col
dist = np.abs(tss[cand.row[m]] - tss[cand.col[m]])
rng = np.random.RandomState(0)
sel = rng.choice(m.sum(), min(2_000_000, m.sum()), replace=False)
stats["hic_cis_vs_linear_distance"] = {
    "edges_hic_cis": len(ea), "edges_linear": len(eb), "edges_shared": len(shared), "jaccard": len(shared) / len(ea | eb),
    "frac_hic_cis_edges_also_in_linear": len(shared) / len(ea),
    "spearman_OE_vs_neg_distance_all_candidate_pairs": float(spearmanr(cand.data[m][sel], -dist[sel])[0]),
    "spearman_OE_vs_neg_distance_kept_edges": float(spearmanr(A.data, -np.abs(tss[A.row] - tss[A.col]))[0]),
    "n_candidate_pairs": int(m.sum())}
if shared:
    ia = {e: v for e, v in zip(zip(A.row.tolist(), A.col.tolist()), A.data)}
    ib = {e: v for e, v in zip(zip(B.row.tolist(), B.col.tolist()), B.data)}
    s = sorted(shared)
    stats["hic_cis_vs_linear_distance"]["spearman_weights_on_shared_edges"] = float(spearmanr([ia[e] for e in s], [ib[e] for e in s])[0])
Path("results/hic/graph_stats.json").write_text(json.dumps(stats, indent=2, default=float))
print(json.dumps(stats["hic_cis_vs_linear_distance"], indent=2))
