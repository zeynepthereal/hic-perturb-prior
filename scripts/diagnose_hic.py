"""python -m scripts.diagnose_hic -> results/hic/diagnostics.md, diagnostics_*.csv, figs/*.png

Per variant: nodes/edges/components/degree distribution, eigenvalue spectrum, how much of each embedding dimension is explained by
chromosome identity and by position along the chromosome, similarity between embeddings, and how many GEARS perturbations
(overall and in each test split) become predictable compared with the PCA embedding."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

from src.baselines import EMBEDDINGS
from src.data import PerturbationDataset, split_from_conditions
from src.hic.genes import HIC_CHROMS

OUT = Path("results/hic"); FIG = OUT / "figs"; FIG.mkdir(parents=True, exist_ok=True)
ART = Path("data/hic/graphs")
VARS = ["hic_cis", "hic_cis_trans", "linear_distance"]
COL = {"hic_cis": "#2a78d6", "hic_cis_trans": "#eb6834", "linear_distance": "#1baf7a"}   # categorical slots 1-3 (validated all-pairs)
INK, GRID = "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.axisbelow": True, "legend.frameon": False, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})

nodes = pd.read_csv(OUT / "gene_table.csv")
stats = json.loads((OUT / "graph_stats.json").read_text())
S = {v: dict(np.load(ART / f"{v}_spectral.npz")) for v in VARS}
W = {v: sp.load_npz(ART / f"{v}.npz") for v in VARS}
lines = ["# Hi-C gene graph diagnostics\n"]

# ---- 1. structure
rows = []
for v in VARS:
    deg = np.diff(W[v].tocsr().indptr); d = deg[deg > 0]
    rows.append(dict(variant=v, nodes_total=len(nodes), nodes_embedded=stats[v]["nodes_embedded"], isolated=stats[v]["isolated_nodes"],
                     edges=stats[v]["edges"], components=stats[v]["n_components"], largest_component=stats[v]["largest_component"],
                     deg_min=int(d.min()), deg_median=float(np.median(d)), deg_mean=round(float(d.mean()), 2), deg_p95=float(np.percentile(d, 95)),
                     deg_max=int(d.max())))
st = pd.DataFrame(rows); st.to_csv(OUT / "diagnostics_structure.csv", index=False)
lines += ["## Graph structure\n", st.to_markdown(index=False), ""]

fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
for v in VARS:
    deg = np.diff(W[v].tocsr().indptr); d = np.sort(deg[deg > 0])
    ax[0].plot(d, 1 - np.arange(len(d)) / len(d), color=COL[v], lw=2, label=v)
ax[0].set_xscale("log"); ax[0].set_xlabel("node degree"); ax[0].set_ylabel("fraction of nodes with degree ≥ x"); ax[0].legend()
ax[0].set_title("Degree distribution", loc="left", color=INK)
for v in VARS:
    ev = S[v]["eigvals"][:50]
    ax[1].plot(np.arange(1, len(ev) + 1), ev, color=COL[v], lw=2, marker="o", ms=4, label=f"{v} ({stats[v]['n_components']} component{'s' if stats[v]['n_components'] > 1 else ''})")
ax[1].set_yscale("log"); ax[1].set_xlabel("non-trivial eigenvector index"); ax[1].set_ylabel("eigenvalue of normalized Laplacian (log)")
ax[1].legend(loc="center left", bbox_to_anchor=(0.02, 0.58), fontsize=8)
ax[1].set_title("Spectrum (trivial zero eigenvalues dropped)", loc="left", color=INK)
fig.tight_layout(); fig.savefig(FIG / "degree_and_spectrum.png", dpi=150); plt.close(fig)

# ---- 2. variance explained by chromosome / position
chr_i = nodes.chr.map({c: i for i, c in enumerate(HIC_CHROMS)}).values
pos = nodes.tss.values.astype(float)
posn = np.zeros(len(nodes))
for c in HIC_CHROMS:
    m = chr_i == HIC_CHROMS.index(c)
    posn[m] = (pos[m] - pos[m].min()) / max(pos[m].max() - pos[m].min(), 1)
onehot = np.eye(len(HIC_CHROMS))[chr_i]
X_chr = onehot
X_lin = np.hstack([onehot, onehot * (posn[:, None] - 0.5)])
X_poly = np.hstack([onehot] + [onehot * ((posn[:, None] - 0.5) ** p) for p in (1, 2, 3)])


def r2(Y, X):
    B, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return 1 - ((Y - X @ B) ** 2).sum(0) / ((Y - Y.mean(0)) ** 2).sum(0)


rows, per_dim = [], {}
for v in VARS:
    E, ok = S[v]["emb"][:, :50], S[v]["embedded"].astype(bool)
    res = {n: r2(E[ok], X[ok]) for n, X in (("chromosome", X_chr), ("chrom+linear_pos", X_lin), ("chrom+poly3_pos", X_poly))}
    per_dim[v] = res
    for kmax in (10, 50):
        rows.append(dict(variant=v, dims=f"first {kmax}", **{f"R2_{n}": round(float(a[:kmax].mean()), 3) for n, a in res.items()}))
# reference: a random Gaussian embedding of the same dims (expected R2 ~ n_params/n_nodes)
Er = np.random.RandomState(0).randn(len(nodes), 10)
rows.append(dict(variant="random_gaussian (reference)", dims="first 10",
                 **{f"R2_{n}": round(float(r2(Er, X).mean()), 3) for n, X in (("chromosome", X_chr), ("chrom+linear_pos", X_lin), ("chrom+poly3_pos", X_poly))}))
r2df = pd.DataFrame(rows); r2df.to_csv(OUT / "diagnostics_r2_chrom_pos.csv", index=False)
pd.concat({v: pd.DataFrame(per_dim[v]).iloc[:50] for v in VARS}).to_csv(OUT / "diagnostics_r2_per_dim.csv")
lines += ["## Variance of embedding dimensions explained by chromosome / position (R², mean over dims)\n",
          r2df.to_markdown(index=False),
          "\n`chromosome` = 23 chromosome indicators; `chrom+linear_pos` adds a per-chromosome linear position term; `chrom+poly3_pos` a per-chromosome cubic. "
          "Low-frequency eigenvectors of a path-like graph are cosines of position, so the cubic term is the fair 'linear genomic position' benchmark.\n"]
flags = [f"**FLAG {v}**: first-10 mean R² of chromosome+cubic position = {r2df[(r2df.variant==v)&(r2df.dims=='first 10')]['R2_chrom+poly3_pos'].iloc[0]:.2f} (> 0.5): the embedding is largely chromosome/position."
         for v in VARS if r2df[(r2df.variant == v) & (r2df.dims == "first 10")]["R2_chrom+poly3_pos"].iloc[0] > 0.5]
lines += flags + [""] if flags else ["No variant exceeds R² 0.5 for chromosome+cubic position on the first 10 dims.\n"]

fig, ax = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
for a, v in zip(ax, VARS):
    d = np.arange(1, 11)
    a.plot(d, per_dim[v]["chromosome"][:10], color="#2a78d6", lw=2, marker="o", ms=4, label="chromosome")
    a.plot(d, per_dim[v]["chrom+poly3_pos"][:10], color="#eb6834", lw=2, marker="o", ms=4, label="chromosome + position (cubic)")
    a.set_title(v, loc="left", color=INK); a.set_xlabel("embedding dimension"); a.set_ylim(0, 1.02)
ax[0].set_ylabel("R²")
h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=2, fontsize=9)
fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(FIG / "r2_chromosome_position.png", dpi=150); plt.close(fig)

# ---- 3. similarity between embeddings (canonical correlations, first 10 dims) and graphs
def cca(A, B):
    A, B = A - A.mean(0), B - B.mean(0)
    qa, _ = np.linalg.qr(A); qb, _ = np.linalg.qr(B)
    return np.linalg.svd(qa.T @ qb, compute_uv=False)


rows = []
for a, b in (("hic_cis", "linear_distance"), ("hic_cis_trans", "linear_distance"), ("hic_cis_trans", "hic_cis")):
    ok = S[a]["embedded"].astype(bool) & S[b]["embedded"].astype(bool)
    cc = cca(S[a]["emb"][ok, :10], S[b]["emb"][ok, :10])
    rows.append(dict(pair=f"{a} vs {b}", n_common_nodes=int(ok.sum()), mean_canonical_corr_10d=round(float(cc.mean()), 3),
                     top3_canonical=np.round(cc[:3], 3).tolist(), min_canonical=round(float(cc.min()), 3)))
cc_df = pd.DataFrame(rows); cc_df.to_csv(OUT / "diagnostics_embedding_similarity.csv", index=False)
g = stats["hic_cis_vs_linear_distance"]
lines += ["## Similarity of Hi-C and linear-distance constructions\n", cc_df.to_markdown(index=False), "",
          "Graph level (after O/E):\n", "```", json.dumps(g, indent=1), "```\n"]

# ---- 4. predictability on the GEARS benchmark
ds = PerturbationDataset.load("data/processed/k562_gears.npz")
names = ["train_pca_expr"] + VARS
ok = {n: EMBEDDINGS[n].predictable(ds) for n in names}
rows = [dict(scope="all 1092 perturbations", **{n: int(ok[n].sum()) for n in names},
             common_pca_and_hic_cis_trans=int((ok["train_pca_expr"] & ok["hic_cis_trans"]).sum()))]
for s in range(1, 6):
    sp_ = split_from_conditions(ds, json.load(open(f"configs/splits/gears_k562_seed{s}.json")))
    te = sp_["test"]
    rows.append(dict(scope=f"test split gears_seed{s} (n={len(te)})", **{n: int(ok[n][te].sum()) for n in names},
                     common_pca_and_hic_cis_trans=int((ok["train_pca_expr"] & ok["hic_cis_trans"])[te].sum())))
pr = pd.DataFrame(rows); pr.to_csv(OUT / "diagnostics_predictable.csv", index=False)
lines += ["## Predictable GEARS perturbations (eval-set sizes)\n", pr.to_markdown(index=False),
          "\n`common` eval set for PCA-vs-Hi-C comparisons = PCA ∩ Hi-C column; `extended:<hic>` = the Hi-C column (vs mean baseline only).\n"]
(OUT / "diagnostics.md").write_text("\n".join(lines))
print("\n".join(lines))
