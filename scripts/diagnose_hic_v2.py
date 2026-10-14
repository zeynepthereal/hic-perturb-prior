"""python -m scripts.diagnose_hic_v2 -> results/hic_v2/diagnostics.md + csv + figs/
All variants: structure (components, degrees), spectrum, R2 of embedding dims on chromosome / chromosome+position (first 10 and 50 dims),
correlation with the A/B compartment track, localisation of the embedding, and eval-set sizes on the GEARS benchmark."""
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
from src.hic.bins import Genome
from src.hic.embeddings import ART, gene_tables, load_spectral
from src.hic.position import position_basis, r2_columns, residualize

OUT = Path("results/hic_v2"); FIG = OUT / "figs"; FIG.mkdir(parents=True, exist_ok=True)
RES = [r for r in (500, 250) if (ART / f"b{r}_hic_trans_spectral.npz").exists()]
CELLS = ["hic_trans", "hic_cis", "lin_trans", "lin_cis"]
VARS = [f"b{r}_{c}" for r in RES for c in CELLS] + ["hic_local25k"]
BLUE, ORANGE, INK, GRID = "#2a78d6", "#eb6834", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.axisbelow": True, "legend.frameon": False, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})
style = lambda v: dict(color=ORANGE if "_lin_" in v else BLUE, ls="-" if v.endswith("trans") else "--") if v != "hic_local25k" else dict(color="#52514e", ls=":")

nodes, sym, node_of_id = gene_tables()
chr_of_gene = nodes.chr.map({c: i for i, c in enumerate(sorted(set(nodes.chr), key=lambda c: (len(c), c)))})
lines = ["# Hi-C graph diagnostics v2 (bin-level, whole genome)\n"]


def node_info(v):
    """returns embedding (n,50), ok, chr_idx, pos, ab(dict name->array or None) for the nodes of variant v"""
    base = v.split(":")[0]
    z = load_spectral(base)
    E, ok = EMBEDDINGS[v]._node_matrix()
    chr_idx, pos = EMBEDDINGS[v]._node_coords(z)
    ab = {}
    if "res" in z:
        t = pd.read_csv(OUT / f"ab_{int(z['res']) // 1000}k.csv")
        ab = {"AB_PC1": t.ab_pc1.values, "AB_best_of_3": t.ab_best_of_3.values}
    else:                                              # gene-level local graph: use the 500 kb track at each gene's TSS bin
        G = Genome(500_000); t = pd.read_csv(OUT / "ab_500k.csv")
        gb = G.gene_bins(nodes.chr.values, nodes.tss.values)
        ab = {"AB_PC1": np.where(gb >= 0, t.ab_pc1.values[np.maximum(gb, 0)], np.nan),
              "AB_best_of_3": np.where(gb >= 0, t.ab_best_of_3.values[np.maximum(gb, 0)], np.nan)}
    return E, ok, chr_idx, pos, ab


# ---------------- structure
rows = []
for v in VARS:
    z = load_spectral(v); W = sp.load_npz(ART / f"{v}.npz"); d = np.diff(W.tocsr().indptr); d = d[d > 0]
    ncomp = int(np.sum(z["comp_sizes"] > 0))
    rows.append(dict(variant=v, nodes_embedded=int(z["embedded"].sum()), isolated=int((~z["embedded"]).sum()), edges=int(W.nnz // 2),
                     components=ncomp, largest_component=int(z["comp_sizes"][0]), deg_median=float(np.median(d)), deg_mean=round(float(d.mean()), 1),
                     deg_p99=float(np.percentile(d, 99)), deg_max=int(d.max()), first_eigval=float(z["eigvals"][0]), eigval_10=float(z["eigvals"][9]),
                     eigval_50=float(z["eigvals"][49])))
st = pd.DataFrame(rows); st.to_csv(OUT / "structure.csv", index=False)
lines += ["## Structure (components, degrees, spectrum)\n", st.round(5).to_markdown(index=False), ""]

fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
for r_i, r in enumerate(RES[:1]):
    for c in CELLS:
        v = f"b{r}_{c}"; W = sp.load_npz(ART / f"{v}.npz"); d = np.sort(np.diff(W.indptr)[np.diff(W.indptr) > 0])
        ax[0].plot(d, 1 - np.arange(len(d)) / len(d), lw=2, label=v, **style(v))
        ax[1].plot(np.arange(1, 51), load_spectral(v)["eigvals"], lw=2, marker="o", ms=3, label=f"{v} ({int(st.set_index('variant').loc[v,'components'])} comp.)", **style(v))
ax[0].set_xscale("log"); ax[0].set_xlabel("node degree"); ax[0].set_ylabel("fraction of nodes with degree ≥ x"); ax[0].set_title(f"Degree distribution (b{RES[0]})", loc="left", color=INK)
ax[0].legend(fontsize=8)
ax[1].set_yscale("log"); ax[1].set_xlabel("non-trivial eigenvector index"); ax[1].set_ylabel("Laplacian eigenvalue (log)")
ax[1].set_title(f"Spectrum, trivial modes dropped (b{RES[0]})", loc="left", color=INK); ax[1].legend(fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig(FIG / "degree_and_spectrum.png", dpi=150); plt.close(fig)

# ---------------- R2 on chromosome / position, localisation, A/B
r2rows, abrows = [], []
ALL = VARS + [f"b{r}_{c}:resid" for r in RES for c in ("hic_trans", "hic_cis")] + [f"b{RES[0]}_hic_trans:random:0"]
for v in ALL:
    E, ok, ci, pos, ab = node_info(v)
    Xc, Xs, Xp = (position_basis(ci, pos, kind=k) for k in ("chromosome", "spline", "cubic"))
    for kmax in (10, 50):
        Y = E[ok][:, :kmax]
        support = float(np.mean((np.abs(E[ok][:, :kmax]) > 1e-8).any(1)))
        # localisation: for each dim the fraction of embedded nodes with non-negligible loading (|x| > 0.05 RMS)
        loc = float(np.mean((np.abs(Y) > 0.05).mean(0)))
        r2rows.append(dict(variant=v, dims=kmax, R2_chromosome=float(r2_columns(Y, Xc[ok]).mean()), R2_chrom_cubic=float(r2_columns(Y, Xp[ok]).mean()),
                           R2_chrom_spline=float(r2_columns(Y, Xs[ok]).mean()), frac_nodes_nonzero_in_a_dim=support, mean_frac_nodes_loaded_per_dim=loc))
    for name, a in ab.items():
        m = ok & np.isfinite(a)
        pos_only = float(r2_columns(a[m][:, None], Xs[m])[0])
        for kmax in (10, 50):
            Y = E[m][:, :kmax]
            rs = [abs(np.corrcoef(Y[:, k], a[m])[0, 1]) if Y[:, k].std() > 0 else 0 for k in range(kmax)]
            R2 = float(r2_columns(a[m][:, None], np.hstack([Y, np.ones((m.sum(), 1))]))[0])
            # beyond position: residualize BOTH on chromosome + spline position, then regress
            ar = a[m] - Xs[m] @ np.linalg.lstsq(Xs[m], a[m], rcond=None)[0]
            Yr = Y - Xs[m] @ np.linalg.lstsq(Xs[m], Y, rcond=None)[0]
            R2b = float(r2_columns(ar[:, None], np.hstack([Yr, np.ones((m.sum(), 1))]))[0])
            abrows.append(dict(variant=v, ab=name, dims=kmax, max_abs_r_single_dim=max(rs), best_dim=int(np.argmax(rs)) + 1, multiple_R=np.sqrt(max(R2, 0)),
                               R2_beyond_position=R2b, R2_AB_from_position_only=pos_only, n_nodes=int(m.sum())))
r2df, abdf = pd.DataFrame(r2rows), pd.DataFrame(abrows)
r2df.to_csv(OUT / "r2_chrom_position.csv", index=False); abdf.to_csv(OUT / "ab_correlation.csv", index=False)
lines += ["## R² of embedding dimensions explained by chromosome / position (mean over dims)\n",
          "`R2_chrom_spline` = chromosome-specific cubic B-spline of position (8 df each); `frac_nodes_nonzero_in_a_dim` and `mean_frac_nodes_loaded_per_dim` measure "
          "how localized the eigenvectors are (cis-only graphs have one localized slow mode per chromosome).\n",
          r2df.round(3).to_markdown(index=False), ""]
lines += ["## Correlation with A/B compartment track (bin level)\n",
          "AB_PC1 = PC1 of the per-chromosome O/E correlation matrix oriented by gene density (the standard definition, as specified); AB_best_of_3 = whichever of PC1-3 "
          "correlates best with gene density (PC1 does not track gene density in this rearranged K562 map; see ab_*k_chrom.csv). `multiple_R` = sqrt(R²) of AB regressed on the first "
          "`dims` embedding dimensions; `R2_beyond_position` repeats the regression after removing chromosome+position from both sides; `R2_AB_from_position_only` is the reference.\n",
          abdf.round(3).to_markdown(index=False), ""]

# ---------------- figures
fig, ax = plt.subplots(1, 2, figsize=(11, 0.45 * len(ALL) + 1.2), sharey=True)
ys = np.arange(len(ALL))[::-1]
for a, col, title in ((ax[0], "R2_chromosome", "R² on chromosome"), (ax[1], "R2_chrom_spline", "R² on chromosome + position (spline)")):
    for k, (kmax, filled) in enumerate(((10, True), (50, False))):
        s = r2df[r2df.dims == kmax].set_index("variant").loc[ALL]
        a.scatter(s[col].values, ys, s=38, facecolors=[("#2a78d6" if "_hic_" in v and "random" not in v else "#eb6834" if "_lin_" in v else "#52514e") if filled else "none" for v in ALL],
                  edgecolors=[("#2a78d6" if "_hic_" in v and "random" not in v else "#eb6834" if "_lin_" in v else "#52514e") for v in ALL], linewidths=1.5, zorder=3,
                  label=f"first {kmax} dims")
    a.set_title(title, loc="left", color=INK); a.set_xlim(-0.02, 1.02); a.set_xlabel("mean R² over dimensions")
ax[0].set_yticks(ys); ax[0].set_yticklabels(ALL, fontsize=8)
h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=2, fontsize=9)
fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(FIG / "r2_chromosome_position.png", dpi=150); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(11, 0.45 * len(ALL) + 1.2), sharey=True)
for a, name in zip(ax, ("AB_PC1", "AB_best_of_3")):
    for kmax, filled in ((10, True), (50, False)):
        s = abdf[(abdf.ab == name) & (abdf.dims == kmax)].set_index("variant").loc[ALL]
        cols = [("#2a78d6" if "_hic_" in v and "random" not in v else "#eb6834" if "_lin_" in v else "#52514e") for v in ALL]
        a.scatter(s.multiple_R.values, ys, s=38, facecolors=cols if filled else "none", edgecolors=cols, linewidths=1.5, zorder=3, label=f"first {kmax} dims")
    a.set_title(f"Multiple R with {name}", loc="left", color=INK); a.set_xlim(-0.02, 1.02); a.set_xlabel("multiple correlation")
ax[0].set_yticks(ys); ax[0].set_yticklabels(ALL, fontsize=8)
h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=2, fontsize=9)
fig.tight_layout(rect=(0, 0.04, 1, 1)); fig.savefig(FIG / "ab_correlation.png", dpi=150); plt.close(fig)

# ---------------- predictable perturbations
ds = PerturbationDataset.load("data/processed/k562_gears.npz")
names = ["train_pca_expr"] + [v for v in VARS]
ok = {n: EMBEDDINGS[n].predictable(ds) for n in names}
prow = [dict(scope="all 1092 perturbations", **{n: int(ok[n].sum()) for n in names})]
for s in range(1, 6):
    te = split_from_conditions(ds, json.load(open(f"configs/splits/gears_k562_seed{s}.json")))["test"]
    prow.append(dict(scope=f"test seed{s} (n={len(te)})", **{n: int(ok[n][te].sum()) for n in names}))
pr = pd.DataFrame(prow); pr.to_csv(OUT / "predictable.csv", index=False)
common_all = np.logical_and.reduce([ok[n] for n in names])
lines += ["## Predictable GEARS perturbations\n", pr.T.to_markdown(), f"\n`common` set of PCA ∩ all Hi-C variants: {int(common_all.sum())} perturbations; "
          f"`common` among bin-level cells only: {int(np.logical_and.reduce([ok[n] for n in VARS if n.startswith('b')]).sum())}.\n"]
(OUT / "diagnostics.md").write_text("\n".join(lines))
print("\n".join(lines))
