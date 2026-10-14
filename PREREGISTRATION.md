# Pre-registration — Hi-C as a prior for unseen-perturbation prediction

**Status: DRAFT v0.1 (2026-10-14). Not frozen.** To be finalized with the project owner and committed **before any Hi-C model is run**.
Items marked **[OPEN]** need a decision. At freeze we will fill in: git commit, `results/hic_v2/MANIFEST.json` (SHA-256 of every graph/embedding/split artifact), and the date.

## 1. Question and hypotheses

Does 3D genome architecture (Hi-C) help predict transcriptional responses to *unseen* CRISPRi perturbations in K562, beyond what 1D genomic position and chromosome identity provide?

**H1 (primary, confirmatory).** Replacing the perturbation embedding P of the linear baseline by the spectral embedding of the bin-level whole-genome Hi-C graph (cis + trans), `b500_hic_trans`, **reduces test L2** relative to the control `b500_lin_cis` (linear-genomic-distance chain per chromosome, no trans), on the `common` eval set.
H0: the mean paired difference is zero. Test is two-sided; support for H1 additionally requires the difference to be in the stated direction.

## 2. Frozen design

| Item | Specification |
|---|---|
| Benchmark | GEARS-processed Replogle K562 essential (zip md5 `84eed779531de88f83d3eb16a773a261`), 1,092 single perturbations, 5,000 read-out genes |
| Splits | GEARS `simulation` splits regenerated with gears 0.1.2, seeds **1–5** (`configs/splits/`); train 738 (incl. ctrl) / val 82 / test 273 per split. 837 unique test perturbations across the five splits (438 appear once, 285 twice, 99 three times, 15 four times) |
| Model | `Y ≈ G W Pᵀ + b`, ridge, closed form (`src/baselines.py`). **G** = top-10 PCs of *training* expression pseudobulks (identical in all arms); **P** = rows of the embedding at the perturbed gene |
| Fixed hyper-parameters | **K = 10, λ = 0.1** (paper settings). No tuning in the primary analysis. Model fit on the training split only (`refit_trainval: false`) |
| Embeddings | first 10 non-trivial eigenvectors of the normalized Laplacian, columns scaled to unit RMS over embedded nodes; a gene takes the embedding of its TSS bin (genes in one bin share it) |
| Arms | **A = `b500_hic_trans`** (Hi-C cis + trans), **C = `b500_lin_cis`** (control). Same node set, identical predictable perturbations |
| Eval set | `common` = test perturbations predictable by both arms. Expected 273 per split (all 1,092 targets map to valid bins) |
| Metric | per-perturbation **L2** between predicted and observed delta (perturbation − control pseudobulk) on the **1,000 genes most expressed in control cells**, **perturbed gene itself excluded** (before selecting the top 1,000) |
| Graphs | `src/hic/bins.py` (`BIN_CFG`): 500 kb bins, genome-wide ICE, cis O/E stratified by diagonal, trans O/E by global mean; kNN (k=10) ∪ per-chromosome maximum spanning forest (cis), kNN (k=10) trans; degree cap 20; Hi-C ENCODE ENCFF616PUW (GRCh38), GENCODE v44 |
| Code | analysis code `src/stats.py` and its tests are committed before the run; the run uses `configs/k562_gears_hic_primary.yaml` |

## 3. Primary analysis

For each (perturbation *i*, split *s*) test case, `d_is = L2_C(i,s) − L2_A(i,s)` (positive = Hi-C arm better).
Effect size: **D = (1/5) Σ_s mean_{i ∈ test_s} d_is**, i.e. equal weight per split.

**Test.** Two-sided **sign-flip permutation test**, 100,000 permutations, **α = 0.05**. All cases of one perturbation share one random sign (a perturbation recurs across splits and its cases are not independent), which is equivalent to testing `D = Σ_i a_i` with `a_i = (1/5) Σ_s d_is / n_s`.
Reported with it: 95 % percentile cluster-bootstrap CI over perturbations (10,000 draws); per-split mean differences with within-split paired t and Wilcoxon p-values (descriptive only); fraction of perturbations with `a_i > 0`.
A single confirmatory test, so no multiplicity correction. Implementation: `src/stats.py::paired_table`, `signflip_test` (tested on synthetic data for calibration and power).

**Interpretation rule.** H1 is *supported* only if `D > 0`, `p < 0.05`, the CI excludes 0, **and** the null controls of §5 (S2) do not reproduce an equal improvement. Otherwise H1 is not supported (including the case `D < 0`).

## 4. What the pre-run diagnostics already tell us (so they are not "discovered" later)

Full tables in `results/hic_v2/diagnostics.md` (graphs and embeddings were fixed before any perturbation model was run).

- `b500_lin_cis` is 99.7 % explained by chromosome + position (R², first 10 dims) — by design a chain per chromosome. `b500_hic_trans` is 69 % (first 10) / 55 % (first 50) explained; `b500_lin_trans` 74 % / 65 %; `b500_hic_cis` 49 % / 42 %; random embedding 3 %.
- **Degenerate control.** With one component per chromosome, the 10 smallest eigenvalues sit on only a few chromosomes: only 62 % of nodes are non-zero in any of the first 10 `b500_lin_cis` dimensions (each dimension loads ~6 % of nodes). Genes on the other chromosomes receive an all-zero P row. A win of `b500_hic_trans` over this control can therefore arise simply because it informs every gene. **[OPEN 1]**
- Compartment signal *beyond position* (R² of the best-of-3 A/B track after removing chromosome + spline position from both sides): `b500_hic_trans` 0.16 with 10 dims / 0.63 with 50; `b500_lin_trans` 0.09 / 0.29; `b500_hic_cis` 0.03 / 0.30; `b500_lin_cis` 0.004 / 0.03. Most Hi-C-specific compartment information sits beyond dimension 10.
- PC1 of the O/E correlation matrix does not track gene density in this rearranged K562 map (median r ≈ 0.02; PC1 is the best of PC1–3 for 2/23 chromosomes), so "A/B" is reported both as PC1 (the specified definition) and best-of-3.
- Trans hubs (bins at chr21p, chr22, the BCR-ABL region) are limited by the degree cap but a few bins retain degree ≈ 240 via the spanning forest. At 250 kb the Hi-C graphs are sparser than their controls (median degree 4 vs 10 in the cis-only cells).

## 5. Secondary and exploratory analyses (all pre-specified; none confirmatory)

Each is analysed with the same paired sign-flip machinery. **Holm correction within each family at α = 0.05**; no correction across families; results labelled exploratory.

- **S1 Factorial (cis source × trans).** Cells `b500_{hic,lin}_{trans,cis}`. Contrasts: hic_trans vs lin_trans (Hi-C cis given trans); hic_cis vs lin_cis (Hi-C cis without trans); hic_trans vs hic_cis and lin_trans vs lin_cis (effect of trans); interaction contrast. S1a, `hic_trans` vs `lin_trans`, is the cleanest isolation of Hi-C-specific cis information. **[OPEN 1]**
- **S2 Nulls.** For `b500_hic_trans`: `:perm:<seed>` (gene-label permutation, 10 seeds), `:rewire:<seed>` (degree-preserving rewire, 5 seeds), `:random:<seed>` (Gaussian, matched dimensions, 10 seeds). Observed D is compared with the null distribution of D (empirical one-sided p from the null draws).
- **S3 Position-residualized** embeddings (`:resid`; chromosome + per-chromosome cubic B-spline of position, 8 df, columns re-scaled): `b500_hic_trans:resid` and `b500_hic_cis:resid` vs the raw versions and vs `b500_lin_cis`.
- **S4 Extended eval set.** Hi-C arms vs the **mean baseline only** on all 273 test perturbations per split, and the PCA embedding on its own predictable subset (≈ 100–105 per split).
- **S5 Metrics.** Pearson delta; top-20 DE genes; all genes; target gene included; for each primary contrast.
- **S6 Cis vs trans read-out genes.** Metrics computed separately on read-out genes on the same chromosome as the perturbed gene vs all others.
- **S7 K.** K ∈ {10, 20, 50} (λ = 0.1), plus validation-tuned K/λ.
- **S8 Resolution/locality.** 250 kb cells (`b250_*`) and the gene-level 25 kb local graph (`hic_local25k`).
- **S9 Reference.** All arms vs the training-PCA embedding on the `common` set with PCA (≈ 100–105 per split) and vs the mean baseline.

## 6. Pre-specified sanity checks (must hold before interpreting anything)

1. `train_pca_expr` and the mean baseline reproduce `results/*k562_gears*` within Monte-Carlo noise when rerun with the primary code path.
2. Random-embedding arms do not beat the mean baseline by more than expected by chance (null distribution of S2).
3. Predictable sets are identical across the four bin-level cells (verified: 1,092/1,092 for every cell).

## 7. Known limitations (stated in advance)

K562 is a rearranged, aneuploid line (BCR-ABL amplicon; translocations create real trans contacts and hubs); Hi-C is bulk, from Rao et al. 2014, reprocessed on GRCh38 by ENCODE; the file has no KR vectors, so balancing is our own ICE; the paper's own GEARS splits could not be recovered (ours are regenerated); absolute L2 on this benchmark is ~4–6 % above the paper's (Pearson and ratios agree); chromosome identity itself (R² ≈ 0.4 of `b500_hic_trans`) may predict expression changes through copy-number or batch-like structure; a gene's perturbation response and its Hi-C neighbours share the chromosome, so "Hi-C helps" and "chromosome helps" are not separated without S1/S3.

## 8. Open decisions for the project owner

1. **[OPEN 1] Control and contrast.** Keep the specified primary (`hic_trans` vs `lin_cis`) even though `lin_cis` is degenerate at K = 10? Options: (a) keep as specified, and report S1a as the headline robustness check; (b) make `hic_trans` vs `lin_trans` a **co-primary** (α = 0.025 each, Holm); (c) change the disconnected-graph rule to a per-component embedding (first K modes of *each* chromosome, sign-fixed) — this changes the spec "first K non-trivial eigenvectors" and must be decided before freezing. *Recommendation: (b), keeping the spec'd eigenvector rule.*
2. **[OPEN 2] Smallest effect of interest.** Propose reporting an equivalence/SESOI bound (e.g. |D| < 0.5 % of the control's mean L2) so a null result can be called a null. Value to be chosen.
3. **[OPEN 3] Resolution.** 500 kb primary, 250 kb secondary (as drafted).
4. **[OPEN 4] Null seed counts** (10/5/10) and α = 0.05 two-sided.
5. **[OPEN 5]** Whether S3 (residualized) should be promoted to co-primary: it is the only analysis that removes the positional confound by construction.
