# Pre-registration — Hi-C as a prior for unseen-perturbation prediction

**Status: DRAFT v0.2 (2026-10-14). Not frozen. No Hi-C perturbation model has been run.**
To be approved by the project owner; at freeze we regenerate `results/hic_v2/MANIFEST.json` on a clean tree and commit (and tag) before any analysis in §5–§6 is run.

### Changes since v0.1
- The primary control `b500_lin_cis` is **replaced** by a connected position-only control, `b500_lin_randtrans` (10 seeds, averaged). `b500_lin_cis` fails the new all-zero-row gate and is demoted to a secondary cell.
- Single primary test (no co-primary). Equivalence testing (TOST) with a pre-specified smallest effect of interest and a three-way decision rule added.
- K = 50 designated the **key secondary** analysis (justification from diagnostics, §5).
- Secondary additions: full factorial with matched random trans, per-component embeddings, K562 rearranged-region exclusion (S10), RPE1 replication plan (R1). Residualized embeddings stay secondary.

## 1. Question and hypothesis

Does 3D genome architecture (Hi-C) help predict transcriptional responses to *unseen* CRISPRi perturbations in K562, beyond what chromosome identity and 1D genomic position provide?

**H1 (single primary hypothesis).** Replacing the perturbation embedding P of the linear baseline by the spectral embedding of the bin-level whole-genome Hi-C graph (cis + trans), `b500_hic_trans`, **reduces test L2** relative to the connected position-only control `b500_lin_randtrans` (linear-distance chain per chromosome + random trans edges matched to the Hi-C trans graph), on the `common` eval set.

## 2. Frozen design

| Item | Specification |
|---|---|
| Benchmark | GEARS-processed Replogle K562 essential (zip md5 `84eed779531de88f83d3eb16a773a261`): 1,092 single perturbations, 5,000 read-out genes |
| Splits | GEARS `simulation` splits regenerated with gears 0.1.2, seeds **1–5** (`configs/splits/`): train 738 (incl. ctrl) / val 82 / test 273. 837 unique test perturbations over the five splits (438 appear in 1 split, 285 in 2, 99 in 3, 15 in 4) |
| Model | `Y ≈ G W Pᵀ + b`, ridge, closed form (`src/baselines.py`). **G** = top-K PCs of *training* expression pseudobulks (identical in all arms); **P** = rows of the embedding at the perturbed gene. Fit on the training split only (`refit_trainval: false`) |
| Primary hyper-parameters | **K = 10, λ = 0.1** (paper settings), no tuning. K sets both the number of PCs in G and of embedding dimensions in P |
| Embedding | first K non-trivial eigenvectors of the normalized Laplacian of the graph, columns scaled to unit RMS over embedded nodes; a gene takes the embedding of its TSS bin |
| Hi-C arm | **A = `b500_hic_trans`** (one deterministic graph) |
| Control arm | **C = `b500_lin_randtrans_s0…s9`**: same node set, same per-chromosome linear chain (weight exp(−d/5 Mb), kNN + spanning forest) as `b500_lin_cis`, plus trans edges obtained by degree-preserving randomization of the Hi-C trans graph, constrained to stay inter-chromosomal (per-node trans degrees, number of edges and weight multiset identical to Hi-C; verified for all 20 random graphs). Seeds 0–9, same sparsification and spectral procedure |
| Eval set | `common` = test perturbations predictable by every arm. Expected 273 / split (all 1,092 targets map to valid 500 kb bins; verified) |
| Metric | per-perturbation **L2** between predicted and observed delta (perturbation − control pseudobulk) on the **1,000 genes most expressed in control cells**, **perturbed gene itself excluded** (before selecting the top 1,000) |
| Graph recipe | `src/hic/bins.py` (`BIN_CFG`): 500 kb bins, genome-wide ICE, cis O/E stratified by diagonal, trans O/E by the global trans mean; kNN k=10 (degree cap 20) ∪ per-chromosome maximum spanning forest (cis), kNN k=10 (trans). Hi-C: ENCODE ENCFF616PUW (GRCh38, md5 `d4b999471d4c17a74f1e66370c0012e3`, read remotely), annotation GENCODE v44 |
| Code / configs | `src/stats.py` (analysis), `src/evaluate.py`, `configs/prereg/*.yaml` (one config per analysis, all written before any is run) |

**Eligibility gate for primary comparators.** Any embedding with **> 5 % all-zero P rows** (over the 1,092 perturbation targets, at the analysed K, or over graph nodes) cannot be a primary comparator. Computed for every embedding in `results/hic_v2/zero_rows.csv`. At K = 10: `b500_hic_trans` 0 %, all `b500_lin_randtrans` seeds 0 %, `b500_lin_trans` 0 % (**eligible**); `b500_lin_cis` 51 % (38 % of nodes), `b500_hic_cis` 56 %, per-component variants 48 %, `hic_local25k` 63 % (**ineligible**; secondary only).

## 3. Primary analysis

For each (perturbation *i*, split *s*) test case, `d_is = L2_C(i,s) − L2_A(i,s)` where `L2_C` is the **average over the 10 control seeds** of the case-level L2 (positive = Hi-C arm better). Effect size **D = (1/5) Σ_s mean_{i ∈ test_s} d_is** (equal weight per split).

- **Test (single, confirmatory).** Two-sided **sign-flip permutation test**, 100,000 permutations, **α = 0.05**. All cases of one perturbation share one random sign (a perturbation recurs across splits); equivalent to testing `D = Σ_i a_i`, `a_i = (1/5) Σ_s d_is / n_s`. Reported with it: 95 % cluster-bootstrap CI over perturbations (10,000 draws), per-split mean differences with within-split paired t/Wilcoxon p (descriptive), fraction of perturbations with `a_i > 0`, and D for each individual control seed (descriptive; averaging over seeds does not model between-seed variance).
- **Smallest effect size of interest (SESOI).** `δ = 0.5 % × (mean control L2 over the primary test cases)`. Indicatively ≈ 0.021 L2 units (0.5 % of 4.256, the non-Hi-C mean-baseline L2 on this metric; for scale, the paper's PCA embedding gains 0.106 over the mean baseline). The exact δ is computed from the data by the fixed formula.
- **Equivalence (TOST).** Two one-sided sign-flip tests against the bounds ±δ (`H0: D ≤ −δ` and `H0: D ≥ +δ`), α = 0.05 each; equivalent if both reject (`src/stats.py::tost_signflip`).
- **Decision rule (pre-specified).**
  1. **Significant** — two-sided p < 0.05: *Hi-C better* if D > 0, *Hi-C worse* if D < 0. If TOST also rejects, report "significant but within the equivalence bounds (negligible)".
  2. **Equivalent** — p ≥ 0.05 and TOST p < 0.05 (D is within ±δ with confidence).
  3. **Inconclusive** — neither.
- **H1 is supported** only if the outcome is "Significant: Hi-C better" **and** the pre-specified nulls (S2) do not reproduce an equal improvement (observed D above the 95th percentile of the null D distributions). Implementation: `src/stats.py::primary_analysis` (tested for calibration, power, TOST and the decision rule on synthetic data).

## 4. Pre-run diagnostics that fix the design (full tables: `results/hic_v2/`)

500 kb; mean over the first 10 / 50 dimensions.

| embedding | R² chromosome + position | A/B signal beyond position (best-of-3 track) | all-zero P rows (K=10) |
|---|---|---|---|
| `b500_hic_trans` (A) | 0.69 / 0.55 | **0.16 / 0.63** | 0 % |
| `b500_lin_randtrans` (C; mean±sd over 10 seeds) | 0.55±0.02 / 0.51±0.01 | **0.02 / 0.07** | 0 % |
| `b500_hic_randtrans` | 0.41±0.08 / 0.35±0.01 | 0.02 / 0.10 | 0 % |
| `b500_lin_trans` | 0.74 / 0.65 | 0.09 / 0.29 | 0 % |
| `b500_hic_cis` | 0.49 / 0.42 | 0.03 / 0.30 | 56 % |
| `b500_lin_cis` (old control) | 1.00 / 0.99 | 0.00 / 0.03 | 51 % |
| random Gaussian | 0.03 | ≈ 0.00 | 0 % |

(Position = chromosome-specific cubic B-spline, 8 df; "beyond position" = R² of the A/B track after removing chromosome + position from track and embedding. PC1 of the O/E correlation matrix does not track gene density in this rearranged K562 map — median r ≈ 0.02, PC1 best of PC1–3 on 2/23 chromosomes — so both PC1 and best-of-3 are reported; PC1 numbers are in the diagnostics.)

Known characteristics, stated in advance: both arms remain substantially explained by chromosome/position (the control by construction); trans hubs (bins at chr21p, chr22, the BCR-ABL region) are limited by a degree cap but a few bins keep degree ≈ 240 through the spanning forest (edges: Hi-C 47.8k vs control 52.3k, i.e. not exactly matched); at 250 kb the Hi-C graphs are sparser than their controls (median degree 4 vs 10 in the cis-only cells).

## 5. Key secondary analysis — K = 50

Same contrast, test and decision rule, with **fixed K = 50, λ = 0.1** (G and P both 50-dimensional). *Justification, from diagnostics only:* most Hi-C-specific compartment signal lies beyond dimension 10 — the A/B track beyond position is explained at R² 0.16 with 10 dims but 0.63 with 50 dims for `b500_hic_trans`, versus 0.02 → 0.07 for the matched control. A K = 10 test may therefore underestimate any Hi-C-specific gain. K = 10 remains primary for fidelity to the paper; a conclusion about H1 rests on the primary only. (K = 20 is also computed, descriptive.)

## 6. Other secondary and exploratory analyses (all pre-specified; configs in `configs/prereg/`)

Each uses the same paired sign-flip machinery; **Holm correction within each family at α = 0.05**; no correction across families; reported as exploratory.

- **S1 Factorial.** cis ∈ {hic, lin} × trans ∈ {hic, random-matched (10 seeds), none}: cells `b500_{hic,lin}_{trans,randtrans_s*,cis}`. Contrasts: Hi-C cis vs linear cis within each trans level; effect of Hi-C trans vs random trans within each cis level; interaction. Cis-only cells also in a **per-component** variant (`:percomp`): the K eigenvectors are drawn per connected component, allocated proportionally to component size (largest remainder, ≥ 1 per component when K ≥ 23 components), each scaled to unit RMS within its component. At K = 10 this cannot cover 23 chromosomes (48 % all-zero rows) so it is also run at K = 50 (0 %).
- **S2 Nulls** for `b500_hic_trans`: gene-label permutation (10 seeds), degree-preserving rewire (5), random Gaussian with matched dimensions (10). Observed D vs the null D distribution (empirical one-sided p).
- **S3 Position-residualized** (secondary, not co-primary): `:resid` (chromosome + per-chromosome cubic B-spline, 8 df regressed out, columns re-scaled) for `b500_hic_trans` and `b500_hic_cis` vs the control.
- **S4 Extended eval set.** Hi-C and control arms vs the **mean baseline only** on all 273 test perturbations per split; PCA embedding on its own predictable subset (≈ 100–105 per split) for reference (S9).
- **S5 Metrics.** Pearson delta; top-20 DE; all genes; target gene included.
- **S6 Cis vs trans read-out genes.** Metrics separately on read-out genes on the perturbed gene's chromosome vs others.
- **S7 K and tuning.** K ∈ {10, 20, 50} (λ = 0.1) and validation-tuned K.
- **S8 Resolution.** 250 kb analogue of the primary contrast (`b250_hic_trans` vs `b250_lin_randtrans`, 10 seeds) and the gene-level 25 kb local graph. *Caveat:* at 250 kb the Hi-C graph has a lower median degree than the controls (4 vs 10 for cis-only cells), so a difference there is partly a sparsity difference.
- **S10 K562 rearrangements.** K562 is highly rearranged (BCR-ABL t(9;22), amplifications) and translocations appear as trans contacts. Re-run the primary contrast with bins in rearranged/amplified regions removed **before balancing** (`b500x_*` graphs) and the perturbations whose target falls in those bins dropped from all arms. **Region list** (`configs/k562_rearranged_regions_hg38.bed`, 94 regions, built by `scripts/make_k562_regions.py`): (i) all 70 segments called `amp` (196.4 Mb) in ENCODE **ENCFF486MJU** (ENCSR000BBO; HAIB Illumina genotyping array copy-number segments, hg19; md5 `b4a78fa77789291cb5e3d528b32aa55a`); (ii) ±1 Mb windows around every inter-chromosomal breakpoint (24) in ENCODE **ENCFF863MPP** (ENCSR053AXS; 10X linked-read SV calls, grocsvs, hg19; md5 `04f1c1f3aad184c0e05c863da7a7d3e0`); lifted hg19→hg38 with UCSC `hg19ToHg38.over.chain.gz` (md5 `35887f73fe5e2231656504d1f6430900`; fetched from `hgdownload.soe.ucsc.edu`, because the TLS certificate of `hgdownload.cse.ucsc.edu` did not verify). Effect: 513 of 6,073 bins excluded (8.4 %), 56 of 1,092 perturbations dropped (256–266 test perturbations per split remain). The thresholds are the sources' own calls (no tuning); the window size (±1 Mb) was fixed a priori.
- **R1 Replication in RPE1 (planned; availability check only — nothing downloaded).** Same primary analysis on hTERT RPE-1. *Hi-C:* ENCODE experiment **ENCSR846DUS** (in situ Hi-C, Aiden lab, 4 replicates merged), GRCh38 `.hic` files `ENCFF876OWE` (MAPQ-thresholded, 2.24 GB, md5 `b6828cb642d1fce2b64c5383f22ec4f9`) and `ENCFF274QIY` (2.24 GB, md5 `68a2ad184e63fc5e1fe22e263a673aec`); remote reading is possible as for K562. This map is ~10× smaller than the K562 file, so a depth/ICE QC gate precedes any analysis (R1 is declared *not evaluable* if the valid-bin fraction or the zero-row gate fails). *Perturb-seq:* GEARS `replogle_rpe1_essential.zip`, 665,914,266 B, md5 `169370a8da6093470f122172ab30f7c5`. Success criterion: same direction and two-sided p < 0.05 under the identical primary test and decision rule. RPE1 is near-diploid, so no rearrangement exclusion is planned.

## 7. Pre-specified sanity checks (before interpreting anything)

1. `train_pca_expr` and the mean baseline reproduce the earlier `k562_gears` results within noise under the primary code path.
2. Random Gaussian arms do not beat the mean baseline by more than expected under the S2 null distribution.
3. Predictable sets are identical across all bin-level arms of a contrast (verified: 1,092/1,092 at 500 kb; 1,036 for `b500x_*`).
4. The eligibility gate (§2) holds for every arm of the primary contrast.

## 8. Known limitations

K562 is aneuploid and rearranged; Hi-C is bulk (Rao et al. 2014, reprocessed on GRCh38 by ENCODE) with no KR vectors (own ICE); the paper's own GEARS splits could not be recovered (ours are regenerated with the GEARS procedure); absolute L2 is ~4–6 % above the paper's (Pearson and ratios agree); chromosome identity (R² ≈ 0.4 of `b500_hic_trans`) may predict expression changes through copy-number or batch-like structure, which the control only partly removes (control: R² 0.55 on chromosome + position, 0.11 on chromosome alone); the A/B track in this map is itself largely positional; control seeds are averaged and their between-seed variance is not part of the test; a perturbed gene's response and its Hi-C neighbours share the chromosome, so "Hi-C helps" and "chromosome helps" are separated only by the control design and S1/S3.

## 9. Housekeeping

- Git history is **not** rewritten; run logs reference existing hashes. `user.email` will be set for *future* commits once the address is provided (pending).
- At freeze: regenerate `results/hic_v2/MANIFEST.json` on a clean tree (SHA-256 of every graph, embedding, split, config, region list and source file; sources' MD5s) and tag the commit.
- Deviations after freeze are logged in a `DEVIATIONS.md`, never edited into this file.
