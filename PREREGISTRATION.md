# Pre-registration — Hi-C as a prior for unseen-perturbation prediction

**Status: FROZEN v1.0 (2026-10-14), tag prereg-v1. No Hi-C perturbation model had been run at freeze.**
The manifest (`results/hic_v2/MANIFEST.json`) was regenerated on a clean tree at freeze, and the commit was tagged `prereg-v1`.

### Changes since v0.4
- **Tier 2 made apples-to-apples (§3, S3):** the Tier-2 test now compares `b500_hic_trans:resid` against **`b500_lin_randtrans:resid`** (all ten control seeds residualized with the identical basis: chromosome + per-chromosome cubic B-spline, 8 df, same column re-scaling), control value = mean of the ten per-seed case-level L2 values exactly as in the primary; same test, α = 0.05 two-sided, K = 10, λ = 0.1, serial gatekeeping after Tier 1. The v0.4 comparison against the *un*residualized control would have been biased against Hi-C whenever position is predictive (that control still carries chromosome/position information, R² 0.55). The residualized control passes the eligibility gate (0 % all-zero rows at K = 10 and 50). `hic_cis:resid` stays exploratory. Nothing else changed.

### Changes v0.3 → v0.4
- **Null distributions enlarged and the H1 support condition made precise (§3, S2):** gene-label permutation 1,000 seeds, random Gaussian 1,000 seeds, degree-preserving rewire 20 seeds; null `D_j` defined with the same pairing and split weighting as `D`; H1 requires `D` above the 95th percentile of the permutation null, above the 95th percentile of the random-Gaussian null, and above the maximum of all 20 rewire seeds; empirical p reported per null.
- **Two-tier interpretation rule added (§3):** Tier 1 (H1 supported) and Tier 2 (mechanistic claim, additionally requires S3 `b500_hic_trans:resid` to be significant in the same direction). Rationale in §8.
- **§7 sanity check 2 replaced** (random-Gaussian arms vs the mean baseline within δ; reported, not gating); check 5 added (lean null evaluator reproduces the standard path).
- R1 formatting fixed; §9 freeze procedure extended (annotated tag `prereg-v1`, push of commit and tag, independent OSF timestamp recorded in `DEVIATIONS.md`).

### Changes v0.2 → v0.3
- **Control wording made explicit (§3):** the control value of a case is the *mean of the ten per-seed case-level L2 values*, **not** the L2 of an averaged prediction (which would give the control an ensemble advantage).
- **Descriptive companion added (§3):** rank of the Hi-C arm's mean test L2 among the ten control seeds, per split and pooled (no p-value).
- **Owner decisions recorded:** S10 region definition approved as written; RPE1 replication kept in the plan, **nothing downloaded**.
- Housekeeping: git history was rewritten once (see §9) before any freeze.

### Changes v0.1 → v0.2
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
| Code / configs | `src/stats.py` (analysis), `src/evaluate.py`, `src/nulls.py` (lean evaluator and driver for the S2 nulls; unit-tested for exact agreement with `evaluate.py` on synthetic data), `configs/prereg/*.yaml` (one config per analysis, all written before any is run) |

**Eligibility gate for primary comparators.** Any embedding with **> 5 % all-zero P rows** (over the 1,092 perturbation targets, at the analysed K, or over graph nodes) cannot be a primary comparator. Computed for every embedding in `results/hic_v2/zero_rows.csv`. At K = 10: `b500_hic_trans` 0 %, all `b500_lin_randtrans` seeds 0 %, `b500_lin_trans` 0 % (**eligible**); `b500_lin_cis` 51 % (38 % of nodes), `b500_hic_cis` 56 %, per-component variants 48 %, `hic_local25k` 63 % (**ineligible**; secondary only).

## 3. Primary analysis

For each (perturbation *i*, split *s*) test case, `d_is = L2_C(i,s) − L2_A(i,s)` (positive = Hi-C arm better), where

> **`L2_C(i,s) = (1/10) Σ_{k=0…9} L2_k(i,s)`, the arithmetic mean of ten case-level L2 values, `L2_k(i,s)` being the L2 of control seed *k*'s *own* prediction for perturbation *i* in split *s*.**
> It is **not** the L2 of an averaged (ensembled) prediction: by Jensen's inequality that quantity is never larger, and would give the control an ensemble advantage the single Hi-C arm does not have.

This is how `src/stats.py::_per_case` computes it (it averages the per-case metric across models, never predictions), with a unit test that fails if predictions were averaged instead. Effect size **D = (1/5) Σ_s mean_{i ∈ test_s} d_is** (equal weight per split).

- **Test (single, confirmatory).** Two-sided **sign-flip permutation test**, 100,000 permutations, **α = 0.05**. All cases of one perturbation share one random sign (a perturbation recurs across splits); equivalent to testing `D = Σ_i a_i`, `a_i = (1/5) Σ_s d_is / n_s`. Reported with it: 95 % cluster-bootstrap CI over perturbations (10,000 draws), per-split mean differences with within-split paired t/Wilcoxon p (descriptive), fraction of perturbations with `a_i > 0`, and D for each individual control seed (descriptive; averaging over seeds does not model between-seed variance).
- **Descriptive companion (not a test, no p-value).** *Rank of the Hi-C arm among the control seeds:* per split, the mean test L2 of `b500_hic_trans` is ranked against the mean test L2 of each of the ten `b500_lin_randtrans` seeds (rank 1 = lowest L2 of the 11 values; also the number of control seeds beaten strictly; ties get the best rank); **pooled** = the same on the mean of the five split means (equal weight per split, as for D). Reported next to the primary result so that the position of the Hi-C arm within the control-seed spread is visible; it plays no role in the decision rule. Implementation: `src/stats.py::seed_rank_summary` (tested), returned by `primary_analysis` as `seed_ranks_descriptive`.
- **Smallest effect size of interest (SESOI).** `δ = 0.5 % × (mean control L2 over the primary test cases)`. Indicatively ≈ 0.021 L2 units (0.5 % of 4.256, the non-Hi-C mean-baseline L2 on this metric; for scale, the paper's PCA embedding gains 0.106 over the mean baseline). The exact δ is computed from the data by the fixed formula.
- **Equivalence (TOST).** Two one-sided sign-flip tests against the bounds ±δ (`H0: D ≤ −δ` and `H0: D ≥ +δ`), α = 0.05 each; equivalent if both reject (`src/stats.py::tost_signflip`).
- **Decision rule (pre-specified).**
  1. **Significant** — two-sided p < 0.05: *Hi-C better* if D > 0, *Hi-C worse* if D < 0. If TOST also rejects, report "significant but within the equivalence bounds (negligible)".
  2. **Equivalent** — p ≥ 0.05 and TOST p < 0.05 (D is within ±δ with confidence).
  3. **Inconclusive** — neither.
- **Two-tier interpretation rule (pre-specified).**
  - **Tier 1 — H1 supported** (all conditions below must hold): *"A spectral embedding of the Hi-C graph improves unseen-perturbation prediction relative to a connected, position-matched control."*
  - **Tier 2 — mechanistic claim.** Tier 1 holds **and** S3 is also **significant in the same direction** (*Hi-C better*), where S3 is the **apples-to-apples contrast `b500_hic_trans:resid` vs `b500_lin_randtrans:resid`**. Both arms are residualized with the *identical* procedure — chromosome + per-chromosome cubic B-spline of position (8 df per chromosome), regressed out on each embedding's own embedded nodes, residual columns re-scaled to unit RMS (`src/hic/position.py::residualize`) — and this includes **all ten control seeds**, so neither arm retains linear chromosome/position information. As in the primary, the control value of a case is the **mean of the ten per-seed case-level L2 values** (each from that seed's own residualized embedding and its own prediction, never an averaged prediction). Same test (two-sided sign-flip permutation, 100,000 permutations), **α = 0.05**, K = 10, λ = 0.1; **serial gatekeeping**: evaluated only after Tier 1 holds, so no multiplicity adjustment. Both residualized arms must pass the eligibility gate (verified: 0 % all-zero rows for `b500_hic_trans:resid` and for each of the ten `b500_lin_randtrans:resid` seeds, at K = 10 and 50). *Why residualize the control too:* the raw control still carries chromosome/position information (R² 0.55), so contrasting a residualized Hi-C arm with it would bias Tier 2 against Hi-C whenever position is predictive. Only if Tier 2 holds may we claim that the gain comes from 3D information **beyond chromosome identity and 1D position**. Tier 1 without Tier 2 is reported as "an improvement that cannot be attributed to 3D organisation beyond chromosome/position". Rationale: §8.
- **Condition for Tier 1 (H1 supported).** All of:
  1. the decision rule outcome is **"Significant: Hi-C better"** (D > 0, p < 0.05); and
  2. **D exceeds all three null distributions**, where for each null draw *j* the null effect is `D_j = (1/5) Σ_s mean_{i ∈ test_s} [ L2_C(i,s) − L2_null_j(i,s) ]` — the same pairing, the same control `L2_C` and the same equal weight per split as `D` (`src/stats.py::null_effect`), with `L2_null_j` the case-level L2 of the *same* linear model (K = 10, λ = 0.1) using null embedding *j* in place of the Hi-C embedding:
     - (a) `D` > the **95th percentile** (`numpy.quantile`, linear interpolation) of the **gene-label-permutation** null, 1,000 seeds (which gene sits on which embedded node is permuted; rows of the same graph's embedding are permuted, **no new eigendecomposition**);
     - (b) `D` > the **95th percentile** of the **random-Gaussian** null, 1,000 seeds (i.i.d. N(0,1) rows with matched dimensions on the same embedded nodes; genes sharing a bin share a row);
     - (c) `D` > the **maximum over all 20 seeds** of the **degree-preserving-rewire** null (the Hi-C graph rewired by double-edge swaps preserving every node degree and the weight multiset, then re-embedded).
  For each null the **empirical p** is reported: `p = (1 + #{D_j ≥ D}) / (1 + n_null)` (smallest attainable: 1/1001 for the permutation and random nulls, 1/21 for the rewire null). The same procedure at fixed K = 50 is reported for the key secondary analysis. Implementation: `src/stats.py::h1_support`, `src/nulls.py` (the lean evaluator fits the identical model with everything seed-independent cached; estimated runtime on this machine ≈ 2 min for both 1,000-seed nulls at K = 10 and 50 together (timed on synthetic data and on P-matrix construction only), ≈ 12 s per rewire seed ≈ 4 min for 20; the standard evaluator would need ≈ 4.5 h per 1,000 seeds). `python -m src.nulls --config configs/prereg/s2_nulls.yaml`.


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
- **S2 Nulls** for `b500_hic_trans` (part of the Tier-1 condition, §3): gene-label permutation (**1,000** seeds), random Gaussian with matched dimensions (**1,000**), degree-preserving rewire (**20**); null effects `D_j` and empirical p as defined in §3; run with `src/nulls.py` at K = 10 (decisive) and K = 50 (key secondary).
- **S3 Position-residualized** (secondary, not co-primary; **its `b500_hic_trans:resid` vs `b500_lin_randtrans:resid` contrast is the Tier-2 test of §3**): `:resid` (chromosome + per-chromosome cubic B-spline, 8 df regressed out, columns re-scaled — the identical procedure for every arm) for `b500_hic_trans`, `b500_hic_cis` and the ten control seeds `b500_lin_randtrans_s0…s9`; control value = mean of the ten per-seed case-level L2 values. The `b500_hic_cis:resid` contrast (vs the residualized control) stays **exploratory** (`b500_hic_cis` fails the zero-row gate, 56 % at K = 10). Also reported descriptively: `b500_hic_trans:resid` vs the *un*residualized control (biased against Hi-C when position is predictive; not used for any claim).
- **S4 Extended eval set.** Hi-C and control arms vs the **mean baseline only** on all 273 test perturbations per split; PCA embedding on its own predictable subset (≈ 100–105 per split) for reference (S9).
- **S5 Metrics.** Pearson delta; top-20 DE; all genes; target gene included.
- **S6 Cis vs trans read-out genes.** Metrics separately on read-out genes on the perturbed gene's chromosome vs others.
- **S7 K and tuning.** K ∈ {10, 20, 50} (λ = 0.1) and validation-tuned K.
- **S8 Resolution.** 250 kb analogue of the primary contrast (`b250_hic_trans` vs `b250_lin_randtrans`, 10 seeds) and the gene-level 25 kb local graph. *Caveat:* at 250 kb the Hi-C graph has a lower median degree than the controls (4 vs 10 for cis-only cells), so a difference there is partly a sparsity difference.
- **S10 K562 rearrangements** *(region definition approved by the project owner, 2026-10-14)*. K562 is highly rearranged (BCR-ABL t(9;22), amplifications) and translocations appear as trans contacts. Re-run the primary contrast with bins in rearranged/amplified regions removed **before balancing** (`b500x_*` graphs) and the perturbations whose target falls in those bins dropped from all arms. **Region list** (`configs/k562_rearranged_regions_hg38.bed`, 94 regions, built by `scripts/make_k562_regions.py`): (i) all 70 segments called `amp` (196.4 Mb) in ENCODE **ENCFF486MJU** (ENCSR000BBO; HAIB Illumina genotyping array copy-number segments, hg19; md5 `b4a78fa77789291cb5e3d528b32aa55a`); (ii) ±1 Mb windows around every inter-chromosomal breakpoint (24) in ENCODE **ENCFF863MPP** (ENCSR053AXS; 10X linked-read SV calls, grocsvs, hg19; md5 `04f1c1f3aad184c0e05c863da7a7d3e0`); lifted hg19→hg38 with UCSC `hg19ToHg38.over.chain.gz` (md5 `35887f73fe5e2231656504d1f6430900`; fetched from `hgdownload.soe.ucsc.edu`, because the TLS certificate of `hgdownload.cse.ucsc.edu` did not verify). Effect: 513 of 6,073 bins excluded (8.4 %), 56 of 1,092 perturbations dropped (256–266 test perturbations per split remain). The thresholds are the sources' own calls (no tuning); the window size (±1 Mb) was fixed a priori.
- **R1 Replication in RPE1** (planned; availability check only — **nothing has been downloaded and nothing will be until the project owner approves**). Same primary analysis on hTERT RPE-1. *Hi-C:* ENCODE experiment **ENCSR846DUS** (in situ Hi-C, Aiden lab, 4 replicates merged), GRCh38 `.hic` files `ENCFF876OWE` (MAPQ-thresholded, 2.24 GB, md5 `b6828cb642d1fce2b64c5383f22ec4f9`) and `ENCFF274QIY` (2.24 GB, md5 `68a2ad184e63fc5e1fe22e263a673aec`); remote reading is possible as for K562. This map is ~10× smaller than the K562 file, so a depth/ICE QC gate precedes any analysis (R1 is declared *not evaluable* if the valid-bin fraction or the zero-row gate fails). *Perturb-seq:* GEARS `replogle_rpe1_essential.zip`, 665,914,266 B, md5 `169370a8da6093470f122172ab30f7c5`. Success criterion: same direction and two-sided p < 0.05 under the identical primary test and decision rule. RPE1 is near-diploid, so no rearrangement exclusion is planned.

## 7. Pre-specified sanity checks (before interpreting anything)

1. `train_pca_expr` and the mean baseline reproduce the earlier `k562_gears` results within noise under the primary code path.
2. **Each random-Gaussian arm's mean L2 differs from the mean baseline by less than δ (the SESOI)** (mean over the primary test cases, equal weight per split, on the same cases; δ as defined in §3). **Reported, not a gate:** the report states the fraction of the 1,000 random arms satisfying it and the range of differences.
3. Predictable sets are identical across all bin-level arms of a contrast (verified: 1,092/1,092 at 500 kb; 1,036 for `b500x_*`).
4. The eligibility gate (§2) holds for every arm of the primary contrast.
5. The lean evaluator (`src/nulls.py`) reproduces the case-level L2 of the primary arms obtained from the standard evaluator (`src/evaluate.py`) to within 1e-8 (the equivalence is also unit-tested on synthetic data). If it does not, the null analysis is not interpreted until the discrepancy is resolved and logged in `DEVIATIONS.md`.

## 8. Known limitations

K562 is aneuploid and rearranged; Hi-C is bulk (Rao et al. 2014, reprocessed on GRCh38 by ENCODE) with no KR vectors (own ICE); the paper's own GEARS splits could not be recovered (ours are regenerated with the GEARS procedure); absolute L2 is ~4–6 % above the paper's (Pearson and ratios agree); **why Tier 1 does not license a mechanistic claim:** the Hi-C embedding encodes chromosome identity far more than the control does (R² on chromosome alone ≈ 0.4 for `b500_hic_trans` vs 0.11 for the control; with position 0.69 vs 0.55), and K562 is aneuploid, so chromosome-level copy number could drive a Tier-1 gain without any 3D mechanism (we do not assert that it does; the design cannot exclude it without S3). Hence the claim that the gain reflects 3D information beyond chromosome identity and 1D position requires Tier 2 (S3); the A/B track in this map is itself largely positional; control seeds are averaged and their between-seed variance is not part of the test; a perturbed gene's response and its Hi-C neighbours share the chromosome, so "Hi-C helps" and "chromosome helps" are separated only by the control design and S1/S3.

## 9. Housekeeping

- **Git history was rewritten once, on 2026-10-14, before this document was frozen or any model was run**, to give every commit the author and committer identity `Zeynep <61983919+zeynepthereal@users.noreply.github.com>` (`git filter-repo --mailmap`; trees, timestamps, messages and parent structure unchanged; HEAD tree hash identical before and after). The 4 commit-hash references in `results/` (3 run logs, the manifest) were remapped through filter-repo's commit-map. No commit hash recorded in this document predates the rewrite. The remote repository is recreated empty and the rewritten history pushed afresh.
- **Freeze procedure (only after the project owner's explicit approval):** (1) regenerate `results/hic_v2/MANIFEST.json` on a clean tree (SHA-256 of every graph, embedding, split, config, region list and source file; sources' MD5s); (2) commit; (3) create the **annotated tag `prereg-v1`**; (4) push the commit **and** the tag to GitHub; (5) the project owner uploads the frozen `PREREGISTRATION.md` to **OSF** as an independent timestamp — the OSF link is recorded in `DEVIATIONS.md` (not in this file, which is frozen at that point).
- Deviations after freeze are logged in a `DEVIATIONS.md`, never edited into this file.
