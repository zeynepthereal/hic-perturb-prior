# hic-perturb-prior

Does 3D genome architecture (Hi-C) help predict transcriptional responses to *unseen* genetic perturbations?
**Week 1: reproduce the linear baseline of Ahlmann-Eltze, Huber & Anders (Nat Methods 2025) on Replogle 2022 K562 essential CRISPRi. No Hi-C yet.**

## Run
```bash
mamba env create -f environment.yml && mamba activate hic-perturb-prior
python -m src.run_baseline --config configs/k562_gears.yaml --preprocess            # PRIMARY: GEARS-processed data, GEARS-procedure splits 1-5
python -m src.run_baseline --config configs/k562_scperturb_5000hvg.yaml --preprocess # scPerturb raw counts (SECONDARY)
# GEARS splits (separate env, needs data/gears/...): mamba env create -f envs_gears.yml; mamba run -n hic-perturb-gears python scripts/make_gears_splits.py --seeds 1 2 3 4 5
python -m src.compare_paper results/<run_dir>
pytest src/test_synthetic.py
```
Each run writes `results/<timestamp>_<name>/` with `config.yaml`, `git.json` (commit + dirty flag; `NO_COMMITS` until the first commit),
metrics, hyper-parameter grid, QC and tables. QC for all perturbations: `results/qc_perturbations.csv`.

## Data
`ReplogleWeissman2022_K562_essential.h5ad` from scPerturb (Zenodo record 10044268, v1.3),
URL `https://zenodo.org/api/records/10044268/files/ReplogleWeissman2022_K562_essential.h5ad/content`,
1,546,729,675 bytes, MD5 `d8cba17576d1a8afc0f7d71b79cad0f7` (verified). 310,385 cells x 8,563 genes, raw UMI counts (dense float32),
perturbation column `perturbation` (= target gene symbol), control label `control` (10,691 cells), 2,057 targeted genes.

### GEARS-processed data (PRIMARY benchmark, `configs/k562_gears.yaml`)
`replogle_k562_essential.zip` from Harvard Dataverse, `https://dataverse.harvard.edu/api/access/datafile/7458695`, 669,703,049 bytes,
MD5 `84eed779531de88f83d3eb16a773a261` (= Dataverse's published MD5), SHA-256 `e769c89aa876061fdb3ac02a8de274f19a741a4f95decf14a52fd621d1eea319`.
Contains only `perturb_processed.h5ad` (1,092 single perturbations incl. 680 whose target is not among its 5000 genes; 412 evaluable by the PCA model).
Splits: `scripts/make_gears_splits.py` calls GEARS 0.1.2's `DataSplitter` exactly as `prepare_split('simulation', seed)` does, without `PertData.__init__`
(which would download/unpickle an unrelated `gene2go_all.pkl`). Output: `configs/splits/gears_k562_seed{1..5}.json` (train 738 incl. ctrl / val 82 / test 273 conditions).
**These do NOT coincide with the paper's own splits**: the paper's evaluated test perturbations per seed overlap ours at chance level for every seed in {0..10, 42}.
GEARS splitting is deterministic (verified across PYTHONHASHSEED), so the paper must have used a different input file version or code path; its split JSONs are not published.

## Method (src/)
- `data.py`: streaming `normalize_total(1e4)` -> `log1p` -> top-N HVG -> pseudobulk mean per perturbation -> `delta = pseudobulk(pert) - pseudobulk(control)`.
- `baselines.py`: no-change, mean, `Y_hat = G W P~' + b` (ridge, closed form as in the paper's Methods / `solve_y_axb`). Embeddings are plug-ins in `EMBEDDINGS`
  `(Y_train, train_idx, ds, K) -> (G, P_all)`; add GO / Hi-C / shuffled / random P there, model and evaluation stay unchanged.
- `evaluate.py`: L2 and Pearson (pred vs observed delta) per perturbation; subsets `all`, `top20_de` (oracle: largest |observed delta|), `top1000_expr` (paper's metric:
  most expressed genes in control). Each perturbation's own target gene is excluded by default (`eval.exclude_target_primary`); "included" metrics are always reported too.
  Perturbations are split (never cells) 70/10/20, 5 seeds; K, lambda tuned on val (grid in config); final models refit on train+val.

## eval_set (every results table states which one it uses and its size)
- `common` (primary): test perturbations predictable by ALL compared embeddings (mean/no-change predict everything).
- `extended:<emb>`: all test perturbations predictable by `<emb>`, compared against the mean baseline only (emitted only when it differs from `common`).
- `<name>` (optional, `eval.extra_eval_sets`): `common` restricted to a perturbation list (diagnostics).
Embeddings declare `predictable(ds)`; G is fit on all training perturbations, W and b only on the training perturbations the embedding can embed (as in the paper's code).

## Reproduction status (see results/*/comparison_vs_paper.md)
Paper (Replogle K562, top-1000 expressed genes, target included): mean L2 4.119 / Pearson 0.404, linear 3.996 / 0.433 (ratios 0.970 / 1.072).
- GEARS data + GEARS-procedure splits 1,2 (`k562_gears_repro`): mean 4.267 / 0.404, linear 4.153 / 0.435 (ratios 0.973 / 1.075). 5 splits (`k562_gears`): mean 4.350 / 0.411, linear 4.251 / 0.438.
  Pearson and ratios reproduce; absolute L2 is ~3.6% (2 splits) to ~6% (5 splits) higher than the paper (within split-to-split variation, but systematic in sign).
- scPerturb raw counts, own preprocessing, all 797 evaluable perturbations: Pearson 0.28 / 0.31 (mean / linear). Restricted to the perturbations GEARS kept (`gears_pertset`): 0.41 / 0.45-0.47.
  => the Pearson gap is explained by GEARS' perturbation filtering; the remaining L2 offset (3.4-3.6 vs 4.15-4.35) is a normalisation/scale difference.

## Leakage audit
- **G (PCA) and W** are fit on training perturbations only (`select_hparams`: train; final fit: train+val). Test deltas never enter PCA, W or b.
- **HVG selection** (and cell normalisation) use *all* cells, including cells of test perturbations, but no perturbation labels (unsupervised, same as the paper/GEARS
  which preprocess the full dataset first). A 2nd-order effect: the HVG set decides which perturbations are evaluable (target must be a read-out gene).
- `top1000_expr` uses control cells only; `top20_de` uses the observed test delta by design (oracle subset, evaluation only, never fed to a model).
- Hyper-parameters are chosen on val perturbations only. The "paper" column uses fixed K=10, lambda=0.1.

## Week 2 part 1: Hi-C gene graphs (no perturbation model run with them yet)
Pipeline: `scripts/build_gene_table.py` -> `scripts/fetch_hic.py cis|trans` -> `scripts/build_graphs.py` -> `scripts/diagnose_hic.py`
(code in `src/hic/`; outputs `results/hic/` incl. `diagnostics.md`, `mapping_report.txt`, `graph_stats.json`, `figs/`; heavy artifacts in gitignored `data/hic/`).
- **Hi-C**: ENCODE `ENCFF616PUW` (experiment ENCSR545YBD, K562 in situ Hi-C of Rao et al. 2014, **GRCh38**, 4 replicates merged, MAPQ-thresholded),
  22.64 GB, ENCODE md5 `d4b999471d4c17a74f1e66370c0012e3`. **Read remotely** with hicstraw (HTTP range requests); the file was not downloaded
  (checksum therefore not verifiable locally; ~0.9 GB total network traffic for all reads). Raw counts are fetched; balancing is our own ICE (the file has no KR vectors; VC/SCALE exist).
  Other candidates (not used): GRCh38 ENCFF080DPJ 20.7 GB / ENCFF662ZGL 22.5 GB / ENCFF539KNF 24.3 GB; hg19 per-replicate files 2.3-4.5 GB (ENCFF013TGD etc.).
- **Annotation**: GENCODE v44 primary assembly basic (GRCh38), MD5 `fc4731226d2d315e6b2415fd1690403d`; TSS = 5' end of the gene; chr1-22,X only (K562 is female; no chrY/chrM signal).
- **Nodes**: 19,970 protein-coding genes + 293 extra GEARS genes = 20,263. GEARS genes: 4,986/5,000 read-outs map (14 fail: 10 chrM, 4 unplaced contigs), 1,092/1,092 targets map.
- **Graph parameters** (`src/hic/graphs.py: GRAPH_CFG`): cis 25 kb, ICE, O/E by per-diagonal mean, edges within 5 Mb with >= 3 raw reads; trans 500 kb genome-wide ICE, O/E by chromosome-pair mean;
  kNN rule identical for all variants: every node keeps k_cis=10 strongest cis edges (+ k_trans=3 trans partner bins for `hic_cis_trans`, one random gene per partner bin), union-symmetrized, weights clipped at the 99.9th percentile.
  `linear_distance`: same per-node degree as the hic_cis directed selection, nearest genes by TSS distance, weight exp(-d/d0).
- **Spectral**: normalized Laplacian, each connected component decomposed separately and its trivial zero mode dropped (ARPACK cannot resolve the degenerate zero eigenspace: it returned 3 of 39),
  first 50 non-trivial eigenvectors saved, columns scaled to unit RMS. Isolated nodes get no embedding.
- **Embeddings** (`src/hic/embeddings.py`, names usable in `eval.embeddings`): `hic_cis`, `hic_cis_trans`, `linear_distance`, plus nulls `<variant>:perm:<seed>` (gene-label permutation),
  `hic_cis_trans:rewire:<seed>` (degree-preserving rewire, built lazily), `<variant>:random:<seed>` (Gaussian, matched dims). G stays the training-data PCA; only P changes.
- **Diagnostics summary** (full tables in `results/hic/diagnostics.md`):
  - components / isolated nodes / edges: hic_cis 39 / 98 / 135k; hic_cis_trans 1 / 7 / 195k; linear_distance 111 / 31 / 119k. Median degree 12 / 16 / 11.
  - R² of the first 10 embedding dims explained by chromosome + per-chromosome cubic position: hic_cis 0.81, hic_cis_trans 0.51, linear_distance 0.62 (random reference 0.005);
    chromosome alone: 0.005 / 0.22 / 0.00. **Flag: all three exceed 0.5; the cis-only Hi-C embedding is more positional than the linear-distance control.**
    Only hic_cis_trans has a substantial non-positional component (first 50 dims: 0.32 explained).
  - cis-only graphs have tiny Laplacian eigenvalues (1e-4..5e-3, near-degenerate, localized on components) so their dimensions are unstable; hic_cis_trans eigenvalues 0.016..0.22.
  - hic_cis vs linear_distance edges: Jaccard 0.14 (23% of hic_cis edges also in the control); Spearman(O/E, -distance) = 0.06 over all 1.06M candidate pairs (decay removed), 0.016 on shared edge weights.
  - predictable GEARS perturbations: PCA 412/1092 (about 103/273 per test split); every Hi-C variant 1092/1092 (273/273 per split) -> extended eval sets ~2.6x larger.

## Week 2 part 2: bin-level whole-genome graphs, factorial controls, pre-registration draft
Supersedes the gene-level `hic_cis_trans` / `linear_distance` graphs of part 1 (random trans gene picks, 10 Mb cis band); `hic_cis` is kept as `hic_local25k`.
Pipeline: `scripts/fetch_hic.py trans 500000|250000` (all chromosome pairs, full cis range; remote read) -> `scripts/build_bin_graphs.py 500000|250000`
-> `scripts/diagnose_hic_v2.py`; outputs in `results/hic_v2/` (`diagnostics.md`, `structure.csv`, `r2_chrom_position.csv`, `ab_correlation.csv`, `ab_*k.csv`, `figs/`, `MANIFEST.json`).
- **Traffic**: 500 kb all-pairs ~16.3 M records (read earlier); 250 kb all-pairs 55.4 M records, **192 MB** (estimate before fetching: 0.3-0.75 GB). No local copy of the .hic.
- **Bin-level graph** (`src/hic/bins.py`): genome-wide ICE (cooler-style low-coverage filter), cis O/E by diagonal (per chromosome, far diagonals with < 20 pairs use the plateau),
  trans O/E by the global trans mean; nodes = valid bins; a gene takes the embedding of its TSS bin. Sparsification (same for all cells): kNN k=10 (union, each node keeps at most 20 edges)
  U per-chromosome maximum spanning forest (cis); kNN k=10 (trans). Union/mutual kNN were tried first: plain union gives hubs of degree 2,428, mutual kNN is far too sparse for Hi-C.
- **Factorial cells**: `b500_hic_trans`, `b500_hic_cis`, `b500_lin_trans`, `b500_lin_cis` (linear control = chain per chromosome, weight exp(-d/5 Mb)); also `b250_*`. All cells share the node set;
  the cis-only cells have exactly 23 components, the trans cells 1.
- **Nulls / variants** (names for `eval.embeddings`): `<cell>:perm:<seed>`, `<cell>:rewire:<seed>`, `<cell>:random:<seed>`, `<cell>:resid` (chromosome + spline position regressed out).
- **A/B compartments**: PC1 of the per-chromosome O/E correlation matrix oriented by gene density. In this rearranged K562 map PC1 does not track gene density (median r 0.02; best of PC1-3: |r| 0.46 at 500 kb),
  so both PC1 and best-of-3 are reported; robustness variants (clipping, log, amplified-bin exclusion) did not change this.
- **Diagnostic headline** (500 kb, first 10 dims, R² on chromosome + position spline): hic_trans 0.69, hic_cis 0.49, lin_trans 0.74, lin_cis 0.997, local25k 0.91, random 0.03.
  A/B signal beyond position (best-of-3, 10 / 50 dims): hic_trans 0.16 / 0.63, lin_trans 0.09 / 0.29, hic_cis 0.03 / 0.30, lin_cis 0.004 / 0.03.
  The control `lin_cis` is localized (62 % of nodes non-zero in its first 10 dims). Details and open decisions in `PREREGISTRATION.md` (DRAFT, not frozen; no Hi-C model has been run).
- `src/stats.py`: paired per-perturbation comparison + sign-flip permutation test for the pre-registered primary analysis (tested). `configs/prereg/*.yaml`: one config per pre-registered analysis, **not to be run before freeze**.

## Week 2 part 3: connected control, equivalence testing, pre-registration v0.4 (DRAFT, not frozen)
- **Primary control replaced**: `b500_lin_cis` has 51 % all-zero P rows at K=10 (one component per chromosome => the slowest modes sit on a few chromosomes). New control
  `b500_lin_randtrans_s0..s9`: linear chain per chromosome + trans edges from a degree-preserving, inter-chromosomal-constrained randomization of the Hi-C trans graph
  (`src/hic/bins.py::rewire_trans`; per-node trans degrees, edge count and weight multiset identical; connected; 0 % zero rows). Control prediction = per-case L2 averaged over the 10 seeds.
- **Gate**: `results/hic_v2/zero_rows.csv`: fraction of all-zero P rows for every embedding (> 5 % => not eligible as primary comparator).
- **Statistics** (`src/stats.py`): paired sign-flip test (perturbation-level; control = mean of per-seed case-level L2, never an averaged prediction), SESOI = 0.5 % of the control mean L2, TOST equivalence, decision rule significant / equivalent / inconclusive, descriptive rank of the Hi-C arm among the control seeds, `primary_analysis`.
- **Extra variants**: `:percomp` (per-component eigenvectors, proportional allocation) for cis-only cells; `b500x_*` = K562 amplified/translocated regions removed before balancing
  (`configs/k562_rearranged_regions_hg38.bed`, sources and md5s in its header; `scripts/make_k562_regions.py`); `eval.extra_fixed_hp` for fixed K=20/50 models.
- **Nulls** (`src/nulls.py`, `configs/prereg/s2_nulls.yaml`, not to be run before freeze): lean evaluator for 1,000 permutation + 1,000 random-Gaussian + 20 rewire seeds; H1 support = D above the 95th percentile of the first two and above all rewire seeds; two-tier interpretation (Tier 2 needs S3 `:resid`).
- **Diagnostics**: `results/hic_v2/diagnostics_randtrans.md`, `zero_rows*.csv`. RPE1 availability (not downloaded) is documented in `PREREGISTRATION.md` (R1).
- `PREREGISTRATION.md` v0.4 (DRAFT) and `configs/prereg/` are the documents to approve; nothing has been run on perturbation data.

## Primary config and caveats
`configs/k562_gears.yaml` is the primary benchmark (directly comparable to the literature). Caveats to carry into every comparison:
(1) absolute L2 is ~4-6% above the paper's (Pearson and linear/mean ratios match); (2) the paper's own train/test splits are not recoverable
(ours are regenerated with the GEARS procedure, seeds 1-5); (3) the PCA model can only predict perturbations whose target is a read-out gene (412 of 1092).
The scPerturb pipeline (`configs/k562_scperturb_5000hvg.yaml`) is kept as a secondary config.

## Known deviations from the paper
- scPerturb pipeline: own scanpy-seurat HVG replica (5000 genes), **no perturbation QC filtering** (QC only logged); not expected to match absolute numbers (see reproduction status).
- GEARS pipeline: paper's exact splits unavailable (see above); linear model and mean baseline fit on the training split only (`refit_trainval: false`), as in the paper.
- Early scPerturb results from commit-less / dirty-flag runs were superseded and deleted; all kept runs have a clean `git.json`.
- Paper's 'val' perturbations are unused by its linear model; ours use them to tune K/lambda (secondary column only).
