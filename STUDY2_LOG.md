# STUDY 2 LOG

Append-only working log of **Study 2** (K562 discovery → RPE1 confirmation). Study 1 (`prereg-v1`, `PREREGISTRATION.md`, its configs, manifest, `DEVIATIONS.md`) is **closed and untouched**; nothing here modifies it.
Everything in this file is **exploratory** until a new pre-registration (`prereg-v2`) is frozen. Times are 2026-10-14 local unless stated.

## Rules for Study 2
1. **K562 = discovery data** (exploratory, may be analysed repeatedly). **RPE1 = confirmation data**, to be analysed **once**, under a new frozen pre-registration (`prereg-v2`).
2. **Until `prereg-v2` is frozen: no RPE1 Perturb-seq data is loaded, inspected or summarised in any way** (the GEARS RPE1 zip is not downloaded). Only the RPE1 *Hi-C* QC gate may be computed (Phase 1, item 3).
3. **Phase 1 = calibration only: no Hi-C vs control comparison of any kind on real data.** (The graph-neighbour model has only been run on synthetic data; GEARS has been run only with its standard GO graph.)
4. Model families: (1) the existing linear model; (2) a graph-neighbour model (new, `src/graphnn.py`); (3) GEARS (cell-gears 0.1.2) with its GO graph replaced/combined with the Hi-C graph.
   Information sources per family: standard (PCA / GO), standard + Hi-C, standard + matched control graph, mean baseline.

## Data-access ledger
| data | accessed? | note |
|---|---|---|
| K562 GEARS-processed `perturb_processed.h5ad` (Study 1 file, md5 of the zip `84eed779…`) | yes | copied locally from the Study 1 scratch workspace (no network); GEARS cell-graph cache built next to it (6.6 GB, git-ignored) |
| GEARS resources `gene2go.pkl`, `essential_all_data_pert_genes.pkl`, `go_essential_all.tar.gz` (Harvard Dataverse, same source as the K562 zip) | yes | see P1.1 |
| RPE1 Hi-C ENCODE ENCFF876OWE (GRCh38, 2,240,397,369 B, md5 `b6828cb642d1fce2b64c5383f22ec4f9`) | yes, **remote read, QC only** | 122 MB network traffic; the .hic was never downloaded |
| **RPE1 Perturb-seq (GEARS `replogle_rpe1_essential.zip`)** | **NO — not downloaded, not loaded, not inspected** | |

## P1.1 GEARS reproduction on K562 (standard GO graph, default benchmark hyper-parameters)
- **Resources.** GEARS' standard GO path needs three small files from Harvard Dataverse (same host as the approved K562 file). Downloaded manually (not by GEARS) into `data/gears_resources/`, sizes and MD5 verified against Dataverse:
  `gene2go.pkl` (id 6153417, 9,462,558 B, md5 `77c9af0c61c30ea4d7a85680f4d122dc`), `essential_all_data_pert_genes.pkl` (6934320, 558,811 B, md5 `b7bc2a91ca513b86f27f090d963711a6`), `go_essential_all.tar.gz` (6934319, 60,654,049 B, md5 `b8bff0d53407f26648330d264df7fe16`).
  Safety inspection before use: `gene2go.pkl` pickle has no code-capable opcodes; the other pickle only reconstructs a numpy array; the tar has two members with safe paths (extracted by us). GEARS then found local copies and made **no network access**.
- **Hyper-parameters** = the Ahlmann-Eltze et al. benchmark script `run_gears.py`: `epochs=20`, `hidden_size=64`, batch 32 / test batch 128, custom split from our GEARS split JSON, `np.random.seed(1)` (+ `torch.manual_seed(1)`, our addition). Driver: `scripts/study2_gears_run.py`; evaluation: `scripts/study2_gears_eval.py`.
- **Environment patch (one line).** cell-gears 0.1.2 indexes a scipy sparse matrix with a pandas Series, which scipy 1.15 rejects; patched `X[obs.condition == 'ctrl']` → `X[(obs.condition == 'ctrl').values]` (`scripts/study2_patch_gears.py`, diff in `results/study2/gears_env_patch.diff`). Semantics unchanged.
- **Split adaptation.** Five perturbations are not in GEARS' GO graph (C7orf26, C14orf178, RPS10-NUDT3, SEM1, FAU) and are dropped by `PertData.load`; the paper made its splits after that filtering, ours were made on the unfiltered file. They were removed from the split lists before training. Split 1: dropped from test: C14orf178, RPS10-NUDT3, FAU (3); from train: C7orf26, SEM1 (2).
- **Data preparation (one-off).** 162,751 cells × 5,000 genes (10,691 control cells; 1,093 conditions). Building the per-cell PyG graphs took ≈ 3 min and ≈ 6 GiB RSS; the cache is 6.6 GB (`data_pyg/cell_graphs.pkl`). Later loads take 15 s. Loader sizes (split 1): train 3,508 / val 379 / test 1,184 batches.
- **Device.** CUDA: not available. **MPS: available and usable** (GEARS trains and tests on it with `PYTORCH_ENABLE_MPS_FALLBACK=1`). A first in-process benchmark (40 steps + limited evaluation through GEARS' own loop) gave CPU 435.9 s vs MPS 443.5 s, which only measured GEARS' single-threaded evaluation code and could not separate the devices. Real training rates, measured by sampling the step counter over 120 s windows:
  **CPU 3.33 steps/s (0.30 s/step, 17.5 min per epoch of 3,508 steps → 5.8 h for 20 epochs, training only); MPS 5.00 steps/s (0.20 s/step, 11.7 min per epoch → 3.9 h for 20 epochs, training only).** MPS was used for the real run; a CPU attempt (≈ 9 min) was aborted and is kept in `results/study2/_aborted_cpu_run_gears_seed1`.
- **Single-split end-to-end run (split 1, MPS):** STATUS — running (projected > 2 h, so per the instruction only this one split is run). Result, wall time and peak memory: see the entry below once finished.
- Concurrency note: unit tests, the RPE1 Hi-C read/QC (network + a few minutes of CPU) and writing code ran concurrently with the GEARS run; the GPU-bound training rate was sampled while the RPE1 read was already finished.

## P1.2 Graph-neighbour model (`src/graphnn.py`, `src/test_graphnn.py`) — synthetic data only
- Prediction for a test perturbation = weighted average of the deltas of **training** perturbations whose targets are graph neighbours; weights `S = w0·I + Σ_(h=1..hops) γ^(h-1) P^h`, `P = D⁻¹W` (edge-weighted, k-hop); fallback to the mean training delta if no training perturbation is reachable; targets without a node → NaN row (not predictable, as in the eval_set mechanism). No fitted parameters; defaults `hops=2, γ=0.5, w0=1`. Interface `fit_predict_graph_neighbour(ds, train_idx, test_idx, GraphSource)` mirrors `fit_predict_linear`; `GraphSource.from_hic_graph(name, ds)` reuses the Study 1 bin graphs.
- Tests: planted modular graph with module-specific responses → L2 < 0.45 × mean baseline and correlation with the module signature > 0.95; random graph (5 seeds) → never better than the mean baseline; exact edge-weight arithmetic on a path graph (1 and 2 hops), isolated-node fallback = mean baseline exactly, NaN rows for unmapped targets, same-node weight; dataset interface. **Not run on real data.**

## P1.3 RPE1 Hi-C QC gate only (ENCODE ENCFF876OWE, GRCh38, read remotely; no Perturb-seq data)
Header: hg38, 26 chromosome records (incl. `All`), resolutions 100 bp … 2.5 Mb (same as K562). All 276 chromosome-pair blocks read at 500 kb (8,110,035 records) and 250 kb (12,427,384 records); traffic **122 MB** (limit 1 GB, not reached). Same code as Study 1 (`src/hic/bins.py`, same ICE and coverage filter); K562 reference from the cached Study 1 blocks.
| | RPE1 500 kb | RPE1 250 kb | K562 500 kb | K562 250 kb |
|---|---|---|---|---|
| total contacts (sum of raw blocks) | 62,479,963 | 62,479,963 | 1,908,747,765 | 1,908,747,769 |
| bins | 6,073 | 12,137 | 6,073 | 12,137 |
| valid bins (ICE coverage filter) | 5,509 | 11,005 | 5,861 | 11,775 |
| **valid-bin fraction** | **90.71 %** | **90.67 %** | 96.51 % | 97.02 % |
| bins with zero coverage | 235 | 498 | 0 | 0 |
| median coverage per valid bin | 17,742 | 9,183 | 484,506 | 246,692 |
| coefficient of variation of coverage | 0.154 | 0.159 | 0.528 | 0.587 |

Graph-level zero-row gate on RPE1 at 500 kb (node level; Study 1 recipe; the perturbation-target level cannot be evaluated without RPE1 target lists, which are Perturb-seq-derived and off limits): `hic_trans`: 1 component, 564 isolated (invalid) bins, all-zero embedding rows at K=10 / 50: 0.0 % / 0.0 %; matched random-trans control (seed 0): 1 component, 0.0 % / 0.0 %.
**Gate (PREREGISTRATION.md §6 R1):** the text says R1 is "not evaluable if the valid-bin fraction or the zero-row gate fails" but **gives no numeric threshold for the valid-bin fraction**, so a formal pass/fail cannot be declared from the registration. The zero-row part passes at node level (0 %, ≤ 5 % rule of §2). The valid-bin fraction is 90.7 % at both resolutions (K562: 96.5 % / 97.0 %) with ≈ 30× fewer contacts (62.5 M vs 1.91 G; median 17.7 k vs 485 k per 500 kb bin). A numeric threshold must be fixed in `prereg-v2` before RPE1 is analysed.
