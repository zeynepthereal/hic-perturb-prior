# hic-perturb-prior

Does 3D genome architecture (Hi-C) help predict transcriptional responses to *unseen* genetic perturbations?
**Week 1: reproduce the linear baseline of Ahlmann-Eltze, Huber & Anders (Nat Methods 2025) on Replogle 2022 K562 essential CRISPRi. No Hi-C yet.**

## Run
```bash
mamba env create -f environment.yml && mamba activate hic-perturb-prior
python -m src.run_baseline --config configs/k562_gears.yaml --preprocess            # GEARS-processed data, splits 1-5 (candidate primary)
python -m src.run_baseline --config configs/k562_scperturb_5000hvg.yaml --preprocess # scPerturb raw counts (secondary)
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

### GEARS-processed data (candidate primary benchmark)
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

## Known deviations from the paper
- scPerturb pipeline: own scanpy-seurat HVG replica (5000 genes), **no perturbation QC filtering** (QC only logged); not expected to match absolute numbers (see reproduction status).
- GEARS pipeline: paper's exact splits unavailable (see above); linear model and mean baseline fit on the training split only (`refit_trainval: false`), as in the paper.
- Early scPerturb results from commit-less / dirty-flag runs were superseded and deleted; all kept runs have a clean `git.json`.
- Paper's 'val' perturbations are unused by its linear model; ours use them to tune K/lambda (secondary column only).
