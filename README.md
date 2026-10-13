# hic-perturb-prior

Does 3D genome architecture (Hi-C) help predict transcriptional responses to *unseen* genetic perturbations?
**Week 1: reproduce the linear baseline of Ahlmann-Eltze, Huber & Anders (Nat Methods 2025) on Replogle 2022 K562 essential CRISPRi. No Hi-C yet.**

## Run
```bash
mamba env create -f environment.yml && mamba activate hic-perturb-prior
python -m src.run_baseline --config configs/k562_repro_5000hvg.yaml --preprocess   # primary
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

## Method (src/)
- `data.py`: streaming `normalize_total(1e4)` -> `log1p` -> top-N HVG -> pseudobulk mean per perturbation -> `delta = pseudobulk(pert) - pseudobulk(control)`.
- `baselines.py`: no-change, mean, `Y_hat = G W P~' + b` (ridge, closed form as in the paper's Methods / `solve_y_axb`). Embeddings are plug-ins in `EMBEDDINGS`
  `(Y_train, train_idx, ds, K) -> (G, P_all)`; add GO / Hi-C / shuffled / random P there, model and evaluation stay unchanged.
- `evaluate.py`: L2 and Pearson (pred vs observed delta) per perturbation; subsets `all`, `top20_de` (oracle: largest |observed delta|), `top1000_expr` (paper's metric:
  most expressed genes in control). Each perturbation's own target gene is excluded by default (`eval.exclude_target_primary`); "included" metrics are always reported too.
  Perturbations are split (never cells) 70/10/20, 5 seeds; K, lambda tuned on val (grid in config); final models refit on train+val.

## Leakage audit
- **G (PCA) and W** are fit on training perturbations only (`select_hparams`: train; final fit: train+val). Test deltas never enter PCA, W or b.
- **HVG selection** (and cell normalisation) use *all* cells, including cells of test perturbations, but no perturbation labels (unsupervised, same as the paper/GEARS
  which preprocess the full dataset first). A 2nd-order effect: the HVG set decides which perturbations are evaluable (target must be a read-out gene).
- `top1000_expr` uses control cells only; `top20_de` uses the observed test delta by design (oracle subset, evaluation only, never fed to a model).
- Hyper-parameters are chosen on val perturbations only. The "paper" column uses fixed K=10, lambda=0.1.

## Known deviations from the paper
- Paper data = GEARS-processed Replogle K562 (5000 genes, filtered perturbations; exact GEARS preprocessing not reproduced). Ours: scPerturb raw counts, own scanpy-seurat HVG replica
  (5000 genes), **no perturbation QC filtering** (QC only logged). Absolute L2/Pearson are therefore not expected to match.
- Paper splits: GEARS `prepare_split('simulation', seed=1|2)` on the GEARS file (2 splits). Not stored anywhere, regenerable only with the GEARS file + package (not done yet).
- Paper's 'val' perturbations are unused by its linear model; ours use them to tune K/lambda (secondary column only).
