# Comparison with Ahlmann-Eltze et al. 2025 (Replogle K562)

Paper rows: top-1000-expressed-genes L2 and Pearson delta, mean ± sd over their 2 splits (GEARS-processed data).
Our rows: same metric definitions, mean ± sd over our splits. Every table states its eval_set and size.

## eval_set = common, target gene included (like-for-like with the paper)

104 test perturbations per split, 2 splits.

|                                                       | L2 (top1000 expr)   | Pearson delta (top1000 expr)   |
|:------------------------------------------------------|:--------------------|:-------------------------------|
| paper: mean                                           | 4.119 ± 0.026       | 0.404 ± 0.009                  |
| paper: linear model                                   | 3.996 ± 0.091       | 0.433 ± 0.017                  |
| paper: GEARS                                          | 3.991 ± 0.029       | 0.387 ± 0.019                  |
| ours: no_change                                       | 4.650 ± 0.057       | n/a                            |
| ours: mean                                            | 4.267 ± 0.054       | 0.404 ± 0.039                  |
| ours: linear[train_pca_expr] (paper hp K=10, lam=0.1) | 4.153 ± 0.121       | 0.435 ± 0.048                  |
| ours: linear[train_pca_expr] (val-tuned)              | 4.125 ± 0.151       | 0.446 ± 0.053                  |

Ratio to the mean baseline (linear / mean; L2 < 1 and Pearson > 1 mean the linear model wins):

|                                                       |   L2 ratio |   Pearson ratio |
|:------------------------------------------------------|-----------:|----------------:|
| ours: linear[train_pca_expr] (paper hp K=10, lam=0.1) |      0.973 |           1.075 |
| ours: linear[train_pca_expr] (val-tuned)              |      0.967 |           1.103 |
| paper: linear model                                   |      0.97  |           1.072 |

## eval_set = common, target gene excluded

104 test perturbations per split, 2 splits.

|                                                       | L2 (top1000 expr)   | Pearson delta (top1000 expr)   |
|:------------------------------------------------------|:--------------------|:-------------------------------|
| paper: mean                                           | 4.119 ± 0.026       | 0.404 ± 0.009                  |
| paper: linear model                                   | 3.996 ± 0.091       | 0.433 ± 0.017                  |
| paper: GEARS                                          | 3.991 ± 0.029       | 0.387 ± 0.019                  |
| ours: no_change                                       | 4.555 ± 0.026       | n/a                            |
| ours: mean                                            | 4.172 ± 0.089       | 0.412 ± 0.042                  |
| ours: linear[train_pca_expr] (paper hp K=10, lam=0.1) | 4.055 ± 0.158       | 0.443 ± 0.052                  |
| ours: linear[train_pca_expr] (val-tuned)              | 4.026 ± 0.188       | 0.454 ± 0.058                  |

Ratio to the mean baseline (linear / mean; L2 < 1 and Pearson > 1 mean the linear model wins):

|                                                       |   L2 ratio |   Pearson ratio |
|:------------------------------------------------------|-----------:|----------------:|
| ours: linear[train_pca_expr] (paper hp K=10, lam=0.1) |      0.972 |           1.075 |
| ours: linear[train_pca_expr] (val-tuned)              |      0.965 |           1.104 |
| paper: linear model                                   |      0.97  |           1.072 |
