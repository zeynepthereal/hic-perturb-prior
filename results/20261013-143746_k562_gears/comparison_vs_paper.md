# Comparison with Ahlmann-Eltze et al. 2025 (Replogle K562)

Paper rows: top-1000-expressed-genes L2 and Pearson delta, mean ± sd over their 2 splits (GEARS-processed data).
Our rows: same metric definitions, mean ± sd over our splits. Every table states its eval_set and size.

## eval_set = common, target gene included (like-for-like with the paper)

103 test perturbations per split, 5 splits.

|                                                       | L2 (top1000 expr)   | Pearson delta (top1000 expr)   |
|:------------------------------------------------------|:--------------------|:-------------------------------|
| paper: mean                                           | 4.119 ± 0.026       | 0.404 ± 0.009                  |
| paper: linear model                                   | 3.996 ± 0.091       | 0.433 ± 0.017                  |
| paper: GEARS                                          | 3.991 ± 0.029       | 0.387 ± 0.019                  |
| ours: no_change                                       | 4.762 ± 0.217       | n/a                            |
| ours: mean                                            | 4.350 ± 0.171       | 0.411 ± 0.024                  |
| ours: linear[train_pca_expr] (paper hp K=10, lam=0.1) | 4.251 ± 0.163       | 0.438 ± 0.025                  |
| ours: linear[train_pca_expr] (val-tuned)              | 4.223 ± 0.168       | 0.442 ± 0.033                  |

Ratio to the mean baseline (linear / mean; L2 < 1 and Pearson > 1 mean the linear model wins):

|                                                       |   L2 ratio |   Pearson ratio |
|:------------------------------------------------------|-----------:|----------------:|
| ours: linear[train_pca_expr] (paper hp K=10, lam=0.1) |      0.977 |           1.064 |
| ours: linear[train_pca_expr] (val-tuned)              |      0.971 |           1.074 |
| paper: linear model                                   |      0.97  |           1.072 |

## eval_set = common, target gene excluded

103 test perturbations per split, 5 splits.

|                                                       | L2 (top1000 expr)   | Pearson delta (top1000 expr)   |
|:------------------------------------------------------|:--------------------|:-------------------------------|
| paper: mean                                           | 4.119 ± 0.026       | 0.404 ± 0.009                  |
| paper: linear model                                   | 3.996 ± 0.091       | 0.433 ± 0.017                  |
| paper: GEARS                                          | 3.991 ± 0.029       | 0.387 ± 0.019                  |
| ours: no_change                                       | 4.668 ± 0.215       | n/a                            |
| ours: mean                                            | 4.256 ± 0.174       | 0.418 ± 0.026                  |
| ours: linear[train_pca_expr] (paper hp K=10, lam=0.1) | 4.150 ± 0.167       | 0.445 ± 0.027                  |
| ours: linear[train_pca_expr] (val-tuned)              | 4.123 ± 0.172       | 0.449 ± 0.035                  |

Ratio to the mean baseline (linear / mean; L2 < 1 and Pearson > 1 mean the linear model wins):

|                                                       |   L2 ratio |   Pearson ratio |
|:------------------------------------------------------|-----------:|----------------:|
| ours: linear[train_pca_expr] (paper hp K=10, lam=0.1) |      0.975 |           1.065 |
| ours: linear[train_pca_expr] (val-tuned)              |      0.969 |           1.074 |
| paper: linear model                                   |      0.97  |           1.072 |
