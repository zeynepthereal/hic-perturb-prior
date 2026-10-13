# Comparison with Ahlmann-Eltze et al. 2025 (Replogle K562)

Paper rows: top-1000-expressed-genes L2 and Pearson delta, mean ± sd over 2 splits (GEARS-processed data, 5000 genes).
Our rows: same metric definitions on our preprocessing, mean ± sd over our splits.

## target gene included (like-for-like with the paper)

|                                    | L2 (top1000 expr)   | Pearson delta (top1000 expr)   |
|:-----------------------------------|:--------------------|:-------------------------------|
| paper: mean                        | 4.119 ± 0.026       | 0.404 ± 0.009                  |
| paper: linear model                | 3.996 ± 0.091       | 0.433 ± 0.017                  |
| paper: GEARS                       | 3.991 ± 0.029       | 0.387 ± 0.019                  |
| ours: no_change                    | 3.087 ± 0.165       | n/a                            |
| ours: mean                         | 3.008 ± 0.150       | 0.273 ± 0.006                  |
| ours: linear (paper K=10, lam=0.1) | 2.939 ± 0.109       | 0.294 ± 0.012                  |
| ours: linear (val-tuned)           | 2.929 ± 0.119       | 0.300 ± 0.010                  |

Ratio to the mean baseline (linear / mean; L2 < 1 and Pearson > 1 mean the linear model wins):

|                                    |   L2 ratio |   Pearson ratio |
|:-----------------------------------|-----------:|----------------:|
| ours: linear (paper K=10, lam=0.1) |      0.977 |           1.079 |
| ours: linear (val-tuned)           |      0.974 |           1.101 |
| paper: linear model                |      0.97  |           1.072 |

## target gene excluded

|                                    | L2 (top1000 expr)   | Pearson delta (top1000 expr)   |
|:-----------------------------------|:--------------------|:-------------------------------|
| paper: mean                        | 4.119 ± 0.026       | 0.404 ± 0.009                  |
| paper: linear model                | 3.996 ± 0.091       | 0.433 ± 0.017                  |
| paper: GEARS                       | 3.991 ± 0.029       | 0.387 ± 0.019                  |
| ours: no_change                    | 3.012 ± 0.160       | n/a                            |
| ours: mean                         | 2.937 ± 0.144       | 0.277 ± 0.007                  |
| ours: linear (paper K=10, lam=0.1) | 2.864 ± 0.101       | 0.299 ± 0.013                  |
| ours: linear (val-tuned)           | 2.855 ± 0.110       | 0.305 ± 0.011                  |

Ratio to the mean baseline (linear / mean; L2 < 1 and Pearson > 1 mean the linear model wins):

|                                    |   L2 ratio |   Pearson ratio |
|:-----------------------------------|-----------:|----------------:|
| ours: linear (paper K=10, lam=0.1) |      0.975 |           1.08  |
| ours: linear (val-tuned)           |      0.972 |           1.101 |
| paper: linear model                |      0.97  |           1.072 |
