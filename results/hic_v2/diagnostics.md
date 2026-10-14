# Hi-C graph diagnostics v2 (bin-level, whole genome)

## Structure (components, degrees, spectrum)

| variant        |   nodes_embedded |   isolated |   edges |   components |   largest_component |   deg_median |   deg_mean |   deg_p99 |   deg_max |   first_eigval |   eigval_10 |   eigval_50 |
|:---------------|-----------------:|-----------:|--------:|-------------:|--------------------:|-------------:|-----------:|----------:|----------:|---------------:|------------:|------------:|
| b500_hic_trans |             5861 |        212 |   47832 |            1 |                5861 |           14 |       16.3 |     41    |       242 |        0.06016 |     0.13985 |     0.3022  |
| b500_hic_cis   |             5861 |        212 |   25384 |           23 |                 493 |            7 |        8.7 |     30.4  |       240 |        0.03773 |     0.0867  |     0.22035 |
| b500_lin_trans |             5861 |        212 |   52261 |            1 |                5861 |           16 |       17.8 |     31    |        37 |        0.02665 |     0.06366 |     0.13089 |
| b500_lin_cis   |             5861 |        212 |   29813 |           23 |                 493 |           10 |       10.2 |     14    |        17 |        0.00017 |     0.00057 |     0.00512 |
| b250_hic_trans |            11775 |        362 |   80129 |            1 |               11775 |           11 |       13.6 |     43    |       408 |        0.05916 |     0.14902 |     0.29888 |
| b250_hic_cis   |            11775 |        362 |   38424 |           23 |                 989 |            4 |        6.5 |     29.26 |       405 |        0.01421 |     0.05917 |     0.1772  |
| b250_lin_trans |            11775 |        362 |  101128 |            1 |               11775 |           16 |       17.2 |     30    |        37 |        0.00467 |     0.05107 |     0.12244 |
| b250_lin_cis   |            11775 |        362 |   59423 |           23 |                 989 |           10 |       10.1 |     14    |        17 |        5e-05   |     0.00015 |     0.00141 |
| hic_local25k   |            20165 |         98 |  135400 |           39 |                1490 |           12 |       13.4 |     35    |        92 |        0.00016 |     0.00046 |     0.00538 |

## R² of embedding dimensions explained by chromosome / position (mean over dims)

`R2_chrom_spline` = chromosome-specific cubic B-spline of position (8 df each); `frac_nodes_nonzero_in_a_dim` and `mean_frac_nodes_loaded_per_dim` measure how localized the eigenvectors are (cis-only graphs have one localized slow mode per chromosome).

| variant                 |   dims |   R2_chromosome |   R2_chrom_cubic |   R2_chrom_spline |   frac_nodes_nonzero_in_a_dim |   mean_frac_nodes_loaded_per_dim |
|:------------------------|-------:|----------------:|-----------------:|------------------:|------------------------------:|---------------------------------:|
| b500_hic_trans          |     10 |           0.42  |            0.605 |             0.691 |                         1     |                            0.887 |
| b500_hic_trans          |     50 |           0.248 |            0.441 |             0.55  |                         1     |                            0.905 |
| b500_hic_cis            |     10 |           0.004 |            0.399 |             0.488 |                         0.462 |                            0.054 |
| b500_hic_cis            |     50 |           0.005 |            0.283 |             0.417 |                         0.981 |                            0.049 |
| b500_lin_trans          |     10 |           0.236 |            0.551 |             0.74  |                         1     |                            0.901 |
| b500_lin_trans          |     50 |           0.137 |            0.398 |             0.653 |                         1     |                            0.885 |
| b500_lin_cis            |     10 |           0     |            0.991 |             0.997 |                         0.616 |                            0.061 |
| b500_lin_cis            |     50 |           0     |            0.853 |             0.994 |                         0.973 |                            0.054 |
| b250_hic_trans          |     10 |           0.261 |            0.436 |             0.536 |                         1     |                            0.862 |
| b250_hic_trans          |     50 |           0.178 |            0.327 |             0.408 |                         1     |                            0.882 |
| b250_hic_cis            |     10 |           0.003 |            0.374 |             0.443 |                         0.45  |                            0.045 |
| b250_hic_cis            |     50 |           0.005 |            0.223 |             0.328 |                         1     |                            0.049 |
| b250_lin_trans          |     10 |           0.11  |            0.336 |             0.489 |                         1     |                            0.635 |
| b250_lin_trans          |     50 |           0.086 |            0.267 |             0.448 |                         1     |                            0.769 |
| b250_lin_cis            |     10 |           0     |            0.997 |             0.999 |                         0.616 |                            0.061 |
| b250_lin_cis            |     50 |           0     |            0.832 |             0.996 |                         0.974 |                            0.056 |
| hic_local25k            |     10 |           0.005 |            0.814 |             0.914 |                         0.388 |                            0.043 |
| hic_local25k            |     50 |           0.004 |            0.63  |             0.855 |                         0.892 |                            0.036 |
| b500_hic_trans:resid    |     10 |           0     |            0     |             0     |                         1     |                            0.835 |
| b500_hic_trans:resid    |     50 |           0     |            0     |             0     |                         1     |                            0.885 |
| b500_hic_cis:resid      |     10 |           0     |           -0     |             0     |                         0.462 |                            0.053 |
| b500_hic_cis:resid      |     50 |           0     |            0     |             0     |                         0.981 |                            0.049 |
| b250_hic_trans:resid    |     10 |          -0     |           -0     |            -0     |                         1     |                            0.775 |
| b250_hic_trans:resid    |     50 |          -0     |           -0     |            -0     |                         1     |                            0.855 |
| b250_hic_cis:resid      |     10 |          -0     |           -0     |            -0     |                         0.45  |                            0.044 |
| b250_hic_cis:resid      |     50 |           0     |           -0     |            -0     |                         1     |                            0.049 |
| b500_hic_trans:random:0 |     10 |           0.004 |            0.017 |             0.032 |                         1     |                            0.961 |
| b500_hic_trans:random:0 |     50 |           0.004 |            0.016 |             0.032 |                         1     |                            0.96  |

## Correlation with A/B compartment track (bin level)

AB_PC1 = PC1 of the per-chromosome O/E correlation matrix oriented by gene density (the standard definition, as specified); AB_best_of_3 = whichever of PC1-3 correlates best with gene density (PC1 does not track gene density in this rearranged K562 map; see ab_*k_chrom.csv). `multiple_R` = sqrt(R²) of AB regressed on the first `dims` embedding dimensions; `R2_beyond_position` repeats the regression after removing chromosome+position from both sides; `R2_AB_from_position_only` is the reference.

| variant                 | ab           |   dims |   max_abs_r_single_dim |   best_dim |   multiple_R |   R2_beyond_position |   R2_AB_from_position_only |   n_nodes |
|:------------------------|:-------------|-------:|-----------------------:|-----------:|-------------:|---------------------:|---------------------------:|----------:|
| b500_hic_trans          | AB_PC1       |     10 |                  0.078 |          3 |        0.118 |                0.078 |                      0.512 |      5861 |
| b500_hic_trans          | AB_PC1       |     50 |                  0.117 |         12 |        0.378 |                0.213 |                      0.512 |      5861 |
| b500_hic_trans          | AB_best_of_3 |     10 |                  0.181 |         10 |        0.283 |                0.163 |                      0.319 |      5861 |
| b500_hic_trans          | AB_best_of_3 |     50 |                  0.309 |         17 |        0.75  |                0.634 |                      0.319 |      5861 |
| b500_hic_cis            | AB_PC1       |     10 |                  0.043 |          7 |        0.07  |                0.036 |                      0.512 |      5861 |
| b500_hic_cis            | AB_PC1       |     50 |                  0.183 |         11 |        0.389 |                0.172 |                      0.512 |      5861 |
| b500_hic_cis            | AB_best_of_3 |     10 |                  0.146 |          9 |        0.199 |                0.032 |                      0.319 |      5861 |
| b500_hic_cis            | AB_best_of_3 |     50 |                  0.194 |         41 |        0.561 |                0.298 |                      0.319 |      5861 |
| b500_lin_trans          | AB_PC1       |     10 |                  0.218 |          2 |        0.305 |                0.074 |                      0.512 |      5861 |
| b500_lin_trans          | AB_PC1       |     50 |                  0.218 |          2 |        0.495 |                0.19  |                      0.512 |      5861 |
| b500_lin_trans          | AB_best_of_3 |     10 |                  0.224 |         10 |        0.32  |                0.093 |                      0.319 |      5861 |
| b500_lin_trans          | AB_best_of_3 |     50 |                  0.224 |         10 |        0.502 |                0.286 |                      0.319 |      5861 |
| b500_lin_cis            | AB_PC1       |     10 |                  0.153 |          4 |        0.278 |                0.004 |                      0.512 |      5861 |
| b500_lin_cis            | AB_PC1       |     50 |                  0.171 |         14 |        0.548 |                0.066 |                      0.512 |      5861 |
| b500_lin_cis            | AB_best_of_3 |     10 |                  0.13  |          2 |        0.259 |                0.004 |                      0.319 |      5861 |
| b500_lin_cis            | AB_best_of_3 |     50 |                  0.13  |          2 |        0.431 |                0.03  |                      0.319 |      5861 |
| b250_hic_trans          | AB_PC1       |     10 |                  0.096 |          4 |        0.152 |                0.07  |                      0.534 |     11775 |
| b250_hic_trans          | AB_PC1       |     50 |                  0.113 |         29 |        0.363 |                0.147 |                      0.534 |     11775 |
| b250_hic_trans          | AB_best_of_3 |     10 |                  0.207 |          8 |        0.299 |                0.151 |                      0.334 |     11775 |
| b250_hic_trans          | AB_best_of_3 |     50 |                  0.229 |         24 |        0.646 |                0.455 |                      0.334 |     11775 |
| b250_hic_cis            | AB_PC1       |     10 |                  0.071 |         10 |        0.113 |                0.024 |                      0.534 |     11775 |
| b250_hic_cis            | AB_PC1       |     50 |                  0.118 |         15 |        0.286 |                0.095 |                      0.534 |     11775 |
| b250_hic_cis            | AB_best_of_3 |     10 |                  0.151 |          4 |        0.26  |                0.037 |                      0.334 |     11775 |
| b250_hic_cis            | AB_best_of_3 |     50 |                  0.16  |         22 |        0.475 |                0.185 |                      0.334 |     11775 |
| b250_lin_trans          | AB_PC1       |     10 |                  0.209 |          2 |        0.328 |                0.063 |                      0.534 |     11775 |
| b250_lin_trans          | AB_PC1       |     50 |                  0.209 |          2 |        0.522 |                0.171 |                      0.534 |     11775 |
| b250_lin_trans          | AB_best_of_3 |     10 |                  0.197 |         10 |        0.268 |                0.046 |                      0.334 |     11775 |
| b250_lin_trans          | AB_best_of_3 |     50 |                  0.258 |         13 |        0.579 |                0.376 |                      0.334 |     11775 |
| b250_lin_cis            | AB_PC1       |     10 |                  0.137 |          4 |        0.253 |                0.001 |                      0.534 |     11775 |
| b250_lin_cis            | AB_PC1       |     50 |                  0.172 |         16 |        0.586 |                0.071 |                      0.534 |     11775 |
| b250_lin_cis            | AB_best_of_3 |     10 |                  0.159 |          4 |        0.32  |                0.002 |                      0.334 |     11775 |
| b250_lin_cis            | AB_best_of_3 |     50 |                  0.159 |          4 |        0.449 |                0.02  |                      0.334 |     11775 |
| hic_local25k            | AB_PC1       |     10 |                  0.156 |          5 |        0.322 |                0.036 |                      0.61  |     20165 |
| hic_local25k            | AB_PC1       |     50 |                  0.171 |         25 |        0.523 |                0.086 |                      0.61  |     20165 |
| hic_local25k            | AB_best_of_3 |     10 |                  0.191 |          7 |        0.32  |                0.005 |                      0.392 |     20165 |
| hic_local25k            | AB_best_of_3 |     50 |                  0.191 |          7 |        0.396 |                0.039 |                      0.392 |     20165 |
| b500_hic_trans:resid    | AB_PC1       |     10 |                  0.098 |          4 |        0.195 |                0.078 |                      0.512 |      5861 |
| b500_hic_trans:resid    | AB_PC1       |     50 |                  0.098 |          4 |        0.322 |                0.213 |                      0.512 |      5861 |
| b500_hic_trans:resid    | AB_best_of_3 |     10 |                  0.215 |         10 |        0.333 |                0.163 |                      0.319 |      5861 |
| b500_hic_trans:resid    | AB_best_of_3 |     50 |                  0.326 |         17 |        0.657 |                0.634 |                      0.319 |      5861 |
| b500_hic_cis:resid      | AB_PC1       |     10 |                  0.071 |          5 |        0.132 |                0.036 |                      0.512 |      5861 |
| b500_hic_cis:resid      | AB_PC1       |     50 |                  0.124 |         50 |        0.289 |                0.172 |                      0.512 |      5861 |
| b500_hic_cis:resid      | AB_best_of_3 |     10 |                  0.109 |          9 |        0.148 |                0.032 |                      0.319 |      5861 |
| b500_hic_cis:resid      | AB_best_of_3 |     50 |                  0.162 |         41 |        0.45  |                0.298 |                      0.319 |      5861 |
| b250_hic_trans:resid    | AB_PC1       |     10 |                  0.086 |          4 |        0.18  |                0.07  |                      0.534 |     11775 |
| b250_hic_trans:resid    | AB_PC1       |     50 |                  0.086 |          4 |        0.262 |                0.147 |                      0.534 |     11775 |
| b250_hic_trans:resid    | AB_best_of_3 |     10 |                  0.186 |          8 |        0.317 |                0.151 |                      0.334 |     11775 |
| b250_hic_trans:resid    | AB_best_of_3 |     50 |                  0.186 |          8 |        0.551 |                0.455 |                      0.334 |     11775 |
| b250_hic_cis:resid      | AB_PC1       |     10 |                  0.054 |          5 |        0.105 |                0.024 |                      0.534 |     11775 |
| b250_hic_cis:resid      | AB_PC1       |     50 |                  0.062 |         28 |        0.211 |                0.095 |                      0.534 |     11775 |
| b250_hic_cis:resid      | AB_best_of_3 |     10 |                  0.095 |          4 |        0.158 |                0.037 |                      0.334 |     11775 |
| b250_hic_cis:resid      | AB_best_of_3 |     50 |                  0.119 |         34 |        0.351 |                0.185 |                      0.334 |     11775 |
| b500_hic_trans:random:0 | AB_PC1       |     10 |                  0.017 |          6 |        0.034 |                0.001 |                      0.512 |      5861 |
| b500_hic_trans:random:0 | AB_PC1       |     50 |                  0.033 |         47 |        0.091 |                0.006 |                      0.512 |      5861 |
| b500_hic_trans:random:0 | AB_best_of_3 |     10 |                  0.019 |          2 |        0.037 |                0.002 |                      0.319 |      5861 |
| b500_hic_trans:random:0 | AB_best_of_3 |     50 |                  0.036 |         32 |        0.084 |                0.006 |                      0.319 |      5861 |

## Predictable GEARS perturbations

|                | 0                      | 1                  | 2                  | 3                  | 4                  | 5                  |
|:---------------|:-----------------------|:-------------------|:-------------------|:-------------------|:-------------------|:-------------------|
| scope          | all 1092 perturbations | test seed1 (n=273) | test seed2 (n=273) | test seed3 (n=273) | test seed4 (n=273) | test seed5 (n=273) |
| train_pca_expr | 412                    | 104                | 105                | 105                | 100                | 103                |
| b500_hic_trans | 1092                   | 273                | 273                | 273                | 273                | 273                |
| b500_hic_cis   | 1092                   | 273                | 273                | 273                | 273                | 273                |
| b500_lin_trans | 1092                   | 273                | 273                | 273                | 273                | 273                |
| b500_lin_cis   | 1092                   | 273                | 273                | 273                | 273                | 273                |
| b250_hic_trans | 1092                   | 273                | 273                | 273                | 273                | 273                |
| b250_hic_cis   | 1092                   | 273                | 273                | 273                | 273                | 273                |
| b250_lin_trans | 1092                   | 273                | 273                | 273                | 273                | 273                |
| b250_lin_cis   | 1092                   | 273                | 273                | 273                | 273                | 273                |
| hic_local25k   | 1092                   | 273                | 273                | 273                | 273                | 273                |

`common` set of PCA ∩ all Hi-C variants: 412 perturbations; `common` among bin-level cells only: 1092.
