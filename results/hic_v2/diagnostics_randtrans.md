# Matched random-trans controls (500 kb)

seeds found: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

trans-degree sequence identical to hic_trans in all 20 random graphs: True

## Structure (mean, sd over seeds)

|                            |   b500_hic_randtrans |   b500_lin_randtrans |
|:---------------------------|---------------------:|---------------------:|
| ('nodes_embedded', 'mean') |            5861      |            5861      |
| ('nodes_embedded', 'std')  |               0      |               0      |
| ('edges', 'mean')          |           47832      |           52261      |
| ('edges', 'std')           |               0      |               0      |
| ('components', 'mean')     |               1      |               1      |
| ('components', 'std')      |               0      |               0      |
| ('deg_median', 'mean')     |              14      |              16      |
| ('deg_median', 'std')      |               0      |               0      |
| ('deg_mean', 'mean')       |              16.3221 |              17.8335 |
| ('deg_mean', 'std')        |               0      |               0      |
| ('deg_p99', 'mean')        |              41      |              31      |
| ('deg_p99', 'std')         |               0      |               0      |
| ('deg_max', 'mean')        |             242      |              37      |
| ('deg_max', 'std')         |               0      |               0      |
| ('eigval_1', 'mean')       |               0.1741 |               0.0289 |
| ('eigval_1', 'std')        |               0.0039 |               0.0004 |
| ('eigval_10', 'mean')      |               0.2519 |               0.1423 |
| ('eigval_10', 'std')       |               0.0106 |               0.0036 |
| ('eigval_50', 'mean')      |               0.3293 |               0.2179 |
| ('eigval_50', 'std')       |               0.0041 |               0.0013 |

## R² on chromosome / position (mean, sd over seeds)

|                            |   ('R2_chromosome', 'mean') |   ('R2_chromosome', 'std') |   ('R2_chrom_cubic', 'mean') |   ('R2_chrom_cubic', 'std') |   ('R2_chrom_spline', 'mean') |   ('R2_chrom_spline', 'std') |   ('frac_nodes_nonzero_in_a_dim', 'mean') |   ('frac_nodes_nonzero_in_a_dim', 'std') |
|:---------------------------|----------------------------:|---------------------------:|-----------------------------:|----------------------------:|------------------------------:|-----------------------------:|------------------------------------------:|-----------------------------------------:|
| ('b500_hic_randtrans', 10) |                       0.305 |                      0.065 |                        0.36  |                       0.076 |                         0.414 |                        0.08  |                                         1 |                                        0 |
| ('b500_hic_randtrans', 50) |                       0.202 |                      0.006 |                        0.293 |                       0.01  |                         0.352 |                        0.011 |                                         1 |                                        0 |
| ('b500_lin_randtrans', 10) |                       0.114 |                      0.004 |                        0.298 |                       0.011 |                         0.548 |                        0.018 |                                         1 |                                        0 |
| ('b500_lin_randtrans', 50) |                       0.085 |                      0.002 |                        0.277 |                       0.006 |                         0.512 |                        0.008 |                                         1 |                                        0 |

Reference cells (same metrics):

| variant        |   dims |   R2_chromosome |   R2_chrom_cubic |   R2_chrom_spline |   frac_nodes_nonzero_in_a_dim |   mean_frac_nodes_loaded_per_dim |
|:---------------|-------:|----------------:|-----------------:|------------------:|------------------------------:|---------------------------------:|
| b500_hic_trans |     10 |           0.42  |            0.605 |             0.691 |                         1     |                            0.887 |
| b500_hic_trans |     50 |           0.248 |            0.441 |             0.55  |                         1     |                            0.905 |
| b500_hic_cis   |     10 |           0.004 |            0.399 |             0.488 |                         0.462 |                            0.054 |
| b500_hic_cis   |     50 |           0.005 |            0.283 |             0.417 |                         0.981 |                            0.049 |
| b500_lin_trans |     10 |           0.236 |            0.551 |             0.74  |                         1     |                            0.901 |
| b500_lin_trans |     50 |           0.137 |            0.398 |             0.653 |                         1     |                            0.885 |
| b500_lin_cis   |     10 |           0     |            0.991 |             0.997 |                         0.616 |                            0.061 |
| b500_lin_cis   |     50 |           0     |            0.853 |             0.994 |                         0.973 |                            0.054 |

## Correlation with A/B track (mean, sd over seeds)

|                                            |   ('multiple_R', 'mean') |   ('multiple_R', 'std') |   ('R2_beyond_position', 'mean') |   ('R2_beyond_position', 'std') |
|:-------------------------------------------|-------------------------:|------------------------:|---------------------------------:|--------------------------------:|
| ('b500_hic_randtrans', 'AB_PC1', 10)       |                    0.087 |                   0.008 |                            0.017 |                           0.006 |
| ('b500_hic_randtrans', 'AB_PC1', 50)       |                    0.232 |                   0.016 |                            0.08  |                           0.006 |
| ('b500_hic_randtrans', 'AB_best_of_3', 10) |                    0.104 |                   0.018 |                            0.017 |                           0.007 |
| ('b500_hic_randtrans', 'AB_best_of_3', 50) |                    0.307 |                   0.022 |                            0.096 |                           0.01  |
| ('b500_lin_randtrans', 'AB_PC1', 10)       |                    0.296 |                   0.009 |                            0.054 |                           0.001 |
| ('b500_lin_randtrans', 'AB_PC1', 50)       |                    0.367 |                   0.009 |                            0.118 |                           0.007 |
| ('b500_lin_randtrans', 'AB_best_of_3', 10) |                    0.155 |                   0.013 |                            0.018 |                           0.002 |
| ('b500_lin_randtrans', 'AB_best_of_3', 50) |                    0.3   |                   0.006 |                            0.066 |                           0.005 |

Reference cells:

| variant        | ab           |   dims |   multiple_R |   R2_beyond_position |   R2_AB_from_position_only |
|:---------------|:-------------|-------:|-------------:|---------------------:|---------------------------:|
| b500_hic_trans | AB_PC1       |     10 |        0.118 |                0.078 |                      0.512 |
| b500_hic_trans | AB_PC1       |     50 |        0.378 |                0.213 |                      0.512 |
| b500_hic_trans | AB_best_of_3 |     10 |        0.283 |                0.163 |                      0.319 |
| b500_hic_trans | AB_best_of_3 |     50 |        0.75  |                0.634 |                      0.319 |
| b500_hic_cis   | AB_PC1       |     10 |        0.07  |                0.036 |                      0.512 |
| b500_hic_cis   | AB_PC1       |     50 |        0.389 |                0.172 |                      0.512 |
| b500_hic_cis   | AB_best_of_3 |     10 |        0.199 |                0.032 |                      0.319 |
| b500_hic_cis   | AB_best_of_3 |     50 |        0.561 |                0.298 |                      0.319 |
| b500_lin_trans | AB_PC1       |     10 |        0.305 |                0.074 |                      0.512 |
| b500_lin_trans | AB_PC1       |     50 |        0.495 |                0.19  |                      0.512 |
| b500_lin_trans | AB_best_of_3 |     10 |        0.32  |                0.093 |                      0.319 |
| b500_lin_trans | AB_best_of_3 |     50 |        0.502 |                0.286 |                      0.319 |
| b500_lin_cis   | AB_PC1       |     10 |        0.278 |                0.004 |                      0.512 |
| b500_lin_cis   | AB_PC1       |     50 |        0.548 |                0.066 |                      0.512 |
| b500_lin_cis   | AB_best_of_3 |     10 |        0.259 |                0.004 |                      0.319 |
| b500_lin_cis   | AB_best_of_3 |     50 |        0.431 |                0.03  |                      0.319 |
