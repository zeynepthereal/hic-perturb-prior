# S4 / S9 — arms vs the MEAN BASELINE (pre-registered secondary; numbers only)

Primary test cases (common set, GEARS splits 1-5, L2 on the top-1000 expressed genes, target excluded). D = (1/5) Σ_s mean[L2_mean − L2_arm]; **positive = the arm is better than the mean baseline**.
Primary run `results/20261014-134117_prereg_primary`; PCA reference run `results/20261014-133938_sanity1_pca_only`. Named arms: frozen sign-flip test (100,000 permutations, 10,000 bootstrap draws). Random-Gaussian arms: 10,000 permutations / 1,000 bootstrap draws per arm.

## K = 10

| arm | mean L2 | mean L2 (mean baseline) | D vs mean baseline | p (2-sided) | 95 % CI |
|---|---|---|---|---|---|
| b500_hic_trans | 4.24097 | 4.22982 | -0.01115 | 0.072659 | [-0.02396, 0.00041] |
| b500_lin_randtrans (10-seed mean of case-level L2) | 4.28209 | 4.22982 | -0.05227 | 0.000010 | [-0.09409, -0.02091] |
| random_gaussian (1,000 arms; summary): median [5 %..95 %] | 4.24945 [4.24204..4.25606] | 4.22982 | -0.01963 [-0.02624..-0.01223]; range -0.03113..-0.00497 | median 0.0001; D>0 in 0/1000; p<0.05 & arm worse in 993/1000; p<0.05 & arm better in 0/1000 | median CI [-0.02845, -0.01080] |

## K = 50

| arm | mean L2 | mean L2 (mean baseline) | D vs mean baseline | p (2-sided) | 95 % CI |
|---|---|---|---|---|---|
| b500_hic_trans | 4.45330 | 4.22982 | -0.22348 | 0.000010 | [-0.33450, -0.14048] |
| b500_lin_randtrans (10-seed mean of case-level L2) | 4.52083 | 4.22982 | -0.29101 | 0.000010 | [-0.40467, -0.19733] |
| random_gaussian (1,000 arms; summary): median [5 %..95 %] | 4.36423 [4.34818..4.37961] | 4.22982 | -0.13441 [-0.14979..-0.11836]; range -0.16545..-0.10556 | median 0.0001; D>0 in 0/1000; p<0.05 & arm worse in 1000/1000; p<0.05 & arm better in 0/1000 | median CI [-0.15560, -0.11389] |

## S4/S9 reference — PCA embedding vs mean baseline on its own predictable subset (517 cases, 319 unique perturbations)

| K | mean L2 (PCA) | mean L2 (mean baseline) | D | p (2-sided) | 95 % CI |
|---|---|---|---|---|---|
| 10 | 4.15013 | 4.25574 | 0.10561 | 0.005120 | [0.03553, 0.18347] |
| 50 | 4.40093 | 4.25574 | -0.14519 | 0.007010 | [-0.24732, -0.03970] |
