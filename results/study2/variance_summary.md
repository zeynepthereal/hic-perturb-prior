# Study 2: variance summary (finished runs: 9 of 9; missing: none)

Top-1000 control-expressed genes; cells are `median [min, max]` over runs (GEARS: up to 9 runs; deterministic baselines: one value per split, up to 3). Pearson_centred = Pearson after subtracting the mean perturbation effect across the evaluated test perturbations (NaN for constant predictors). Amplitude ratios: mean-delta = ||mean predicted delta|| / ||mean observed delta||; per-pert norm = mean ||predicted delta|| / mean ||observed delta||; collapse ratio = collapse(pred)/collapse(obs).


## A: GEARS & PCA predictable - target excluded

| model | n runs | L2 | Pearson | Pearson (centred) | amp. ratio, mean delta | amp. ratio, per-pert norm | collapse ratio pred/obs |
|---|---|---|---|---|---|---|---|
| GEARS | 9 | 4.1365 [3.8722, 4.7229] | 0.4115 [0.3530, 0.4316] | 0.2371 [0.1523, 0.3304] | 0.8615 [0.7251, 1.2655] | 0.5694 [0.4707, 0.7162] | 1.5627 [1.4351, 1.8104] |
| GEARS [shuffled across perts] | 9 | 4.7107 [4.6669, 5.1127] | 0.2894 [0.2322, 0.3676] | 0.0109 [0.0010, 0.0131] | n.a. | n.a. | n.a. |
| linear[PCA] (K=10, lam=0.1) | 3 | 4.1665 [3.9425, 4.3857] | 0.4560 [0.4054, 0.4796] | 0.2181 [0.1416, 0.2645] | 1.0294 [0.8385, 1.0678] | 0.5413 [0.4882, 0.5961] | 1.7267 [1.7174, 1.9728] |
| linear[PCA] (K=10, lam=0.1) [shuffled across perts] | 3 | 4.3999 [4.3922, 4.7622] | 0.3879 [0.3565, 0.3955] | 0.0156 [-0.0017, 0.0194] | n.a. | n.a. | n.a. |
| mean | 3 | 4.2350 [4.1091, 4.5359] | 0.4416 [0.3816, 0.4436] | n.a. | 0.8940 [0.7800, 0.9975] | 0.4598 [0.4094, 0.4682] | 1.9441 [1.9053, 2.1304] |
| no_change | 3 | 4.5728 [4.5364, 5.0358] | n.a. | n.a. | n.a. | n.a. | n.a. |

## A: GEARS & PCA predictable - target included

| model | n runs | L2 | Pearson | Pearson (centred) | amp. ratio, mean delta | amp. ratio, per-pert norm | collapse ratio pred/obs |
|---|---|---|---|---|---|---|---|
| GEARS | 9 | 4.2148 [4.0005, 4.8173] | 0.4019 [0.3477, 0.4248] | 0.2324 [0.1460, 0.3209] | 0.8548 [0.7217, 1.2593] | 0.5554 [0.4618, 0.7027] | 1.5788 [1.4597, 1.8292] |
| GEARS [shuffled across perts] | 9 | 4.8039 [4.7425, 5.1778] | 0.2857 [0.2375, 0.3633] | 0.0103 [0.0052, 0.0186] | n.a. | n.a. | n.a. |
| linear[PCA] (K=10, lam=0.1) | 3 | 4.2392 [4.0675, 4.4971] | 0.4472 [0.4007, 0.4685] | 0.2111 [0.1388, 0.2551] | 1.0213 [0.8344, 1.0620] | 0.5329 [0.4790, 0.5816] | 1.7559 [1.7419, 1.9931] |
| linear[PCA] (K=10, lam=0.1) [shuffled across perts] | 3 | 4.4868 [4.4524, 4.8665] | 0.3818 [0.3552, 0.3899] | 0.0109 [0.0096, 0.0299] | n.a. | n.a. | n.a. |
| mean | 3 | 4.3051 [4.2282, 4.6348] | 0.4319 [0.3770, 0.4357] | n.a. | 0.8866 [0.7761, 0.9920] | 0.4485 [0.4016, 0.4609] | 1.9771 [1.9324, 2.1525] |
| no_change | 3 | 4.6903 [4.6098, 5.1345] | n.a. | n.a. | n.a. | n.a. | n.a. |

## B: GEARS predictable - target excluded

| model | n runs | L2 | Pearson | Pearson (centred) | amp. ratio, mean delta | amp. ratio, per-pert norm | collapse ratio pred/obs |
|---|---|---|---|---|---|---|---|
| GEARS | 9 | 4.0887 [3.9272, 4.6970] | 0.3781 [0.3626, 0.4351] | 0.2561 [0.1296, 0.3511] | 0.9110 [0.8773, 1.4904] | 0.5674 [0.5045, 0.7561] | 1.6685 [1.5699, 1.9711] |
| GEARS [shuffled across perts] | 9 | 4.7538 [4.6541, 5.0053] | 0.2848 [0.2515, 0.3254] | 0.0175 [0.0064, 0.0190] | n.a. | n.a. | n.a. |
| mean | 3 | 4.2374 [4.1893, 4.3507] | 0.4009 [0.3848, 0.4109] | n.a. | 0.9458 [0.9203, 1.0253] | 0.4536 [0.4362, 0.4713] | 2.1098 [2.0852, 2.1753] |
| no_change | 3 | 4.6369 [4.5072, 4.7273] | n.a. | n.a. | n.a. | n.a. | n.a. |

## B: GEARS predictable - target included

| model | n runs | L2 | Pearson | Pearson (centred) | amp. ratio, mean delta | amp. ratio, per-pert norm | collapse ratio pred/obs |
|---|---|---|---|---|---|---|---|
| GEARS | 9 | 4.1387 [3.9769, 4.7334] | 0.3755 [0.3601, 0.4320] | 0.2542 [0.1273, 0.3473] | 0.9084 [0.8755, 1.4873] | 0.5638 [0.5005, 0.7502] | 1.6756 [1.5805, 1.9825] |
| GEARS [shuffled across perts] | 9 | 4.8143 [4.6914, 5.0458] | 0.2800 [0.2468, 0.3227] | 0.0128 [0.0028, 0.0164] | n.a. | n.a. | n.a. |
| mean | 3 | 4.2836 [4.2163, 4.3888] | 0.3979 [0.3830, 0.4071] | n.a. | 0.9430 [0.9183, 1.0234] | 0.4492 [0.4327, 0.4684] | 2.1220 [2.0992, 2.1847] |
| no_change | 3 | 4.6824 [4.5355, 4.7654] | n.a. | n.a. | n.a. | n.a. | n.a. |

## collapse(obs) per split (a property of the data split, not of the model; set A / B, target included)

| set | split | collapse(obs) |
|---|---|---|
| A: GEARS & PCA predictable | seed1 | 0.4646 |
| A: GEARS & PCA predictable | seed2 | 0.5058 |
| A: GEARS & PCA predictable | seed3 | 0.5175 |
| B: GEARS predictable | seed1 | 0.4577 |
| B: GEARS predictable | seed2 | 0.4764 |
| B: GEARS predictable | seed3 | 0.4712 |

## Per split (set A, target included): median [min, max] over the 3 training seeds for GEARS; deterministic baselines one value

| split | model | L2 | Pearson | Pearson (centred) | amp. ratio, per-pert norm |
|---|---|---|---|---|---|
| seed1 | GEARS | 4.2148 [4.1823, 4.3029] | 0.3620 [0.3477, 0.3838] | 0.2324 [0.1972, 0.2338] | 0.5479 [0.5115, 0.5778] |
| seed1 | linear[PCA] (K=10, lam=0.1) | 4.2392 [4.2392, 4.2392] | 0.4007 [0.4007, 0.4007] | 0.1388 [0.1388, 0.1388] | 0.5329 [0.5329, 0.5329] |
| seed1 | mean | 4.3051 [4.3051, 4.3051] | 0.3770 [0.3770, 0.3770] | n.a. | 0.4609 [0.4609, 0.4609] |
| seed1 | no_change | 4.6098 [4.6098, 4.6098] | n.a. | n.a. | n.a. |
| seed1 | GEARS [shuffled across perts] | 4.7448 [4.7425, 4.7612] | 0.2732 [0.2375, 0.2857] | 0.0177 [0.0075, 0.0186] | n.a. |
| seed2 | GEARS | 4.0408 [4.0005, 4.1337] | 0.4047 [0.4019, 0.4233] | 0.3027 [0.2519, 0.3209] | 0.5554 [0.5473, 0.5621] |
| seed2 | linear[PCA] (K=10, lam=0.1) | 4.0675 [4.0675, 4.0675] | 0.4685 [0.4685, 0.4685] | 0.2551 [0.2551, 0.2551] | 0.5816 [0.5816, 0.5816] |
| seed2 | mean | 4.2282 [4.2282, 4.2282] | 0.4319 [0.4319, 0.4319] | n.a. | 0.4485 [0.4485, 0.4485] |
| seed2 | no_change | 4.6903 [4.6903, 4.6903] | n.a. | n.a. | n.a. |
| seed2 | GEARS [shuffled across perts] | 4.8039 [4.8008, 4.8381] | 0.2693 [0.2662, 0.2998] | 0.0103 [0.0101, 0.0115] | n.a. |
| seed3 | GEARS | 4.6479 [4.5383, 4.8173] | 0.4155 [0.3930, 0.4248] | 0.1665 [0.1460, 0.1928] | 0.6058 [0.4618, 0.7027] |
| seed3 | linear[PCA] (K=10, lam=0.1) | 4.4971 [4.4971, 4.4971] | 0.4472 [0.4472, 0.4472] | 0.2111 [0.2111, 0.2111] | 0.4790 [0.4790, 0.4790] |
| seed3 | mean | 4.6348 [4.6348, 4.6348] | 0.4357 [0.4357, 0.4357] | n.a. | 0.4016 [0.4016, 0.4016] |
| seed3 | no_change | 5.1345 [5.1345, 5.1345] | n.a. | n.a. | n.a. |
| seed3 | GEARS [shuffled across perts] | 5.0482 [5.0282, 5.1778] | 0.3601 [0.3197, 0.3633] | 0.0058 [0.0052, 0.0147] | n.a. |

## GEARS spread decomposition (set A, target included)

| metric | IQR over runs | between-split range (of per-split medians) | mean within-split range (between seeds) |
|---|---|---|---|
| L2 | 0.4046 | 0.6072 | 0.1776 |
| Pearson | 0.0317 | 0.0535 | 0.0297 |
| Pearson_centred | 0.0590 | 0.1363 | 0.0508 |
| mean_delta_ratio | 0.2493 | 0.2490 | 0.2793 |
| pert_norm_ratio | 0.0305 | 0.0579 | 0.1073 |
| collapse_ratio | 0.2287 | 0.2804 | 0.1802 |

Collapse index = ||mean_p X_p|| / mean_p ||X_p|| (1 = all perturbations along one shared direction); GEARS collapse ratio pred/obs > 1 in 9 of 9 runs.

GEARS mean-delta amplitude ratio > 1 in 3 of 9 runs; per-perturbation norm ratio < 1 in 9 of 9.
