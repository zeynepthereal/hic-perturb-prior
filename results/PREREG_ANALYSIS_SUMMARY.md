# Pre-registered analysis — results record (numbers and pre-registered labels only)

Pre-registration: `PREREGISTRATION.md` at tag `prereg-v1` (commit 7bb6c399992a06a236b8698c28e51eea2934a8be). All runs below logged code commit **fd3da67621befe4fc8ead49ea4c963ed6a98acc2** (clean tree).
Run directories: primary `results/20261014-134117_prereg_primary`, nulls `results/20261014-134435_prereg_s2_nulls`, S3 `results/20261014-135031_prereg_s3_resid`; sanity checks `results/20261014-133938_sanity1_pca_only`, `results/20261014-134015_sanity5_lean_vs_standard`.
Raw per-perturbation metrics of the primary and S3 runs are committed gzip-compressed (`per_perturbation_metrics.csv.gz`; `gunzip -k` to reuse with `scripts/run_primary_analysis.py`).

## §7 sanity checks
| check | result |
|---|---|
| 1 PCA + mean baseline reproduce earlier k562_gears | PASS (max \|ΔL2\| = 0 over 18 combinations; mean L2 4.2557 / 4.1501 in both) |
| 2 random-Gaussian arms vs mean baseline, threshold δ (reported) | K=10: 679/1000 arms within δ, differences +0.0050…+0.0311 (δ = 0.02141); K=50: 0/1000 within δ, differences +0.1056…+0.1655 |
| 3 predictable sets identical across arms | PASS (primary 1,092; S3 1,092; S10 1,036) |
| 4 eligibility gate | PASS (0.0 % all-zero rows, all 11 primary arms, K=10 and 50) |
| 5 lean evaluator = standard evaluator | PASS (worst max \|ΔL2\| = 3.6e-15; 22 arm×K combinations on the primary run + 4 on the sanity run) |

## Primary analysis (§3), K = 10, λ = 0.1 — b500_hic_trans vs b500_lin_randtrans (10 seeds)
| quantity | value |
|---|---|
| **decision label** | **SIGNIFICANT: Hi-C better** |
| D (control − Hi-C; positive = Hi-C better) | 0.041119 |
| two-sided p (sign-flip, 100,000 perm.) | 0.021720 |
| 95 % cluster-bootstrap CI | [0.007542, 0.084656] |
| SESOI δ (0.5 % × mean control L2 4.2821) | 0.021410 |
| TOST p (equivalence) | 0.826412 (equivalent: False) |
| fraction of perturbations with a_i > 0 | 0.5305 (837 unique perturbations, 1365 cases) |
| per-split mean difference | gears_seed1: 0.15049, gears_seed2: 0.05147, gears_seed3: 0.01974, gears_seed4: -0.01878, gears_seed5: 0.00268 |
| descriptive seed rank (per split; pooled) | 1: 1/11, 2: 1/11, 3: 1/11, 4: 11/11, 5: 6/11; pooled 1/11 (controls beaten 10/10) |
| per-control-seed D | 0.04489, 0.03296, 0.04182, 0.03690, 0.03190, 0.06870, 0.04445, 0.01805, 0.05167, 0.03985 (seeds 0–9) |

## S2 nulls (Tier-1 condition)
| K | D | permutation null (1,000): q95 / p | random-Gaussian null (1,000): q95 / p | rewire null (20): max / p | all conditions |
|---|---|---|---|---|---|
| 10 | 0.041119 | 0.040332 / 0.036963 | 0.040042 / 0.026973 | 0.031735 / 0.047619 | True |
| 50 | 0.067531 | 0.158780 / 0.976024 | 0.172652 / 1.000000 | 0.051934 / 0.047619 | False |

K = 10 is decisive; K = 50 is the key secondary (null part only; the K = 50 contrast itself has not been run).

## Tier outcomes (§3)
| tier | outcome |
|---|---|
| Tier 1 (H1 supported): decision "SIGNIFICANT: Hi-C better" and D above all three K=10 nulls | **True** |
| Tier 2 (mechanistic): Tier 1 and S3 (`b500_hic_trans:resid` vs `b500_lin_randtrans:resid`) significant, Hi-C better | **False** — S3 decision label: **INCONCLUSIVE** (D = -0.009257, p = 0.226298, 95 % CI [-0.024779, 0.004643], TOST p = 0.054909, δ = 0.021254) |

No other secondary analysis has been run.
