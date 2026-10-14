"""python -m scripts.make_prereg_configs -> configs/prereg/*.yaml : every analysis named in PREREGISTRATION.md, written BEFORE any is run.
All share the frozen design: GEARS splits 1-5, K=10 / lambda=0.1 paper settings (+ fixed K=20/50 at lambda 0.1 as secondary), training-split fit only."""
from pathlib import Path

import yaml

OUT = Path("configs/prereg"); OUT.mkdir(parents=True, exist_ok=True)
CTRL = [f"b500_lin_randtrans_s{k}" for k in range(10)]


def cfg(name, embeddings, note):
    d = {"name": f"prereg_{name}", "data": {"source": "gears", "raw_path": "data/gears/replogle_k562_essential/perturb_processed.h5ad",
                                            "processed_path": "data/processed/k562_gears.npz", "qc_path": "results/qc_perturbations_gears.csv",
                                            "pert_col": "condition", "control_label": "ctrl"},
         "eval": {"seeds": [], "split_fracs": [0.7, 0.1, 0.2],
                  "split_files": {f"gears_seed{s}": f"configs/splits/gears_k562_seed{s}.json" for s in range(1, 6)},
                  "embeddings": embeddings, "refit_trainval": False, "paper_K": 10, "paper_lambda": 0.1, "K_grid": [10, 20, 50], "lambda_grid": [0.1],
                  "extra_fixed_hp": [[20, 0.1], [50, 0.1]], "n_top_de": 20, "n_top_expr": 1000, "exclude_target_primary": True}}
    (OUT / f"{name}.yaml").write_text(f"# PRE-REGISTERED ({note}). DO NOT RUN before the pre-registration is frozen.\n" + yaml.safe_dump(d, sort_keys=False))


cfg("primary", ["b500_hic_trans"] + CTRL, "H1 primary: b500_hic_trans vs the 10 matched random-trans control seeds, K=10 lambda=0.1; K=50 fixed = KEY secondary")
cfg("s1_factorial", ["b500_hic_trans", "b500_hic_cis", "b500_lin_trans", "b500_lin_cis"] + [f"b500_hic_randtrans_s{k}" for k in range(10)] + CTRL +
    ["b500_hic_cis:percomp", "b500_lin_cis:percomp"], "S1 factorial cis{hic,lin} x trans{hic,randtrans,none} + per-component cis-only variants")
# S2 nulls are evaluated by the lean driver `python -m src.nulls --config configs/prereg/s2_nulls.yaml` (fixed K and lambda only):
# permutation x1000 and random Gaussian x1000 (no new eigendecomposition), rewire x20; K=10 decides H1 support, K=50 is the key-secondary analogue.
(OUT / "s2_nulls.yaml").write_text("# PRE-REGISTERED (S2 nulls + H1 support condition). DO NOT RUN before the pre-registration is frozen.\n" + yaml.safe_dump({
    "name": "prereg_s2_nulls", "data": {"processed_path": "data/processed/k562_gears.npz"}, "hic": "b500_hic_trans", "control": CTRL, "K": [10, 50], "lambda": 0.1,
    "nulls": {"permutation": {"n": 1000}, "random_gaussian": {"n": 1000}, "rewire": {"n": 20}},
    "eval": {"split_files": {f"gears_seed{s}": f"configs/splits/gears_k562_seed{s}.json" for s in range(1, 6)}, "metric": "L2 top1000 expressed genes, target excluded",
             "refit_trainval": False}}, sort_keys=False))
cfg("s3_resid", ["b500_hic_trans", "b500_hic_trans:resid", "b500_hic_cis:resid"] + CTRL + [f"{c}:resid" for c in CTRL],
    "S3 Tier-2 test: b500_hic_trans:resid vs b500_lin_randtrans:resid (identical residualization, 10 control seeds); hic_cis:resid exploratory")
cfg("s4_pca_extended", ["train_pca_expr", "b500_hic_trans"] + CTRL, "S4/S9 PCA reference on its predictable subset + extended eval set vs mean baseline")
cfg("s8_resolution", ["b250_hic_trans"] + [f"b250_lin_randtrans_s{k}" for k in range(10)] + ["b250_hic_cis", "b250_lin_cis", "b250_lin_trans", "hic_local25k"],
    "S8 250 kb analogue of the primary contrast + local 25 kb graph")
cfg("s10_exclusion", ["b500x_hic_trans"] + [f"b500x_lin_randtrans_s{k}" for k in range(10)], "S10 K562 rearranged/amplified regions excluded (configs/k562_rearranged_regions_hg38.bed)")
print(sorted(p.name for p in OUT.glob("*.yaml")))
