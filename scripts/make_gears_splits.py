"""Regenerate GEARS 'simulation' splits exactly as `PertData.prepare_split(split='simulation', seed=s)` does
(gears==0.1.2, defaults train_gene_set_size=0.75, combo_seen2_train_frac=0.75), WITHOUT instantiating PertData:
its constructor downloads and unpickles gene2go_all.pkl from Dataverse, which splitting does not need.

Run in the separate env:  mamba run -n hic-perturb-gears python scripts/make_gears_splits.py --seeds 1 2 3 4 5
Writes configs/splits/gears_k562_seed{s}.json = {"train": [...], "val": [...], "test": [...]} (condition names 'GENE+ctrl')."""
import argparse
import json
from pathlib import Path

import scanpy as sc
from gears.data_utils import DataSplitter

ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", default="data/gears/replogle_k562_essential/perturb_processed.h5ad")
ap.add_argument("--seeds", type=int, nargs="+", required=True)
ap.add_argument("--out", default="configs/splits")
a = ap.parse_args()

adata = sc.read_h5ad(a.h5ad)
Path(a.out).mkdir(parents=True, exist_ok=True)
for seed in a.seeds:
    DS = DataSplitter(adata, split_type="simulation")
    ad2, subgroup = DS.split_data(train_gene_set_size=0.75, combo_seen2_train_frac=0.75, seed=seed,
                                  test_perts=None, only_test_set_perts=False)
    s2c = dict(ad2.obs.groupby("split").agg({"condition": lambda x: x}).condition)   # verbatim from prepare_split
    s2c = {k: v.unique().tolist() for k, v in s2c.items()}
    with open(Path(a.out) / f"gears_k562_seed{seed}.json", "w") as f:
        json.dump(s2c, f)
    print(seed, {k: len(v) for k, v in s2c.items()})
