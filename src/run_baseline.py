"""python -m src.run_baseline --config configs/k562_repro_5000hvg.yaml [--preprocess]"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from .data import PerturbationDataset, predictable_subset, preprocess
from .evaluate import format_table, get_splits, run_split, summarize
from .utils import REPO_ROOT, load_config, make_run_dir


def qc_summary(qc: pd.DataFrame) -> str:
    t = qc[qc.target_in_matrix]
    pct = 100 * (1 - 2 ** t.target_log2fc)
    lines = [f"perturbations: {len(qc)}",
             f"  <25 cells: {(qc.n_cells < 25).sum()}   (<10: {(qc.n_cells < 10).sum()}, <50: {(qc.n_cells < 50).sum()})",
             f"  target gene not in 8563-gene matrix: {(~qc.target_in_matrix).sum()}",
             f"  target in matrix but not in HVG read-out set (cannot be predicted by G-rows P): "
             f"{(qc.target_in_matrix & ~qc.target_in_readout).sum()}",
             f"  target in read-out set (evaluable): {qc.target_in_readout.sum()}",
             f"knockdown (log2FC of mean normalised target expr vs control; {len(t)} targets in matrix):",
             f"  no knockdown (log2FC >= -0.1, <~7% reduction): {(t.target_log2fc >= -0.1).sum()}",
             f"  weak (-0.58 < log2FC < -0.1, 7-33% reduction): {((t.target_log2fc < -0.1) & (t.target_log2fc > -0.58)).sum()}",
             f"  >=33% reduction (log2FC <= -0.58): {(t.target_log2fc <= -0.58).sum()}",
             f"  median reduction among targets: {np.nanmedian(pct):.0f}%",
             f"  low-expressed targets (ctrl mean norm < 0.5, knockdown poorly measurable): {(t.target_ctrl_mean_norm < 0.5).sum()}"]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--preprocess", action="store_true", help="(re)build the processed dataset")
    a = ap.parse_args()
    cfg = load_config(a.config)
    d = cfg["data"]
    proc, qc_path = REPO_ROOT / d["processed_path"], REPO_ROOT / d["qc_path"]
    if a.preprocess or not proc.exists():
        ds, qc = preprocess(REPO_ROOT / d["raw_path"], cfg)
        ds.save(proc)
        qc_path.parent.mkdir(parents=True, exist_ok=True)
        qc.to_csv(qc_path, index=False)
    ds_all = PerturbationDataset.load(proc)
    qc = pd.read_csv(qc_path)
    ds = predictable_subset(ds_all)

    run_dir = make_run_dir(cfg, cfg["name"])
    qc.to_csv(run_dir / "qc_perturbations.csv", index=False)
    (run_dir / "qc_summary.txt").write_text(qc_summary(qc) + "\n")
    print(qc_summary(qc))

    per_pert, infos, grids = [], [], []
    for label, seed, split in get_splits(ds, cfg):
        pp, info, grid = run_split(ds, cfg, split, label, seed)
        per_pert.append(pp); infos.append(info); grids.append(grid.assign(split=label))
    per_pert = pd.concat(per_pert)
    summ = summarize(per_pert)
    per_pert.to_csv(run_dir / "per_perturbation_metrics.csv", index=False)
    summ.to_csv(run_dir / "summary.csv", index=False)
    pd.concat(grids).to_csv(run_dir / "val_grid.csv", index=False)
    (run_dir / "chosen_hparams.json").write_text(json.dumps(infos, indent=2))
    (run_dir / "dataset_meta.json").write_text(json.dumps({**ds.meta, "n_perts_evaluated": len(ds.perts)}, indent=2))
    primary = "excluded" if cfg["eval"]["exclude_target_primary"] else "included"
    other = "included" if primary == "excluded" else "excluded"
    md = f"## target gene {primary} (primary)\n\n{format_table(summ, primary)}\n\n## target gene {other}\n\n{format_table(summ, other)}\n"
    (run_dir / "summary.md").write_text(md)
    print(md, "\n->", run_dir)


if __name__ == "__main__":
    main()
