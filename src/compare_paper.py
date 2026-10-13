"""python -m src.compare_paper results/<run_dir>  -> writes <run_dir>/comparison_vs_paper.md

Paper numbers: configs/paper_reference_replogle_k562.csv (from the authors' source_data/single_perturbation_prediction.xlsx,
Panel A, replogle_k562_essential, train=='test'; mean over perturbations per split, then mean/sd over their 2 splits;
L2 over the 1000 most expressed genes; 'r2_delta' column = Pearson delta). The paper does not exclude target genes, so the
like-for-like comparison is our target='included', subset='top1000_expr'."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .evaluate import lm_name
from .utils import REPO_ROOT


def load_paper():
    r = pd.read_csv(REPO_ROOT / "configs/paper_reference_replogle_k562.csv", header=[0, 1], index_col=[0, 1])
    r.columns = ["l2_mean", "l2_std", "pearson_mean", "pearson_std"]
    return r.droplevel(0)


def main(run_dir: str):
    run = Path(run_dir)
    ours = pd.read_csv(run / "summary.csv")
    paper = load_paper()
    f = lambda m, s: "n/a" if np.isnan(m) else f"{m:.3f} ± {0 if np.isnan(s) else s:.3f}"
    pmap = {"paper: mean": "mean", "paper: linear model": "lpm_selftrained", "paper: GEARS": "gears"}
    rows = []
    for k, m in pmap.items():
        r = paper.loc[m]
        rows.append((k, f(r.l2_mean, r.l2_std), f(r.pearson_mean, r.pearson_std)))
    ref_ratio = (paper.loc["lpm_selftrained"].l2_mean / paper.loc["mean"].l2_mean,
                 paper.loc["lpm_selftrained"].pearson_mean / paper.loc["mean"].pearson_mean)
    emb = [m for m in ours.model.unique() if m.startswith("linear[")][0].split("]")[0][len("linear["):]
    PAPER_LM, TUNED_LM = lm_name(emb, False), lm_name(emb, True)
    out = ["# Comparison with Ahlmann-Eltze et al. 2025 (Replogle K562)\n",
           "Paper rows: top-1000-expressed-genes L2 and Pearson delta, mean ± sd over their 2 splits (GEARS-processed data).",
           "Our rows: same metric definitions, mean ± sd over our splits. Every table states its eval_set and size.\n"]
    for es in [e for e in ("common", "paper_pertset") if e in set(ours.eval_set)]:
      for target in ("included", "excluded"):
        t = ours[(ours.eval_set == es) & (ours.target == target) & (ours.subset == "top1000_expr")].set_index("model")
        n, k = t.n_perts.mean(), t.n_splits.max()
        tbl = list(rows)
        for m in ("no_change", "mean", PAPER_LM, TUNED_LM):
            r = t.loc[m]
            tbl.append((f"ours: {m}", f(r.l2_mean, r.l2_std), f(r.pearson_mean, r.pearson_std)))
        df = pd.DataFrame(tbl, columns=["", "L2 (top1000 expr)", "Pearson delta (top1000 expr)"]).set_index("")
        out.append(f"## eval_set = {es}, target gene {target}" + (" (like-for-like with the paper)" if target == "included" else "")
                   + f"\n\n{n:.0f} test perturbations per split, {k} splits.\n")
        out.append(df.to_markdown() + "\n")
        rt = []
        for m in (PAPER_LM, TUNED_LM):
            rt.append((f"ours: {m}", t.loc[m].l2_mean / t.loc["mean"].l2_mean, t.loc[m].pearson_mean / t.loc["mean"].pearson_mean))
        rt.append(("paper: linear model", *ref_ratio))
        out.append("Ratio to the mean baseline (linear / mean; L2 < 1 and Pearson > 1 mean the linear model wins):\n")
        out.append(pd.DataFrame(rt, columns=["", "L2 ratio", "Pearson ratio"]).set_index("").round(3).to_markdown() + "\n")
    (run / "comparison_vs_paper.md").write_text("\n".join(out))
    print("\n".join(out))


if __name__ == "__main__":
    main(sys.argv[1])
