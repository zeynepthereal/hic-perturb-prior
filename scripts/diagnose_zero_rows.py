"""python -m scripts.diagnose_zero_rows -> results/hic_v2/zero_rows.csv
Fraction of all-zero P rows for EVERY embedding (pre-registered gate: > 5 % all-zero rows => cannot be a primary comparator).
Two views: (a) over the P rows of the 1,092 GEARS perturbation targets (what the model sees), (b) over all embedded graph nodes.
Also reports the fraction of perturbations that are not predictable at all (NaN rows)."""
import re
from pathlib import Path

import numpy as np
import pandas as pd

from src.baselines import EMBEDDINGS
from src.data import PerturbationDataset
from src.hic.embeddings import ART

ds = PerturbationDataset.load("data/processed/k562_gears.npz")
names = sorted({f.name[: -len("_spectral.npz")] for f in ART.glob("*_spectral.npz") if "__rewire" not in f.name})
extra = [f"{b}:{s}" for b in ("b500_hic_trans", "b500_hic_cis", "b500_lin_cis", "b500_lin_trans") for s in ("resid", "perm:0", "random:0", "rewire:0") if
         not (s == "resid" and "lin" in b)] + [f"{b}:percomp" for b in ("b500_hic_cis", "b500_lin_cis", "b250_hic_cis", "b250_lin_cis")]
extra += [f"b500_lin_randtrans_s{k}:resid" for k in range(10)]   # Tier-2 control (v0.5): residualized with the same basis as b500_hic_trans:resid
names = [n for n in names if re.match(r"^(hic_local25k|b\d+x?_(hic|lin)_(trans|cis|randtrans_s\d+))$", n)] + [n for n in extra if (ART / (n.split(":")[0] + "_spectral.npz")).exists()]
rows = []
for n in names:
    e = EMBEDDINGS[n]
    for K in (10, 50):
        P = e.perturbation_matrix(ds, K) if hasattr(e, "perturbation_matrix") else e.fit(ds.delta[:50], np.arange(50), ds, K)[1]
        nan = np.isnan(P).any(1)
        zero = (~nan) & (np.abs(P) < 1e-8).all(1)
        E, ok = e._node_matrix(K) if hasattr(e, "_node_matrix") else (None, None)
        zn = float(((np.abs(E[ok][:, :K]) < 1e-8).all(1)).mean()) if E is not None else np.nan
        rows.append(dict(embedding=n, K=K, perts_predictable=int((~nan).sum()), frac_pert_rows_allzero=float(zero.sum() / max(1, (~nan).sum())),
                         frac_node_rows_allzero=zn))
df = pd.DataFrame(rows)
df["primary_comparator_eligible"] = (df.frac_pert_rows_allzero <= 0.05) & (df.frac_node_rows_allzero <= 0.05)
df.to_csv("results/hic_v2/zero_rows.csv", index=False)
agg = df.assign(group=df.embedding.str.replace(r"_s\d+", "_sXX", regex=True)).groupby(["group", "K"]).agg(
    n=("embedding", "count"), perts_predictable=("perts_predictable", "mean"), pert_rows_allzero_mean=("frac_pert_rows_allzero", "mean"),
    pert_rows_allzero_max=("frac_pert_rows_allzero", "max"), node_rows_allzero_mean=("frac_node_rows_allzero", "mean"),
    eligible_all=("primary_comparator_eligible", "all")).reset_index()
print(agg.round(4).to_markdown(index=False))
agg.to_csv("results/hic_v2/zero_rows_summary.csv", index=False)
