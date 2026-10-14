"""python -m scripts.diagnose_randtrans -> results/hic_v2/randtrans/*.csv + diagnostics_randtrans.md
Diagnostics for the matched-random-trans cells only (lin_randtrans = primary control, hic_randtrans = factorial secondary), mean +- sd over seeds,
next to the already-computed reference cells."""
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from src.hic.diagnostics import ab_rows, r2_rows, structure
from src.hic.embeddings import ART

OUT = Path("results/hic_v2"); (OUT / "randtrans").mkdir(exist_ok=True)
seeds = sorted(int(f.name.split("_s")[1].split("_")[0]) for f in ART.glob("b500_lin_randtrans_s*_spectral.npz"))
lines = ["# Matched random-trans controls (500 kb)\n", f"seeds found: {seeds}\n"]

# sanity: trans degrees identical to hic_trans, trans edges inter-chromosomal, cis part untouched
T_deg = np.diff(sp.load_npz(ART / "b500_hic_trans.npz").indptr) - np.diff(sp.load_npz(ART / "b500_hic_cis.npz").indptr)
chk = []
for s in seeds:
    for cis in ("hic", "lin"):
        W = sp.load_npz(ART / f"b500_{cis}_randtrans_s{s}.npz"); C = sp.load_npz(ART / f"b500_{cis}_cis.npz")
        chk.append(bool((np.diff(W.indptr) - np.diff(C.indptr) == T_deg).all()))
lines += [f"trans-degree sequence identical to hic_trans in all {len(chk)} random graphs: {all(chk)}\n"]

st, r2, ab = [], [], []
for cis in ("lin", "hic"):
    for s in seeds:
        v = f"b500_{cis}_randtrans_s{s}"
        st.append({**structure(v), "group": f"b500_{cis}_randtrans"}); r2 += [{**r, "group": f"b500_{cis}_randtrans"} for r in r2_rows(v, OUT)]
        ab += [{**r, "group": f"b500_{cis}_randtrans"} for r in ab_rows(v, OUT)]
st, r2, ab = pd.DataFrame(st), pd.DataFrame(r2), pd.DataFrame(ab)
for n, d in (("structure", st), ("r2", r2), ("ab", ab)):
    d.to_csv(OUT / "randtrans" / f"{n}_per_seed.csv", index=False)
agg = lambda d, keys: d.groupby(keys).agg(["mean", "std"]).round(4)
lines += ["## Structure (mean, sd over seeds)\n", agg(st.drop(columns="variant"), ["group"]).T.to_markdown(), ""]
ref_r2 = pd.read_csv(OUT / "r2_chrom_position.csv"); ref_ab = pd.read_csv(OUT / "ab_correlation.csv")
refs = ["b500_hic_trans", "b500_lin_trans", "b500_hic_cis", "b500_lin_cis"]
r2s = r2.drop(columns="variant").groupby(["group", "dims"]).agg(["mean", "std"]).round(3)
lines += ["## R² on chromosome / position (mean, sd over seeds)\n", r2s.to_markdown(), "",
          "Reference cells (same metrics):\n", ref_r2[ref_r2.variant.isin(refs)].round(3).to_markdown(index=False), ""]
abs_ = ab.drop(columns="variant").groupby(["group", "ab", "dims"]).agg(["mean", "std"]).round(3)
lines += ["## Correlation with A/B track (mean, sd over seeds)\n", abs_.to_markdown(), "",
          "Reference cells:\n", ref_ab[ref_ab.variant.isin(refs)][["variant", "ab", "dims", "multiple_R", "R2_beyond_position", "R2_AB_from_position_only"]].round(3).to_markdown(index=False), ""]
(OUT / "diagnostics_randtrans.md").write_text("\n".join(lines))
print("\n".join(lines))
