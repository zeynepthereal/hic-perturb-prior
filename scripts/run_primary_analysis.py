"""Thin command-line wrapper around the frozen analysis code (src/stats.py as of tag prereg-v1). It makes NO analytical choice of its own:
every number comes from `src.stats.primary_analysis` (paired sign-flip test, TOST at the SESOI, decision rule, descriptive seed ranks).

  python -m scripts.run_primary_analysis RUN_DIR                              # primary: b500_hic_trans vs b500_lin_randtrans_s0..9, K=10, lambda=0.1
  python -m scripts.run_primary_analysis RUN_DIR --kind fixed --K 50 --lam 0.1  # key secondary (K = 50)
  python -m scripts.run_primary_analysis S3_RUN --hic b500_hic_trans:resid --control-template 'b500_lin_randtrans_s{k}:resid'   # Tier-2 test (S3)
  python -m scripts.run_primary_analysis RUN_DIR --nulls results/<null_run>/nulls.json [--tier2-from S3_RUN/primary_analysis.json]

RUN_DIR holds `per_perturbation_metrics.csv` written by `python -m src.run_baseline --config configs/prereg/<analysis>.yaml`.
Optional `--nulls` / `--tier2-from` only TRANSCRIBE the Tier-1 / Tier-2 rule of PREREGISTRATION.md section 3 from results computed by the frozen
code (`src/nulls.py` -> nulls.json; a second run of this script); if this wrapper and the prose ever disagree, the prose governs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from src.evaluate import fixed_name, lm_name
from src.stats import primary_analysis
from src.utils import git_state

BETTER = "SIGNIFICANT: Hi-C better"


def model_name(emb: str, kind: str, K: int, lam: float) -> str:
    return lm_name(emb, False) if kind == "paper" else fixed_name(emb, K, lam)      # 'paper hp K=10, lam=0.1' | 'fixed K=.., lam=..'


def build_parser() -> argparse.ArgumentParser:
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("run_dir")
    a.add_argument("--hic", default="b500_hic_trans", help="embedding name of the Hi-C arm")
    a.add_argument("--control", nargs="+", default=None, help="explicit control embedding names (default: --control-template with --n-control seeds)")
    a.add_argument("--control-template", default="b500_lin_randtrans_s{k}"); a.add_argument("--n-control", type=int, default=10)
    a.add_argument("--kind", choices=["paper", "fixed"], default="paper", help="paper = K=10, lambda=0.1 model; fixed = --K/--lam model")
    a.add_argument("--K", type=int, default=10); a.add_argument("--lam", type=float, default=0.1)
    a.add_argument("--eval-set", default="common"); a.add_argument("--target", default="excluded")
    a.add_argument("--subset", default="top1000_expr"); a.add_argument("--metric", default="l2")
    a.add_argument("--n-perm", type=int, default=100_000); a.add_argument("--seed", type=int, default=0)
    a.add_argument("--nulls", default=None, help="nulls.json written by `python -m src.nulls` (adds the Tier-1 null condition)")
    a.add_argument("--tier2-from", default=None, help="primary_analysis.json of the S3 contrast (resid vs resid); requires --nulls")
    a.add_argument("--out", default=None, help="output JSON (default RUN_DIR/primary_analysis.json)")
    return a


def run(argv=None) -> dict:
    args = build_parser().parse_args(argv)
    run_dir = Path(args.run_dir)
    pp = pd.read_csv(run_dir / "per_perturbation_metrics.csv")
    controls = args.control or [args.control_template.format(k=k) for k in range(args.n_control)]
    new = model_name(args.hic, args.kind, args.K, args.lam)
    ctrl = [model_name(c, args.kind, args.K, args.lam) for c in controls]
    r = primary_analysis(pp, new, ctrl, n_perm=args.n_perm, seed=args.seed, eval_set=args.eval_set, target=args.target, subset=args.subset, metric=args.metric)
    r["cli"] = {"args": vars(args), "hic_model": new, "control_models": ctrl, "git": git_state()}
    if args.nulls:
        K = 10 if args.kind == "paper" else args.K
        nk = json.load(open(args.nulls))[str(K)]
        r["null_condition"] = {k: v for k, v in nk.items() if k != "null_D"}
        r["D_diff_nulls_vs_primary"] = abs(nk["D_obs"] - r["D"])           # sanity check 5: lean evaluator vs standard evaluator
        r["tier1_supported"] = bool(r["decision"] == BETTER and nk["support"])
        if args.tier2_from:
            r["tier2_supported"] = bool(r["tier1_supported"] and json.load(open(args.tier2_from))["decision"] == BETTER)
    elif args.tier2_from:
        raise SystemExit("--tier2-from requires --nulls (Tier 2 is only evaluated after Tier 1)")
    out = Path(args.out) if args.out else run_dir / "primary_analysis.json"
    out.write_text(json.dumps(r, indent=1, default=float))
    print(r["decision"], "| D =", round(r["D"], 4), "| p =", r["p_signflip"], *(["| tier1:", r["tier1_supported"]] if args.nulls else []), "->", out)
    return r


if __name__ == "__main__":
    run()
