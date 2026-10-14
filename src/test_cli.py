"""scripts/run_primary_analysis.py must reproduce, exactly, what the one-liner reported after the freeze produces (synthetic metrics file only)."""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.run_primary_analysis import run

ROOT = Path(__file__).resolve().parents[1]
ONE_LINER = ("import sys, json, pandas as pd; from src.stats import primary_analysis; run = sys.argv[1]; pp = pd.read_csv(run + '/per_perturbation_metrics.csv'); "
             "ctrl = [f'linear[b500_lin_randtrans_s{k}] (paper hp K=10, lam=0.1)' for k in range(10)]; "
             "r = primary_analysis(pp, 'linear[b500_hic_trans] (paper hp K=10, lam=0.1)', ctrl); json.dump(r, open(run + '/primary_analysis.json', 'w'), indent=1, default=float)")


def _synthetic(path, names_hic, names_ctrl, n_split=3, n_pert=80, seed=0):
    rng = np.random.RandomState(seed); rows = []
    for s in range(1, n_split + 1):
        for i in rng.choice(400, n_pert, replace=False):
            base = 4.2 + 0.3 * rng.randn()
            for m in [names_hic] + names_ctrl:
                rows.append((i, f"gears_seed{s}", m, "common", "excluded", "top1000_expr", base + (0.0 if m == names_hic else 0.03) + 0.02 * rng.randn(), 0.3))
    path.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["pert_i", "split", "model", "eval_set", "target", "subset", "l2", "pearson"]).to_csv(path / "per_perturbation_metrics.csv", index=False)


def test_cli_reproduces_the_one_liner_exactly(tmp_path):
    hic = "linear[b500_hic_trans] (paper hp K=10, lam=0.1)"
    _synthetic(tmp_path / "run", hic, [f"linear[b500_lin_randtrans_s{k}] (paper hp K=10, lam=0.1)" for k in range(10)])
    subprocess.run([sys.executable, "-c", ONE_LINER, str(tmp_path / "run")], check=True, cwd=ROOT, capture_output=True)      # default 100,000 permutations
    ref = json.load(open(tmp_path / "run" / "primary_analysis.json"))
    run([str(tmp_path / "run"), "--out", str(tmp_path / "cli.json")])                                                       # same defaults
    got = json.load(open(tmp_path / "cli.json"))
    assert {k: v for k, v in got.items() if k != "cli"} == ref                    # every number identical (deterministic seeds)
    assert got["cli"]["hic_model"] == hic and len(got["cli"]["control_models"]) == 10


def test_cli_fixed_k_resid_names_and_tier_rules(tmp_path):
    hic = "linear[b500_hic_trans:resid] (fixed K=50, lam=0.1)"
    _synthetic(tmp_path / "s3", hic, [f"linear[b500_lin_randtrans_s{k}:resid] (fixed K=50, lam=0.1)" for k in range(10)])
    r3 = run([str(tmp_path / "s3"), "--hic", "b500_hic_trans:resid", "--control-template", "b500_lin_randtrans_s{k}:resid", "--kind", "fixed", "--K", "50", "--n-perm", "2000"])
    assert r3["decision"] == "SIGNIFICANT: Hi-C better" and r3["cli"]["hic_model"] == hic
    # tier rules: transcription of section 3 from nulls.json (+ the S3 json)
    _synthetic(tmp_path / "p", "linear[b500_hic_trans] (paper hp K=10, lam=0.1)", [f"linear[b500_lin_randtrans_s{k}] (paper hp K=10, lam=0.1)" for k in range(10)])
    def nulls(support, D):
        (tmp_path / "nulls.json").write_text(json.dumps({"10": {"K": 10, "D_obs": D, "support": support, "p": {}, "thresholds": {}, "null_D": {"permutation": [0.0]}}}))
    nulls(True, 0.03)
    r = run([str(tmp_path / "p"), "--n-perm", "2000", "--nulls", str(tmp_path / "nulls.json"), "--tier2-from", str(tmp_path / "s3" / "primary_analysis.json"), "--out", str(tmp_path / "o1.json")])
    assert r["tier1_supported"] is True and r["tier2_supported"] is True and "null_D" not in r["null_condition"]
    nulls(False, 0.03)
    r = run([str(tmp_path / "p"), "--n-perm", "2000", "--nulls", str(tmp_path / "nulls.json"), "--tier2-from", str(tmp_path / "s3" / "primary_analysis.json"), "--out", str(tmp_path / "o2.json")])
    assert r["tier1_supported"] is False and r["tier2_supported"] is False        # Tier 2 only after Tier 1
