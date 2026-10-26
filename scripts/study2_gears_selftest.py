"""Self-test of scripts/study2_gears_train.py against GEARS.train itself, on a handful of real batches (CPU, fixed seeds). Run in env hic-perturb-gears:
   python -u scripts/study2_gears_selftest.py [--batches 6]
T1  new loop (eval_train=True) == GEARS.train(epochs=2): final model and best model identical
T2  crash + resume: 1 epoch -> checkpoint -> a FRESH differently-initialised object resumes to epoch 2 == uninterrupted 2-epoch run (bit-identical on CPU)
T3  early stopping bookkeeping is consistent with the history (eval_train=False path), and a finished checkpoint resumes without training
No research result is produced; no Hi-C; no RPE1 data."""
import argparse
import itertools
import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gears import GEARS, PertData                                                      # noqa: E402
from study2_gears_train import train_gears                                             # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--batches", type=int, default=6); ap.add_argument("--split", default="seed1"); a = ap.parse_args()
pert_data = PertData(Path("data/gears")); pert_data.load(data_path="data/gears/replogle_k562_essential")
set2c = json.load(open(f"configs/splits/gears_k562_{a.split}.json"))
valid = set(pert_data.dataset_processed.keys())
set2c = {k: [c for c in v if c == "ctrl" or c in valid] for k, v in set2c.items()}
pert_data.set2conditions = set2c; pert_data.split = "custom"; pert_data.subgroup = None; pert_data.seed = 1; pert_data.train_gene_set_size = 0.75
pert_data.get_dataloader(batch_size=32, test_batch_size=128)


class Limited:
    def __init__(self, loader, n): self.loader, self.n = loader, n
    def __iter__(self): return itertools.islice(iter(self.loader), self.n)
    def __len__(self): return min(self.n, len(self.loader))


def make(seed):
    np.random.seed(seed); torch.manual_seed(seed)
    g = GEARS(pert_data, device="cpu"); g.model_initialize(hidden_size=64)
    g.dataloader = {"train_loader": Limited(pert_data.dataloader["train_loader"], a.batches), "val_loader": Limited(pert_data.dataloader["val_loader"], 3)}   # no test loader: GEARS.train returns after validation
    return g


def maxdiff(m1, m2):
    s1, s2 = m1.state_dict(), m2.state_dict()
    return max(float((s1[k].float().cpu() - s2[k].float().cpu()).abs().max()) for k in s1)


quiet = lambda *x: None
res = {}
tmp = Path(tempfile.mkdtemp(prefix="gears_selftest_"))
try:
    # T1: GEARS.train vs ours
    g1 = make(7); t = time.time(); g1.train(epochs=2); res["gears_train_seconds"] = time.time() - t
    g2 = make(7); s2 = train_gears(g2, epochs=2, patience=99, monitor="mse", eval_train=True, ckpt_dir=tmp / "t1", meta={"seed": 7}, log=quiet)
    res["T1_final_model_max_abs_diff"] = maxdiff(g1.model, g2.model); res["T1_best_model_max_abs_diff"] = maxdiff(g1.best_model, g2.best_model)
    # T2: crash + resume
    g3 = make(7); train_gears(g3, epochs=1, patience=99, monitor="mse", eval_train=True, ckpt_dir=tmp / "t2", meta={"seed": 7}, log=quiet)
    g4 = make(123)                                                                       # different initial weights and RNG: everything must come from the checkpoint
    s4 = train_gears(g4, epochs=2, patience=99, monitor="mse", eval_train=True, ckpt_dir=tmp / "t2", meta={"seed": 7}, log=quiet)
    res["T2_resumed_vs_uninterrupted_final_max_abs_diff"] = maxdiff(g2.model, g4.model); res["T2_resumed_vs_uninterrupted_best_max_abs_diff"] = maxdiff(g2.best_model, g4.best_model)
    res["T2_epochs_in_history_after_resume"] = len(s4["history"])
    # T3: early-stopping bookkeeping, eval_train=False
    g5 = make(11); s5 = train_gears(g5, epochs=8, patience=1, monitor="mse", eval_train=False, ckpt_dir=tmp / "t3", meta={"seed": 11}, log=quiet)
    h = s5["history"]; vals = [r["val_mse"] for r in h]
    best, bad, stop_at = float("inf"), 0, None
    for i, v in enumerate(vals):
        if v < best: best, bad = v, 0
        else:
            bad += 1
            if bad >= 1: stop_at = i + 1; break
    res["T3_epochs_run"] = s5["epochs_run"]; res["T3_expected_stop_epoch_from_history"] = stop_at; res["T3_stopped_early_flag"] = s5["stopped_early"]
    g6 = make(5); s6 = train_gears(g6, epochs=8, patience=1, monitor="mse", eval_train=False, ckpt_dir=tmp / "t3", meta={"seed": 11}, log=quiet)   # same checkpoint, finished -> no training
    res["T3_finished_checkpoint_resumes_without_training"] = bool(s6["epochs_run"] == s5["epochs_run"] and maxdiff(g5.best_model, g6.best_model) == 0.0)
    ok = (res["T1_final_model_max_abs_diff"] < 1e-5 and res["T1_best_model_max_abs_diff"] < 1e-5 and res["T2_resumed_vs_uninterrupted_final_max_abs_diff"] < 1e-5
          and res["T2_resumed_vs_uninterrupted_best_max_abs_diff"] < 1e-5 and res["T2_epochs_in_history_after_resume"] == 2
          and (stop_at is None or res["T3_epochs_run"] == stop_at) and res["T3_finished_checkpoint_resumes_without_training"])
    res["ALL_PASS"] = bool(ok)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
print(json.dumps(res, indent=1))
sys.exit(0 if res.get("ALL_PASS") else 1)
