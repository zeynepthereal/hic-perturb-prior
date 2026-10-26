"""STUDY 2 / Phase 1 (calibration). Run in the SEPARATE env:  mamba run -n hic-perturb-gears python -u scripts/study2_gears_run.py --label gears_seed1
Trains GEARS (cell-gears 0.1.2) on GEARS-processed K562 with the STANDARD GO perturbation graph and the default hyper-parameters of the Ahlmann-Eltze et al. benchmark script
(run_gears.py: epochs=20 -> here max 12 + early stopping with checkpoints, hidden_size=64, batch 32 / test batch 128, custom split = our GEARS split JSON, np.random.seed(seed)). Writes predictions for the TEST perturbations
(no Hi-C of any kind; no RPE1 data).  Instrumented: wall time per phase, RSS sampled every 20 s, optional time-limited CPU-vs-MPS benchmark through GEARS' own training loop."""
import argparse
import itertools
import json
import os
import resource
import sys
import threading
import time
from pathlib import Path

import numpy as np
import psutil
import torch

ap = argparse.ArgumentParser()
ap.add_argument("--label", default="gears_seed1"); ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--epochs", type=int, default=12, help="MAX epochs (early stopping may end earlier)"); ap.add_argument("--hidden", type=int, default=64)
ap.add_argument("--patience", type=int, default=3); ap.add_argument("--monitor", default="mse", choices=["mse", "mse_de"], help="validation metric for early stopping (GEARS keeps its best model by top-20-DE MSE regardless)")
ap.add_argument("--eval-train", action="store_true", help="also evaluate the full training set every epoch (GEARS default; only for printing)")
ap.add_argument("--ckpt-root", default="data/study2_ckpt"); ap.add_argument("--no-resume", action="store_true")
ap.add_argument("--batch", type=int, default=32); ap.add_argument("--test-batch", type=int, default=128)
ap.add_argument("--device", default="mps", help="cpu | mps | auto (benchmark both, use the faster that works)")
ap.add_argument("--bench-steps", type=int, default=40); ap.add_argument("--bench-only", action="store_true")
ap.add_argument("--out", default=None)
a = ap.parse_args()
out = Path(a.out or f"results/study2/gears_k562_{a.label}_s{a.seed}"); out.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
T0 = time.time(); proc = psutil.Process(); peak = [0.0]; log_f = open(out / "run.log", "a")


def log(*m):
    s = f"[{time.time()-T0:8.1f}s rss {proc.memory_info().rss/2**30:5.1f} GiB] " + " ".join(map(str, m)); print(s, flush=True); log_f.write(s + "\n"); log_f.flush()


def sampler():
    while True:
        r = proc.memory_info().rss / 2**30; peak[0] = max(peak[0], r)
        with open(out / "memory_samples.tsv", "a") as f: f.write(f"{time.time()-T0:.0f}\t{r:.2f}\n")
        time.sleep(20)


threading.Thread(target=sampler, daemon=True).start()


class _TS:                                                                          # timestamp every stdout line (GEARS prints step counters without times)
    def __init__(self, f): self.f, self.buf = f, ""
    def write(self, x):
        self.buf += x
        while "\n" in self.buf:
            line, self.buf = self.buf.split("\n", 1); self.f.write(f"[t={time.time()-T0:8.1f}s] {line}\n")
    def flush(self): self.f.flush()


sys.stdout = _TS(sys.__stdout__)
from gears import GEARS, PertData                                                     # noqa: E402

timing = {"label": a.label, "args": vars(a), "torch": torch.__version__, "cuda": torch.cuda.is_available(), "mps_available": torch.backends.mps.is_available()}
np.random.seed(a.seed); torch.manual_seed(a.seed)
t = time.time(); pert_data = PertData(Path("data/gears")); pert_data.load(data_path="data/gears/replogle_k562_essential"); timing["t_load_pertdata_s"] = time.time() - t
log("PertData loaded (cell graphs built or read from cache):", f"{timing['t_load_pertdata_s']:.0f}s")
set2c = json.load(open(f"configs/splits/gears_k562_{a.label.replace('gears_', '')}.json"))
valid = set(pert_data.dataset_processed.keys())      # conditions GEARS can use: those in its GO graph (the paper made its splits AFTER this filtering; ours were made on the raw file)
dropped = {k: [c for c in v if c != "ctrl" and c not in valid] for k, v in set2c.items()}
set2c = {k: [c for c in v if c == "ctrl" or c in valid] for k, v in set2c.items()}
timing["conditions_dropped_not_in_go_graph"] = dropped; log("conditions dropped (not in GEARS GO graph):", {k: v for k, v in dropped.items() if v})
pert_data.set2conditions = set2c; pert_data.split = "custom"; pert_data.subgroup = None; pert_data.seed = 1; pert_data.train_gene_set_size = 0.75
t = time.time(); pert_data.get_dataloader(batch_size=a.batch, test_batch_size=a.test_batch); timing["t_dataloader_s"] = time.time() - t
log("dataloaders ready", f"{timing['t_dataloader_s']:.0f}s", {k: len(v) for k, v in pert_data.dataloader.items()})


class Limited:                                                              # first n batches of a loader (benchmark only)
    def __init__(self, loader, n): self.loader, self.n = loader, n
    def __iter__(self): return itertools.islice(iter(self.loader), self.n)
    def __len__(self): return min(self.n, len(self.loader))


def make(device):
    g = GEARS(pert_data, device=device); g.model_initialize(hidden_size=a.hidden); return g


device = a.device
if device == "auto":
    bench = {}
    for d in ("cpu", "mps"):
        try:
            g = make(d); g.dataloader = {k: Limited(v, a.bench_steps) for k, v in pert_data.dataloader.items()}
            t = time.time(); g.train(epochs=1); dt = time.time() - t
            bench[d] = {"ok": True, "seconds_for_%d_train_steps_plus_limited_eval" % a.bench_steps: dt}; log("benchmark", d, f"{dt:.1f}s")
            del g
        except Exception as e:                                                  # noqa: BLE001
            bench[d] = {"ok": False, "error": repr(e)[:300]}; log("benchmark", d, "FAILED", repr(e)[:200])
    timing["device_benchmark"] = bench
    ok = {d: v for d, v in bench.items() if v["ok"]}
    device = min(ok, key=lambda d: list(ok[d].values())[1]) if ok else "cpu"
    log("device chosen:", device)
timing["device_used"] = device
json.dump(timing, open(out / "timing.json", "w"), indent=1)
if a.bench_only:
    sys.exit(0)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from study2_gears_train import train_gears                                            # noqa: E402
g = make(device)
ckpt = Path(a.ckpt_root) / f"{a.label}_s{a.seed}"
t = time.time()
summary = train_gears(g, epochs=a.epochs, patience=a.patience, monitor=a.monitor, ckpt_dir=ckpt, resume=not a.no_resume, eval_train=a.eval_train, log=log,
                      meta={"label": a.label, "seed": a.seed, "hidden": a.hidden, "batch": a.batch, "device": device, "max_epochs": a.epochs})
timing["t_train_s_this_process"] = time.time() - t; timing["training_summary"] = summary
log("training finished", f"{timing['t_train_s_this_process']:.0f}s this process;", {k: v for k, v in summary.items() if k != "history"})
test_conds = [c for c in set2c["test"] if c != "ctrl"]
genes = [c.split("+") for c in test_conds]; genes = [[x for x in gl if x != "ctrl"] for gl in genes]
ok_idx = [i for i, gl in enumerate(genes) if all(x in g.pert_list for x in gl)]
t = time.time(); preds = g.predict([genes[i] for i in ok_idx]); timing["t_predict_s"] = time.time() - t
keys = ["_".join(genes[i]) for i in ok_idx]
P = np.stack([np.asarray(preds[k]) for k in keys]).astype(np.float32)
np.savez_compressed(out / "predictions_test.npz", conditions=np.array([test_conds[i] for i in ok_idx]), pred_expression=P,
                    gene_names=np.array(pert_data.adata.var["gene_name"].values, dtype=str),
                    n_test_total=len(test_conds), n_test_in_go_graph=len(ok_idx))
timing.update({"n_test": len(test_conds), "n_test_predicted": len(ok_idx), "t_total_s": time.time() - T0,
               "peak_rss_sampled_GiB": peak[0], "ru_maxrss_GiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30})
json.dump(timing, open(out / "timing.json", "w"), indent=1)
log("DONE", json.dumps({k: (round(v, 1) if isinstance(v, float) else v) for k, v in timing.items() if k.startswith(("t_", "peak", "ru_", "n_", "device_used"))}))
