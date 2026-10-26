"""Study 2: time the pieces of one GEARS epoch that are INFERENCE-only (no training): validation pass, a slice of the train-set pass, checkpoint write. Run in env hic-perturb-gears.
python -u scripts/study2_gears_time_eval.py --device mps"""
import argparse
import itertools
import json
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1])); sys.path.insert(0, str(Path(__file__).resolve().parent))
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
from gears import GEARS, PertData                                                      # noqa: E402
import gears.gears as gg                                                               # noqa: E402
from src.study2_train_utils import atomic_torch_save                                   # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--device", default="mps"); ap.add_argument("--train-batches", type=int, default=400); a = ap.parse_args()
pert_data = PertData(Path("data/gears")); pert_data.load(data_path="data/gears/replogle_k562_essential")
set2c = json.load(open("configs/splits/gears_k562_seed1.json")); valid = set(pert_data.dataset_processed.keys())
pert_data.set2conditions = {k: [c for c in v if c == "ctrl" or c in valid] for k, v in set2c.items()}; pert_data.split = "custom"; pert_data.subgroup = None; pert_data.seed = 1; pert_data.train_gene_set_size = 0.75
pert_data.get_dataloader(batch_size=32, test_batch_size=128)
np.random.seed(1); torch.manual_seed(1)
g = GEARS(pert_data, device=a.device); g.model_initialize(hidden_size=64); g.model = g.model.to(a.device)
res = {"device": a.device, "n_val_batches": len(pert_data.dataloader["val_loader"]), "n_train_batches": len(pert_data.dataloader["train_loader"])}
torch.mps.synchronize() if a.device == "mps" else None
t = time.time(); vr = gg.evaluate(pert_data.dataloader["val_loader"], g.model, g.config["uncertainty"], a.device); res["val_evaluate_s"] = time.time() - t
t = time.time(); gg.compute_metrics(vr); res["val_compute_metrics_s"] = time.time() - t
class Lim:
    def __init__(self, l, n): self.l, self.n = l, n
    def __iter__(self): return itertools.islice(iter(self.l), self.n)
    def __len__(self): return min(self.n, len(self.l))
t = time.time(); tr = gg.evaluate(Lim(pert_data.dataloader["train_loader"], a.train_batches), g.model, g.config["uncertainty"], a.device); dt = time.time() - t
res[f"train_evaluate_{a.train_batches}_batches_s"] = dt
res["train_evaluate_full_set_extrapolated_s"] = dt * len(pert_data.dataloader["train_loader"]) / a.train_batches
t = time.time(); gg.compute_metrics(tr); res["train_compute_metrics_slice_s"] = time.time() - t
import torch.optim as optim
opt = optim.Adam(g.model.parameters()); d = Path(tempfile.mkdtemp())
t = time.time(); atomic_torch_save({"model": g.model.state_dict(), "optimizer": opt.state_dict(), "best": g.model.state_dict(), "meta": {}}, d / "last.pt"); res["checkpoint_write_s"] = time.time() - t
res["checkpoint_size_MB"] = os.path.getsize(d / "last.pt") / 1e6
print(json.dumps(res, indent=1))
