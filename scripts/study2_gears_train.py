"""Checkpointed, early-stopping replica of `GEARS.train` (cell-gears 0.1.2). Imported by scripts/study2_gears_run.py; run in the env `hic-perturb-gears`.

IDENTICAL to GEARS.train: Adam(lr=1e-3, weight_decay=5e-4), StepLR(step_size=1, gamma=0.5) stepped once per epoch, loss_fct / uncertainty_loss_fct, gradient value clipping 1.0,
best model = the epoch with the lowest VALIDATION TOP-20-DE MSE (strict '<'), validation computed with gears.inference.evaluate / compute_metrics.
ADDED (Study 2):
  * per-epoch atomic checkpoint `last.pt` (model, optimizer, scheduler, best model, early-stopper, history, all RNG states) and automatic resume;
  * early stopping: stop after `patience` epochs without a strict improvement of the monitored validation metric (`monitor` = 'mse' overall or 'mse_de' top-20 DE); max `epochs`;
  * both validation metrics are logged every epoch (history.json);
  * the full-train-set evaluation that GEARS runs every epoch purely for printing is skipped unless eval_train=True (it consumes DataLoader RNG, so with eval_train=False the data order
    after epoch 1 differs from GEARS' stream; training dynamics and model selection are otherwise unchanged), and the final test-set analysis of GEARS.train is skipped
    (predictions are produced separately with GEARS.predict)."""
import random
import sys
import time
from copy import deepcopy
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.optim.lr_scheduler import StepLR

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.study2_train_utils import EarlyStopper, atomic_json_save, atomic_torch_save  # noqa: E402


def _rng(device):
    s = {"torch": torch.get_rng_state(), "numpy": np.random.get_state(), "python": random.getstate()}
    if device == "mps":
        try:
            s["mps"] = torch.mps.get_rng_state()
        except Exception:                                                              # noqa: BLE001
            pass
    return s


def _set_rng(s, device):
    torch.set_rng_state(s["torch"].cpu()); np.random.set_state(s["numpy"]); random.setstate(s["python"])
    if "mps" in s and device == "mps":
        torch.mps.set_rng_state(s["mps"].cpu())


def _cpu_state(model):
    return {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}


def train_gears(g, *, epochs=12, patience=3, monitor="mse", lr=1e-3, weight_decay=5e-4, ckpt_dir=None, meta=None, eval_train=False, log=print, resume=True, step_log=50):
    """Train GEARS object `g` (after model_initialize). Sets g.best_model. Returns a summary dict."""
    import gears.gears as gg                                                          # evaluate, compute_metrics, loss_fct, uncertainty_loss_fct
    assert monitor in ("mse", "mse_de")
    device = g.device
    train_loader, val_loader = g.dataloader["train_loader"], g.dataloader["val_loader"]
    g.model = g.model.to(device)
    optimizer = optim.Adam(g.model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = StepLR(optimizer, step_size=1, gamma=0.5)
    stopper = EarlyStopper(patience)
    best_state, best_de, best_epoch = _cpu_state(g.model), float("inf"), 0               # GEARS starts from best_model = deepcopy(untrained model)
    history, start, stopped_early, finished = [], 0, False, False
    meta = dict(meta or {}, patience=patience, monitor=monitor, lr=lr, weight_decay=weight_decay, eval_train=eval_train)
    ck = Path(ckpt_dir) / "last.pt" if ckpt_dir else None
    if ck:
        ck.parent.mkdir(parents=True, exist_ok=True)
    if resume and ck is not None and ck.exists():
        c = torch.load(ck, map_location="cpu", weights_only=False)
        if c["meta"] != meta:
            raise RuntimeError(f"checkpoint {ck} was written with different settings: {c['meta']} != {meta}")
        g.model.load_state_dict(c["model"]); optimizer.load_state_dict(c["optimizer"]); scheduler.load_state_dict(c["scheduler"])
        stopper.load_state_dict(c["stopper"]); best_state, best_de, best_epoch = c["best_state"], c["best_de"], c["best_epoch"]
        history, start, stopped_early, finished = c["history"], c["epoch"], c["stopped_early"], c["finished"]
        _set_rng(c["rng"], device)
        log(f"RESUMED from {ck}: {start} epochs already done, finished={finished}, best (top-20-DE) epoch {best_epoch}")

    for epoch in range(start, epochs):
        if finished:
            break
        t0 = time.time(); g.model.train()
        for step, batch in enumerate(train_loader):
            batch.to(device); optimizer.zero_grad(); y = batch.y
            if g.config["uncertainty"]:
                pred, logvar = g.model(batch)
                loss = gg.uncertainty_loss_fct(pred, logvar, y, batch.pert, reg=g.config["uncertainty_reg"], ctrl=g.ctrl_expression, dict_filter=g.dict_filter,
                                               direction_lambda=g.config["direction_lambda"])
            else:
                pred = g.model(batch)
                loss = gg.loss_fct(pred, y, batch.pert, ctrl=g.ctrl_expression, dict_filter=g.dict_filter, direction_lambda=g.config["direction_lambda"])
            loss.backward(); nn.utils.clip_grad_value_(g.model.parameters(), clip_value=1.0); optimizer.step()
            if step % step_log == 0:
                log(f"Epoch {epoch + 1} Step {step + 1} Train Loss: {loss.item():.4f}")
        scheduler.step()
        if eval_train:                                                                  # GEARS evaluates the training set first (this consumes DataLoader RNG)
            gg.compute_metrics(gg.evaluate(train_loader, g.model, g.config["uncertainty"], device))
        val_metrics, _ = gg.compute_metrics(gg.evaluate(val_loader, g.model, g.config["uncertainty"], device))
        if val_metrics["mse_de"] < best_de:                                             # GEARS' own best-model rule
            best_de, best_epoch, best_state = float(val_metrics["mse_de"]), epoch + 1, _cpu_state(g.model)
        stop = stopper.update(float(val_metrics[monitor]), epoch + 1)
        rec = {"epoch": epoch + 1, "val_mse": float(val_metrics["mse"]), "val_mse_de": float(val_metrics["mse_de"]), "val_pearson_de": float(val_metrics.get("pearson_de", float("nan"))),
               "seconds": time.time() - t0, "lr_next": optimizer.param_groups[0]["lr"], "monitor_best_epoch": stopper.best_epoch, "bad_epochs": stopper.bad, "best_de_epoch": best_epoch}
        history.append(rec)
        log(f"Epoch {epoch + 1}: val overall MSE {rec['val_mse']:.4f} | val top-20-DE MSE {rec['val_mse_de']:.4f} | monitor={monitor} best epoch {stopper.best_epoch}, "
            f"bad epochs {stopper.bad}/{patience} | {rec['seconds'] / 60:.1f} min")
        stopped_early = bool(stop and epoch + 1 < epochs); finished = bool(stop)                  # 'finished' = early-stopped; the epoch cap is judged from `epochs` at (re)start
        if ck is not None:
            atomic_torch_save({"meta": meta, "epoch": epoch + 1, "model": _cpu_state(g.model), "optimizer": optimizer.state_dict(), "scheduler": scheduler.state_dict(),
                               "stopper": stopper.state_dict(), "best_state": best_state, "best_de": best_de, "best_epoch": best_epoch, "history": history,
                               "stopped_early": stopped_early, "finished": finished, "rng": _rng(device)}, ck)
            atomic_json_save({"history": history, "finished": finished, "stopped_early": stopped_early, "best_top20DE_epoch": best_epoch}, Path(ckpt_dir) / "history.json")
        if finished:
            break
    bm = deepcopy(g.model); bm.load_state_dict(best_state); g.best_model = bm.to(device)
    return {"history": history, "epochs_run": len(history), "best_top20DE_epoch": best_epoch, "monitor": monitor, "monitor_best_epoch": stopper.best_epoch, "stopped_early": stopped_early,
            "max_epochs": epochs, "patience": patience}
