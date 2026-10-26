"""Dependency-free helpers for the Study 2 GEARS training loop (importable from both conda envs): early stopping and atomic checkpoint writes."""
from __future__ import annotations

import json
import os


class EarlyStopper:
    """Stop when the monitored value (lower = better) has not improved for `patience` consecutive epochs. Improvement = strictly lower than the best so far
    (the same strict rule GEARS uses to keep its best model). `update` returns True when training should stop."""

    def __init__(self, patience: int):
        self.patience, self.best, self.best_epoch, self.bad = patience, float("inf"), 0, 0

    def update(self, value: float, epoch: int) -> bool:
        if value < self.best:
            self.best, self.best_epoch, self.bad = value, epoch, 0
        else:
            self.bad += 1
        return self.bad >= self.patience

    def state_dict(self) -> dict:
        return dict(patience=self.patience, best=self.best, best_epoch=self.best_epoch, bad=self.bad)

    def load_state_dict(self, d: dict) -> None:
        self.patience, self.best, self.best_epoch, self.bad = d["patience"], d["best"], d["best_epoch"], d["bad"]


def atomic_torch_save(obj, path) -> None:
    """Write to a temporary file, then rename: a crash or power loss mid-write can never corrupt the existing checkpoint."""
    import torch
    tmp = f"{path}.tmp"
    torch.save(obj, tmp)
    os.replace(tmp, path)


def atomic_json_save(obj, path) -> None:
    tmp = f"{path}.tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1, default=float)
    os.replace(tmp, path)
