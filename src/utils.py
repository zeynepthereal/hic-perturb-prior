"""Seeding and run logging (config + git hash) shared by all entry points."""
from __future__ import annotations

import json
import os
import random
import subprocess
import time
from pathlib import Path

import numpy as np
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)


def load_config(path: str | Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def _git(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=REPO_ROOT, stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return ""


def git_state() -> dict:
    """Commit hash + dirty flag ('NO_COMMITS' if the repo has no commit yet)."""
    return {
        "commit": _git("rev-parse", "HEAD") or "NO_COMMITS",
        "dirty": bool(_git("status", "--porcelain")),
    }


def make_run_dir(cfg: dict, name: str) -> Path:
    """Create results/<timestamp>_<name>/ and write config.yaml + git.json into it."""
    run_dir = REPO_ROOT / "results" / f"{time.strftime('%Y%m%d-%H%M%S')}_{name}"
    run_dir.mkdir(parents=True, exist_ok=False)
    with open(run_dir / "config.yaml", "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    with open(run_dir / "git.json", "w") as f:
        json.dump(git_state(), f, indent=2)
    return run_dir
