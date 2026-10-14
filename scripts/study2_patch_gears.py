"""Study 2: one-line compatibility patch of the INSTALLED cell-gears 0.1.2 in the separate env `hic-perturb-gears` (scipy 1.15 rejects a pandas Series as sparse-matrix row index).
Semantics-preserving: `X[obs.condition == 'ctrl']` -> `X[(obs.condition == 'ctrl').values]`.   python scripts/study2_patch_gears.py   (idempotent; writes the diff)"""
import difflib
import importlib.util
from pathlib import Path

spec = importlib.util.find_spec("gears"); f = Path(spec.origin).parent / "gears.py"
old = "np.mean(self.adata.X[self.adata.obs.condition == 'ctrl'],"
new = "np.mean(self.adata.X[(self.adata.obs.condition == 'ctrl').values],"
src = f.read_text()
if new in src:
    print("already patched")
else:
    assert src.count(old) == 1
    f.write_text(src.replace(old, new))
    d = "".join(difflib.unified_diff(src.splitlines(True), f.read_text().splitlines(True), "gears/gears.py (original)", "gears/gears.py (patched)"))
    Path("results/study2").mkdir(parents=True, exist_ok=True); Path("results/study2/gears_env_patch.diff").write_text(d); print(d)
