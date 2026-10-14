"""python -m scripts.build_percomp -> data/hic/graphs/<cell>_percomp.npz for the disconnected (cis-only) cells: per-component modes (secondary analysis)."""
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from src.hic.spectral import component_modes

ART = Path("data/hic/graphs")
for cell in ("b500_hic_cis", "b500_lin_cis", "b250_hic_cis", "b250_lin_cis"):
    W = sp.load_npz(ART / f"{cell}.npz")
    cm = component_modes(W, m_max=50)
    meta = np.load(ART / f"{cell}_spectral.npz")
    np.savez_compressed(ART / f"{cell}_percomp.npz", **cm, node_chr=meta["node_chr"], node_pos=meta["node_pos"], res=meta["res"])
    print(cell, "components", len(cm["sizes"]), "sizes", sorted(cm["sizes"].tolist(), reverse=True)[:6], "...", "modes available min", int(cm["avail"].min()))
