"""Study 2 metrics added on top of the Study 1 evaluator (which is left untouched)."""
from __future__ import annotations

import numpy as np


def centred_pearson(pred: np.ndarray, obs: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Pearson correlation per perturbation AFTER subtracting the mean perturbation effect across the evaluated test perturbations.
    pred, obs: (n_pert, n_gene) deltas (perturbation minus control) of the evaluated perturbations; mask: (n_pert, n_gene) boolean genes used per perturbation.
    For every gene the mean over the n perturbations is removed from BOTH the prediction and the observation (the shared response that a mean-baseline already captures
    is eliminated), then Pearson is computed over the perturbation's masked genes. Constant-over-perturbations predictions (mean baseline, no-change) become all-zero vectors
    -> NaN: they carry no perturbation-specific information. n < 2 -> NaN."""
    pred, obs = np.asarray(pred, float), np.asarray(obs, float)
    out = np.full(pred.shape[0], np.nan)
    if pred.shape[0] < 2:
        return out
    pc, oc = pred - pred.mean(0, keepdims=True), obs - obs.mean(0, keepdims=True)
    for i in range(pred.shape[0]):
        a, b = pc[i, mask[i]], oc[i, mask[i]]
        a, b = a - a.mean(), b - b.mean()
        den = np.sqrt((a * a).sum() * (b * b).sum())
        if den > 1e-10 and np.sqrt((pc[i, mask[i]] ** 2).sum()) > 1e-10:
            out[i] = (a * b).sum() / den
    return out
