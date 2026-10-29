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


def amplitude_ratios(pred: np.ndarray, obs: np.ndarray, mask: np.ndarray) -> dict:
    """Amplitude calibration of a prediction over the evaluated perturbations (genes outside each perturbation's mask are zeroed).
    mean_delta_ratio = ||mean over perts of predicted delta|| / ||mean over perts of observed delta||  (the shared-response amplitude);
    pert_norm_ratio  = mean over perts of ||predicted delta|| / mean over perts of ||observed delta||.   > 1 = overshoot, < 1 = undershoot."""
    P, O = np.where(mask, np.asarray(pred, float), 0.0), np.where(mask, np.asarray(obs, float), 0.0)
    return {"mean_delta_ratio": float(np.linalg.norm(P.mean(0)) / np.linalg.norm(O.mean(0))),
            "pert_norm_ratio": float(np.linalg.norm(P, axis=1).mean() / np.linalg.norm(O, axis=1).mean()),
            "collapse_pred": collapse_index(P, mask), "collapse_obs": collapse_index(O, mask)}


def collapse_index(X: np.ndarray, mask: np.ndarray) -> float:
    """collapse(X) = ||mean_p X_p|| / mean_p ||X_p|| over the evaluated perturbations (masked genes only). 1 = all perturbations point the same way (fully collapsed onto one shared
    direction); ~1/sqrt(n) = directions unrelated to each other. Identity: mean_delta_ratio = collapse(pred)/collapse(obs) * pert_norm_ratio."""
    X = np.where(mask, np.asarray(X, float), 0.0)
    return float(np.linalg.norm(X.mean(0)) / np.linalg.norm(X, axis=1).mean())
