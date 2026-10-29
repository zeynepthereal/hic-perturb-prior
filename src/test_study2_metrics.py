import numpy as np

from .study2_metrics import centred_pearson


def test_centred_pearson_removes_shared_response():
    rng = np.random.RandomState(0)
    n, g = 40, 200
    shared = rng.randn(g) * 3.0                                   # response shared by all perturbations
    specific = rng.randn(n, g)
    obs = shared + specific
    mask = np.ones((n, g), bool)
    # a predictor that knows only the shared response: raw Pearson is high, centred Pearson is undefined (no perturbation-specific information)
    const = np.tile(shared, (n, 1))
    raw = np.array([np.corrcoef(const[i], obs[i])[0, 1] for i in range(n)])
    assert raw.mean() > 0.8 and np.isnan(centred_pearson(const, obs, mask)).all()
    # perfect predictor + an arbitrary extra shared offset in the prediction: centred Pearson = 1
    pred = obs + rng.randn(g) * 5.0
    assert np.allclose(centred_pearson(pred, obs, mask), 1.0)
    # a predictor with half the specific signal: clearly positive, below 1
    half = shared + 0.5 * specific + 0.5 * rng.randn(n, g)
    c = centred_pearson(half, obs, mask)
    assert 0.1 < np.nanmean(c) < 0.9
    # random predictor: ~0
    assert abs(np.nanmean(centred_pearson(rng.randn(n, g), obs, mask))) < 0.1
    # masks are honoured; n < 2 -> NaN
    m2 = mask.copy(); m2[:, 100:] = False
    assert np.allclose(centred_pearson(pred, obs, m2), 1.0) and np.isnan(centred_pearson(pred[:1], obs[:1], mask[:1])).all()


def test_amplitude_ratios():
    from .study2_metrics import amplitude_ratios
    rng = np.random.RandomState(1)
    obs = rng.randn(30, 50) + 2.0
    mask = np.ones_like(obs, bool)
    r = amplitude_ratios(1.5 * obs, obs, mask)
    assert abs(r["mean_delta_ratio"] - 1.5) < 1e-9 and abs(r["pert_norm_ratio"] - 1.5) < 1e-9
    assert abs(amplitude_ratios(0.5 * obs, obs, mask)["mean_delta_ratio"] - 0.5) < 1e-9


def test_collapse_index():
    from .study2_metrics import amplitude_ratios, collapse_index
    rng = np.random.RandomState(2)
    mask = np.ones((60, 80), bool)
    v = rng.randn(80)
    assert abs(collapse_index(np.tile(v, (60, 1)), mask) - 1.0) < 1e-9          # all along one direction
    assert collapse_index(rng.randn(60, 80), mask) < 0.2                          # unrelated directions: ~1/sqrt(n)
    obs = rng.randn(60, 80); pred = np.tile(v, (60, 1)) * 0.5 + 0.2 * rng.randn(60, 80)
    r = amplitude_ratios(pred, obs, mask)
    assert r["collapse_pred"] > 5 * r["collapse_obs"]
    assert abs(r["mean_delta_ratio"] - r["collapse_pred"] / r["collapse_obs"] * r["pert_norm_ratio"]) < 1e-9
