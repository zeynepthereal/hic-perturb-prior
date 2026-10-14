"""scripts/study2_gears_eval.score on synthetic data: a 'GEARS' that predicts the mean baseline must reproduce the mean baseline's metrics exactly; a perfect predictor has L2 = 0."""
import numpy as np

from scripts.study2_gears_eval import score

from .test_synthetic import _toy_ds


def test_score_reproduces_mean_baseline_and_perfect_predictor():
    ds = _toy_ds(n_pert=150, n_gene=80, n_missing=20, seed=5)
    perm = np.random.RandomState(1).permutation(150)
    split = {"train": np.sort(perm[:100]), "val": np.sort(perm[100:110]), "test": np.sort(perm[110:])}
    te = split["test"]
    gd = {int(i): ds.delta[split["train"]].mean(0) for i in te[::2]}                  # "GEARS" predicts the mean baseline for half of the test perturbations
    M = score(ds, split, gd, n_top_expr=30, n_top_de=5)
    for (s, t), g in M[M.model.isin(["GEARS", "mean"])].groupby(["set", "target"]):
        a = g[g.model == "GEARS"].set_index(["pert_i", "subset"]).l2.sort_index(); b = g[g.model == "mean"].set_index(["pert_i", "subset"]).l2.sort_index()
        assert len(a) > 0 and np.allclose(a.values, b.loc[a.index].values)
    gd2 = {int(i): ds.delta[i] for i in te}
    M2 = score(ds, split, gd2, n_top_expr=30, n_top_de=5)
    assert np.allclose(M2[M2.model == "GEARS"].l2.values, 0) and set(M2.set) == {"A: GEARS & PCA predictable", "B: GEARS predictable"}
