"""Sanity checks on synthetic data (run: pytest src/test_synthetic.py)."""
import numpy as np

from .baselines import LinearModel, MeanBaseline, pca_gene_embedding
from .evaluate import gene_masks, per_pert_metrics


def test_linear_recovers_planted_model():
    rng = np.random.RandomState(0)
    n_gene, n_pert, K = 200, 150, 5
    G = rng.randn(n_gene, K)
    W = rng.randn(K, K)
    P = G[:n_pert] - G[:100].mean(0)     # centred on train perts so empirical row-mean b equals planted b
    b = rng.randn(n_gene, 1)
    Y = (G @ W @ P.T + b).T
    m = LinearModel(lam=1e-8).fit(Y[:100], G, P[:100])
    pred = m.predict(P[100:])
    assert np.abs(pred - Y[100:]).max() < 1e-3


def test_mean_baseline_and_masks():
    Y = np.random.RandomState(1).randn(10, 50)
    assert np.allclose(MeanBaseline().fit(Y).predict(3)[0], Y.mean(0))
    tgt = np.arange(10)
    m = gene_masks(Y, np.arange(50.0), tgt, False, n_de=5, n_expr=7)
    assert m["top20_de"].sum(1).tolist() == [5] * 10 and m["top1000_expr"].sum(1)[0] == 7
    df = per_pert_metrics(Y, Y, m)
    assert np.allclose(df.l2, 0) and np.allclose(df.pearson, 1)


def test_hvg_replica_matches_scanpy():
    import anndata as ad
    import scanpy as sc
    from .data import seurat_hvg_from_moments
    rng = np.random.RandomState(0)
    counts = rng.negative_binomial(2, 0.05 + 0.5 * rng.rand(3000)[None, :], size=(400, 3000)).astype(float)
    norm = counts / np.maximum(counts.sum(1, keepdims=True), 1) * 1e4
    a = ad.AnnData(np.log1p(norm).astype(np.float64))
    sc.pp.highly_variable_genes(a, n_top_genes=500, flavor="seurat")
    mine = seurat_hvg_from_moments(norm.mean(0), norm.var(0, ddof=1), 500)
    assert (mine == a.var["highly_variable"].values).all()


def test_target_exclusion():
    Y = np.random.RandomState(2).randn(10, 50)
    tgt = np.arange(10)
    Y[np.arange(10), tgt] = 100.0          # target gene is the strongest DE gene -> must be dropped when excluded
    inc = gene_masks(Y, np.arange(50.0), tgt, False, n_de=5, n_expr=7)
    exc = gene_masks(Y, np.arange(50.0), tgt, True, n_de=5, n_expr=7)
    assert inc["top20_de"][np.arange(10), tgt].all() and not exc["top20_de"][np.arange(10), tgt].any()
    assert (exc["top20_de"].sum(1) == 5).all() and (exc["top1000_expr"].sum(1) == 7).all()
    assert (exc["all"].sum(1) == 49).all()
