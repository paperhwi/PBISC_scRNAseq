import numpy as np

from pbisc_metrics import gene_mean as gm
from pbisc_metrics._utils import depth, thin_to_depth


def test_identical_sources_give_ratio_one(real_a):
    t = gm.compare_gene_means(real_a["X"], real_a["genes"], real_a["X"], real_a["genes"], depth_match=False)
    assert np.allclose(t.ratio_cp10k.dropna(), 1.0)
    c = gm.correlation_summary(t)
    assert c["pearson_log1p_cp10k"] > 0.999999
    b = gm.band_table(t, by="rank", edges=(0, 100, 200, np.inf))
    assert np.allclose(b.ratio, 1.0) and np.allclose(b.det_ratio, 1.0)


def test_thinning_hits_target_depth_and_never_adds_counts(real_a):
    X = real_a["X"]
    Y = thin_to_depth(X, 500.0, seed=0)
    assert abs(np.median(depth(Y)) - 500) / 500 < 0.05
    assert (Y - X).max() <= 0  # thinning only removes
    Z = thin_to_depth(X, 1e9, seed=0)  # target above depth: unchanged
    assert (Z != X).nnz == 0


def test_independent_real_replicate_is_close_and_blunted_tail_is_detected(world, real_a, real_b):
    from pbisc_metrics import synthetic as S
    t = gm.compare_gene_means(real_a["X"], real_a["genes"], real_b["X"], real_b["genes"], depth_match=True)
    assert gm.correlation_summary(t)["pearson_log1p_cp10k"] > 0.97
    cliff = S.sample_cohort(world, n_per_group=14, seed=3, lowexpr_floor=1.0, prefix="c")
    tc = gm.compare_gene_means(real_a["X"], real_a["genes"], cliff["X"], cliff["genes"], depth_match=True)
    bands = gm.band_table(tc, by="rank", edges=(0, 100, 200, 400))
    assert bands.det_ratio.iloc[0] > 0.9 and bands.det_ratio.iloc[-1] < 0.5  # only the low-expression tail loses detection
    assert tc.gen_mean_umi.eq(0).sum() > t.gen_mean_umi.eq(0).sum() + 10


def test_gene_alignment_dedupes_and_reorders(real_a):
    g = real_a["genes"].copy()
    g2 = np.concatenate([g[::-1], g[:3]])  # reversed + 3 duplicated symbols
    X2 = np.hstack([real_a["X"].toarray()[:, ::-1], real_a["X"].toarray()[:, :3]])
    t = gm.compare_gene_means(real_a["X"], g, X2, g2, depth_match=False)
    assert len(t) == len(g) and np.allclose(t.ratio_cp10k.dropna(), 1.0)
