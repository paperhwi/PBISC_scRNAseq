import numpy as np

from pbisc_metrics import ilisi as I


def test_lisi_two_perfectly_mixed_and_separated_sources():
    rng = np.random.default_rng(0)
    n = 600
    src = np.array(["a"] * n + ["b"] * n)
    mixed = rng.normal(size=(2 * n, 5))
    sep = mixed.copy(); sep[n:] += 25.0
    assert I.lisi_stats(mixed, src, n_sub=800)[0] > 1.8     # close to the maximum of 2
    assert I.lisi_stats(sep, src, n_sub=800)[0] < 1.05      # no mixing


def test_harmony_removes_a_batch_shift_between_shared_clusters():
    rng = np.random.default_rng(1)
    n = 600
    cent = rng.normal(scale=5, size=(4, 10))
    Z = cent[rng.integers(0, 4, 2 * n)] + rng.normal(size=(2 * n, 10))
    src = np.array(["a"] * n + ["b"] * n)
    Z[n:] += 2.0  # batch shift
    pre = I.lisi_stats(Z, src, n_sub=1200)[0]
    post = I.lisi_stats(I.pairwise_harmony(Z, src), src, n_sub=1200)[0]
    assert pre < 1.1 and post > 1.5


def test_patient_bootstrap_ci_contains_median():
    rng = np.random.default_rng(0)
    v = rng.uniform(1, 2, 3000); pat = rng.integers(0, 15, 3000)
    med, lo, hi = I.patient_bootstrap_median(v, pat, n_boot=300)
    assert lo <= med <= hi and 1 < lo and hi < 2


def test_report_orders_ceiling_above_distorted_source(world, real_a, real_b):
    from pbisc_metrics import synthetic as S
    gen = S.sample_cohort(world, n_per_group=14, seed=5, lowexpr_floor=0.2, within_scale=0.0, prefix="g")  # blunted, distorted generator
    rep = I.ilisi_report(real_a["X"], real_a["genes"],
                         {"real_B": (real_b["X"], real_b["genes"]), "generated": (gen["X"], gen["genes"])},
                         ceiling="real_B", n_per_source=3000, n_hvg=300, n_pcs=20, n_draws=2)
    assert abs(rep.loc["real_B", "fraction_of_ceiling_post"] - 1) < 1e-9
    assert rep.loc["real_B", "ilisi_post"] > 1.7
    assert rep.loc["generated", "ilisi_post"] < rep.loc["real_B", "ilisi_post"]
    assert 0 <= rep.loc["generated", "fraction_of_ceiling_post"] < 1
