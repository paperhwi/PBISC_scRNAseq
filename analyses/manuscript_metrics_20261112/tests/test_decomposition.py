import numpy as np
import pandas as pd

from pbisc_metrics import decomposition as D
from pbisc_metrics import synthetic as S


def _agg(c):
    meta = c["obs"].drop_duplicates("sample").set_index("sample")[["group", "patient"]]
    return D.aggregate(c["X"], c["genes"], c["obs"]["sample"], c["obs"]["parent"], meta)


def test_composition_plus_within_equals_total_exactly(real_a):
    a = _agg(real_a)
    r = D.decompose(a, (a.meta.group == "case").to_numpy(), (a.meta.group == "control").to_numpy())
    e = r["effects"]
    assert np.allclose(e.composition + e.within, e.total, atol=1e-8)
    assert abs(r["dprop_pp"].sum()) < 1e-8 and r["dprop_pp"]["Mono"] > 10  # synthetic truth: monocytes expand


def test_pure_composition_and_pure_within_effects_separate(world):
    comp_only = _agg(S.sample_cohort(world, n_per_group=20, seed=11, within_scale=0.0, prefix="x"))
    with_only = _agg(S.sample_cohort(world, n_per_group=20, seed=12, comp_scale=0.0, prefix="y"))
    f = lambda a: D.decompose(a, (a.meta.group == "case").to_numpy(), (a.meta.group == "control").to_numpy())["effects"]
    ec, ew = f(comp_only), f(with_only)
    assert ec.within.abs().mean() < 0.25 * ec.composition.abs().mean()
    assert ew.composition.abs().mean() < 0.25 * ew.within.abs().mean()


def test_split_half_slope_recovers_known_attenuation(world, real_a):
    real = _agg(real_a)
    half = _agg(S.sample_cohort(world, n_per_group=14, seed=21, within_scale=0.5, comp_scale=0.5, prefix="h"))  # effects halved
    same = _agg(S.sample_cohort(world, n_per_group=14, seed=22, prefix="s"))
    s_half, *_ = D.decomposition_report(real, half, "case", "control", n_splits=40)
    s_same, *_ = D.decomposition_report(real, same, "case", "control", n_splits=40)
    assert 0.3 < s_half.loc["within", "slope"] < 0.7
    assert 0.75 < s_same.loc["within", "slope"] < 1.25
    assert s_same.loc["within", "r_dis"] > s_half.loc["total", "r_dis"] - 0.5
    assert (s_same.n_splits_used > 0).all()


def test_split_is_deterministic(world, real_a):
    real = _agg(real_a)
    gen = _agg(S.sample_cohort(world, n_per_group=14, seed=23, prefix="d"))
    a, *_ = D.decomposition_report(real, gen, "case", "control", n_splits=10, seed=3)
    b, *_ = D.decomposition_report(real, gen, "case", "control", n_splits=10, seed=3)
    pd.testing.assert_frame_equal(a, b)
