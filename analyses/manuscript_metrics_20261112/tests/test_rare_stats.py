import numpy as np

from pbisc_metrics import rare_stats as R
from pbisc_metrics import synthetic as S


def test_power_formula():
    assert abs(R.power_no_cell(0.001, 200) - 0.999 ** 200) < 1e-12
    assert R.power_no_cell(0.0, 100) == 1.0


def test_identical_sources_are_recovered_and_missing_state_is_below(world, real_a, real_b):
    t = R.rare_table(real_a["obs"], real_b["obs"], n_boot=200)
    assert (t.verdict == "recovered").all()
    gen = S.sample_cohort(world, n_per_group=14, seed=7, drop_fine=("NK_bright",), prefix="g")
    t2 = R.rare_table(real_a["obs"], gen["obs"], obs_ctrl=real_b["obs"], n_boot=200, rare_threshold=0.1).set_index("fine")
    assert t2.loc["NK_bright", "verdict"] == "below" and t2.loc["NK_bright", "n_gen"] == 0 and t2.loc["NK_bright", "ratio"] == 0
    assert t2.loc["NK_bright", "is_rare"] and t2.loc["NK_bright", "low_n"]
    assert t2.loc["T_naive", "verdict"] == "recovered"
    assert t2.loc["T_cyc", "ratio_ci_lo"] <= t2.loc["T_cyc", "ratio"] <= t2.loc["T_cyc", "ratio_ci_hi"] + 1e-9
    s = R.verdict_summary(t2.reset_index())
    assert s["below"] >= 1 and s.sum() == t2.is_rare.sum()


def test_reference_absent_state_not_assessable(world, real_a):
    ref = real_a["obs"][real_a["obs"].fine != "Mono_inflam"]
    t = R.rare_table(ref, real_a["obs"], n_boot=0).set_index("fine")
    assert t.loc["Mono_inflam", "verdict"] == "not_assessable"


def test_state_lost_in_the_real_control_is_not_assessable(world, real_a, real_b):
    ctrl = real_b["obs"][real_b["obs"].fine != "T_cyc"]
    t = R.rare_table(real_a["obs"], real_b["obs"], obs_ctrl=ctrl, n_boot=0).set_index("fine")
    assert t.loc["T_cyc", "verdict"] == "not_assessable"
