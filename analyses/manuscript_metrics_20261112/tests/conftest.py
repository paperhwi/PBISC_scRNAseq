import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pbisc_metrics import synthetic as S  # noqa: E402


@pytest.fixture(scope="session")
def world():
    return S.make_world(n_genes=400, seed=0)


@pytest.fixture(scope="session")
def real_a(world):
    return S.sample_cohort(world, n_per_group=14, seed=1, prefix="a")


@pytest.fixture(scope="session")
def real_b(world):
    return S.sample_cohort(world, n_per_group=14, seed=2, prefix="b")
