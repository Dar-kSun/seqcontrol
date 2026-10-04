import numpy as np
import pytest

from seqcontrol import metrics
from seqcontrol.artefact import (
    Binormal,
    check_auroc_invariance,
    origin_fit,
    rank_preserving_surrogate,
    scale_ratio,
    sens_spec,
    share_explained,
    verdict_a3,
)


@pytest.fixture
def scores():
    rng = np.random.default_rng(0)
    y = np.r_[np.ones(40), np.zeros(30)].astype(int)
    dl = np.where(y == 1, rng.normal(-0.006, 0.004, 70), rng.normal(-0.001, 0.002, 70))
    return y, dl


def test_scaling_never_changes_auroc(scores):
    y, dl = scores
    check_auroc_invariance(y, dl, [1.0, 0.7, 0.3, 0.1, 0.01])


def test_class_specific_scaling_does_change_auroc(scores):
    # Why the A3 null must use one alpha: scaling the classes differently reorders.
    y, dl = scores
    per_class = dl * np.where(y == 1, 0.2, 1.0)
    assert metrics.auroc(y, -per_class) != metrics.auroc(y, -dl)


def test_compression_lowers_sensitivity_and_raises_specificity(scores):
    y, dl = scores
    t = -0.003
    s1, p1 = sens_spec(y, dl, t)
    s2, p2 = sens_spec(y, 0.4 * dl, t)
    assert s2 < s1 and p2 >= p1


def test_scale_estimates_recover_a_pure_compression(scores):
    _, dl = scores
    assert scale_ratio(dl, 0.6 * dl) == pytest.approx(0.6)
    alpha, r2 = origin_fit(dl, 0.6 * dl)
    assert alpha == pytest.approx(0.6) and r2 == pytest.approx(1.0)


def test_surrogate_has_control_values_and_native_order(scores):
    y, dl = scores
    control = np.random.default_rng(1).permutation(0.5 * dl)
    sur = rank_preserving_surrogate(dl, control)
    assert sorted(sur) == sorted(control)
    assert metrics.spearman(sur, dl) == pytest.approx(1.0)
    assert metrics.auroc(y, -sur) == pytest.approx(metrics.auroc(y, -dl))


def test_share_and_verdict():
    assert share_explained(0.8, 0.6, 0.5) == pytest.approx(2 / 3)
    assert verdict_a3(0.75, True).startswith("primarily")
    assert verdict_a3(0.75, False) == "indeterminate at this sample size"
    assert verdict_a3(0.3, True).startswith("hypothesis wrong")


def test_binormal_calibration_hits_the_paper_operating_point():
    t = -0.0081
    b = Binormal.calibrated(t, sens=0.658, spec=0.785, sd_p=0.01)
    sens, spec = b.sens_spec(t)
    assert sens == pytest.approx(0.658) and spec == pytest.approx(0.785)
    # Youden-optimal: densities equal at t, so a small move of t changes J by ~0.
    j = lambda tt: sum(b.sens_spec(tt)) - 1  # noqa: E731
    assert j(t) >= j(t + 1e-5) and j(t) >= j(t - 1e-5)
    a = b.alpha_for_sensitivity(t, 0.051)
    assert b.sens_spec(t, a)[0] == pytest.approx(0.051)
    assert b.sens_spec(t, a)[1] > spec  # specificity rises under compression
