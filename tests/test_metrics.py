import numpy as np
import pytest

from seqcontrol.metrics import (
    auprc,
    auroc,
    bootstrap,
    cdi,
    sensitivity_specificity,
    spearman,
    youden_threshold,
)


def test_auroc_perfect_inverted_and_tied():
    assert auroc([0, 0, 1, 1], [1, 2, 3, 4]) == 1.0
    assert auroc([0, 0, 1, 1], [4, 3, 2, 1]) == 0.0
    assert auroc([0, 1], [5, 5]) == 0.5


def test_auroc_counts_pairs():
    # Pairs (pos, neg): (3,1) win, (3,4) loss, (5,1) win, (5,4) win -> 3/4.
    assert auroc([1, 0, 1, 0], [3, 1, 5, 4]) == 0.75


def test_auprc_hand_computed():
    # Ranked: 0.9 pos, 0.8 neg, 0.7 pos, 0.1 neg -> AP = 0.5*1 + 0.5*(2/3).
    assert auprc([1, 0, 1, 0], [0.9, 0.8, 0.7, 0.1]) == pytest.approx(0.5 + 1 / 3)


def test_auprc_tied_scores_are_one_threshold():
    # All tied: one threshold, precision = prevalence.
    assert auprc([1, 0, 0, 0], [1, 1, 1, 1]) == pytest.approx(0.25)


def test_spearman():
    assert spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    assert spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)


def test_cdi():
    assert cdi(0.9, 0.9) == 0.0
    assert cdi(0.9, 0.5) == pytest.approx(1.0)
    assert cdi(0.8, 0.65) == pytest.approx(0.5)
    with pytest.raises(ValueError):
        cdi(0.5, 0.4)


def test_youden_and_sensitivity_specificity():
    y, s = [0, 0, 0, 1, 1], [0.1, 0.2, 0.6, 0.5, 0.9]
    t = youden_threshold(y, s)
    sens, spec = sensitivity_specificity(y, s, t)
    assert (t, sens, spec) == (0.5, 1.0, 2 / 3)


def test_bootstrap_interval_brackets_point_and_is_reproducible():
    rng = np.random.default_rng(1)
    y = np.r_[np.ones(40), np.zeros(60)].astype(int)
    s = y + rng.normal(0, 1, 100)
    a = bootstrap(auroc, y, s, n_boot=500, seed=3)
    b = bootstrap(auroc, y, s, n_boot=500, seed=3)
    assert a == b
    assert a.lo < a.value < a.hi
    assert 0.6 < a.value < 0.9


def test_bad_inputs_are_rejected():
    with pytest.raises(ValueError):
        auroc([1, 1], [0.1, 0.2])
    with pytest.raises(ValueError):
        auroc([0, 2], [0.1, 0.2])
    with pytest.raises(ValueError):
        auroc([0, 1], [0.1, np.nan])


def test_paired_bootstrap_of_a_difference_is_zero_for_identical_scores():
    from seqcontrol.metrics import paired_bootstrap

    y = np.r_[np.ones(20), np.zeros(30)].astype(int)
    s = np.random.default_rng(0).normal(size=50) + y
    est = paired_bootstrap(lambda i: auroc(y[i], s[i]) - auroc(y[i], s[i]), y, n_boot=200)
    assert (est.value, est.lo, est.hi) == (0.0, 0.0, 0.0)
