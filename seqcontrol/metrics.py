"""Metrics, defined once and reused everywhere.

Convention: `labels` are 1 for pathogenic and 0 for benign; `scores` are higher for
"more likely pathogenic". For Evo 2 that means the negated variant score, -(ΔL),
because a damaging variant makes the sequence less probable.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass

import numpy as np


def _check(labels, scores) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    if y.shape != s.shape or y.ndim != 1:
        raise ValueError("labels and scores must be 1-D arrays of the same length")
    if not set(np.unique(y)) <= {0, 1}:
        raise ValueError("labels must be 0 or 1")
    if not np.all(np.isfinite(s)):
        raise ValueError("scores must be finite")
    return y, s


def _ranks(x: np.ndarray) -> np.ndarray:
    """1-based ranks; tied values share their average rank."""
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    ranks = np.empty(len(x))
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and xs[j + 1] == xs[i]:
            j += 1
        ranks[order[i : j + 1]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def auroc(labels, scores) -> float:
    """Area under the ROC curve: P(random pathogenic scores above random benign), ties = 1/2."""
    y, s = _check(labels, scores)
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    if n_pos == 0 or n_neg == 0:
        raise ValueError("AUROC needs both classes")
    r = _ranks(s)
    return float((r[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def auprc(labels, scores) -> float:
    """Average precision: sum over distinct thresholds of precision x recall gained.

    The no-skill baseline is the fraction of pathogenic variants, not 0.5.
    """
    y, s = _check(labels, scores)
    n_pos = int(y.sum())
    if n_pos == 0:
        raise ValueError("AUPRC needs at least one pathogenic variant")
    order = np.argsort(-s, kind="mergesort")
    s, y = s[order], y[order]
    last_of_tie = np.r_[np.diff(s) != 0, True]  # evaluate only once all tied scores are in
    tp = np.cumsum(y)[last_of_tie]
    k = np.arange(1, len(y) + 1)[last_of_tie]
    recall_gain = np.diff(np.r_[0, tp]) / n_pos
    return float(np.sum(recall_gain * tp / k))


def spearman(a, b) -> float:
    """Spearman rank correlation, with average ranks for ties."""
    ra, rb = _ranks(np.asarray(a, float)), _ranks(np.asarray(b, float))
    return float(np.corrcoef(ra, rb)[0, 1])


def cdi(auroc_native: float, auroc_control: float) -> float:
    """Context-dependence index: the share of above-chance AUROC lost under a control.

    CDI = (AUROC_native - AUROC_control) / (AUROC_native - 0.5)
    Near 1: essentially all signal was context. Near 0: the model reads the sequence.
    Undefined when the native AUROC is at chance.
    """
    if auroc_native <= 0.5:
        raise ValueError("CDI is undefined when native AUROC is at or below chance")
    return (auroc_native - auroc_control) / (auroc_native - 0.5)


def youden_threshold(labels, scores) -> float:
    """The score cut-off maximising sensitivity + specificity - 1 (call pathogenic if >=)."""
    y, s = _check(labels, scores)
    best_j, best_t = -np.inf, None
    for t in np.unique(s):
        sens, spec = sensitivity_specificity(y, s, t)
        if sens + spec - 1 > best_j:
            best_j, best_t = sens + spec - 1, t
    return float(best_t)


def sensitivity_specificity(labels, scores, threshold: float) -> tuple[float, float]:
    """Fraction of pathogenic scored >= threshold, and of benign scored below it."""
    y, s = _check(labels, scores)
    called = s >= threshold
    return float(called[y == 1].mean()), float((~called[y == 0]).mean())


@dataclass(frozen=True)
class Estimate:
    value: float
    lo: float
    hi: float
    n_boot: int

    def __str__(self) -> str:
        return f"{self.value:.3f} [{self.lo:.3f}, {self.hi:.3f}]"

    def as_dict(self) -> dict:
        return {"value": self.value, "ci95": [self.lo, self.hi], "n_boot": self.n_boot}


def bootstrap(
    metric: Callable[[np.ndarray, np.ndarray], float],
    labels,
    scores,
    n_boot: int = 2000,
    seed: int = 0,
) -> Estimate:
    """Point estimate plus a 95% percentile interval from resampling variants.

    Resamples are stratified by label, so every resample keeps the original number
    of pathogenic and benign variants (and the metric is always defined).
    """
    y, s = _check(labels, scores)
    values = [metric(y[idx], s[idx]) for idx in stratified_resamples(y, n_boot, seed)]
    lo, hi = np.percentile(values, [2.5, 97.5])
    return Estimate(metric(y, s), float(lo), float(hi), n_boot)


def stratified_resamples(labels, n_boot: int = 2000, seed: int = 0) -> Iterator[np.ndarray]:
    """Index arrays resampling variants with replacement, keeping each class's count.

    For paired comparisons (native vs control on the same variants), compute both
    statistics on the same index array.
    """
    y = np.asarray(labels, dtype=int)
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    for _ in range(n_boot):
        yield np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]


def paired_bootstrap(
    statistic: Callable[[np.ndarray], float], labels, n_boot: int = 2000, seed: int = 0
) -> Estimate:
    """Estimate for `statistic(idx)`, a function of a variant index array.

    Called once with all variants for the point value, then on each stratified
    resample for the 95% percentile interval.
    """
    y = np.asarray(labels, dtype=int)
    values = [statistic(idx) for idx in stratified_resamples(y, n_boot, seed)]
    lo, hi = np.percentile(values, [2.5, 97.5])
    return Estimate(statistic(np.arange(len(y))), float(lo), float(hi), n_boot)


def within_group_auroc(labels, scores, groups) -> float:
    """AUROC over pathogenic-benign pairs from the same group only (ties count 1/2).

    A stratified Mann-Whitney statistic: comparisons between groups never enter, so
    a score that only tells groups apart gets 0.5. Undefined (ValueError) when no
    group holds both classes.
    """
    y, s = _check(labels, scores)
    g = np.asarray(groups)
    wins, pairs = 0.0, 0
    for name in np.unique(g):
        m = g == name
        pos, neg = s[m & (y == 1)], s[m & (y == 0)]
        if len(pos) and len(neg):
            diff = pos[:, None] - neg[None, :]
            wins += float((diff > 0).sum() + 0.5 * (diff == 0).sum())
            pairs += diff.size
    if pairs == 0:
        raise ValueError("no group contains both classes")
    return wins / pairs


def within_group_pairs(labels, groups) -> int:
    y, g = np.asarray(labels, dtype=int), np.asarray(groups)
    return int(sum((y[g == k] == 1).sum() * (y[g == k] == 0).sum() for k in np.unique(g)))


@dataclass(frozen=True)
class ClusterResample:
    idx: np.ndarray  # variant indices; a cluster drawn twice contributes its variants twice
    draws: int  # attempts needed before this resample was accepted
    groups: np.ndarray  # cluster label per index; duplicated draws get distinct labels


def cluster_resamples(
    groups, accept: Callable[[np.ndarray, np.ndarray], bool], n_boot: int = 2000, seed: int = 0
) -> Iterator[ClusterResample]:
    """Resample whole clusters (genes) with replacement.

    Draws the same number of clusters as the data has. A resample is redrawn until
    `accept(idx, groups)` is true (e.g. both classes present); `draws` records how
    many attempts that took. A cluster drawn twice gets two distinct labels, so
    within-cluster statistics treat the copies as separate clusters.
    """
    g = np.asarray(groups)
    names = np.unique(g)
    members = {k: np.flatnonzero(g == k) for k in names}
    rng = np.random.default_rng(seed)
    for _ in range(n_boot):
        draws = 0
        while True:
            draws += 1
            picked = rng.choice(names, len(names))
            idx = np.concatenate([members[k] for k in picked])
            labels = np.concatenate([np.full(len(members[k]), i) for i, k in enumerate(picked)])
            if accept(idx, labels):
                break
        yield ClusterResample(idx, draws, labels)


def cluster_bootstrap(
    statistic: Callable[[np.ndarray, np.ndarray], float],
    groups,
    accept: Callable[[np.ndarray, np.ndarray], bool],
    n_boot: int = 2000,
    seed: int = 0,
) -> tuple[Estimate, float]:
    """Estimate for `statistic(idx, resample_groups)` with a gene-level 95% interval.

    Returns the estimate and the mean number of draws per accepted resample.
    """
    g = np.asarray(groups)
    values, draws = [], []
    for r in cluster_resamples(g, accept, n_boot, seed):
        values.append(statistic(r.idx, r.groups))
        draws.append(r.draws)
    lo, hi = np.percentile(values, [2.5, 97.5])
    point = statistic(np.arange(len(g)), g)
    return Estimate(point, float(lo), float(hi), n_boot), float(np.mean(draws))
