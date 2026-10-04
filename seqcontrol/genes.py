"""Gene-level helpers for the gene-identity analyses (v0.2 B, v0.3)."""

from __future__ import annotations

import numpy as np

from seqcontrol import metrics


def gene_priors(y: np.ndarray, groups: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(leave-one-variant-out, in-sample) pathogenic fraction of each variant's gene.

    Both ignore the variant itself and score only which gene it is in. The
    leave-one-out version is biased low (removing a variant pushes its gene's
    fraction away from its own label); the in-sample version is biased high.
    """
    n_all, p_all = len(y), y.sum()
    lovo, insample = np.empty(len(y)), np.empty(len(y))
    for g in np.unique(groups):
        m = groups == g
        n, p = m.sum(), y[m].sum()
        insample[m] = p / n
        if n > 1:
            lovo[m] = (p - y[m]) / (n - 1)
        else:  # alone in its gene: fall back to every other variant
            lovo[m] = (p_all - y[m]) / (n_all - 1)
    return lovo, insample


def both_classes(y):
    """Accept a cluster resample only if it holds both pathogenic and benign variants."""
    return lambda idx, grp: 0 < y[idx].sum() < len(idx)


def has_pairs(y):
    """Accept a cluster resample only if it holds at least one within-gene pair."""
    return lambda idx, grp: metrics.within_group_pairs(y[idx], grp) > 0


def fmt(d: dict, spec: str = ".3f") -> str:
    """'value [lo, hi]' for an Estimate.as_dict()."""
    return f"{d['value']:{spec}} [{d['ci95'][0]:{spec}}, {d['ci95'][1]:{spec}}]"
