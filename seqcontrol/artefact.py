"""Analysis A: how much of a fixed-threshold sensitivity collapse is score compression?

All functions work on dL (mean log-likelihood, alt minus ref): more negative is more
damaging, a variant is called pathogenic when dL <= t, and AUROC uses -dL.
Compression means dL -> alpha * dL with 0 < alpha <= 1: every score moves toward
zero, the ordering is unchanged (so AUROC is unchanged exactly), but fewer scores
cross a cut-off fixed in absolute units.
"""

from __future__ import annotations

from dataclasses import dataclass
from statistics import NormalDist

import numpy as np

from seqcontrol import metrics

_N = NormalDist()


def sens_spec(labels, dl, t: float) -> tuple[float, float]:
    """Sensitivity and specificity when dL <= t is called pathogenic."""
    y, d = np.asarray(labels, int), np.asarray(dl, float)
    called = d <= t
    return float(called[y == 1].mean()), float((~called[y == 0]).mean())


def check_auroc_invariance(labels, dl, alphas, tol: float = 1e-12) -> None:
    """AUROC must be identical under any positive rescaling; anything else is a bug."""
    base = metrics.auroc(labels, -np.asarray(dl))
    for a in alphas:
        diff = abs(metrics.auroc(labels, -a * np.asarray(dl)) - base)
        if diff >= tol:
            raise AssertionError(f"AUROC changed by {diff} under scaling by {a}")


def scale_ratio(native, control) -> float:
    """Robust compression estimate: median |dL_control| / median |dL_native|."""
    return float(np.median(np.abs(control)) / np.median(np.abs(native)))


def origin_fit(native, control) -> tuple[float, float]:
    """Fit control ~ alpha * native through the origin; returns (alpha, R^2).

    R^2 is 1 - SSR / SS about the mean of `control` (not the uncentred version, which
    is inflated for fits through the origin), so it can be negative for a poor fit.
    """
    x, y = np.asarray(native, float), np.asarray(control, float)
    alpha = float(x @ y / (x @ x))
    ssr = float(((y - alpha * x) ** 2).sum())
    sst = float(((y - y.mean()) ** 2).sum())
    return alpha, 1 - ssr / sst


def rank_preserving_surrogate(native, control) -> np.ndarray:
    """Control's score values, reassigned so their ordering is native's.

    The variant with the k-th lowest native dL gets the k-th lowest control dL. It has
    the control's scale and distribution, and native's ranking.
    """
    native, control = np.asarray(native, float), np.asarray(control, float)
    out = np.empty_like(control)
    out[np.argsort(native, kind="mergesort")] = np.sort(control)
    return out


def share_explained(sens_native: float, sens_null: float, sens_observed: float) -> float:
    """Fraction of the observed sensitivity drop that pure compression reproduces."""
    return (sens_native - sens_null) / (sens_native - sens_observed)


def verdict_a3(share: float, auroc_drop_ci_includes_zero: bool) -> str:
    """The pre-declared interpretation rule (docs/plan-v0.2-threshold-and-gene.md, A3)."""
    if share >= 0.70 and auroc_drop_ci_includes_zero:
        return ("primarily a threshold artefact; genuine ranking loss is small and not "
                "resolvable at n = 67")  # fmt: skip
    if share < 0.40:
        return "hypothesis wrong: compression does not explain the drop"
    return "indeterminate at this sample size"


@dataclass(frozen=True)
class Binormal:
    """Gaussian dL distributions for pathogenic (p) and benign (b) variants."""

    mu_p: float
    sd_p: float
    mu_b: float
    sd_b: float

    def sens_spec(self, t: float, alpha: float = 1.0) -> tuple[float, float]:
        """At cut-off t after compressing every score by alpha."""
        sens = _N.cdf((t / alpha - self.mu_p) / self.sd_p)
        spec = 1 - _N.cdf((t / alpha - self.mu_b) / self.sd_b)
        return sens, spec

    def auroc(self) -> float:
        """P(pathogenic dL < benign dL); unchanged by compression."""
        return _N.cdf((self.mu_b - self.mu_p) / np.hypot(self.sd_p, self.sd_b))

    @classmethod
    def calibrated(cls, t: float, sens: float, spec: float, sd_p: float) -> Binormal:
        """Gaussians with the given sensitivity and specificity at t, t Youden-optimal.

        sens and spec fix each class's z-score at t; Youden-optimality (equal densities
        at t) fixes sd_b / sd_p. That leaves one free scale, sd_p, which the caller
        chooses and should vary to show its effect.
        """
        z_p = _N.inv_cdf(sens)  # (t - mu_p) / sd_p
        z_b = _N.inv_cdf(1 - spec)  # (t - mu_b) / sd_b
        sd_b = sd_p * _N.pdf(z_b) / _N.pdf(z_p)
        return cls(t - z_p * sd_p, sd_p, t - z_b * sd_b, sd_b)

    def alpha_for_sensitivity(self, t: float, target: float) -> float:
        """The compression alpha at which sensitivity at t falls to `target`."""
        z = _N.inv_cdf(target)
        return t / (self.mu_p + z * self.sd_p)
