"""Compare two per-variant score files variant by variant.

Usage:
    python scripts/compare_runs.py results/baseline_fp8-delayed.csv results/baseline_fp8-current.csv

Matches variants by position and alleles, then reports, overall and per region:
Spearman correlation of ΔL (score stability), the fraction of variants whose ΔL keeps
its sign, and the AUROC difference with a paired bootstrap 95% interval (the same
resampled variants are used for both runs). Runs anywhere; needs no GPU.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from seqcontrol import metrics

N_BOOT = 2000


def load(path: Path) -> dict[tuple, dict]:
    with open(path, newline="") as f:
        return {(r["pos"], r["ref"], r["alt"]): r for r in csv.DictReader(f)}


def paired_auroc_diff(y, a, b, seed=0) -> metrics.Estimate:
    """AUROC(a) - AUROC(b), with label-stratified bootstrap over shared variants."""
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    diffs = np.empty(N_BOOT)
    for i in range(N_BOOT):
        idx = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        diffs[i] = metrics.auroc(y[idx], a[idx]) - metrics.auroc(y[idx], b[idx])
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return metrics.Estimate(metrics.auroc(y, a) - metrics.auroc(y, b), lo, hi, N_BOOT)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_a", type=Path)
    parser.add_argument("run_b", type=Path)
    args = parser.parse_args()

    a, b = load(args.run_a), load(args.run_b)
    shared = sorted(set(a) & set(b), key=lambda k: int(k[0]))
    if len(shared) != len(a) or len(shared) != len(b):
        print(f"warning: runs share {len(shared)} variants ({len(a)} in A, {len(b)} in B)")

    print(f"A = {args.run_a}\nB = {args.run_b}\n")
    print(f"{'region':15s} {'n':>4s}  {'Spearman dL':>11s}  {'same sign':>9s}  AUROC(A) - AUROC(B)")
    regions = ["all"] + sorted({a[k]["region"] for k in shared})
    for region in regions:
        keys = [k for k in shared if region == "all" or a[k]["region"] == region]
        y = np.array([int(a[k]["label"]) for k in keys])
        da = np.array([float(a[k]["delta"]) for k in keys])
        db = np.array([float(b[k]["delta"]) for k in keys])
        rho = metrics.spearman(da, db)
        same = np.mean(np.sign(da) == np.sign(db))
        if min(y.sum(), (1 - y).sum()) >= 2:
            diff = str(paired_auroc_diff(y, -da, -db))
        else:
            diff = "(too few of one class)"
        print(f"{region:15s} {len(keys):4d}  {rho:11.3f}  {same:9.1%}  {diff}")


if __name__ == "__main__":
    main()
