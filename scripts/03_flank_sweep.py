"""Flank-shuffle sweep on mtDNA tRNA variants: how far out does Evo 2's context reach?

Usage (inside the WSL evo2 environment, after 01_baseline.py):
    python scripts/03_flank_sweep.py              # score on the GPU (~30 min)
    python scripts/03_flank_sweep.py --from-csv   # recompute summaries, no GPU

For each tRNA variant's 1,025 bp window, the tRNA is held fixed and everything more
than r bp outside it is dinucleotide-shuffled (seqcontrol/controls/flank_shuffle.py),
for r in RADII and 10 random seeds. r = 0 shuffles all context; at r = 400 only the
outer ~40-80 bp on each side are shuffled. Radii beyond ~480 bp would fall outside the
window and change nothing, so the sweep stops there (see docs/findings.md).

Per radius, with 95% paired bootstrap intervals over variants and the spread across
seeds: Spearman of per-variant dL against native (score stability), AUROC, and CDI.
Writes results/flank_sweep.csv and results/flank_sweep.json.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import subprocess
import time

import numpy as np

from seqcontrol import config, metrics
from seqcontrol.controls.base import Interval, Region
from seqcontrol.controls.flank_shuffle import FlankShuffle
from seqcontrol.controls.permute import merge_overlapping
from seqcontrol.data import mtdna
from seqcontrol.models.evo2 import PRECISIONS, Evo2Adapter
from seqcontrol.variants import variant_windows

WINDOW = 1025
HALF = WINDOW // 2
RADII = [0, 25, 50, 100, 200, 300, 400]
N_SEEDS = 10
SEED = 20261004
N_BOOT = 2000


def git_commit() -> str:
    out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=config.ROOT,
                         capture_output=True, text=True)  # fmt: skip
    return out.stdout.strip() or "unknown"


def with_alt(seq: str, i: int, ref: str, alt: str) -> str:
    if seq[i] != ref:
        raise AssertionError(f"expected {ref} at index {i}, found {seq[i]}")
    return seq[:i] + alt + seq[i + 1 :]


def score(args):
    data = mtdna.load()
    trnas = [g for g in data.genes if g.biotype == "Mt_tRNA"]
    units = merge_overlapping([Interval(g.name, g.start - 1, g.end) for g in trnas])
    variants = [v for v in data.variants if any(u.start <= v.pos - 1 < u.end for u in units)]
    n = len(data.sequence)
    shuffler = FlankShuffle(RADII, N_SEEDS)

    pairs: dict[tuple[int, int], list[tuple[str, str]]] = {}  # (radius, seed) -> per variant
    native = []
    for v in variants:
        ref_w, alt_w = variant_windows(data.sequence, v, WINDOW, circular=True)
        native.append((ref_w, alt_w))
        start = (v.pos - 1 - HALF) % n
        unit = next(u for u in units if u.start <= v.pos - 1 < u.end)
        g0 = (unit.start - start) % n
        region = Region(ref_w, genes=(Interval(unit.name, g0, g0 + len(unit)),))
        rng = np.random.default_rng([SEED, v.pos, "ACGT".index(v.alt)])
        for shuffled in shuffler.apply(region, rng):
            key = (shuffled.meta["radius"], shuffled.meta["seed"])
            alt = with_alt(shuffled.sequence, HALF, v.ref, v.alt)
            pairs.setdefault(key, []).append((shuffled.sequence, alt))

    unique = sorted({w for ps in [native, *pairs.values()] for p in ps for w in p})
    print(f"{len(variants)} tRNA variants, {len(RADII)} radii x {N_SEEDS} seeds; "
          f"{len(unique)} distinct windows to score")  # fmt: skip
    model = Evo2Adapter("evo2_1b_base", precision=args.precision)
    model.load()
    t0 = time.perf_counter()
    ll = dict(zip(unique, model.score_sequences(unique), strict=True))
    elapsed = time.perf_counter() - t0
    print(f"scored in {elapsed:.0f} s")

    def dl(ps):
        return np.array([ll[a] - ll[r] for r, a in ps])

    rows = []
    for key, d in [((-1, -1), dl(native))] + [(k, dl(ps)) for k, ps in pairs.items()]:
        for v, x in zip(variants, d, strict=True):
            rows.append([v.pos, v.ref, v.alt, v.label, data.gene_of(v).name, *key, x])
    provenance = {"model": model.label, "scored_at_commit": git_commit(),
                  "score_date": dt.date.today().isoformat(),
                  "scoring_seconds": round(elapsed, 1)}  # fmt: skip
    return rows, provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--precision", choices=PRECISIONS, default="fp8-delayed")
    parser.add_argument(
        "--from-csv", action="store_true", help="recompute summaries from the saved CSV (no GPU)"
    )
    args = parser.parse_args()
    out = config.ROOT / "results"
    header = ["pos", "ref", "alt", "label", "gene", "radius", "seed", "delta"]

    if args.from_csv:
        with open(out / "flank_sweep.csv", newline="") as f:
            rows = [[r[h] for h in header] for r in csv.DictReader(f)]
        old = json.loads((out / "flank_sweep.json").read_text())
        provenance = {k: old[k] for k in ("model", "scored_at_commit", "score_date",
                                          "scoring_seconds")}  # fmt: skip
    else:
        rows, provenance = score(args)
        with open(out / "flank_sweep.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows)

    # Arrange as delta[(radius, seed)] -> per-variant array, in a fixed variant order.
    keys, labels, delta = [], {}, {}
    for pos, ref, alt, label, _gene, radius, seed, d in rows:
        k = (int(pos), ref, alt)
        if k not in labels:
            keys.append(k)
            labels[k] = int(label)
        delta.setdefault((int(radius), int(seed)), []).append(float(d))
    delta = {s: np.array(d) for s, d in delta.items()}
    y = np.array([labels[k] for k in keys])
    native = delta[(-1, -1)]

    with open(out / f"baseline_{args.precision}.csv", newline="") as f:
        base = {(int(r["pos"]), r["ref"], r["alt"]): float(r["delta"]) for r in csv.DictReader(f)}
    if not np.array_equal(native, np.array([base[k] for k in keys])):
        raise AssertionError("native scores differ from the M3 baseline")
    print("native scores reproduce the baseline exactly")

    au_native = metrics.auroc(y, -native)
    per_radius = {}
    print(f"\nflank shuffle, {provenance['model']}, {len(y)} tRNA variants "
          f"({y.sum()} P / {len(y) - y.sum()} B); native AUROC {au_native:.3f}")  # fmt: skip
    print(f"{'radius':>6s}  {'Spearman vs native':>26s}  {'AUROC':>22s}  {'CDI':>22s}  seed SD")
    for r in RADII:
        runs = np.stack([delta[(r, s)] for s in range(N_SEEDS)])  # seeds x variants

        def rho(i, runs=runs):
            return float(np.mean([metrics.spearman(native[i], d[i]) for d in runs]))

        def auc(i, runs=runs):
            return float(np.mean([metrics.auroc(y[i], -d[i]) for d in runs]))

        def cdi(i, auc=auc):
            return metrics.cdi(metrics.auroc(y[i], -native[i]), auc(i))

        est = {
            "spearman": metrics.paired_bootstrap(rho, y, N_BOOT),
            "auroc": metrics.paired_bootstrap(auc, y, N_BOOT),
            "cdi": metrics.paired_bootstrap(cdi, y, N_BOOT),
        }
        seed_rho = [metrics.spearman(native, d) for d in runs]
        seed_auc = [metrics.auroc(y, -d) for d in runs]
        per_radius[r] = {
            **{k: e.as_dict() for k, e in est.items()},
            "spearman_per_seed": seed_rho,
            "auroc_per_seed": seed_auc,
            "spearman_seed_sd": float(np.std(seed_rho, ddof=1)),
            "auroc_seed_sd": float(np.std(seed_auc, ddof=1)),
        }
        sd = f"rho {np.std(seed_rho, ddof=1):.3f}, AUROC {np.std(seed_auc, ddof=1):.3f}"
        print(f"{r:6d}  {est['spearman']!s:>26s}  {est['auroc']!s:>22s}  {est['cdi']!s:>22s}  {sd}")

    results = {
        "control": "flank-shuffle (dinucleotide, Altschul-Erikson)",
        **provenance,
        "window_bp": WINDOW,
        "radii_bp": RADII,
        "radius_meaning": "bp of untouched sequence kept on each side of the tRNA",
        "n_seeds": N_SEEDS,
        "n_pathogenic": int(y.sum()),
        "n_benign": int(len(y) - y.sum()),
        "auroc_native": au_native,
        "per_radius": {str(r): v for r, v in per_radius.items()},
        "bootstrap": f"{N_BOOT} label-stratified resamples of variants; each statistic is "
        "the mean over seeds, recomputed per resample",
        "analysis_commit": git_commit(),
    }
    (out / "flank_sweep.json").write_text(json.dumps(results, indent=2) + "\n")
    wrote = "results/flank_sweep.json" if args.from_csv else "results/flank_sweep.csv and .json"
    print(f"\nwrote {wrote}")


if __name__ == "__main__":
    main()
