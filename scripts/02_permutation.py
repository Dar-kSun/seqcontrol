"""Context-swap control on mtDNA tRNA variants: does Evo 2 read the gene or its address?

Usage (inside the WSL evo2 environment, after 01_baseline.py):
    python scripts/02_permutation.py --control trna-swap        # Mathur & Sachidanandam 2026
    python scripts/02_permutation.py --control window-rotation  # CLAUDE.md 5.1

Both controls keep each tRNA's own bases byte-identical (asserted for every perturbed
sequence) and change only its surroundings; see seqcontrol/controls/permute.py.
Every tRNA variant is scored under every setting (tRNA shift k, or rotation offset r),
plus the unperturbed setting 0, which must reproduce results/baseline_*.csv exactly.

Reports, with 95% paired bootstrap intervals over variants:
  AUROC native vs control (control = mean of per-setting AUROCs), and
  CDI = (AUROC_native - AUROC_control) / (AUROC_native - 0.5);
  sensitivity/specificity at a threshold fixed on native scores (the paper's design);
  Spearman of per-variant dL, native vs each setting;
  the median dL shift, to separate a shifted score distribution from lost ranking.
Writes results/permutation_<control>.csv and .json.
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
from seqcontrol.controls.permute import TrnaSwap, WindowRotation, merge_overlapping
from seqcontrol.data import mtdna
from seqcontrol.models.evo2 import PRECISIONS, Evo2Adapter
from seqcontrol.variants import Variant, variant_windows

WINDOW = 1025
HALF = WINDOW // 2
PAPER_THRESHOLD_DL = -0.0081
N_BOOT = 2000


def git_commit() -> str:
    out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=config.ROOT,
                         capture_output=True, text=True)  # fmt: skip
    return out.stdout.strip() or "unknown"


def with_alt(seq: str, i: int, ref: str, alt: str) -> str:
    if seq[i] != ref:
        raise AssertionError(f"expected {ref} at index {i}, found {seq[i]}")
    return seq[:i] + alt + seq[i + 1 :]


def trna_swap_windows(data, variants, units) -> dict[int, list[tuple[str, str]]]:
    """setting (shift k) -> per-variant (ref, alt) windows on the permuted chromosome."""
    chrom = Region(data.sequence, genes=tuple(units), circular=True)
    out = {}
    for k in range(len(units)):  # k = 0 is the unperturbed chromosome
        region = TrnaSwap.shift(chrom, units, k)
        out[k] = []
        for v in variants:
            j = region.map_index(v.pos - 1)
            moved = Variant(v.chrom, j + 1, v.ref, v.alt)
            out[k].append(variant_windows(region.sequence, moved, WINDOW, circular=True))
    return out


def window_rotation_windows(data, variants, units) -> dict[int, list[tuple[str, str]]]:
    """setting (offset r) -> per-variant (ref, alt) windows, rotated.

    Only offsets that leave every variant's own tRNA unit uncut are used, so every
    setting has every variant.
    """
    n = len(data.sequence)
    per_variant = []
    for v in variants:
        ref_w, _ = variant_windows(data.sequence, v, WINDOW, circular=True)
        start = (v.pos - 1 - HALF) % n  # window index 0 on the chromosome
        unit = next(u for u in units if u.start <= v.pos - 1 < u.end)
        gene = Interval(unit.name, (unit.start - start) % n, (unit.start - start) % n + len(unit))
        if gene.end > WINDOW:
            raise ValueError(f"{unit.name} does not fit inside the window around {v}")
        per_variant.append((v, Region(ref_w, genes=(gene,))))
    rot = WindowRotation(step=64)
    offsets = sorted(set.intersection(*(set(rot.valid_offsets(r)) for _, r in per_variant)))
    out = {}
    for r in [0, *offsets]:
        out[r] = []
        for v, region in per_variant:
            new = region if r == 0 else WindowRotation.rotate(region, r)
            i = new.map_index(HALF)
            out[r].append((new.sequence, with_alt(new.sequence, i, v.ref, v.alt)))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--control", choices=["trna-swap", "window-rotation"], required=True)
    parser.add_argument("--precision", choices=PRECISIONS, default="fp8-delayed")
    args = parser.parse_args()

    data = mtdna.load()
    trnas = [g for g in data.genes if g.biotype == "Mt_tRNA"]
    units = merge_overlapping([Interval(g.name, g.start - 1, g.end) for g in trnas])
    variants = [v for v in data.variants if any(u.start <= v.pos - 1 < u.end for u in units)]
    y = np.array([v.label for v in variants])
    print(f"{len(variants)} tRNA variants ({y.sum()} P / {len(y) - y.sum()} B), "
          f"{len(units)} tRNA units")  # fmt: skip

    build = trna_swap_windows if args.control == "trna-swap" else window_rotation_windows
    windows = build(data, variants, units)
    settings = [s for s in windows if s != 0]
    unique = sorted({w for pairs in windows.values() for pair in pairs for w in pair})
    print(f"{len(settings)} control settings + native; {len(unique)} distinct windows to score")

    model = Evo2Adapter("evo2_1b_base", precision=args.precision)
    model.load()
    t0 = time.perf_counter()
    ll = dict(zip(unique, model.score_sequences(unique), strict=True))
    elapsed = time.perf_counter() - t0
    print(f"scored in {elapsed:.0f} s")
    dl = {s: np.array([ll[a] - ll[r] for r, a in pairs]) for s, pairs in windows.items()}

    # Setting 0 must be exactly the M3 baseline.
    with open(config.ROOT / "results" / f"baseline_{args.precision}.csv", newline="") as f:
        base = {(int(r["pos"]), r["ref"], r["alt"]): float(r["delta"]) for r in csv.DictReader(f)}
    native_from_baseline = np.array([base[(v.pos, v.ref, v.alt)] for v in variants])
    if not np.array_equal(dl[0], native_from_baseline):
        diff = np.abs(dl[0] - native_from_baseline).max()
        raise AssertionError(f"unperturbed setting differs from baseline by up to {diff}")
    print("setting 0 reproduces the baseline exactly")

    native = dl[0]
    control = np.stack([dl[s] for s in settings])  # settings x variants

    def auroc_native(i):
        return metrics.auroc(y[i], -native[i])

    def auroc_control(i):
        return float(np.mean([metrics.auroc(y[i], -c[i]) for c in control]))

    def cdi(i):
        return metrics.cdi(auroc_native(i), auroc_control(i))

    est = {
        "auroc_native": metrics.paired_bootstrap(auroc_native, y, N_BOOT),
        "auroc_control": metrics.paired_bootstrap(auroc_control, y, N_BOOT),
        "cdi": metrics.paired_bootstrap(cdi, y, N_BOOT),
        "auroc_control_of_mean_dL": metrics.bootstrap(
            metrics.auroc, y, -control.mean(axis=0), N_BOOT
        ),
    }
    per_setting_auroc = [metrics.auroc(y, -c) for c in control]
    rho = [metrics.spearman(native, c) for c in control]

    t_native = -metrics.youden_threshold(y, -native)  # as a dL cut-off: call P if dL <= t

    def sens_spec(d, t):
        return metrics.sensitivity_specificity(y, -d, -t)

    thresholds = {}
    for name, t in [("paper", PAPER_THRESHOLD_DL), ("native_youden", t_native)]:
        nat = sens_spec(native, t)
        ctl = np.mean([sens_spec(c, t) for c in control], axis=0)
        thresholds[name] = {
            "threshold_dL": t,
            "native": {"sensitivity": nat[0], "specificity": nat[1]},
            "control_mean": {"sensitivity": float(ctl[0]), "specificity": float(ctl[1])},
        }
    shift = {
        cls: {"native_median_dL": float(np.median(native[y == lab])),
              "control_median_dL": float(np.median(control[:, y == lab]))}
        for cls, lab in [("pathogenic", 1), ("benign", 0)]
    }  # fmt: skip

    results = {
        "control": args.control,
        "model": model.label,
        "window_bp": WINDOW,
        "n_pathogenic": int(y.sum()),
        "n_benign": int(len(y) - y.sum()),
        "settings": settings,
        "setting_meaning": "tRNA shift k" if args.control == "trna-swap" else "rotation bp",
        **{k: v.as_dict() for k, v in est.items()},
        "auroc_per_setting": dict(zip(map(str, settings), per_setting_auroc, strict=True)),
        "spearman_native_vs_setting": {
            "mean": float(np.mean(rho)),
            "min": min(rho),
            "max": max(rho),
        },  # fmt: skip
        "thresholds": thresholds,
        "median_dL_shift": shift,
        "bootstrap": f"{N_BOOT} label-stratified resamples of variants, paired across settings",
        "git_commit": git_commit(),
        "run_date": dt.date.today().isoformat(),
        "scoring_seconds": round(elapsed, 1),
    }

    out = config.ROOT / "results"
    stem = f"permutation_{args.control}"
    with open(out / f"{stem}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pos", "ref", "alt", "label", "gene", "setting", "delta"])
        for s, d in dl.items():
            for v, x in zip(variants, d, strict=True):
                w.writerow([v.pos, v.ref, v.alt, v.label, data.gene_of(v).name, s, x])
    (out / f"{stem}.json").write_text(json.dumps(results, indent=2) + "\n")

    print(f"\n{args.control}, {model.label}, tRNA variants")
    print(f"  AUROC native             {est['auroc_native']}")
    print(f"  AUROC control (mean of {len(settings)}) {est['auroc_control']}")
    print(f"  per-setting AUROC range  {min(per_setting_auroc):.3f} - {max(per_setting_auroc):.3f}")
    print(f"  CDI                      {est['cdi']}")
    print(f"  Spearman dL native vs control: mean {np.mean(rho):.3f} "
          f"(range {min(rho):.3f} - {max(rho):.3f}); FP8-recipe noise floor ~0.95")  # fmt: skip
    for name, t in thresholds.items():
        n, c = t["native"], t["control_mean"]
        sens = f"sens {n['sensitivity']:.3f} -> {c['sensitivity']:.3f}"
        spec = f"spec {n['specificity']:.3f} -> {c['specificity']:.3f}"
        print(f"  at {name} cut-off dL <= {t['threshold_dL']:.4f}: {sens}, {spec}")
    for cls, m in shift.items():
        print(
            f"  median dL {cls:10s} {m['native_median_dL']:+.5f} -> {m['control_median_dL']:+.5f}"
        )
    print(f"\nwrote results/{stem}.csv and results/{stem}.json")


if __name__ == "__main__":
    main()
