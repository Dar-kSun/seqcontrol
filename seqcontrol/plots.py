"""Figures. Needs the `plot` extra: pip install -e ".[plot]"."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from seqcontrol import metrics

# Reference palette (dataviz skill, light mode): one data series plus recessive inks.
SURFACE = "#fcfcfb"
SERIES = "#2a78d6"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
MUTED = "#a3a29c"
GRID = "#e4e3df"


def noise_floor_spearman(results: Path, region: str = "Mt_tRNA") -> float:
    """Spearman of dL between the two FP8 recipes on the baseline variants of `region`."""

    def load(name):
        with open(results / name, newline="") as f:
            rows = [r for r in csv.DictReader(f) if r["region"] == region]
        return {(r["pos"], r["ref"], r["alt"]): float(r["delta"]) for r in rows}

    a, b = load("baseline_fp8-delayed.csv"), load("baseline_fp8-current.csv")
    keys = sorted(a)
    return metrics.spearman([a[k] for k in keys], [b[k] for k in keys])


def _style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=TEXT_2, labelsize=9)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def _reference(ax, y: float, label: str, x_text: float, va: str = "bottom") -> None:
    ax.axhline(y, color=MUTED, linewidth=1.2, linestyle=(0, (4, 3)), zorder=1)
    offset = 0.012 if va == "bottom" else -0.012
    ax.text(x_text, y + offset, label, color=TEXT_2, fontsize=8.5, ha="right", va=va)


def _series(ax, radii, est, per_seed) -> None:
    lo = [e["ci95"][0] for e in est]
    hi = [e["ci95"][1] for e in est]
    val = [e["value"] for e in est]
    ax.fill_between(radii, lo, hi, color=SERIES, alpha=0.14, linewidth=0, zorder=2)
    rng = np.random.default_rng(0)  # small horizontal jitter so seeds don't stack
    for r, seeds in zip(radii, per_seed, strict=True):
        jitter = rng.uniform(-5, 5, len(seeds))
        ax.scatter(r + jitter, seeds, s=9, color=SERIES, alpha=0.35, linewidths=0, zorder=3)
    ax.plot(radii, val, color=SERIES, linewidth=2, zorder=4)
    ax.scatter(radii, val, s=64, color=SERIES, edgecolors=SURFACE, linewidths=2, zorder=5)


def flank_sweep_figure(results: Path, out: Path) -> None:
    """The signature figure: score stability and AUROC against untouched flank radius."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sweep = json.loads((results / "flank_sweep.json").read_text())
    radii = sweep["radii_bp"]
    per = [sweep["per_radius"][str(r)] for r in radii]
    floor = noise_floor_spearman(results)
    swap = json.loads((results / "permutation_trna-swap.json").read_text())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2), facecolor=SURFACE)
    for ax in (ax1, ax2):
        _style(ax)
        ax.set_xlim(-20, max(radii) + 20)
        ax.set_xticks(radii)
        ax.set_xlabel("Untouched flank kept on each side of the tRNA (bp)", color=TEXT_2)
    x_text = max(radii) + 15

    _series(ax1, radii, [p["spearman"] for p in per], [p["spearman_per_seed"] for p in per])
    _reference(ax1, floor, f"FP8 rounding-recipe noise floor {floor:.2f}", x_text)
    swap_rho = swap["spearman_native_vs_setting"]["mean"]
    _reference(ax1, swap_rho, f"tRNA moved to another tRNA's slot {swap_rho:.2f}", x_text)
    ax1.set_ylim(0, 1.02)
    ax1.set_ylabel("Spearman ρ of variant scores vs native", color=TEXT_2)
    ax1.set_title("How much do variant scores change?", color=TEXT, loc="left", fontsize=11)

    _series(ax2, radii, [p["auroc"] for p in per], [p["auroc_per_seed"] for p in per])
    native = sweep["auroc_native"]
    _reference(ax2, native, f"native {native:.3f}", x_text)
    _reference(ax2, 0.5, "chance 0.5", x_text)
    ax2.set_ylim(0.4, 1.0)
    ax2.set_ylabel("AUROC, pathogenic vs benign", color=TEXT_2)
    ax2.set_title("Can it still tell pathogenic from benign?", color=TEXT, loc="left", fontsize=11)

    n_p, n_b = sweep["n_pathogenic"], sweep["n_benign"]
    fig.suptitle(
        f"Flank shuffle on mtDNA tRNA variants · {sweep['model']} · "
        f"{n_p} pathogenic / {n_b} benign · {sweep['window_bp']} bp windows",
        color=TEXT, fontsize=11.5, x=0.01, ha="left",
    )  # fmt: skip
    fig.text(
        0.01, 0.01,
        f"Flanks beyond the radius are dinucleotide-shuffled ({sweep['n_seeds']} seeds; faint "
        "dots). Bands: 95% bootstrap over variants. Reproduce: scripts/03_flank_sweep.py.",
        color=TEXT_2, fontsize=8, ha="left",
    )  # fmt: skip
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(out, dpi=160, facecolor=SURFACE)
    plt.close(fig)
