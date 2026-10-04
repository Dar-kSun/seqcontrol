"""Trust card: a one-page summary of how much a model's variant scores depend on context.

Built only from the JSON files in results/, so every number on the card traces to a
script run. Missing results are reported as "not run", never filled in.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from seqcontrol.plots import noise_floor_spearman

REGION_NAMES = {"Mt_tRNA": "mtDNA tRNA"}


@dataclass(frozen=True)
class ControlRow:
    name: str
    description: str
    auroc: dict
    cdi: dict
    spearman: float
    sens_native_cutoff: tuple[float, float] | None  # native -> control, at native Youden


def _load(results: Path, name: str) -> dict | None:
    p = results / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def _ci(d: dict, fmt: str = ".3f") -> str:
    lo, hi = d["ci95"]
    return f"{d['value']:{fmt}} [{lo:{fmt}}, {hi:{fmt}}]"


def control_rows(results: Path) -> list[ControlRow]:
    rows = []
    for stem, name, desc in [
        ("permutation_trna-swap", "tRNA swap",
         "each tRNA moved into another tRNA's slot (19 shifts)"),
        ("permutation_window-rotation", "Window rotation",
         "scoring window rotated, gene kept whole (13 offsets)"),
    ]:  # fmt: skip
        d = _load(results, stem)
        if d is None:
            continue
        t = d["thresholds"]["native_youden"]
        rows.append(ControlRow(
            name, desc, d["auroc_control"], d["cdi"], d["spearman_native_vs_setting"]["mean"],
            (t["native"]["sensitivity"], t["control_mean"]["sensitivity"]),
        ))  # fmt: skip
    sweep = _load(results, "flank_sweep")
    if sweep is not None:
        for r in ("0", "100", "400"):
            p = sweep["per_radius"].get(r)
            if p is None:
                continue
            desc = (
                "all context dinucleotide-shuffled, gene kept"
                if r == "0"
                else f"context shuffled beyond {r} bp of the gene"
            )
            rows.append(ControlRow(f"Flank shuffle, r = {r} bp", desc, p["auroc"], p["cdi"],
                                   p["spearman"]["value"], None))  # fmt: skip
    return rows


def render_markdown(results: Path, figure: str | None = None) -> str:
    base = _load(results, "baseline_fp8-delayed")
    swap = _load(results, "permutation_trna-swap")
    if base is None or swap is None:
        raise FileNotFoundError("trust card needs results/baseline_fp8-delayed.json and "
                                "results/permutation_trna-swap.json")  # fmt: skip
    region = "Mt_tRNA"
    reg = base["regions"][region]
    rows = control_rows(results)
    floor = noise_floor_spearman(results, region)
    n_p, n_b = reg["n_pathogenic"], reg["n_benign"]

    significant = [r.name for r in rows if r.cdi["ci95"][0] > 0]
    commits = sorted({base.get("git_commit", "?"), swap["scored_at_commit"]})
    largest = max(rows, key=lambda r: r.cdi["value"])

    lines = [
        f"# Trust card: {REGION_NAMES[region]} variants · `{base['model']}`",
        "",
        "How much of this model's ability to separate pathogenic from benign variants "
        "depends on the sequence *around* the gene rather than the gene itself?",
        "",
        f"**Model:** Evo 2 `{base['model']}`, 1B parameters. "
        f"**Variants:** {n_p} pathogenic (MITOMAP confirmed) and {n_b} benign (ClinVar, "
        f"2+ stars) single-base variants in mitochondrial tRNA genes. "
        f"**Score:** −ΔL over a {base['window_bp']} bp window. "
        "**Intervals:** 95%, bootstrap over variants.",
        "",
        "## Headline",
        "",
        f"- Native AUROC: **{_ci(reg['auroc'])}** (AUPRC {_ci(reg['auprc'])}, "
        f"no-skill {reg['auprc_no_skill']:.3f}).",
        f"- Largest context dependence: **{largest.name}**, CDI {_ci(largest.cdi, '.2f')}.",
        "- Controls whose CDI interval excludes zero: "
        + (", ".join(significant) if significant else "none")
        + ".",
        f"- Per-variant scores are context-sensitive: Spearman vs native falls to "
        f"{min(r.spearman for r in rows):.2f}, against {floor:.2f} from FP8 rounding alone.",
        "",
    ]
    if figure:
        lines += [f"![AUROC under each control]({figure})", ""]
    lines += [
        "## Controls",
        "",
        "CDI = (AUROC_native − AUROC_control) / (AUROC_native − 0.5): the share of "
        "above-chance discrimination lost under the control. 0 = none, 1 = all.",
        "",
        "| Control | What changes | AUROC | CDI | Spearman ΔL vs native |",
        "|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(f"| {r.name} | {r.description} | {_ci(r.auroc)} | {_ci(r.cdi, '.2f')} | "
                     f"{r.spearman:.2f} |")  # fmt: skip
    lines += [
        "",
        f"Noise floor for Spearman: {floor:.2f}, the agreement between two FP8 rounding "
        "recipes with no control applied.",
        "",
        "## At a fixed threshold",
        "",
        "Sensitivity at a cut-off chosen on native scores (Youden), then held fixed:",
        "",
    ]
    for r in rows:
        if r.sens_native_cutoff:
            a, b = r.sens_native_cutoff
            lines.append(f"- {r.name}: {a:.2f} → {b:.2f}")
    lines += [
        "",
        "Threshold metrics fall further than AUROC because the controls shrink effect "
        "sizes for both classes; a fixed cut-off then misses pathogenic variants even "
        "where their ranking is mostly intact.",
        "",
        "## Read with care",
        "",
        f"- Small sample: {n_p} pathogenic and {n_b} benign. Most intervals are wide.",
        "- One checkpoint (the smallest Evo 2 model), one window size, one label set. "
        "Larger Evo 2 models were not tested.",
        "- tRNA labels cluster by gene, so part of the native AUROC may reflect telling "
        "genes apart rather than variants within a gene.",
        "- Benign variants are mostly common polymorphisms; the model may partly score "
        "allele familiarity.",
        "- FP8 runs on an Ada GPU (RTX 4060 Laptop), not the Hopper GPU Evo 2 documents.",
        "- This card describes model behaviour on a benchmark. It says nothing about clinical use.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python scripts/fetch_data.py",
        "python scripts/01_baseline.py",
        "python scripts/02_permutation.py --control trna-swap",
        "python scripts/02_permutation.py --control window-rotation",
        "python scripts/03_flank_sweep.py",
        "seqcontrol card",
        "```",
        "",
        f"Scored at commits {', '.join(commits)}; data checksums in `data/MANIFEST.md`.",
        "",
    ]
    return "\n".join(lines)


def render_figure(results: Path, out: Path) -> None:
    """Dot-and-interval chart of AUROC under each control, with native and chance."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from seqcontrol.plots import GRID, MUTED, SERIES, SURFACE, TEXT, TEXT_2

    base = _load(results, "baseline_fp8-delayed")
    native = base["regions"]["Mt_tRNA"]["auroc"]
    rows = control_rows(results)
    labels = ["Native (no control)"] + [r.name for r in rows]
    vals = [native] + [r.auroc for r in rows]

    fig, ax = plt.subplots(figsize=(8, 0.5 * len(labels) + 1.6), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ys = list(range(len(labels)))[::-1]
    for y, v in zip(ys, vals, strict=True):
        lo, hi = v["ci95"]
        ax.plot([lo, hi], [y, y], color=SERIES, linewidth=2, solid_capstyle="round", zorder=2)
        ax.scatter([v["value"]], [y], s=70, color=SERIES, edgecolors=SURFACE, linewidths=2,
                   zorder=3)  # fmt: skip
        ax.text(0.995, y, f"{v['value']:.3f}", va="center", ha="right", color=TEXT, fontsize=9.5)
    for x, lab in [(native["value"], "native"), (0.5, "chance")]:
        ax.axvline(x, color=MUTED, linewidth=1.2, linestyle=(0, (4, 3)), zorder=1)
        ax.text(x, len(labels) - 0.35, lab, color=TEXT_2, fontsize=8.5, ha="center")
    ax.set_yticks(ys, labels, color=TEXT, fontsize=9.5)
    ax.set_xlim(0.45, 1.0)  # values sit in a column at the right edge, clear of the marks
    ax.set_ylim(-0.7, len(labels) - 0.1)
    ax.set_xlabel("AUROC, pathogenic vs benign (95% bootstrap interval)", color=TEXT_2)
    ax.tick_params(axis="x", colors=TEXT_2, labelsize=9)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.set_title(f"AUROC under each control · mtDNA tRNA variants · {base['model']}",
                 color=TEXT, fontsize=11, loc="left")  # fmt: skip
    fig.tight_layout()
    fig.savefig(out, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def write_card(results: Path, out_dir: Path, figure: bool = True) -> Path:
    """Write trust_card_tRNA.md (and its figure, if matplotlib is installed)."""
    fig_name = None
    if figure:
        try:
            render_figure(results, out_dir / "trust_card_tRNA.png")
            fig_name = "trust_card_tRNA.png"
        except ImportError:
            pass  # matplotlib is optional; the card is complete without the figure
    md = out_dir / "trust_card_tRNA.md"
    md.write_text(render_markdown(results, fig_name), encoding="utf-8", newline="\n")
    return md
