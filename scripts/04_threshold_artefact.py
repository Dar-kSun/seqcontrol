"""Analysis A (docs/plan-v0.2-threshold-and-gene.md): is the collapse a threshold artefact?

Usage (no GPU; uses committed results/ CSVs):
    pip install -e ".[plot]"   # for the figure
    python scripts/04_threshold_artefact.py

A0  AUROC is unchanged by any scaling dL -> alpha*dL (asserted, not a result).
A1  sensitivity/specificity at fixed cut-offs as native scores are compressed.
A2  how much the real tRNA swap compressed scores: median-ratio and origin-fit alpha.
A3  pure-compression null vs the observed swap, and the share of the sensitivity
    drop that compression explains, judged by the pre-declared rule.
A4  rank-preserving surrogate: the swap's values in native order.
A5  simulation of the paper's operating point (synthetic Gaussians, NOT the paper's
    data and NOT a reproduction).

Choices fixed before running (not specified by the plan):
- Swap scores are pooled over all 19 shifts (shift x variant), not averaged per
  variant first; averaging would itself shrink scores and inflate compression.
- The A3 rule is applied with both alpha estimates; if they disagree, the verdict is
  "indeterminate".
- A5's free scale (pathogenic SD) is set to the SD of this repo's native pathogenic
  tRNA dL, and the result is also reported at half and double that.

Writes results/threshold_artefact.{csv,json,png}.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import platform
import subprocess
from collections import defaultdict

import numpy as np

from seqcontrol import artefact, config, metrics

RES = config.ROOT / "results"
T_PAPER = -0.0081
SENS_PAPER, SPEC_PAPER, SENS_PAPER_SWAP, SPEC_PAPER_SWAP = 0.658, 0.785, 0.051, 0.938
ALPHAS = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
N_BOOT = 2000
SEED = 0


def git_commit() -> str:
    out = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], cwd=config.ROOT, capture_output=True, text=True
    )
    return out.stdout.strip() or "unknown"


def load() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Labels, native dL (variants) and swap dL (shifts x variants), aligned."""
    with open(RES / "baseline_fp8-delayed.csv", newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["region"] == "Mt_tRNA"]
    keys = [(int(r["pos"]), r["ref"], r["alt"]) for r in rows]
    y = np.array([int(r["label"]) for r in rows])
    native = np.array([float(r["delta"]) for r in rows])
    by_shift = defaultdict(dict)
    with open(RES / "permutation_trna-swap.csv", newline="") as f:
        for r in csv.DictReader(f):
            by_shift[int(r["setting"])][(int(r["pos"]), r["ref"], r["alt"])] = float(r["delta"])
    if not np.array_equal(np.array([by_shift[0][k] for k in keys]), native):
        raise AssertionError("swap file's unperturbed setting differs from the baseline")
    swap = np.array([[by_shift[s][k] for k in keys] for s in sorted(by_shift) if s != 0])
    return y, native, swap


def swap_sens_spec(y, swap, t) -> tuple[float, float]:
    """Mean over shifts, as in scripts/02_permutation.py."""
    ss = np.array([artefact.sens_spec(y, row, t) for row in swap])
    return float(ss[:, 0].mean()), float(ss[:, 1].mean())


def alphas_hat(native, swap) -> tuple[float, float, float]:
    pooled_native = np.tile(native, len(swap))
    ratio = artefact.scale_ratio(native, swap.ravel())
    fit, r2 = artefact.origin_fit(pooled_native, swap.ravel())
    return ratio, fit, r2


def main() -> None:
    y, native, swap = load()
    t_youden = -metrics.youden_threshold(y, -native)
    cutoffs = {"native_youden": t_youden, "paper": T_PAPER}
    print(
        f"{len(y)} tRNA variants ({y.sum()} P / {len(y) - y.sum()} B), {len(swap)} shifts; "
        f"native-Youden cut-off dL <= {t_youden:.4f}"
    )

    # A0
    fine = np.round(np.linspace(1.0, 0.05, 96), 4)
    artefact.check_auroc_invariance(y, native, [*ALPHAS, *fine])
    auroc_native = metrics.auroc(y, -native)
    print(f"A0 AUROC identical (|diff| < 1e-12) under all {len(ALPHAS) + len(fine)} alphas")

    # A1
    a1 = []
    for a in sorted({*ALPHAS, *fine}, reverse=True):
        row = {"alpha": a, "auroc": metrics.auroc(y, -a * native)}
        for name, t in cutoffs.items():
            row[f"sens_{name}"], row[f"spec_{name}"] = artefact.sens_spec(y, a * native, t)
        a1.append(row)
    obs = {name: swap_sens_spec(y, swap, t) for name, t in cutoffs.items()}
    reached = [r["alpha"] for r in a1 if r["sens_native_youden"] <= obs["native_youden"][0]]
    alpha_matching = max(reached) if reached else None

    # A2
    ratio, fit, r2 = alphas_hat(native, swap)
    cls_ratio = {
        c: artefact.scale_ratio(native[y == lab], swap[:, y == lab].ravel())
        for c, lab in [("pathogenic", 1), ("benign", 0)]
    }
    boots = np.array(
        [alphas_hat(native[i], swap[:, i]) for i in metrics.stratified_resamples(y, N_BOOT, SEED)]
    )
    ci = {
        k: np.percentile(boots[:, j], [2.5, 97.5]).tolist()
        for j, k in enumerate(["ratio", "fit", "r2"])
    }
    a2 = {
        "alpha_median_ratio": {"value": ratio, "ci95": ci["ratio"], "n_boot": N_BOOT},
        "alpha_origin_fit": {"value": fit, "ci95": ci["fit"], "n_boot": N_BOOT},
        "origin_fit_r2": {"value": r2, "ci95": ci["r2"], "n_boot": N_BOOT},
        "class_median_ratio_descriptive": cls_ratio,
        "spearman_native_vs_swap_mean_over_shifts": float(
            np.mean([metrics.spearman(native, row) for row in swap])
        ),
    }

    # A3
    t = t_youden
    sens_nat, spec_nat = artefact.sens_spec(y, native, t)
    sens_obs, spec_obs = obs["native_youden"]
    swap_json = json.loads((RES / "permutation_trna-swap.json").read_text())
    drop = swap_json["auroc_drop"]
    drop_includes_zero = drop["ci95"][0] <= 0 <= drop["ci95"][1]
    a3 = {
        "cut_off_dL": t,
        "native": {"sensitivity": sens_nat, "specificity": spec_nat, "auroc": auroc_native},
        "observed_swap": {
            "sensitivity": sens_obs,
            "specificity": spec_obs,
            "auroc": swap_json["auroc_control"]["value"],
        },
        "auroc_drop_variant_bootstrap": drop,
        "nulls": {},
    }
    verdicts = []
    for name, a in [("median_ratio", ratio), ("origin_fit", fit)]:
        sens_null, spec_null = artefact.sens_spec(y, a * native, t)
        share = artefact.share_explained(sens_nat, sens_null, sens_obs)
        shares, skipped = [], 0
        for i in metrics.stratified_resamples(y, N_BOOT, SEED):
            ab = alphas_hat(native[i], swap[:, i])[0 if name == "median_ratio" else 1]
            sn = artefact.sens_spec(y[i], native[i], t)[0]
            so = swap_sens_spec(y[i], swap[:, i], t)[0]
            if sn == so:
                skipped += 1  # no observed drop in this resample: share undefined
                continue
            shares.append(
                artefact.share_explained(sn, artefact.sens_spec(y[i], ab * native[i], t)[0], so)
            )
        verdict = artefact.verdict_a3(share, drop_includes_zero)
        verdicts.append(verdict)
        a3["nulls"][name] = {
            "alpha": a,
            "sensitivity": sens_null,
            "specificity": spec_null,
            "auroc": auroc_native,
            "share_of_sensitivity_drop_explained": {
                "value": share,
                "ci95": np.percentile(shares, [2.5, 97.5]).tolist(),
                "n_boot": N_BOOT - skipped,
                "resamples_skipped_no_drop": skipped,
            },
            "verdict": verdict,
        }
    a3["verdict"] = verdicts[0] if len(set(verdicts)) == 1 else "indeterminate at this sample size"

    # A4
    a4 = {}
    for name, tt in cutoffs.items():
        sur = np.array(
            [
                artefact.sens_spec(y, artefact.rank_preserving_surrogate(native, row), tt)
                for row in swap
            ]
        )
        a4[name] = {
            "surrogate_sensitivity": float(sur[:, 0].mean()),
            "surrogate_specificity": float(sur[:, 1].mean()),
            "observed_sensitivity": obs[name][0],
            "observed_specificity": obs[name][1],
        }

    # A5 (simulation)
    sd_p = float(np.std(native[y == 1], ddof=1))
    a5 = {
        "label": "SIMULATION of the mechanism at the paper's reported operating point, "
        "using synthetic Gaussian scores; not the paper's data, not a reproduction",
        "inputs": {
            "cut_off_dL": T_PAPER,
            "native_sensitivity": SENS_PAPER,
            "native_specificity": SPEC_PAPER,
            "target_sensitivity": SENS_PAPER_SWAP,
            "paper_reported_specificity_after": SPEC_PAPER_SWAP,
        },
        "by_scale": {},
    }
    for mult in (0.5, 1.0, 2.0):
        b = artefact.Binormal.calibrated(T_PAPER, SENS_PAPER, SPEC_PAPER, sd_p * mult)
        alpha = b.alpha_for_sensitivity(T_PAPER, SENS_PAPER_SWAP)
        sens_a, spec_a = b.sens_spec(T_PAPER, alpha)
        a5["by_scale"][f"sd_p x{mult}"] = {
            "sd_pathogenic": sd_p * mult,
            "alpha_for_5.1pct_sensitivity": alpha,
            "specificity_at_that_alpha": spec_a,
            "sensitivity_check": sens_a,
            "auroc_before_and_after": b.auroc(),
        }

    results = {
        "analysis": "A, threshold artefact (docs/plan-v0.2-threshold-and-gene.md)",
        "model": swap_json["model"],
        "n_pathogenic": int(y.sum()),
        "n_benign": int(len(y) - y.sum()),
        "n_shifts": len(swap),
        "cut_offs_dL": cutoffs,
        "A1_alpha_matching_observed_swap_sensitivity": alpha_matching,
        "A1_observed_swap": {k: {"sensitivity": v[0], "specificity": v[1]} for k, v in obs.items()},
        "A2": a2,
        "A3": a3,
        "A4": a4,
        "A5_simulation": a5,
        "bootstrap": f"{N_BOOT} label-stratified resamples of variants; seed {SEED}",
        "analysis_commit": git_commit(),
        "run_date": dt.date.today().isoformat(),
        "versions": {"python": platform.python_version(), "numpy": np.__version__},
    }
    (RES / "threshold_artefact.json").write_text(json.dumps(results, indent=2) + "\n", newline="\n")
    with open(RES / "threshold_artefact.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(a1[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(a1)

    def ci(d, spec=".3f"):
        return f"{d['value']:{spec}} [{d['ci95'][0]:{spec}}, {d['ci95'][1]:{spec}}]"

    print(
        f"A1 alpha at which pure compression reaches the swap's sensitivity "
        f"({sens_obs:.3f}): {alpha_matching}"
    )
    print(
        f"A2 alpha (median ratio) {ci(a2['alpha_median_ratio'])}; (origin fit) "
        f"{ci(a2['alpha_origin_fit'])}; R^2 {ci(a2['origin_fit_r2'])}"
    )
    p, b = cls_ratio["pathogenic"], cls_ratio["benign"]
    print(f"   class ratios: pathogenic {p:.3f}, benign {b:.3f}")
    print(f"A3 at dL <= {t:.4f}:  sensitivity / specificity / AUROC")
    print(f"   native            {sens_nat:.3f} / {spec_nat:.3f} / {auroc_native:.3f}")
    for name, n in a3["nulls"].items():
        row = f"{n['sensitivity']:.3f} / {n['specificity']:.3f} / {n['auroc']:.3f}"
        share = ci(n["share_of_sensitivity_drop_explained"], ".2f")
        print(f"   null x{n['alpha']:.3f} ({name:12s}) {row}   share explained {share}")
    print(
        f"   observed swap     {sens_obs:.3f} / {spec_obs:.3f} / {a3['observed_swap']['auroc']:.3f}"
    )
    print(f"   -> {a3['verdict']}")
    for name, v in a4.items():
        print(
            f"A4 {name:13s} surrogate sens {v['surrogate_sensitivity']:.3f} vs observed "
            f"{v['observed_sensitivity']:.3f}; spec {v['surrogate_specificity']:.3f} vs "
            f"{v['observed_specificity']:.3f}"
        )
    for k, v in a5["by_scale"].items():
        print(
            f"A5 (simulation, {k}) alpha {v['alpha_for_5.1pct_sensitivity']:.3f} gives sens "
            f"{v['sensitivity_check']:.3f}, spec {v['specificity_at_that_alpha']:.3f}; "
            f"AUROC {v['auroc_before_and_after']:.3f} throughout"
        )

    try:
        figure(results, a1, RES / "threshold_artefact.png")
        print("wrote results/threshold_artefact.{csv,json,png}")
    except ImportError:
        print("wrote results/threshold_artefact.{csv,json} (install the plot extra for the figure)")


def figure(results: dict, a1: list[dict], out) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from seqcontrol.plots import GRID, MUTED, SURFACE, TEXT, TEXT_2

    blue, orange, aqua = "#2a78d6", "#eb6834", "#1baf7a"  # reference palette slots 1-3
    alpha = [r["alpha"] for r in a1]
    fig, ax = plt.subplots(figsize=(9, 4.6), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    series = [
        ("sens_native_youden", "sensitivity", blue),
        ("spec_native_youden", "specificity", orange),
        ("auroc", "AUROC", aqua),
    ]
    for key, label, color in series:
        vals = [r[key] for r in a1]
        ax.plot(
            alpha,
            vals,
            color=color,
            linewidth=2,
            drawstyle="steps-post" if key != "auroc" else "default",
        )
        y_lab = max(vals[-1], 0.03)  # keep the label off the x-axis
        ax.text(0.04, y_lab, f" {label}", color=TEXT, fontsize=9, va="center", ha="left")
    a3 = results["A3"]
    a_hat = results["A2"]["alpha_median_ratio"]["value"]
    ax.axvline(a_hat, color=MUTED, linestyle=(0, (4, 3)), linewidth=1.2)
    ax.text(
        a_hat, 1.03, f"swap's compression α̂ = {a_hat:.2f}", color=TEXT_2, fontsize=8.5, ha="center"
    )
    obs = a3["observed_swap"]
    for val, color in [
        (obs["sensitivity"], blue),
        (obs["specificity"], orange),
        (obs["auroc"], aqua),
    ]:
        ax.scatter(
            [a_hat], [val], s=80, facecolors=SURFACE, edgecolors=color, linewidths=2, zorder=5
        )
    ax.text(
        a_hat + 0.02,
        obs["sensitivity"],
        "observed tRNA swap (open circles)",
        color=TEXT_2,
        fontsize=8.5,
        va="center",
        ha="right",
    )
    ax.set_xlim(1.02, 0.0)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel(
        "Compression α (every native score multiplied by α)  →  more compression", color=TEXT_2
    )
    ax.set_ylabel(f"At fixed cut-off ΔL ≤ {a3['cut_off_dL']:.4f}", color=TEXT_2)
    ax.tick_params(colors=TEXT_2, labelsize=9)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.set_title(
        "Lines: pure compression at a fixed cut-off  ·  Circles: the observed tRNA swap",
        color=TEXT,
        fontsize=10.5,
        loc="left",
    )
    fig.text(
        0.01,
        0.01,
        f"{results['model']}, {results['n_pathogenic']} pathogenic / {results['n_benign']} benign "
        "tRNA variants.\nLines: native scores × α (AUROC identical by construction). "
        "Reproduce: scripts/04_threshold_artefact.py.",
        color=TEXT_2,
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out, dpi=160, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    main()
