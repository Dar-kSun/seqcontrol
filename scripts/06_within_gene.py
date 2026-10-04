"""Plan v0.3 (docs/plan-v0.3-within-gene.md): within-gene discrimination and context.

Usage (no GPU; reads committed results/ CSVs):
    python scripts/06_within_gene.py --arm 1   # protein-coding variants, strict labels
    python scripts/06_within_gene.py --arm 2   # tRNA, widened labels (needs the
                                               # run_v03_expanded.sh results)

Arm 1: gene-prior AUROC (leave-one-variant-out and in-sample), whole-set and
within-gene AUROC of the model, and the leave-one-gene-out range of the
within-gene AUROC.
Arm 2: the same on tRNA variants with benign widened to ClinVar 1+ star, plus
within-gene AUROC and within-gene CDI under the tRNA swap and the flank shuffle.

Intervals: 2,000 cluster-bootstrap resamples over genes, redrawn until usable
(mean draws reported). A resample whose native within-gene AUROC is at or below
0.5 leaves CDI undefined; those are excluded and counted.
The pre-declared decision rules are applied in code (verdict_arm1, verdict_arm2).
Writes results/within_gene_arm<N>.json.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import platform
import subprocess
from collections import defaultdict

import numpy as np

from seqcontrol import config, metrics
from seqcontrol.genes import both_classes, fmt, gene_priors, has_pairs

RES = config.ROOT / "results"
N_BOOT = 2000
SEED = 0


def git_commit() -> str:
    out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=config.ROOT,
                         capture_output=True, text=True)  # fmt: skip
    return out.stdout.strip() or "unknown"


def load_native(name: str, region: str):
    with open(RES / name, newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["region"] == region]
    keys = [(int(r["pos"]), r["ref"], r["alt"]) for r in rows]
    y = np.array([int(r["label"]) for r in rows])
    dl = np.array([float(r["delta"]) for r in rows])
    genes = np.array([r["gene"] for r in rows])
    return keys, y, dl, genes


def verdict_arm1(within: dict) -> str:
    if within["ci95"][0] > 0.5:
        return ("the model separates pathogenic from benign variants within the same "
                "protein-coding gene")  # fmt: skip
    return "no evidence of within-gene discrimination on protein-coding variants"


def verdict_arm2(cdi: dict, native_within: dict) -> str:
    if native_within["ci95"][0] <= 0.5:
        return "not interpreted: no within-gene signal to lose (native lower bound <= 0.5)"
    lo, hi = cdi["ci95"]
    if lo > 0:
        return "within genes, the control removes a measurable part of the discrimination"
    if hi < 0.30:
        return "within-gene discrimination is largely independent of this context"
    return "indeterminate at this sample size"


def gene_analysis(y, s, genes) -> dict:
    """Whole-set, gene-prior and within-gene AUROC with gene-cluster intervals."""

    def prior(which):
        def stat(idx, grp):
            lv, ins = gene_priors(y[idx], grp)
            return metrics.auroc(y[idx], lv if which == "lovo" else ins)

        return stat

    whole, _ = metrics.cluster_bootstrap(
        lambda idx, grp: metrics.auroc(y[idx], s[idx]), genes, both_classes(y), N_BOOT, SEED
    )
    lovo, draws_b = metrics.cluster_bootstrap(prior("lovo"), genes, both_classes(y), N_BOOT, SEED)
    ins, _ = metrics.cluster_bootstrap(prior("insample"), genes, both_classes(y), N_BOOT, SEED)
    within, draws_w = metrics.cluster_bootstrap(
        lambda idx, grp: metrics.within_group_auroc(y[idx], s[idx], grp),
        genes, has_pairs(y), N_BOOT, SEED,
    )  # fmt: skip
    pairs = {str(g): int((y[genes == g] == 1).sum() * (y[genes == g] == 0).sum())
             for g in np.unique(genes)}  # fmt: skip
    logo = {}
    for g in np.unique(genes):
        keep = genes != g
        if metrics.within_group_pairs(y[keep], genes[keep]) > 0:
            logo[str(g)] = metrics.within_group_auroc(y[keep], s[keep], genes[keep]) - within.value
    return {
        "n_variants": int(len(y)),
        "n_pathogenic": int(y.sum()),
        "n_genes": int(len(np.unique(genes))),
        "within_gene_pairs": sum(pairs.values()),
        "pairs_by_gene": {g: n for g, n in pairs.items() if n},
        "auroc_whole_set": whole.as_dict(),
        "gene_prior_auroc_leave_one_variant_out": lovo.as_dict(),
        "gene_prior_auroc_in_sample_upper_bound": ins.as_dict(),
        "within_gene_auroc": within.as_dict(),
        "within_gene_auroc_change_when_gene_dropped": logo,
        "mean_draws_per_resample": {"both_classes": draws_b, "within_gene_pairs": draws_w},
    }


def control_dl(csv_name: str, keys: list[tuple], select) -> np.ndarray:
    """Per-variant control dL averaged over the selected settings/seeds."""
    per = defaultdict(list)
    with open(RES / csv_name, newline="") as f:
        for r in csv.DictReader(f):
            if select(r):
                per[(int(r["pos"]), r["ref"], r["alt"])].append(float(r["delta"]))
    missing = [k for k in keys if k not in per]
    if missing:
        raise ValueError(f"{csv_name} lacks {len(missing)} variants, e.g. {missing[0]}")
    return np.array([np.mean(per[k]) for k in keys])


def within_cdi(y, s_native, s_control, genes) -> dict:
    """Within-gene AUROC under a control and the within-gene CDI, paired, gene clusters."""
    values, draws, undefined = [], [], 0
    ctl_vals = []
    for r in metrics.cluster_resamples(genes, has_pairs(y), N_BOOT, SEED):
        nat = metrics.within_group_auroc(y[r.idx], s_native[r.idx], r.groups)
        ctl = metrics.within_group_auroc(y[r.idx], s_control[r.idx], r.groups)
        ctl_vals.append(ctl)
        draws.append(r.draws)
        if nat <= 0.5:
            undefined += 1
            continue
        values.append(metrics.cdi(nat, ctl))
    nat0 = metrics.within_group_auroc(y, s_native, genes)
    ctl0 = metrics.within_group_auroc(y, s_control, genes)
    return {
        "within_gene_auroc_control": {
            "value": ctl0,
            "ci95": np.percentile(ctl_vals, [2.5, 97.5]).tolist(),
            "n_boot": N_BOOT,
        },  # fmt: skip
        "within_gene_cdi": {
            "value": metrics.cdi(nat0, ctl0),
            "ci95": np.percentile(values, [2.5, 97.5]).tolist(),
            "n_boot": len(values),
            "resamples_cdi_undefined": undefined,
        },  # fmt: skip
        "mean_draws_per_resample": float(np.mean(draws)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--arm", type=int, choices=[1, 2], required=True)
    args = parser.parse_args()

    if args.arm == 1:
        keys, y, dl, genes = load_native("baseline_fp8-delayed.csv", "protein_coding")
        res = gene_analysis(y, -dl, genes)
        res["verdict"] = verdict_arm1(res["within_gene_auroc"])
        label = "protein-coding variants, strict labels"
    else:
        keys, y, dl, genes = load_native("baseline_fp8-delayed_expanded.csv", "Mt_tRNA")
        s = -dl
        res = gene_analysis(y, s, genes)
        res["controls"] = {}
        swap = control_dl("permutation_trna-swap_expanded.csv", keys, lambda r: r["setting"] != "0")
        controls = {"tRNA swap": swap}
        for rad in ("0", "100", "400"):
            controls[f"flank shuffle r={rad}"] = control_dl(
                "flank_sweep_expanded.csv", keys, lambda r, rad=rad: r["radius"] == rad
            )
        for name, c_dl in controls.items():
            c = within_cdi(y, s, -c_dl, genes)
            whole_ctl = metrics.auroc(y, -c_dl)
            c["auroc_whole_set_control"] = whole_ctl
            c["cdi_whole_set"] = metrics.cdi(res["auroc_whole_set"]["value"], whole_ctl)
            pre_declared = name in ("tRNA swap", "flank shuffle r=0")
            c["verdict"] = (verdict_arm2(c["within_gene_cdi"], res["within_gene_auroc"])
                            if pre_declared else "descriptive (no pre-declared rule)")  # fmt: skip
            res["controls"][name] = c
        label = "tRNA variants, widened labels (ClinVar benign 1+ star)"

    res.update({
        "plan": "docs/plan-v0.3-within-gene.md",
        "arm": args.arm,
        "description": label,
        "model": "evo2_1b_base (fp8-delayed)",
        "bootstrap": f"{N_BOOT} cluster resamples over genes, redrawn until usable; seed {SEED}",
        "analysis_commit": git_commit(),
        "run_date": dt.date.today().isoformat(),
        "versions": {"python": platform.python_version(), "numpy": np.__version__},
    })  # fmt: skip
    out = RES / f"within_gene_arm{args.arm}.json"
    out.write_text(json.dumps(res, indent=2) + "\n", newline="\n")

    print(f"Arm {args.arm}: {label}")
    print(f"  {res['n_variants']} variants ({res['n_pathogenic']} P) in {res['n_genes']} genes; "
          f"{res['within_gene_pairs']} within-gene pairs")  # fmt: skip
    print(f"  model AUROC, whole set          {fmt(res['auroc_whole_set'])}")
    print(f"  gene prior, leave-one-out       {fmt(res['gene_prior_auroc_leave_one_variant_out'])}")
    print(f"  gene prior, in-sample           {fmt(res['gene_prior_auroc_in_sample_upper_bound'])}")
    print(f"  model AUROC, within-gene pairs  {fmt(res['within_gene_auroc'])}")
    logo = res["within_gene_auroc_change_when_gene_dropped"]
    worst = max(logo, key=lambda g: abs(logo[g]))
    print(f"  one gene dropped: within-gene AUROC changes {min(logo.values()):+.3f} to "
          f"{max(logo.values()):+.3f} (largest: {worst})")  # fmt: skip
    if args.arm == 1:
        print(f"  -> {res['verdict']}")
    for name, c in res.get("controls", {}).items():
        print(f"  {name:20s} within-gene AUROC {fmt(c['within_gene_auroc_control'])}  "
              f"within-gene CDI {fmt(c['within_gene_cdi'], '.2f')}  "
              f"(whole-set AUROC {c['auroc_whole_set_control']:.3f})")  # fmt: skip
        print(f"  {'':20s} -> {c['verdict']}")
    print(f"wrote {out.relative_to(config.ROOT)}")


if __name__ == "__main__":
    main()
