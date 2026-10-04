"""Analysis B (docs/plan-v0.2-threshold-and-gene.md): does gene identity explain the AUROC?

Usage (no GPU; needs the committed results/ CSVs and data/raw/ for allele frequencies):
    python scripts/05_gene_confound.py

B1  gene-prior baseline: score each tRNA variant by its gene's pathogenic fraction,
    leave-one-variant-out (primary, biased low) and in-sample (upper bound).
B2  within-gene AUROC: only pathogenic-benign pairs from the same gene.
B3  leave-one-gene-out: native AUROC with each gene dropped in turn.
B4  within-gene AUROC under each control (underpowered, reported as such).
B5  within benign variants, Spearman of dL with population allele frequency.

Intervals: 2,000 cluster-bootstrap resamples over genes (B1-B4), redrawn until
usable, with the mean number of draws reported; B5 resamples variants.
Pre-declared decision rules are applied in code; see verdict_b1 / verdict_b2.
Writes results/gene_confound.csv and results/gene_confound.json.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import platform
import subprocess
from collections import Counter, defaultdict

import numpy as np

from seqcontrol import config, metrics
from seqcontrol.data import mitomap

N_BOOT = 2000
SEED = 0
RES = config.ROOT / "results"


def git_commit() -> str:
    out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=config.ROOT,
                         capture_output=True, text=True)  # fmt: skip
    return out.stdout.strip() or "unknown"


def read_csv(name: str) -> list[dict]:
    with open(RES / name, newline="") as f:
        return list(csv.DictReader(f))


def gene_priors(y: np.ndarray, groups: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(leave-one-variant-out, in-sample) pathogenic fraction of each variant's gene."""
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


def fmt(d: dict, spec: str = ".3f") -> str:
    return f"{d['value']:{spec}} [{d['ci95'][0]:{spec}}, {d['ci95'][1]:{spec}}]"


def both_classes(y):
    return lambda idx, grp: 0 < y[idx].sum() < len(idx)


def has_pairs(y):
    return lambda idx, grp: metrics.within_group_pairs(y[idx], grp) > 0


def verdict_b1(auroc_lovo: float) -> str:
    if auroc_lovo >= 0.75:
        return "native AUROC is largely gene identity"
    if auroc_lovo > 0.60:
        return "gene identity is a substantial part of the native AUROC"
    return "gene identity is a minor part; the headline stands as a statement about variants"


def verdict_b2(n_pairs: int) -> str:
    if n_pairs < 100:
        return (
            "fewer than 100 within-gene pairs: this dataset cannot separate gene identity "
            "from variant effect"
        )
    return "enough within-gene pairs to interpret the within-gene AUROC"  # pragma: no cover


def control_scores(keys: list[tuple]) -> dict[str, np.ndarray]:
    """Per-variant control dL, averaged over settings/seeds, aligned to `keys`."""
    out = {}
    for name, csv_name, select in [
        ("tRNA swap", "permutation_trna-swap.csv", lambda r: r["setting"] != "0"),
        ("window rotation", "permutation_window-rotation.csv", lambda r: r["setting"] != "0"),
        ("flank shuffle r=0", "flank_sweep.csv", lambda r: r["radius"] == "0"),
        ("flank shuffle r=100", "flank_sweep.csv", lambda r: r["radius"] == "100"),
        ("flank shuffle r=400", "flank_sweep.csv", lambda r: r["radius"] == "400"),
    ]:
        per = defaultdict(list)
        for r in read_csv(csv_name):
            if select(r):
                per[(int(r["pos"]), r["ref"], r["alt"])].append(float(r["delta"]))
        out[name] = np.array([np.mean(per[k]) for k in keys])
    return out


def main() -> None:
    rows = [r for r in read_csv("baseline_fp8-delayed.csv") if r["region"] == "Mt_tRNA"]
    keys = [(int(r["pos"]), r["ref"], r["alt"]) for r in rows]
    y = np.array([int(r["label"]) for r in rows])
    dl = np.array([float(r["delta"]) for r in rows])
    genes = np.array([r["gene"] for r in rows])
    s = -dl  # pathogenicity score
    print(f"{len(y)} tRNA variants ({y.sum()} P / {len(y) - y.sum()} B) in "
          f"{len(np.unique(genes))} genes")  # fmt: skip

    # --- B1: gene-prior baseline
    lovo, insample = gene_priors(y, genes)

    def prior_auroc(which):
        def stat(idx, grp):
            lv, ins = gene_priors(y[idx], grp)
            return metrics.auroc(y[idx], lv if which == "lovo" else ins)

        return stat

    native_auc, draws_native = metrics.cluster_bootstrap(
        lambda idx, grp: metrics.auroc(y[idx], s[idx]), genes, both_classes(y), N_BOOT, SEED
    )
    lovo_auc, draws_b1 = metrics.cluster_bootstrap(
        prior_auroc("lovo"), genes, both_classes(y), N_BOOT, SEED
    )
    ins_auc, _ = metrics.cluster_bootstrap(
        prior_auroc("insample"), genes, both_classes(y), N_BOOT, SEED
    )

    # --- B2: within-gene AUROC
    pairs_by_gene = {g: int((y[genes == g] == 1).sum() * (y[genes == g] == 0).sum())
                     for g in np.unique(genes)}  # fmt: skip
    n_pairs = sum(pairs_by_gene.values())
    within, draws_b2 = metrics.cluster_bootstrap(
        lambda idx, grp: metrics.within_group_auroc(y[idx], s[idx], grp),
        genes, has_pairs(y), N_BOOT, SEED,
    )  # fmt: skip

    # --- B3: leave-one-gene-out
    logo = {}
    for g in np.unique(genes):
        keep = genes != g
        if 0 < y[keep].sum() < keep.sum():
            logo[g] = metrics.auroc(y[keep], s[keep]) - native_auc.value
    most = max(logo, key=lambda g: abs(logo[g]))
    no_tl1 = genes != "MT-TL1"
    tl1_auc, _ = metrics.cluster_bootstrap(
        lambda idx, grp: metrics.auroc(y[no_tl1][idx], s[no_tl1][idx]),
        genes[no_tl1], both_classes(y[no_tl1]), N_BOOT, SEED,
    )  # fmt: skip

    # --- B4: within-gene AUROC under controls (underpowered)
    b4 = {}
    for name, c_dl in control_scores(keys).items():
        c = -c_dl
        ctl, _ = metrics.cluster_bootstrap(
            lambda idx, grp, c=c: metrics.within_group_auroc(y[idx], c[idx], grp),
            genes, has_pairs(y), N_BOOT, SEED,
        )  # fmt: skip
        drop, _ = metrics.cluster_bootstrap(
            lambda idx, grp, c=c: metrics.within_group_auroc(y[idx], s[idx], grp)
            - metrics.within_group_auroc(y[idx], c[idx], grp),
            genes, has_pairs(y), N_BOOT, SEED,
        )  # fmt: skip
        b4[name] = {"within_gene_auroc": ctl.as_dict(), "drop_vs_native": drop.as_dict()}

    # --- B5: familiarity, within benign variants
    af = mitomap.allele_frequencies(config.DATA_RAW / "mitomap_polymorphisms.vcf")
    all_rows = read_csv("baseline_fp8-delayed.csv")
    b5 = {}
    for scope, sel in [
        ("all benign", lambda r: True),
        ("tRNA benign", lambda r: r["region"] == "Mt_tRNA"),
    ]:
        ben = [r for r in all_rows if r["label"] == "0" and sel(r)]
        have = [r for r in ben if (int(r["pos"]), r["ref"], r["alt"]) in af]
        d = np.array([float(r["delta"]) for r in have])
        f = np.array([af[(int(r["pos"]), r["ref"], r["alt"])] for r in have])
        rng = np.random.default_rng(SEED)
        boots = [metrics.spearman(d[i], f[i])
                 for i in (rng.integers(0, len(d), len(d)) for _ in range(N_BOOT))]  # fmt: skip
        lo, hi = np.nanpercentile(boots, [2.5, 97.5])
        b5[scope] = {
            "n_benign": len(ben),
            "n_with_allele_frequency": len(have),
            "spearman_dL_vs_af": metrics.Estimate(metrics.spearman(d, f), lo, hi, N_BOOT).as_dict(),
            "median_af": float(np.median(f)),
            "fraction_af_ge_1pct": float(np.mean(f >= 0.01)),
        }

    results = {
        "analysis": "B, gene-identity confound (docs/plan-v0.2-threshold-and-gene.md)",
        "model": "evo2_1b_base (fp8-delayed)",
        "variants": {
            "n": int(len(y)),
            "pathogenic": int(y.sum()),
            "genes": int(len(np.unique(genes))),
            "per_gene": {g: dict(Counter(int(v) for v in y[genes == g])) for g in np.unique(genes)},
        },  # fmt: skip
        "native_auroc_gene_cluster_ci": native_auc.as_dict(),
        "B1": {
            "gene_prior_auroc_leave_one_variant_out": lovo_auc.as_dict(),
            "gene_prior_auroc_in_sample_upper_bound": ins_auc.as_dict(),
            "mean_draws_per_resample": draws_b1,
            "verdict": verdict_b1(lovo_auc.value),
        },
        "B2": {
            "within_gene_pairs": n_pairs,
            "total_pairs": int(y.sum() * (len(y) - y.sum())),
            "pairs_by_gene": {g: n for g, n in pairs_by_gene.items() if n},
            "within_gene_auroc": within.as_dict(),
            "mean_draws_per_resample": draws_b2,
            "verdict": verdict_b2(n_pairs),
        },
        "B3": {
            "auroc_change_when_gene_dropped": logo,
            "range": [min(logo.values()), max(logo.values())],
            "most_influential_gene": most,
            "auroc_without_MT-TL1": tl1_auc.as_dict(),
            "MT-TL1_moves_auroc_by_more_than_0.05": abs(logo["MT-TL1"]) > 0.05,
        },
        "B4_underpowered": b4,
        "B5": b5,
        "bootstrap": f"{N_BOOT} cluster resamples over genes (B1-B4), redrawn until usable; "
        f"B5 over variants; seed {SEED}",
        "analysis_commit": git_commit(),
        "run_date": dt.date.today().isoformat(),
        "versions": {"python": platform.python_version(), "numpy": np.__version__},
    }
    (RES / "gene_confound.json").write_text(json.dumps(results, indent=2) + "\n", newline="\n")
    with open(RES / "gene_confound.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(
            ["pos", "ref", "alt", "label", "gene", "delta", "prior_lovo", "prior_in_sample", "af"]
        )
        for k, lab, g, x, a, b in zip(keys, y, genes, dl, lovo, insample, strict=True):
            w.writerow([*k, lab, g, x, a, b, af.get(k, "")])

    print(f"\nnative AUROC (gene-cluster CI)       {native_auc}")
    print(f"B1 gene-prior AUROC, leave-one-out   {lovo_auc}   (draws/resample {draws_b1:.2f})")
    print(f"B1 gene-prior AUROC, in-sample       {ins_auc}")
    print(f"   -> {results['B1']['verdict']}")
    print(f"B2 within-gene pairs {n_pairs} of {results['B2']['total_pairs']}; "
          f"by gene {results['B2']['pairs_by_gene']}")  # fmt: skip
    print(f"B2 within-gene AUROC                 {within}   (draws/resample {draws_b2:.2f})")
    print(f"   -> {results['B2']['verdict']}")
    lo3, hi3 = results["B3"]["range"]
    print(f"B3 AUROC change with one gene dropped: {lo3:+.3f} to {hi3:+.3f}")
    print(f"   most influential {most} ({logo[most]:+.3f}); MT-TL1 {logo['MT-TL1']:+.3f}")
    print(f"   AUROC without MT-TL1               {tl1_auc}")
    for name, v in b4.items():
        a, d = fmt(v["within_gene_auroc"]), fmt(v["drop_vs_native"], "+.3f")
        print(f"B4 {name:20s} within-gene AUROC {a}  drop {d}")
    for scope, v in b5.items():
        n = f"{v['n_with_allele_frequency']}/{v['n_benign']}"
        print(f"B5 {scope}: n {n} with AF; Spearman(dL, AF) {fmt(v['spearman_dL_vs_af'])}; "
              f"median AF {v['median_af']:.4f}")  # fmt: skip
    print("\nwrote results/gene_confound.csv and results/gene_confound.json")


if __name__ == "__main__":
    main()
