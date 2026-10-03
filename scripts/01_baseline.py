"""Native baseline: how well does Evo 2 separate pathogenic from benign mtDNA variants?

Usage (inside the WSL evo2 environment, after scripts/fetch_data.py):
    python scripts/01_baseline.py                          # Evo 2's own FP8 recipe
    python scripts/01_baseline.py --precision fp8-current  # robustness check

Scores every labelled variant in data/MANIFEST.md on a 1,025 bp window centred on it
(512 bp either side, wrapping around the circular chromosome, as in Mathur &
Sachidanandam 2026). Variant score ΔL = mean log-likelihood(alt) - (ref); the
pathogenicity score is -ΔL. Writes per-variant scores to
results/baseline_<precision>.csv and metrics with 95% bootstrap intervals to
results/baseline_<precision>.json.
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
from seqcontrol.data import mtdna
from seqcontrol.data.cache import sha256
from seqcontrol.models.evo2 import PRECISIONS, Evo2Adapter
from seqcontrol.variants import variant_windows

WINDOW = 1025
PAPER_THRESHOLD_DL = -0.0081  # Mathur & Sachidanandam 2026, Youden-optimal ΔL on their set
REGIONS = ["all", "Mt_tRNA", "protein_coding", "Mt_rRNA"]
N_BOOT = 2000


def git_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], cwd=config.ROOT, capture_output=True, text=True
        )
        return out.stdout.strip() or "unknown"
    except OSError:
        return "unknown"


def region_metrics(labels: np.ndarray, delta: np.ndarray) -> dict:
    path_score = -delta
    n_pos, n_neg = int(labels.sum()), int((1 - labels).sum())
    out = {"n_pathogenic": n_pos, "n_benign": n_neg}
    if n_pos < 2 or n_neg < 2:
        out["note"] = "too few variants of one class for metrics"
        return out
    out["auroc"] = metrics.bootstrap(metrics.auroc, labels, path_score, N_BOOT).as_dict()
    out["auprc"] = metrics.bootstrap(metrics.auprc, labels, path_score, N_BOOT).as_dict()
    out["auprc_no_skill"] = n_pos / (n_pos + n_neg)
    t = metrics.youden_threshold(labels, path_score)
    sens, spec = metrics.sensitivity_specificity(labels, path_score, t)
    out["youden_in_sample"] = {"threshold_dL": -t, "sensitivity": sens, "specificity": spec}
    sens, spec = metrics.sensitivity_specificity(labels, path_score, -PAPER_THRESHOLD_DL)
    out["at_paper_threshold"] = {
        "threshold_dL": PAPER_THRESHOLD_DL,
        "sensitivity": sens,
        "specificity": spec,
    }
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--precision", choices=PRECISIONS, default="fp8-delayed")
    args = parser.parse_args()

    data = mtdna.load()
    windows = [variant_windows(data.sequence, v, WINDOW, circular=True) for v in data.variants]
    unique = sorted({w for pair in windows for w in pair})
    print(f"{len(data.variants)} variants, {len(unique)} distinct {WINDOW} bp windows to score")

    model = Evo2Adapter("evo2_1b_base", precision=args.precision)
    model.load()
    t0 = time.perf_counter()
    ll = dict(zip(unique, model.score_sequences(unique), strict=True))
    elapsed = time.perf_counter() - t0
    print(f"scored in {elapsed:.0f} s")

    rows = []
    for v, (ref_w, alt_w) in zip(data.variants, windows, strict=True):
        gene = data.gene_of(v)
        rows.append(
            {
                "chrom": v.chrom,
                "pos": v.pos,
                "ref": v.ref,
                "alt": v.alt,
                "label": v.label,
                "gene": gene.name if gene else "",
                "region": gene.biotype if gene else "non-coding",
                "ll_ref": ll[ref_w],
                "ll_alt": ll[alt_w],
                "delta": ll[alt_w] - ll[ref_w],
            }
        )

    labels = np.array([r["label"] for r in rows])
    delta = np.array([r["delta"] for r in rows])
    region = np.array([r["region"] for r in rows])
    results = {
        "model": model.label,
        "window_bp": WINDOW,
        "score": "pathogenicity = -(mean LL(alt window) - mean LL(ref window))",
        "bootstrap": f"{N_BOOT} resamples, stratified by label, 95% percentile interval",
        "note": "youden_in_sample picks its threshold on the same variants it is scored on, "
        "so it is optimistic; at_paper_threshold uses the paper's fixed cut-off",
        "regions": {
            name: region_metrics(
                labels[region == name if name != "all" else slice(None)],
                delta[region == name if name != "all" else slice(None)],
            )
            for name in REGIONS
        },
        "data_sha256": {
            p.name: sha256(p)
            for p in sorted(config.DATA_RAW.iterdir())
            if p.is_file() and not p.name.endswith(".part")
        },
        "git_commit": git_commit(),
        "run_date": dt.date.today().isoformat(),
        "scoring_seconds": round(elapsed, 1),
    }

    out = config.ROOT / "results"
    out.mkdir(exist_ok=True)
    stem = f"baseline_{args.precision}"
    with open(out / f"{stem}.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (out / f"{stem}.json").write_text(json.dumps(results, indent=2) + "\n")

    print(f"\n{model.label}, {WINDOW} bp windows")
    for name in REGIONS:
        r = results["regions"][name]
        head = f"{name:15s} {r['n_pathogenic']:3d} P / {r['n_benign']:3d} B"
        if "auroc" not in r:
            print(f"{head}  ({r['note']})")
            continue
        au = metrics.Estimate(r["auroc"]["value"], *r["auroc"]["ci95"], N_BOOT)
        ap = metrics.Estimate(r["auprc"]["value"], *r["auprc"]["ci95"], N_BOOT)
        p = r["at_paper_threshold"]
        print(
            f"{head}  AUROC {au}  AUPRC {ap} (no-skill {r['auprc_no_skill']:.3f})  "
            f"at paper cut-off: sens {p['sensitivity']:.3f}, spec {p['specificity']:.3f}"
        )
    print(f"\nwrote results/{stem}.csv and results/{stem}.json")


if __name__ == "__main__":
    main()
