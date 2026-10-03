"""Smoke test: load one genomic LM, score 1 kb of real mtDNA, report time and peak memory.

Usage:
    python scripts/00_smoke_test.py --model evo2_1b_base

Scores the first 1 kb of human mtDNA (tests/fixtures/chrM_1_1000.fa) plus a few
shuffled copies of it, and prints wall-clock load time, scoring time, peak GPU
memory and each sequence's mean log-likelihood per base.

Sanity check: a working model should find the real sequence more probable than
its shuffles. If it does not, the model loaded but its numerics are suspect
(e.g. FP8 on an unsupported GPU), and the script exits non-zero. Nothing here is
a result; it only confirms the model runs correctly on this hardware.
See docs/model-choice.md.
"""

from __future__ import annotations

import argparse
import platform
import time
from pathlib import Path

import numpy as np

FIXTURE = Path(__file__).resolve().parent.parent / "tests" / "fixtures" / "chrM_1_1000.fa"


def read_fasta(path: Path) -> str:
    lines = path.read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">")).upper()


def shuffled(seq: str, seed: int) -> str:
    # Mononucleotide shuffle: fine for a sanity check, not for the controls (see CLAUDE.md 5.2).
    rng = np.random.default_rng(seed)
    return "".join(rng.permutation(list(seq)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--model", default="evo2_1b_base")
    parser.add_argument("--n-shuffles", type=int, default=3)
    args = parser.parse_args()

    import torch

    from seqcontrol.models.evo2 import Evo2Adapter

    print(
        f"python {platform.python_version()}  torch {torch.__version__}  cuda {torch.version.cuda}"
    )
    print(f"gpu {torch.cuda.get_device_name(0)}")
    torch.cuda.reset_peak_memory_stats()

    t0 = time.perf_counter()
    model = Evo2Adapter(args.model)
    model.load()
    torch.cuda.synchronize()
    load_s = time.perf_counter() - t0
    load_mem = torch.cuda.max_memory_allocated() / 2**30

    real = read_fasta(FIXTURE)
    shuffles = [shuffled(real, seed) for seed in range(args.n_shuffles)]

    model.score_sequences([real])  # warm-up, so the timing below excludes one-off kernel setup
    torch.cuda.synchronize()
    t0 = time.perf_counter()
    real_score = float(model.score_sequences([real])[0])
    torch.cuda.synchronize()
    score_s = time.perf_counter() - t0
    shuffle_scores = [float(x) for x in model.score_sequences(shuffles)]
    peak_mem = torch.cuda.max_memory_allocated() / 2**30

    print(f"model {args.model}")
    print(f"load_time_s {load_s:.1f}")
    print(f"score_time_s {score_s:.3f}  (one {len(real)} bp sequence)")
    print(f"peak_gpu_mem_gib after_load {load_mem:.2f}  after_score {peak_mem:.2f}")
    print(f"mean_loglik real      {real_score:.4f}")
    print(f"mean_loglik shuffled  {' '.join(f'{x:.4f}' for x in shuffle_scores)}")
    print(f"uniform baseline      {np.log(0.25):.4f}")

    if not all(np.isfinite([real_score, *shuffle_scores])):
        raise SystemExit("FAIL: non-finite score; the model ran but its output is invalid")
    if real_score <= max(shuffle_scores):
        raise SystemExit("FAIL: real sequence not more probable than its shuffles")
    print("PASS: real sequence scores above all shuffles")


if __name__ == "__main__":
    main()
