# Model choice

**Model used: Evo 2, `evo2_1b_base` (1 billion parameters, 8 kb context).**
Every result in this repo comes from this checkpoint unless a result says
otherwise. It is the smallest Evo 2 model, and much smaller than the 7B and 40B
models used in the Evo 2 paper. Results here should not be read as results for
those larger models.

Decided 2026-10-03.

## Hardware

| | |
|---|---|
| GPU | NVIDIA GeForce RTX 4060 Laptop GPU, 8 GB, compute capability 8.9 (Ada) |
| Driver | 595.97 |
| Host | Windows 11, 24 GB RAM |
| Runtime | WSL2, Ubuntu 24.04, 11 GB RAM visible to Linux |

## Why this checkpoint

| Option | Verdict |
|---|---|
| Evo 2 natively on Windows | Not supported: Evo 2 needs Linux, flash-attn and (for 1B) Transformer Engine. Not attempted. |
| `evo2_7b` (bf16, no FP8 needed) | Weights alone are ~14 GB in bf16; does not fit in 8 GB VRAM. Not attempted. |
| **`evo2_1b_base` under WSL2** | **Works.** See measurements below. |
| `evo2_20b`, `evo2_40b` | Far beyond 8 GB. Not attempted. |
| Hosted API (NVIDIA NIM) or HyenaDNA | Held in reserve; not needed for v0.1. |

## Measurements

From `scripts/00_smoke_test.py --model evo2_1b_base`, run twice with
identical scores both times:

| | |
|---|---|
| Weights on disk | 2.6 GB |
| Load time (weights cached) | 3.4 s |
| Load time (first run, includes download) | 784 s |
| Scoring one 1 kb sequence | 0.18 s (after a warm-up call) |
| Peak GPU memory | 2.16 GiB after load, 2.35 GiB after scoring |
| Peak host RAM | 3.1 GB |

Sanity check on the first 1 kb of human mtDNA (`tests/fixtures/chrM_1_1000.fa`),
mean log-likelihood per base:

| Sequence | Score |
|---|---|
| Real | −0.9649 |
| 3 mononucleotide shuffles | −1.3769, −1.3732, −1.3732 |
| Uniform guessing (ln 0.25) | −1.3863 |

The model finds real mtDNA clearly more probable than shuffles of it, which sit
near chance. Batched and one-at-a-time scoring agree to 6 decimal places.

## Caveat: FP8 on a non-Hopper GPU

The Evo 2 README says the 1B model needs FP8 via Transformer Engine "and a
Nvidia Hopper GPU". This GPU is Ada (sm_89), which supports FP8 in hardware,
and the model runs here with FP8 input projections enabled
(`use_fp8_input_projections: True`) and no errors. The sanity check above shows
the output is meaningful, but it is a weak check: it does **not** show the
scores match what the same checkpoint produces on Hopper. If a reference score
from Hopper hardware becomes available, compare against it and record the
difference here.

Other notes from the load:
- Transformer Engine warns that it supports flash-attn ≤ 2.7.4.post1; we use
  2.8.0.post2 because the Evo 2 README pins it.
- The checkpoint loader reports unused `_extra_state` keys (FP8 metadata) and an
  unused `unembed.weight` key. Scores look sane, so these appear harmless, but
  they are noted in case a later discrepancy needs explaining.

## Software versions

Python 3.12.14, PyTorch 2.6.0 (CUDA 12.6, conda-forge),
transformer-engine-torch 2.3.0, flash-attn 2.8.0.post2, evo2 0.6.0.

Reproduce the environment from inside WSL with `bash scripts/setup_evo2_wsl.sh`.
