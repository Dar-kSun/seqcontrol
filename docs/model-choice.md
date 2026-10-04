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
near chance. Repeat runs at the same batch size give bit-identical scores.

(Correction: an earlier version of this file, in commit 96ccdc1, said batched
and one-at-a-time scoring agree to 6 decimal places. Both runs in that check
actually used batch size 1, so it tested nothing about batching. The real
comparison is below.)

## Batch size changes scores: always score one sequence at a time

From `python scripts/check_batch_size.py`, run 2026-10-04: the ref and alt
501 bp windows for 8 single-base variants in `tests/fixtures/chrM_1_1000.fa`.

| Comparison | Max absolute score difference |
|---|---|
| Batch size 1, run twice | 0.0 |
| Batch size 1 vs 8 | 0.0033 |
| Batch size 1 vs 16 | 0.0030 |
| *For scale: median absolute variant effect (alt − ref), batch size 1* | *0.00087* |

Batching moves scores by more than a typical variant effect. Between batch
size 1 and 16, 3 of the 8 variant scores changed sign. A likely cause, not
verified, is that FP8 scaling factors are computed over the whole batch, so a
sequence's score depends on its batch-mates.

**Decision:** `Evo2Adapter` always uses batch size 1 and does not expose a
batch-size option. This is slower but exactly reproducible.

## Precision: FP8 recipe, history, and bf16

`Evo2Adapter(precision=...)` offers three modes: `fp8-delayed` (Evo 2's shipped
recipe, Transformer Engine `DelayedScaling` with a 16-step amax history; the
default and the one all headline results use), `fp8-current` (TE current
scaling) and `bf16` (no FP8).

From `python scripts/check_precision.py`, run 2026-10-04:

- **No dependence on scoring history.** A fixed target window scored after 16
  random, poly-A or GC-rich sequences gets an identical score in every mode,
  despite the delayed recipe's amax history. So scoring order does not matter
  at batch size 1.
- **bf16 is not usable.** With FP8 off, real mtDNA scores −1.356 per base,
  close to uniform guessing (−1.386), against −1.095 with FP8. This matches
  the Evo 2 README's statement that the 1B model needs FP8 for accuracy. It is
  also possible that disabling FP8 after loading is not equivalent to a
  native bf16 configuration; this was not investigated further.
- **The two FP8 recipes give AUROC within 0.007 of each other on the full
  variant set, but per-variant scores differ noticeably** (Spearman 0.944, 10%
  of variants change sign). Details in `docs/findings.md`.

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
