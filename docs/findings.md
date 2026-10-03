# Findings log

The running results log. Each entry records the date, the model checkpoint, the
command that produced the number, and the number with its confidence interval.

No variant-effect results yet.

## 2026-10-04: the labelled mtDNA set, and how it differs from the paper

`python scripts/fetch_data.py` builds 94 pathogenic (MITOMAP `Cfrm-[P]`/`Cfrm-[LP]`)
and 228 benign (ClinVar, 2+ stars) single-base chrM variants. Counts per region
are in `data/MANIFEST.md`. The tRNA subset, which the permutation control (M4)
targets, is small: **44 pathogenic, 23 benign**. Confidence intervals there will
be wide.

Compared with the benchmark this repo builds on (Mathur & Sachidanandam,
bioRxiv 2026.03.10.710786):

| | Paper | This repo |
|---|---|---|
| Pathogenic | MITOMAP confirmed, n = 130 | MITOMAP `Cfrm-[P]`/`[LP]`, SNVs only, n = 94 |
| Benign | ClinVar 2+ stars **plus ClinGen frequency data**, n = 623 | ClinVar 2+ stars only, n = 228 |
| tRNA pathogenic | 79 | 44 |
| Model | not stated for the mtDNA experiment | `evo2_1b_base` |

Three points from the paper's methods that M4 must take into account:

1. **The 65.8% to 5.1% collapse is sensitivity at a fixed threshold**
   (Youden-optimal ΔL = −0.0081), not AUROC, and specificity rose to 93.8%
   at the same time. That pattern fits a shift of the whole score
   distribution, which a fixed threshold detects but which need not reduce
   ranking (AUROC). The paper reports no AUROC after permutation, so this
   repo's AUROC-based CDI is a new measurement, not a reproduction. M4 will
   report both sensitivity at a threshold and AUROC.
2. **What "cyclic permutation" means.** The paper "cyclically permuted the
   positions of all mitochondrial tRNAs while leaving their internal sequences
   intact", which reads as moving tRNA genes between one another's positions.
   CLAUDE.md §5.1 describes rotating the sequence window instead. These are
   different perturbations; M4 must pick one, or run both, and say which.
3. The paper's window is 512 bp either side of the variant (1,025 bp).

## 2026-10-04: batch size changes Evo 2 1B scores more than a variant does

Model: `evo2_1b_base` (1B, FP8, RTX 4060 Laptop GPU). On 8 single-base mtDNA
variants (501 bp windows), changing batch size from 1 to 8 or 16 moved
per-sequence scores by up to 0.0033, against a median variant effect of
0.00087, and flipped 3 of 8 variant signs. All scoring in this repo therefore
uses batch size 1. Reproduce: `python scripts/check_batch_size.py`. Details:
`docs/model-choice.md`. Single run on 16 sequences.
