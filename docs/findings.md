# Findings log

The running results log. Each entry records the date, the model checkpoint, the
command that produced the number, and the number with its confidence interval.

No variant-effect results yet.

## 2026-10-04: batch size changes Evo 2 1B scores more than a variant does

Model: `evo2_1b_base` (1B, FP8, RTX 4060 Laptop GPU). On 8 single-base mtDNA
variants (501 bp windows), changing batch size from 1 to 8 or 16 moved
per-sequence scores by up to 0.0033, against a median variant effect of
0.00087, and flipped 3 of 8 variant signs. All scoring in this repo therefore
uses batch size 1. Reproduce: `python scripts/check_batch_size.py`. Details:
`docs/model-choice.md`. Single run on 16 sequences.
