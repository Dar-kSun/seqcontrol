# Findings log

The running results log. Each entry records the date, the model checkpoint, the
command that produced the number, and the number with its confidence interval.

## 2026-10-04: native baseline (M3)

**Model: Evo 2 `evo2_1b_base` (1B parameters), Evo 2's shipped FP8 recipe.**
1,025 bp windows centred on each variant (wrapping around the circular
chromosome). Pathogenicity score = −ΔL, where ΔL = mean log-likelihood of the
alt window minus that of the ref window. 95% intervals are from 2,000
label-stratified bootstrap resamples over variants.

Reproduce: `python scripts/01_baseline.py` (per-variant scores in
`results/baseline_fp8-delayed.csv`, metrics in `.json`). Run twice; the
per-variant scores were bit-identical.

| Region | Pathogenic / benign | AUROC | AUPRC (no-skill) |
|---|---|---|---|
| All | 94 / 228 | **0.856** [0.805, 0.903] | 0.760 [0.683, 0.836] (0.292) |
| tRNA | 44 / 23 | **0.824** [0.713, 0.915] | 0.910 [0.847, 0.960] (0.657) |
| Protein-coding | 48 / 193 | 0.914 [0.864, 0.955] | 0.776 [0.671, 0.877] (0.199) |
| rRNA | 2 / 12 | not interpretable (2 pathogenic) | |

The tRNA AUROC of 0.824 is the number the permutation control (M4) will be
measured against.

How to read this:

- **It is one model, one window size and one label set.** The paper reports an
  overall AUROC of 0.896 with an unstated checkpoint and a larger benign set
  (623, including ClinGen frequency data). The two numbers are not directly
  comparable.
- **The paper's threshold does not transfer.** At its cut-off (ΔL ≤ −0.0081),
  tRNA sensitivity here is 22.7% (specificity 100%), against the paper's 65.8%
  native. A fixed threshold depends on the model and the dataset, which is why
  this repo leads with AUROC.
- **tRNA AUPRC looks high because most tRNA variants in the set are pathogenic**
  (no-skill 0.657). Compare it with its no-skill value, not with 0.5.
- **Possible confound, not tested:** most ClinVar-benign mtDNA variants are
  common population polymorphisms. The model may partly score how *familiar*
  an allele is from training data rather than its effect on function.

### FP8 recipe sensitivity

Re-running with Transformer Engine's current-scaling FP8 recipe
(`--precision fp8-current`) gives overall AUROC 0.849 and tRNA 0.813.
`python scripts/compare_runs.py results/baseline_fp8-delayed.csv results/baseline_fp8-current.csv`:

| Region | Spearman of ΔL | Same sign | AUROC difference (paired bootstrap) |
|---|---|---|---|
| All | 0.944 | 90.1% | 0.007 [−0.013, 0.027] |
| tRNA | 0.953 | 92.5% | 0.011 [−0.035, 0.057] |
| Protein-coding | 0.942 | 88.8% | 0.009 [−0.012, 0.029] |

**The rounding recipe matters for single variants but not for AUROC.**
Per-variant scores move by a median of 31% of a typical |ΔL|, and 1 in 10
variants changes sign, but the AUROC difference is small and not significant.

**Consequence for the controls:** per-variant score stability (Spearman) under
a control has to be read against this ~0.94 floor. A control that leaves
Spearman near 0.94 has changed scores no more than switching FP8 recipes does.

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
