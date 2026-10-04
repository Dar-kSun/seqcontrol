# Findings log

The running results log. Each entry records the date, the model checkpoint, the
command that produced the number, and the number with its confidence interval.

## 2026-10-04: flank-shuffle sweep (M5)

**Model: Evo 2 `evo2_1b_base` (1B), Evo 2's shipped FP8 recipe. 67 tRNA
variants (44 pathogenic, 23 benign), 1,025 bp windows, 10 shuffle seeds per
radius.**

![Flank-shuffle sweep](../results/flank_sweep.png)

Each tRNA is held fixed. Everything more than *r* bp outside it, on either
side, is replaced by a dinucleotide-preserving (Altschul–Erikson) shuffle of
itself. At r = 0 all context is scrambled; at r = 400 only the outer
~40–80 bp at each end of the window are. Base composition and dinucleotide
counts never change, so only the *order* of the context is tested.

Reproduce: `python scripts/03_flank_sweep.py` (GPU, ~40 min) or `--from-csv`
(seconds), then `python scripts/plot_flank_sweep.py`.

| Untouched flank *r* | Spearman of ΔL vs native | AUROC (native 0.824) | CDI | Seed SD (ρ / AUROC) |
|---|---|---|---|---|
| 0 bp | 0.47 [0.30, 0.62] | 0.710 [0.615, 0.803] | **0.35 [0.01, 0.63]** | 0.054 / 0.045 |
| 25 bp | 0.56 [0.39, 0.71] | 0.694 [0.581, 0.802] | 0.40 [0.10, 0.72] | 0.063 / 0.027 |
| 50 bp | 0.67 [0.52, 0.79] | 0.739 [0.627, 0.837] | 0.26 [0.01, 0.53] | 0.043 / 0.046 |
| 100 bp | 0.77 [0.66, 0.86] | 0.765 [0.652, 0.861] | 0.18 [−0.12, 0.44] | 0.028 / 0.029 |
| 200 bp | 0.86 [0.77, 0.92] | 0.799 [0.698, 0.886] | 0.08 [−0.20, 0.27] | 0.026 / 0.022 |
| 300 bp | 0.91 [0.85, 0.94] | 0.800 [0.700, 0.887] | 0.07 [−0.12, 0.23] | 0.014 / 0.020 |
| 400 bp | 0.93 [0.89, 0.95] | 0.797 [0.688, 0.888] | 0.09 [−0.06, 0.22] | 0.013 / 0.019 |

Each statistic is the mean over the 10 seeds; intervals are 95% from 2,000
label-stratified bootstrap resamples of variants (paired with native); seed
SD is the spread across seeds on all variants.

What this shows:

1. **Scrambling all context costs about a third of the above-chance
   signal.** At r = 0, AUROC falls to 0.710 and CDI is 0.35, the only
   setting in this repo so far whose CDI interval excludes zero (and only
   just: lower bound 0.006). Scrambling the context removes more signal
   than moving the tRNA to another tRNA's real neighbourhood (M4: AUROC
   0.750, CDI 0.23).
2. **The model's dependence on context is local.** Score stability rises
   steadily with the untouched radius and reaches 0.93 at 400 bp, close to
   the 0.95 floor set by FP8 rounding noise alone. Discrimination recovers
   faster: from r = 200 bp on, AUROC is ~0.80 and CDI intervals include
   zero.
3. **Seeds agree.** Across 10 shuffles, AUROC varies by SD ≤ 0.046 at any
   radius, so the variant sample (n = 67), not the shuffle, dominates the
   uncertainty.
4. **Even the gene alone keeps most of the signal.** At r = 0, AUROC 0.710 is
   still well above 0.5: most of the discrimination survives with nothing
   but the tRNA's own sequence in its true order.

### Robustness: FP8 recipe

All three controls were re-run with Transformer Engine's current-scaling FP8
recipe (`--precision fp8-current`; results files with that suffix). The
conclusions do not change:

| | Shipped recipe (headline) | Current scaling |
|---|---|---|
| Native tRNA AUROC | 0.824 | 0.813 |
| tRNA swap: AUROC under control / CDI | 0.750 / 0.23 [−0.12, 0.45] | 0.747 / 0.21 [−0.12, 0.44] |
| Window rotation: CDI | 0.09 [−0.05, 0.23] | 0.08 [−0.06, 0.20] |
| Flank shuffle r = 0: AUROC / CDI | 0.710 / 0.35 [0.01, 0.63] | 0.701 / 0.36 [0.03, 0.63] |
| Flank shuffle r = 400: Spearman | 0.93 | 0.94 |

Provenance note: `results/flank_sweep.json` was scored by code at commit
05b672a; the run read HEAD only when it finished and first stamped 91eb3c5.
The stamp was corrected by hand. The scoring code for default arguments is
the same at both commits.

## 2026-10-04: context-swap controls on tRNA variants (M4)

**Model: Evo 2 `evo2_1b_base` (1B), Evo 2's shipped FP8 recipe. 67 tRNA
variants (44 pathogenic, 23 benign). Single scoring run per control.**

Two versions of the control, both keeping every tRNA's own bases
byte-identical (asserted on every perturbed sequence):

- **tRNA swap** (the source paper's design): every tRNA moves into the slot
  of the tRNA *k* places along the chromosome, carrying its own sequence;
  all 19 shifts. Overlapping tRNAs (MT-TI+MT-TQ, MT-TC+MT-TY) move as one unit.
- **Window rotation:** each 1,025 bp scoring window is rotated by 64, 128, …
  bp, skipping offsets that would cut the tRNA; 13 offsets.

The unperturbed setting reproduces the M3 baseline scores exactly.

Reproduce: `python scripts/02_permutation.py --control trna-swap` (GPU,
~10 min), or `--from-csv` to recompute every number below from
`results/permutation_*.csv` in seconds without a GPU.

| | tRNA swap | Window rotation |
|---|---|---|
| AUROC native | 0.824 [0.713, 0.915] | 0.824 [0.713, 0.915] |
| AUROC under control (mean over settings) | 0.750 [0.661, 0.826] | 0.793 [0.692, 0.877] |
| Range across settings | 0.657 – 0.835 | 0.745 – 0.843 |
| **AUROC drop** | **0.074 [−0.030, 0.162]** | 0.031 [−0.015, 0.072] |
| **CDI** | **0.23 [−0.12, 0.45]** | 0.09 [−0.05, 0.23] |
| Spearman of ΔL, native vs control (mean) | **0.56** (0.40 – 0.71) | 0.90 (0.77 – 0.96) |
| Sensitivity / specificity at native-Youden cut-off (ΔL ≤ −0.0030) | 0.82 / 0.83 → **0.56 / 0.83** | 0.82 / 0.83 → 0.69 / 0.80 |
| Sensitivity / specificity at the paper's cut-off (ΔL ≤ −0.0081) | 0.23 / 1.00 → 0.16 / 1.00 | 0.23 / 1.00 → 0.22 / 1.00 |
| Median ΔL, pathogenic | −0.0056 → −0.0036 | −0.0056 → −0.0048 |
| Median ΔL, benign | −0.0013 → −0.0007 | −0.0013 → −0.0013 |

95% intervals: 2,000 label-stratified bootstrap resamples of variants, paired
across native and control.

What this shows:

1. **Moving a tRNA to another tRNA's address changes its variant scores a
   lot.** Per-variant Spearman falls to 0.56, far below the ~0.95 that
   changing only the FP8 recipe produces. Context is a large part of each
   individual score.
2. **But ranking survives better than the scores do.** AUROC falls from
   0.824 to 0.750, still well above chance, and the drop (0.074) has an
   interval that includes zero at this sample size. CDI is 0.23: on the
   point estimate, about a quarter of the above-chance signal depends on
   the tRNA's genomic address, but the data are also consistent with none
   and with nearly half.
3. **Threshold metrics exaggerate the effect.** Under the swap, effect sizes
   shrink for both classes (median pathogenic ΔL −0.0056 → −0.0036), so a
   cut-off fixed on native scores loses 26 points of sensitivity while
   specificity is unchanged. This is the mechanism suggested in the M2 entry
   for the paper's 65.8% → 5.1% collapse: a compressed score distribution
   crossing a fixed threshold, rather than (only) lost discrimination.
4. **Window rotation is a milder control** (Spearman 0.90, CDI 0.09). It
   keeps the same bases in the window and only rearranges them, so it tests
   sensitivity to arrangement and to an artificial junction, not to a
   genuinely different neighbourhood.

**Relation to the paper.** The direction is reproduced: moving tRNAs lowers
sensitivity at a fixed threshold. The magnitude is not: 0.82 → 0.56 here
against 0.658 → 0.051 in the paper, and specificity does not rise as it did
there. The differences in setup are large (model size unstated in the paper,
a different benign set, a different threshold and possibly a different
permutation), so this is not evidence the paper is wrong. It does show that,
for `evo2_1b_base` on this variant set, the tRNA signal is far from purely
contextual.

Limitations specific to this result:
- n = 67, with 23 benign. Intervals are wide; the AUROC drop is not
  significant.
- tRNA labels cluster by gene (MT-TL1 has 13 pathogenic, 0 benign). Part of
  the native AUROC may be the model telling genes apart rather than variants
  within them; this was not separated out.
- The swap leaves genes in reference-strand orientation. A minus-strand tRNA
  placed in a plus-strand neighbourhood is part of the perturbation.
- Single run per control; FP8-recipe sensitivity was checked for the
  baseline only.

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
