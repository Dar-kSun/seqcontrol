# Findings log

A running log of results, newest first. Each entry gives the date, the model
checkpoint, the command that produced the numbers, and the numbers with their
confidence intervals.

## 2026-10-04: v0.3 Arm 1, within-gene discrimination on protein-coding variants

Pre-declared in `docs/plan-v0.3-within-gene.md`, which was committed (96cb968)
before anything was run. Model: Evo 2 `evo2_1b_base` (1B), shipped FP8 recipe.
Strict labels. This reuses the M3 baseline scores; nothing new was scored.
Reproduce with `python scripts/06_within_gene.py --arm 1`. Intervals come from
2,000 cluster-bootstrap resamples over genes.

There are 241 protein-coding variants (48 pathogenic) in 13 genes, and every
gene has both classes. That gives 893 within-gene pathogenic–benign pairs, most
of them in MT-ND1 (264), MT-ATP6 (225) and MT-ND5 (138).

| Score | AUROC |
|---|---|
| Model, whole set | 0.914 [0.853, 0.965] |
| Gene prior, leave-one-variant-out (ignores the variant) | 0.542 [0.279, 0.629] |
| Gene prior, in-sample (upper bound) | 0.697 [0.598, 0.759] |
| Model, within-gene pairs only | **0.934 [0.847, 0.978]** |

Dropping any one gene changes the within-gene AUROC by -0.012 to +0.018.

The pre-declared verdict is that the model separates pathogenic from benign
variants within the same protein-coding gene. Gene identity can't explain this
one. The gene prior is weak (0.697 at most, and that's when it is allowed to
see each variant's own label), and the model does as well within genes as
across them. Whatever is going on with the tRNAs, on protein-coding variants
`evo2_1b_base` is reading the variant.

No context control was run here, because these genes are longer than the
scoring window. Single scoring run.

## 2026-10-04: v0.2 analyses B and A (pre-declared)

Both analyses follow `docs/plan-v0.2-threshold-and-gene.md`, which was
committed before either was run (the draft is d685ed8 and the amendments are
d76533a; the amendments are also listed at the top of the plan). The decision
rules are applied in code. Neither needs a GPU, since both work from the
committed per-variant scores. Model: Evo 2 `evo2_1b_base` (1B), shipped FP8
recipe; 67 tRNA variants (44 pathogenic, 23 benign) in 20 genes; single
scoring run.

### B: does gene identity explain the native AUROC?

Reproduce with `python scripts/05_gene_confound.py`, which writes
`results/gene_confound.*`. The intervals here come from 2,000 cluster-bootstrap
resamples over genes, because 13 variants in MT-TL1 are not 13 independent
observations. That makes them wider than the variant-level intervals elsewhere
in this log; on this basis the native AUROC is 0.824 [0.724, 0.918].

| | AUROC |
|---|---|
| Model (−ΔL) | 0.824 [0.724, 0.918] |
| B1: gene prior, leave-one-variant-out (ignores the variant; biased low) | 0.627 [0.306, 0.802] |
| B1: gene prior, in-sample (uses the variant's own label; biased high) | 0.901 [0.770, 0.969] |
| B2: model, within-gene pairs only (45 pairs) | 0.867 [0.704, 0.938] |
| B3: model, MT-TL1 removed | 0.776 [0.704, 0.878] |

- B1. The leave-one-out gene prior lands in the plan's 0.60–0.75 band, so the
  verdict is that gene identity is a substantial part of the native AUROC. A
  score that only knows which tRNA a variant is in reaches somewhere between
  0.627 and 0.901, and the model's 0.824 sits inside that range.
- B2. Only 45 of the 1,012 pathogenic–benign pairs share a gene (25 of them in
  MT-TS1), short of the 100 the plan asked for. So the verdict is that this
  dataset can't separate gene identity from variant effect. On those 45 pairs
  the model ranks the pathogenic variant higher 87% of the time, which is
  encouraging but only descriptive.
- B3. Dropping any one gene changes AUROC by −0.049 to +0.034. MT-TL1 (13
  pathogenic, 0 benign) matters most, at −0.0485, just under the plan's 0.05
  threshold.
- B4. Scored on within-gene pairs only, every control's AUROC drop has an
  interval that includes zero (tRNA swap 0.044 [−0.071, 0.182]; flank shuffle
  at r = 0, 0.111 [−0.021, 0.417]). As the plan expected, this is
  underpowered, so no CDI conclusions are drawn from it.
- B5. Among all 228 benign variants, ΔL correlates with population allele
  frequency at Spearman 0.109 [−0.022, 0.240]; among the 23 benign tRNA
  variants, −0.245 [−0.618, 0.225]. There's no clear sign that the model is
  scoring familiarity, but this is a weak check (the median benign allele
  frequency is only 0.5%).

### A: is the published collapse a threshold artefact?

Reproduce with `python scripts/04_threshold_artefact.py`, which writes
`results/threshold_artefact.*`.

![Threshold artefact](../results/threshold_artefact.png)

A0. Multiplying every score by some α > 0 leaves AUROC exactly the same (this
is asserted for 106 values of α) but moves scores across any cut-off fixed in
absolute units. So shrinking all scores toward zero can wipe out sensitivity
without touching discrimination at all.

A2. How much did the tRNA swap shrink the scores? Either α̂ = 0.772
[0.552, 0.922] (ratio of median |ΔL|) or 0.744 [0.517, 0.892] (a fit through
the origin). A pure change of scale only explains part of what the swap did:
the fit has R² = 0.49 [0.15, 0.71], which fits with the per-variant Spearman
of 0.56. The two classes also shrink by different amounts. Pathogenic |ΔL|
drops to 0.65 of its original size, benign only to 0.89. (The plan's "~0.64
and ~0.54" came from medians of signed ΔL. Benign ΔL sits on both sides of
zero, so signed and absolute ratios come out differently.) When the classes
shrink differently, the ranking changes too, not just the scale.

A3. The decomposition, at the cut-off chosen on the original scores
(Youden), ΔL ≤ −0.0030:

| | Sensitivity | Specificity | AUROC |
|---|---|---|---|
| Native | 0.818 | 0.826 | 0.824 |
| Pure-compression null (native × α̂, either estimate) | 0.636 | 0.870 | 0.824 (identical by construction) |
| Observed tRNA swap | 0.557 | 0.828 | 0.750 |

- Pure compression reproduces 0.697 [0.23, 1.15] of the sensitivity drop,
  with either estimate of α̂.
- What's left over as lost ranking is an AUROC drop of 0.074
  [−0.030, 0.162], which isn't significant.
- The plan's rule was: 0.70 or more with a non-significant AUROC drop means a
  threshold artefact, under 0.40 means the idea is wrong, and anything in
  between is indeterminate. So the verdict is indeterminate at this sample
  size, because the share falls just short of 0.70. With 44 pathogenic
  variants, sensitivity moves in steps of 1/44 = 0.023, and one more variant
  crossing the cut-off would have given 0.78. The rule is applied as written
  anyway.
- One thing argues against pure compression. If all scores just shrank,
  specificity should rise (0.826 to 0.870) as benign scores also move away
  from the cut-off. Under the real swap it stayed at 0.828.

A4. A rank-preserving surrogate takes the swap's score values and hands them
out in the original order. At the paper's cut-off it matches the swap's
sensitivity almost exactly (0.160 vs 0.159). At the native cut-off it covers
most of the drop (0.592 vs 0.557) but not the specificity (0.895 vs 0.828). So
the change in the *distribution* of scores accounts for most of the lost
sensitivity, and reordering accounts for the rest, and for specificity
staying flat.

A5. A simulation of the paper's operating point. These are synthetic scores,
not the paper's data, and this is not a reproduction. Two Gaussians are
calibrated to the paper's native 65.8% sensitivity and 78.5% specificity at a
Youden-optimal cut-off of ΔL ≤ −0.0081, using this repo's pathogenic SD for
the one free scale. Shrinking every score by α = 0.42 (0.27–0.59 if the scale
is halved or doubled) brings sensitivity down to the paper's 5.1% with AUROC
unchanged at 0.791. Specificity rises, as the paper reports, but to ~100%
rather than 93.8%. So compression alone is enough to produce a collapse of
that size and shape, and it overshoots the specificity rise, which fits
compression plus some reordering.

### What A and B mean together

- The repo's headline AUROC is partly about which tRNA a variant is in. A
  gene prior that ignores the variant gets close to the model, and the data
  can't tell how much of the model's signal comes from within genes.
- Under the tRNA swap, about 70% of the lost sensitivity at a fixed threshold
  is what pure compression would produce. But that falls just short of the
  bar set in advance, compression alone doesn't explain why specificity
  stayed flat, and pathogenic scores shrink more than benign ones. Compression
  is a large part of the story, not all of it.
- None of this means the source paper is wrong. It measured a real effect at
  a real operating point; this shows one mechanism that produces the same
  pattern.

## 2026-10-04: wider windows (4,097 bp)

Model: Evo 2 `evo2_1b_base` (1B); one run per configuration. With 1,025 bp
windows the flank sweep can't go past r = 400 bp, so the baseline and sweep
were repeated with 4,097 bp windows and radii up to 1,900 bp. Reproduce with
`python scripts/01_baseline.py --window 4097` (adding `--precision
fp8-current` for the second recipe) and
`python scripts/03_flank_sweep.py --window 4097 --radii 0,50,100,250,500,1000,1500,1900`
(~3.5 h on the GPU). The comparison below comes from
`python scripts/compare_runs.py results/baseline_fp8-delayed_w4097.csv results/baseline_fp8-current_w4097.csv`.

At 4,097 bp the rounding noise is about four times worse. A variant's effect
is averaged over four times as many positions, so the median |ΔL| falls from
0.0028 to 0.00075, but the FP8 rounding noise doesn't shrink with it:

| | 1,025 bp | 4,097 bp |
|---|---|---|
| Spearman of ΔL between the two FP8 recipes, tRNA | 0.95 | **0.84** |
| Same sign between recipes, all variants | 90.1% | 86.0% |
| Native AUROC, all variants (shipped / current-scaling recipe) | 0.856 / 0.849 | 0.815 / 0.834 |
| Native AUROC, tRNA (shipped / current-scaling) | 0.824 / 0.813 | 0.779 / 0.760 |
| AUROC difference between recipes, protein-coding | 0.009 [−0.012, 0.029] | −0.031 [−0.063, −0.004] |

At 4,097 bp the choice of rounding recipe alone shifts the protein-coding
AUROC by a detectable amount. Longer windows don't help this model
discriminate (all-variant AUROC 0.815–0.834, against 0.849–0.856 at 1,025 bp)
and make its scores noisier; per-variant ΔL at the two window sizes
correlates at Spearman 0.86. The 1,025 bp results stay as the headline.

The flank sweep at 4,097 bp (shipped recipe, native tRNA AUROC 0.779):

| Untouched flank *r* | Spearman of ΔL vs native | AUROC | CDI |
|---|---|---|---|
| 0 bp | 0.49 [0.30, 0.64] | 0.733 [0.637, 0.830] | 0.16 [−0.36, 0.48] |
| 50 bp | 0.63 [0.45, 0.76] | 0.742 [0.635, 0.845] | 0.13 [−0.36, 0.45] |
| 100 bp | 0.68 [0.50, 0.81] | 0.752 [0.642, 0.851] | 0.09 [−0.38, 0.41] |
| 250 bp | 0.70 [0.54, 0.82] | 0.782 [0.679, 0.868] | −0.01 [−0.52, 0.28] |
| 500 bp | 0.80 [0.65, 0.90] | 0.792 [0.689, 0.878] | −0.05 [−0.46, 0.20] |
| 1,000 bp | 0.79 [0.64, 0.90] | 0.800 [0.699, 0.887] | −0.08 [−0.49, 0.15] |
| 1,500 bp | 0.86 [0.78, 0.91] | 0.775 [0.669, 0.863] | 0.01 [−0.29, 0.21] |
| 1,900 bp | 0.86 [0.77, 0.91] | 0.773 [0.670, 0.863] | 0.02 [−0.28, 0.22] |

Measured against the 0.84 noise floor, shuffling context more than ~500 bp
from the tRNA has no detectable effect on per-variant scores, and none on
AUROC beyond ~250 bp. That matches the 1,025 bp sweep: the model uses nearby
context. The CDI intervals here are much wider than at 1,025 bp, because the
native AUROC is closer to chance and the scores are noisier, so this run can't
say how much signal the full scramble takes away.

## 2026-10-04: flank-shuffle sweep (M5)

Model: Evo 2 `evo2_1b_base` (1B) with Evo 2's shipped FP8 recipe. 67 tRNA
variants (44 pathogenic, 23 benign), 1,025 bp windows, 10 shuffle seeds per
radius.

![Flank-shuffle sweep](../results/flank_sweep.png)

Each tRNA is held fixed, and everything more than *r* bp outside it, on either
side, is replaced by a dinucleotide-preserving (Altschul–Erikson) shuffle of
itself. At r = 0 all the context is scrambled; at r = 400 only the outer
~40–80 bp at each end of the window are. Base composition and dinucleotide
counts never change, so what's being tested is the *order* of the context.

Reproduce with `python scripts/03_flank_sweep.py` (GPU, ~40 min) or with
`--from-csv` (seconds), then `python scripts/plot_flank_sweep.py`.

| Untouched flank *r* | Spearman of ΔL vs native | AUROC (native 0.824) | CDI | Seed SD (ρ / AUROC) |
|---|---|---|---|---|
| 0 bp | 0.47 [0.30, 0.62] | 0.710 [0.615, 0.803] | **0.35 [0.01, 0.63]** | 0.054 / 0.045 |
| 25 bp | 0.56 [0.39, 0.71] | 0.694 [0.581, 0.802] | 0.40 [0.10, 0.72] | 0.063 / 0.027 |
| 50 bp | 0.67 [0.52, 0.79] | 0.739 [0.627, 0.837] | 0.26 [0.01, 0.53] | 0.043 / 0.046 |
| 100 bp | 0.77 [0.66, 0.86] | 0.765 [0.652, 0.861] | 0.18 [−0.12, 0.44] | 0.028 / 0.029 |
| 200 bp | 0.86 [0.77, 0.92] | 0.799 [0.698, 0.886] | 0.08 [−0.20, 0.27] | 0.026 / 0.022 |
| 300 bp | 0.91 [0.85, 0.94] | 0.800 [0.700, 0.887] | 0.07 [−0.12, 0.23] | 0.014 / 0.020 |
| 400 bp | 0.93 [0.89, 0.95] | 0.797 [0.688, 0.888] | 0.09 [−0.06, 0.22] | 0.013 / 0.019 |

Each statistic is the mean over the 10 seeds. The intervals are 95%, from
2,000 label-stratified bootstrap resamples of variants (paired with native).
Seed SD is the spread across seeds, computed on all variants.

What this shows:

1. Scrambling all the context costs about a third of the above-chance signal.
   At r = 0, AUROC drops to 0.710 and CDI is 0.35, the only setting in the
   repo so far whose CDI interval excludes zero (and only just, with a lower
   bound of 0.006). Scrambling the context takes away more signal than moving
   the tRNA into another tRNA's real neighbourhood (M4: AUROC 0.750, CDI 0.23).
2. The context the model depends on is local. Score stability rises steadily
   with the untouched radius and reaches 0.93 at 400 bp, close to the 0.95
   floor from FP8 rounding noise alone. Discrimination comes back faster: from
   r = 200 bp on, AUROC is ~0.80 and the CDI intervals include zero.
3. The seeds agree with each other. Across 10 shuffles, AUROC varies by an SD
   of at most 0.046 at any radius, so the uncertainty comes from the small
   variant sample (n = 67), not from the shuffling.
4. The gene alone keeps most of the signal. At r = 0, an AUROC of 0.710 is
   still well above 0.5, so most of the discrimination survives with nothing
   but the tRNA's own sequence in its real order.

### Robustness: FP8 recipe

All three controls were re-run with Transformer Engine's current-scaling FP8
recipe (`--precision fp8-current`; the results files carry that suffix). None
of the conclusions change:

| | Shipped recipe (headline) | Current scaling |
|---|---|---|
| Native tRNA AUROC | 0.824 | 0.813 |
| tRNA swap: AUROC under control / CDI | 0.750 / 0.23 [−0.12, 0.45] | 0.747 / 0.21 [−0.12, 0.44] |
| Window rotation: CDI | 0.09 [−0.05, 0.23] | 0.08 [−0.06, 0.20] |
| Flank shuffle r = 0: AUROC / CDI | 0.710 / 0.35 [0.01, 0.63] | 0.701 / 0.36 [0.03, 0.63] |
| Flank shuffle r = 400: Spearman | 0.93 | 0.94 |

A note on provenance: `results/flank_sweep.json` was scored by the code at
commit e25d6d0, but the run only read HEAD when it finished, so it first
stamped 774e095. The stamp was corrected by hand. The scoring code for default
arguments is the same at both commits.

## 2026-10-04: context-swap controls on tRNA variants (M4)

Model: Evo 2 `evo2_1b_base` (1B) with Evo 2's shipped FP8 recipe. 67 tRNA
variants (44 pathogenic, 23 benign). One scoring run per control.

There are two versions of the control. Both keep every tRNA's own bases
exactly as they are, and this is checked on every perturbed sequence.

- tRNA swap (the source paper's design): every tRNA moves into the slot of the
  tRNA *k* places along the chromosome, carrying its own sequence, for all 19
  shifts. Overlapping tRNAs (MT-TI+MT-TQ, MT-TC+MT-TY) move together as one
  unit.
- Window rotation: each 1,025 bp scoring window is rotated by 64, 128, … bp,
  skipping any offset that would cut the tRNA, which leaves 13 offsets.

The unperturbed setting reproduces the M3 baseline scores exactly.

Reproduce with `python scripts/02_permutation.py --control trna-swap` (GPU,
~10 min), or add `--from-csv` to recompute every number below from
`results/permutation_*.csv` in seconds without a GPU.

| | tRNA swap | Window rotation |
|---|---|---|
| AUROC native | 0.824 [0.713, 0.915] | 0.824 [0.713, 0.915] |
| AUROC under control (mean over settings) | 0.750 [0.661, 0.826] | 0.793 [0.692, 0.877] |
| Range across settings | 0.657 – 0.835 | 0.745 – 0.843 |
| AUROC drop | 0.074 [−0.030, 0.162] | 0.031 [−0.015, 0.072] |
| CDI | 0.23 [−0.12, 0.45] | 0.09 [−0.05, 0.23] |
| Spearman of ΔL, native vs control (mean) | 0.56 (0.40 – 0.71) | 0.90 (0.77 – 0.96) |
| Sensitivity / specificity at native-Youden cut-off (ΔL ≤ −0.0030) | 0.82 / 0.83 → 0.56 / 0.83 | 0.82 / 0.83 → 0.69 / 0.80 |
| Sensitivity / specificity at the paper's cut-off (ΔL ≤ −0.0081) | 0.23 / 1.00 → 0.16 / 1.00 | 0.23 / 1.00 → 0.22 / 1.00 |
| Median ΔL, pathogenic | −0.0056 → −0.0036 | −0.0056 → −0.0048 |
| Median ΔL, benign | −0.0013 → −0.0007 | −0.0013 → −0.0013 |

The 95% intervals come from 2,000 label-stratified bootstrap resamples of
variants, paired between native and control.

What this shows:

1. Moving a tRNA to another tRNA's address changes its variant scores a lot.
   Per-variant Spearman drops to 0.56, far below the ~0.95 you get from
   switching the FP8 recipe and nothing else. Context is a big part of each
   individual score.
2. The ranking holds up better than the scores do. AUROC falls from 0.824 to
   0.750, which is still well above chance, and the interval on the drop
   (0.074) includes zero at this sample size. CDI is 0.23: taking the point
   estimate, about a quarter of the above-chance signal depends on where the
   tRNA sits in the genome, but the data are also consistent with none and
   with nearly half.
3. Fixed thresholds make it look worse. Under the swap, scores shrink toward
   zero for both classes (median pathogenic ΔL −0.0056 to −0.0036), so a
   cut-off fixed on the original scores loses 26 points of sensitivity while
   specificity doesn't move. This is the mechanism suggested in the M2 entry
   for the paper's 65.8% to 5.1% collapse: a compressed score distribution
   crossing a fixed threshold, rather than (only) a loss of discrimination.
4. Window rotation is a milder control (Spearman 0.90, CDI 0.09). The window
   keeps the same bases and only rearranges them, so it tests sensitivity to
   arrangement and to an artificial junction, not to a genuinely different
   neighbourhood.

How this relates to the paper: the direction holds, since moving tRNAs lowers
sensitivity at a fixed threshold. The size doesn't: 0.82 to 0.56 here,
against 0.658 to 0.051 in the paper, and specificity doesn't rise the way it
did there. The setups differ a lot (the paper doesn't say which model size it
used, and the benign set, the threshold and possibly the permutation are
different), so this isn't evidence that the paper is wrong. What it does show
is that for `evo2_1b_base` on this variant set, the tRNA signal is far from
purely contextual.

Limitations specific to this result:
- n = 67, with 23 benign. The intervals are wide and the AUROC drop isn't
  significant.
- The tRNA labels cluster by gene (MT-TL1 has 13 pathogenic and 0 benign).
  Part of the native AUROC may be the model telling genes apart rather than
  variants within them. This wasn't separated out here (see the v0.2 entry
  above).
- The swap leaves genes in their reference-strand orientation, so a
  minus-strand tRNA landing in a plus-strand neighbourhood is part of the
  perturbation.
- One run per control. Sensitivity to the FP8 recipe was checked for the
  baseline only.

## 2026-10-04: native baseline (M3)

Model: Evo 2 `evo2_1b_base` (1B parameters) with Evo 2's shipped FP8 recipe.
Each variant gets a 1,025 bp window centred on it, wrapping around the
circular chromosome. The pathogenicity score is −ΔL, where ΔL is the mean
log-likelihood of the alt window minus that of the ref window. The 95%
intervals come from 2,000 label-stratified bootstrap resamples over variants.

Reproduce with `python scripts/01_baseline.py` (per-variant scores go to
`results/baseline_fp8-delayed.csv`, metrics to the `.json`). It was run twice
and the per-variant scores were bit-identical.

| Region | Pathogenic / benign | AUROC | AUPRC (no-skill) |
|---|---|---|---|
| All | 94 / 228 | **0.856** [0.805, 0.903] | 0.760 [0.683, 0.836] (0.292) |
| tRNA | 44 / 23 | **0.824** [0.713, 0.915] | 0.910 [0.847, 0.960] (0.657) |
| Protein-coding | 48 / 193 | 0.914 [0.864, 0.955] | 0.776 [0.671, 0.877] (0.199) |
| rRNA | 2 / 12 | not interpretable (2 pathogenic) | |

The tRNA AUROC of 0.824 is what the permutation control (M4) is measured
against.

Some things to keep in mind when reading this:

- It's one model, one window size and one label set. The paper reports an
  overall AUROC of 0.896 with an unstated checkpoint and a bigger benign set
  (623, including ClinGen frequency data), so the two numbers aren't directly
  comparable.
- The paper's threshold doesn't carry over. At its cut-off (ΔL ≤ −0.0081),
  tRNA sensitivity here is 22.7% (specificity 100%), against the paper's 65.8%
  native. A fixed threshold depends on both the model and the dataset, which
  is why this repo leads with AUROC.
- The tRNA AUPRC looks high only because most tRNA variants in the set are
  pathogenic (no-skill 0.657). Compare it with that, not with 0.5.
- A possible confound, not tested here: most ClinVar-benign mtDNA variants are
  common population polymorphisms, so the model may partly be scoring how
  *familiar* an allele looks from its training data rather than what it does.

### FP8 recipe sensitivity

Re-running with Transformer Engine's current-scaling FP8 recipe
(`--precision fp8-current`) gives an overall AUROC of 0.849 and 0.813 for tRNA.
From `python scripts/compare_runs.py results/baseline_fp8-delayed.csv results/baseline_fp8-current.csv`:

| Region | Spearman of ΔL | Same sign | AUROC difference (paired bootstrap) |
|---|---|---|---|
| All | 0.944 | 90.1% | 0.007 [−0.013, 0.027] |
| tRNA | 0.953 | 92.5% | 0.011 [−0.035, 0.057] |
| Protein-coding | 0.942 | 88.8% | 0.009 [−0.012, 0.029] |

The rounding recipe matters for single variants but not for AUROC. Per-variant
scores move by a median of 31% of a typical |ΔL|, and 1 in 10 variants changes
sign, but the AUROC difference is small and not significant.

This matters for the controls: any per-variant stability (Spearman) under a
control has to be read against this ~0.94 floor. A control that leaves
Spearman near 0.94 has changed the scores no more than switching FP8 recipes
does.

## 2026-10-04: the labelled mtDNA set, and how it differs from the paper

`python scripts/fetch_data.py` builds 94 pathogenic (MITOMAP `Cfrm-[P]`/`Cfrm-[LP]`)
and 228 benign (ClinVar, 2+ stars) single-base chrM variants; the counts per
region are in `data/MANIFEST.md`. The tRNA subset, which the permutation
control (M4) targets, is small: 44 pathogenic and 23 benign. Confidence
intervals there will be wide.

Compared with the benchmark this repo builds on (Mathur & Sachidanandam,
bioRxiv 2026.03.10.710786):

| | Paper | This repo |
|---|---|---|
| Pathogenic | MITOMAP confirmed, n = 130 | MITOMAP `Cfrm-[P]`/`[LP]`, SNVs only, n = 94 |
| Benign | ClinVar 2+ stars plus ClinGen frequency data, n = 623 | ClinVar 2+ stars only, n = 228 |
| tRNA pathogenic | 79 | 44 |
| Model | not stated for the mtDNA experiment | `evo2_1b_base` |

Three things from the paper's methods that M4 has to take into account:

1. The 65.8% to 5.1% collapse is sensitivity at a fixed threshold
   (Youden-optimal ΔL = −0.0081), not AUROC, and specificity rose to 93.8% at
   the same time. That pattern fits the whole score distribution shifting,
   which a fixed threshold picks up but which needn't hurt the ranking
   (AUROC). The paper reports no AUROC after permutation, so this repo's
   AUROC-based CDI is a new measurement rather than a reproduction. M4 will
   report both sensitivity at a threshold and AUROC.
2. What "cyclic permutation" means. The paper "cyclically permuted the
   positions of all mitochondrial tRNAs while leaving their internal sequences
   intact", which reads as moving tRNA genes into one another's positions.
   CLAUDE.md §5.1 describes rotating the sequence window instead. Those are
   different perturbations, so M4 has to pick one, or run both, and say which.
3. The paper's window is 512 bp either side of the variant (1,025 bp).

## 2026-10-04: batch size changes Evo 2 1B scores more than a variant does

Model: `evo2_1b_base` (1B, FP8, RTX 4060 Laptop GPU). On 8 single-base mtDNA
variants (501 bp windows), changing the batch size from 1 to 8 or 16 moved
per-sequence scores by up to 0.0033, against a median variant effect of
0.00087, and flipped the sign of 3 of the 8 variants. So all scoring in this
repo uses a batch size of 1. Reproduce with `python scripts/check_batch_size.py`;
details are in `docs/model-choice.md`. Single run on 16 sequences.
