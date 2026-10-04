# Plan v0.2 — pre-declared: the threshold artefact, and the gene-identity confound

**Status: pre-declared. Written before any of these numbers were computed.**
Commit this file *before* running anything in it, as `docs/plan-v0.2-*.md` was
done for the earlier milestones. Pre-declaring the analyses and the decision
rules is what stops this turning into a search for a flattering result.

Two analyses. **A** explains the published collapse. **B** tests whether this
repo's own headline number survives its biggest confound. B can invalidate A's
framing, so if time is short, **run B first.**

---

## Background: what we already know

From `docs/findings.md`:

- Native tRNA AUROC: **0.824** [0.713, 0.915], n = 67 (44 pathogenic, 23 benign).
- tRNA swap: AUROC **0.750**, drop 0.074 [−0.030, 0.162] (not significant).
- Sensitivity at the native-Youden cut-off (ΔL ≤ −0.0030): **0.82 → 0.56**;
  specificity **0.83 → 0.83** (unchanged).
- Median ΔL, pathogenic: **−0.0056 → −0.0036**. Benign: −0.0013 → −0.0007.
- The source paper reports sensitivity **65.8% → 5.1%** at a fixed cut-off
  (ΔL ≤ −0.0081) with **specificity rising to 93.8%**. It reports no AUROC
  after permutation.
- Noise floor: changing only the FP8 rounding recipe gives per-variant
  Spearman ≈ 0.95. Nothing below that is a real effect.

**The hypothesis both analyses serve:** the published collapse is largely a
*threshold artefact*. Moving a gene compresses every score toward zero, so
fewer variants cross a fixed cut-off. A fixed threshold registers that as
catastrophic; the ranking need not change much.

---

# Analysis A — Is the collapse a threshold artefact?

No GPU. Operates on committed per-variant scores. Should run in seconds.

New script: `scripts/04_threshold_artefact.py`. New module:
`seqcontrol/artefact.py`.

## A0. The mathematical fact to exploit (state it in the write-up)

AUROC depends only on the *ordering* of scores. Multiplying every score by a
constant α > 0 is strictly monotone, so **AUROC is unchanged exactly** — not
approximately. But a cut-off fixed in absolute units does change: scoring
ΔL ≤ t after scaling by α is identical to scoring ΔL ≤ t/α on the native
scores, i.e. a moving threshold.

So compression is a scenario where **sensitivity can collapse while
discrimination is untouched, by construction.**

**Implementation check:** assert `|AUROC(α) − AUROC(1)| < 1e-12` for every α.
If that assertion fails, the metric code has a bug. This is a test, not a result.

## A1. The compression sweep (the headline figure)

For α in `[1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]`, applied to
native ΔL:

Report, per α:
- Sensitivity and specificity at the **paper's** cut-off (ΔL ≤ −0.0081).
- Sensitivity and specificity at **this repo's native-Youden** cut-off (−0.0030).
- AUROC (flat by construction — plot it to make the point visually).

**Figure (`results/threshold_artefact.png`):** α on x. Sensitivity and
specificity as lines; AUROC as a flat line across the top. Mark the α that
reproduces the observed swap sensitivity. One glance should say: *the cliff is
the threshold, not the model.*

## A2. How much compression did the real swap apply?

Estimate the scale factor α̂ relating swap scores to native scores, **two ways**,
and report both:

1. **Robust scale ratio:** `median(|ΔL_swap|) / median(|ΔL_native|)`.
2. **Regression through the origin:** fit `ΔL_swap ≈ α·ΔL_native`, report α̂ with
   a bootstrap CI and the R² of that fit.

The R² is important and must be reported prominently: it says how well a
*pure scale change* describes what the swap did. Per-variant Spearman under the
swap is 0.56, well below the 0.95 noise floor, so a pure-scale model is known to
be incomplete. **Quantify the shortfall rather than hiding it.**

## A3. The decomposition — the actual scientific claim

Three quantities, reported as one table:

| | Sensitivity at fixed cut-off | AUROC |
|---|---|---|
| Native | 0.82 | 0.824 |
| **Pure-compression null** (native scores × α̂) | *compute* | 0.824 *(identical by construction)* |
| Observed swap | 0.56 | 0.750 |

Then state plainly:

- **Share of the sensitivity drop explained by compression alone:**
  `(sens_native − sens_null) / (sens_native − sens_swap)`.
- **Residual ranking loss beyond compression:** `AUROC_native − AUROC_swap`
  = 0.074 [−0.030, 0.162], which at n = 67 is **not significant**.

### Pre-declared interpretation rules — fix these now

- If compression explains **≥ 70%** of the sensitivity drop and the AUROC drop
  stays non-significant → conclude: *on this setup, the collapse is primarily a
  threshold artefact; genuine ranking loss is small and not resolvable at
  n = 67.*
- If compression explains **< 40%** → the hypothesis is wrong. Say so, and
  report that the swap genuinely degrades discrimination.
- Anything between → report as indeterminate at this sample size. **Do not
  pick whichever framing reads better.**

## A4. Rank-preserving surrogate (isolating scale from rank)

Build a surrogate: take the swap's score *distribution* but the native
*ordering* (sort native scores, assign swapped-score values by rank). By
construction it has the swap's scale and the native ranking.

- If the surrogate's sensitivity at the fixed cut-off ≈ the real swap's →
  the sensitivity change is entirely a scale effect.
- Any gap is the part attributable to reordering.

## A5. Speaking to the paper's 65.8% → 5.1% specifically

**This cannot be reproduced on this dataset** — native sensitivity here at the
paper's cut-off is 22.7%, not 65.8%. The operating points differ.

So do this as an explicitly-labelled **illustration, not a reproduction**:

1. Construct two Gaussian score distributions calibrated so that at the
   Youden-optimal cut-off, sensitivity = 65.8% and specificity matches the
   paper's reported native value.
2. Apply compression α and recompute sensitivity and specificity at the
   **fixed** cut-off.
3. Report the α that lands at 5.1% sensitivity, and confirm specificity *rises*
   (the paper reports 93.8%) — the rising-specificity signature is the strongest
   evidence for compression, because real discrimination loss does not push both
   metrics the same direction.
4. Confirm AUROC is unchanged throughout.

**Labelling requirement:** every mention of this says it is a simulation of the
*mechanism* under the paper's reported operating point, using synthetic scores.
It is **not** a claim about the paper's data, and it is not evidence the paper is
wrong. The paper measured a real effect; this illustrates what produces that
pattern.

---

# Analysis B — Does gene identity explain the native AUROC?

**Run this first.** If native AUROC is mostly gene identity, A's framing needs
restating, because then the controls are disrupting gene recognition rather than
variant scoring.

The problem: MT-TL1 contributes 13 pathogenic and 0 benign variants. A model
that only recognised *which gene it is in* could score well without reading the
variant at all.

New script: `scripts/05_gene_confound.py`.

## B1. The gene-prior baseline (do this first — it's one function)

A "model" that ignores the variant entirely and scores each variant by its
gene's pathogenic fraction in this dataset. Evaluate with **leave-one-gene-out**
(a gene's own variants never inform its own prior).

This single number reframes everything:

- Gene-prior AUROC **≈ 0.80** → native 0.824 is largely gene identity. Say so
  prominently, in the README, not only here.
- Gene-prior AUROC **≈ 0.5–0.6** → gene identity is a minor part, and the
  headline stands as a statement about variants.

## B2. Within-gene (stratified) AUROC

Compute AUROC using **only pathogenic–benign pairs from the same gene**
(a stratified Mann–Whitney statistic). This removes gene identity by
construction.

- **Report the number of usable within-gene pairs first.** With 67 variants
  across ~22 genes, many genes carry one class only. If usable pairs are few
  (say < 100), the honest conclusion is *this dataset cannot separate gene
  identity from variant effect*, and that is a legitimate, publishable-shaped
  finding about the benchmark — not a failure.
- CIs by **cluster bootstrap over genes**, not over variants.

## B3. Leave-one-gene-out sensitivity

Recompute native AUROC with each gene dropped in turn. Report the full range
and name the most influential gene. If dropping MT-TL1 moves AUROC by more than
~0.05, the headline is substantially one gene's doing and the README must say so.

## B4. Re-run the controls under the gene-robust metric

Recompute the swap, rotation and flank-shuffle results using the within-gene
AUROC from B2 (if B2 has enough pairs). State whether the CDI conclusions hold.

## B5. The familiarity confound — scope it, don't solve it

Most benign mtDNA variants here are common population polymorphisms, so the
model may be scoring *how familiar this spelling looks* rather than functional
effect.

Cheap partial check: if allele frequency is available for the benign set,
report the correlation between ΔL and allele frequency within benign variants
alone. A strong correlation means part of the signal is familiarity.

**Do not attempt to fully resolve this.** Measure it, state it as a limitation,
move on.

---

# Statistical rules (both analyses)

- **Bootstrap over genes, not variants**, anywhere gene clustering is in play
  (B2, B3, B4). Variant-level bootstrapping treats 13 MT-TL1 variants as 13
  independent facts, which is the exact error this analysis exists to catch.
- 2,000 resamples, paired against native, matching the existing convention.
- Every reported difference carries a CI. No bare point estimates.
- Fix seeds; record commit, seed and package versions in `run.json` as the
  existing scripts do.

# Deliverables

1. `seqcontrol/artefact.py`, `scripts/04_threshold_artefact.py`,
   `scripts/05_gene_confound.py`.
2. `results/threshold_artefact.{csv,json,png}`,
   `results/gene_confound.{csv,json}`.
3. New `docs/findings.md` sections for A and B, in the existing format
   (command, model, numbers with CIs, what it shows, limitations).
4. README updated **only after both complete**, with:
   - the decomposition table from A3,
   - the gene-prior AUROC from B1 stated near the headline, whatever it says,
   - the headline rewritten if B1 or B3 undermines it.

# Honesty requirements specific to this work

- The repo's position is **"here is the mechanism behind the published
  number,"** never "the paper is wrong." The paper measured a real effect at a
  real operating point. Setups differ (model size unstated there, different
  benign set, different threshold, possibly a different permutation).
- A5 is labelled a simulation everywhere it appears.
- If B1 shows the native AUROC is mostly gene identity, **that becomes the
  headline finding** and the context result is demoted. Agreeing to that now,
  before seeing the number, is the point of pre-declaring.
- Report A3's decomposition whichever way it falls, per the rules in A3.

# Definition of done

A reader can see, in one table, how much of a published sensitivity collapse is
explained by score compression rather than lost discrimination — and can see
whether this dataset is even able to tell variant effects from gene identity.
