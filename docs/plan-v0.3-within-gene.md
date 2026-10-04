# Plan v0.3 — pre-declared: does the model read variants *within* genes?

**Status: pre-declared. Committed before anything in it is run.**

Written by Claude (the coding assistant) on 2026-10-04, after v0.2's Analysis B
found that gene identity is a substantial part of the tRNA AUROC and that the
strict tRNA set has too few within-gene pairs (45) to separate gene identity
from variant effect.

**Known before commit, from labels and gene membership only (no model
scores).** These counts chose the design:

| Variant set | Within-gene pathogenic–benign pairs | Genes with both classes |
|---|---|---|
| Strict tRNA set (v0.1–v0.2) | 45 | 9 |
| Protein-coding, strict labels (already scored in M3) | 893 | 13 of 13 |
| tRNA, benign widened to ClinVar ≥ 1 star | 544 | 16 |
| tRNA, benign widened by population frequency (AF ≥ 0.1–1%) | 48–86 | 11–15 |

The frequency-based option is dropped: too few pairs.

---

## Arm 1 — protein-coding variants, within genes (primary; no GPU)

Uses the committed M3 scores (`results/baseline_fp8-delayed.csv`, region
`protein_coding`: 48 pathogenic, 193 benign, 13 genes). New script
`scripts/06_within_gene.py`, reusing the v0.2 metrics.

1. Gene-prior AUROC, leave-one-variant-out and in-sample (as v0.2 B1).
2. **Within-gene AUROC** of the model (as v0.2 B2), with a cluster bootstrap
   over genes (2,000 resamples).
3. Leave-one-gene-out range of the within-gene AUROC.

**Decision rule for Arm 1 (fixed now):**
- Within-gene AUROC lower 95% bound **> 0.5** → *the model separates
  pathogenic from benign variants within the same protein-coding gene*; gene
  identity cannot be the whole story there.
- Lower bound ≤ 0.5 → *no evidence of within-gene discrimination on
  protein-coding variants.*

No context control is run on protein-coding variants in this plan: the genes
are longer than the scoring window, so "hold the gene, change its context"
does not apply as built.

## Arm 2 — tRNA variants with a widened benign set (sensitivity analysis; GPU)

Benign labels widened from ClinVar ≥ 2 stars to **≥ 1 star** (adds single-
submitter Benign / Likely benign). Pathogenic labels unchanged (MITOMAP
confirmed). These labels are noisier; results here are a sensitivity analysis
and **do not replace** the strict-label headline.

Implementation: a `--labels expanded` option for `01_baseline.py`,
`02_permutation.py` (`--control trna-swap`) and `03_flank_sweep.py`, writing
files with an `_expanded` suffix. Same model, recipe (fp8-delayed) and 1,025 bp
windows. Flank shuffle at r = 0, 100 and 400 bp, 10 seeds.

Report, with gene-cluster bootstrap intervals (2,000 resamples):

1. Native whole-set AUROC and within-gene AUROC; gene-prior AUROC (both
   versions).
2. For the tRNA swap and the flank shuffle at r = 0: **within-gene AUROC under
   the control**, and the **within-gene CDI**
   = (within-gene AUROC_native − within-gene AUROC_control) /
   (within-gene AUROC_native − 0.5), paired against native.
3. The same whole-set metrics as v0.1, for comparison with the strict set.

**Decision rules for Arm 2 (fixed now), applied to the within-gene CDI of
each control:**
- Lower 95% bound **> 0** → *within genes, the control removes a measurable
  part of the model's discrimination.*
- Upper 95% bound **< 0.30** and lower bound ≤ 0 → *within-gene
  discrimination is largely independent of this context.*
- Otherwise → *indeterminate at this sample size.*

The within-gene rules apply only if native within-gene AUROC's lower bound is
above 0.5. If it is not, there is no within-gene signal to lose, and the CDI
is reported but not interpreted.

## Rules for both arms

- Gene-cluster bootstrap wherever gene clustering is in play; redraw resamples
  without a usable within-gene pair or without both classes, and report the
  mean number of draws.
- Provenance (commit, seed, package versions) in each results JSON.
- Every number in the README is checked against the results JSON by script.
- README changes only after both arms finish. If Arm 1 finds within-gene
  discrimination, the README says so next to the gene-identity caveat. If Arm
  2's within-gene CDI is interpretable, it is reported as a sensitivity
  analysis on widened labels, never as the headline.
