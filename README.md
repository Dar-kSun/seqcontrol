# seqcontrol

Control experiments for genomic language models. The question: when a DNA
model scores a mutation, how much of that score comes from the gene it is
supposed to be reading, and how much from the sequence around it?

On 67 mitochondrial tRNA variants, Evo 2 1B separates pathogenic from benign
with an AUROC of 0.824. A good part of that turns out to be *which tRNA* the
variant is in. A score that never looks at the variant, and only uses how
often variants in its gene are pathogenic, reaches between 0.627
(leave-one-out, biased low) and 0.901 (in-sample, biased high). Only 45
pathogenic–benign pairs share a gene, which is too few to pull gene identity
apart from the effect of the variant itself.

Where there are enough pairs, the model does read the variant. On
protein-coding variants (893 within-gene pairs) it ranks pathogenic above
benign within the same gene with an AUROC of 0.934 [0.847, 0.978]. With a
wider, noisier set of tRNA benign labels (544 pairs) it does the same inside
tRNA genes, at 0.789 [0.698, 0.865].

With that caveat in mind: scrambling everything around the gene, while keeping
the gene and the base composition of its surroundings, lowers the AUROC to
0.710. That is a context-dependence index of 0.35 [0.01, 0.63]. Individual
variant scores move a lot more than the ranking does (Spearman 0.47 against
the unscrambled scores).

Every number here is for Evo 2 `evo2_1b_base` (1B parameters), the smallest
Evo 2 checkpoint, run on a single laptop GPU. The larger Evo 2 models were not
tested.

## What it's for

DNA language models like Evo 2 are increasingly used to rank mutations by how
damaging they are likely to be, including the millions of "variants of
uncertain significance" that sit unresolved in genetic databases. These models
are usually judged by accuracy on a benchmark. A good benchmark score doesn't
tell you *why* the model got there, and if the score really comes from the
surrounding sequence, or from recognising which gene a variant is in, it can
look good on the benchmark and fail on the cases that matter.

seqcontrol is the set of control experiments a wet-lab biologist would run
before trusting an assay, applied to a model:

- **Before relying on a model's variant scores**, run the controls on a
  labelled set where you know the answers, and see how much of the signal
  survives when the context changes but the gene doesn't.
- **When comparing models or checkpoints**, compare how much of each one's
  accuracy depends on context, not only the accuracy itself. The one-page
  [trust card](results/trust_card_tRNA.md) summarises this per region.
- **When reading a benchmark result**, check whether it holds up within genes
  and at more than one threshold. This repo found, for example, that a
  dramatic published drop in sensitivity is partly a property of the fixed
  threshold rather than of the model.

It measures models; it doesn't try to beat them, and nothing here is meant for
clinical use. The controls and metrics don't depend on the model. Adding
another model means writing a small adapter class (a name, a context length,
`load()` and `score_sequences()`; see `seqcontrol/models/base.py`) and
pointing the scripts at it, since they currently build the Evo 2 adapter
directly.

![Flank-shuffle sweep](results/flank_sweep.png)

*Left: how closely per-variant scores track the original scores as more of the
real flanking sequence is kept around each tRNA (the rest is
dinucleotide-shuffled). Right: pathogenic-vs-benign AUROC at the same radii.
Bands are 95% bootstrap intervals over variants; the faint dots are the 10
shuffle seeds.*

## Quickstart

The committed per-variant scores are enough to reproduce the main numbers, with
no GPU and nothing to download:

```bash
pip install -e ".[dev]"
seqcontrol run sweep --from-csv    # flank shuffle: the r = 0 row gives AUROC 0.710, CDI 0.351
```

`seqcontrol run swap --from-csv` and `seqcontrol run rotation --from-csv` do the
same for the other two controls, and `seqcontrol card` rebuilds the one-page
[trust card](results/trust_card_tRNA.md).

Re-scoring from scratch needs a Linux GPU environment with Evo 2
(`scripts/setup_evo2_wsl.sh`, explained in `docs/model-choice.md`):

```bash
seqcontrol run data       # ~200 MB of MITOMAP, ClinVar and Ensembl data; see data/MANIFEST.md
seqcontrol run baseline   # ~2 min on an RTX 4060 Laptop GPU
seqcontrol run swap       # ~10 min
seqcontrol run rotation   # ~7 min
seqcontrol run sweep      # ~40 min
```

## What is measured

The variants are 94 pathogenic (MITOMAP, confirmed pathogenic or likely
pathogenic) and 228 benign (ClinVar, two or more review stars) single-base
changes in human mtDNA (GRCh38). The controls use the 67 that fall in tRNA
genes: 44 pathogenic and 23 benign.

Each variant is scored as ΔL: the mean log-likelihood of a 1,025 bp window
carrying the alternate base, minus that of the same window with the reference
base, centred on the variant. More negative means the model thinks the change
is more damaging.

The context-dependence index (CDI) is the share of above-chance discrimination
that a control takes away:

> CDI = (AUROC_native − AUROC_control) / (AUROC_native − 0.5)

A CDI near 0 means the model is reading the gene; near 1 means nearly all of
the signal came from context. Every estimate has a 95% interval from 2,000
bootstrap resamples of variants (stratified by label, and paired between the
original and the control).

## The controls

All three keep the tRNA's own bases exactly as they are. The code checks this
on every perturbed sequence, and the tests check the check.

- **tRNA swap.** Each tRNA is moved into another tRNA's place on the
  chromosome, taking its own sequence with it, for all 19 cyclic shifts. This
  follows the benchmark that motivated the repo (Mathur & Sachidanandam 2026).
- **Window rotation.** The scoring window is rotated, which rearranges the
  context on each side and joins the two ends of the window together. Thirteen
  offsets, none of which cut the gene. The content stays the same; only its
  arrangement changes.
- **Flank shuffle.** Everything more than *r* bp from the tRNA is replaced by a
  dinucleotide-preserving (Altschul–Erikson) shuffle of itself, for r from 0
  to 400 bp and 10 seeds each. Composition doesn't change, so this tests
  whether the order of the context matters.

## Results

tRNA variants, `evo2_1b_base` with Evo 2's shipped FP8 recipe, 1,025 bp windows.
The unperturbed AUROC is 0.824 [0.713, 0.915].

| Control | AUROC under control | CDI | Spearman ΔL vs native |
|---|---|---|---|
| Flank shuffle, r = 0 (all context scrambled) | 0.710 [0.615, 0.803] | 0.35 [0.01, 0.63] | 0.47 |
| Flank shuffle, r = 100 bp | 0.765 [0.652, 0.861] | 0.18 [−0.12, 0.44] | 0.77 |
| Flank shuffle, r = 400 bp | 0.797 [0.688, 0.888] | 0.09 [−0.06, 0.22] | 0.93 |
| tRNA swap | 0.750 [0.661, 0.826] | 0.23 [−0.12, 0.45] | 0.56 |
| Window rotation | 0.793 [0.692, 0.877] | 0.09 [−0.05, 0.23] | 0.90 |

Read the Spearman column against 0.95, which is how well two FP8 rounding
recipes agree with no control applied at all. Across all 322 variants the
unperturbed AUROC is 0.856 [0.805, 0.903].

A few things stand out.

1. Context moves the scores much more than it moves the ranking. The tRNA swap
   and the stronger flank shuffles push per-variant scores well below the
   rounding-noise floor (only r = 400 bp gets close to it), but AUROC stays
   well above chance throughout. The full scramble at r = 0 is the only
   control whose CDI interval excludes zero, and only just.
2. The context that matters is nearby. Score stability climbs steadily as more
   real flank is kept and reaches 0.93 at 400 bp.
3. Fixed thresholds make the effect look bigger than it is. Under the tRNA
   swap, sensitivity at a cut-off set on the original scores falls from 0.82
   to 0.56, while AUROC only falls from 0.824 to 0.750. The source benchmark
   reports its collapse (65.8% to 5.1%) as sensitivity at a fixed threshold.
   This repo sees the same direction but nothing like the same size, and a
   pre-declared analysis (below) finds that score compression explains much of
   the drop but not all of it.
4. The FP8 rounding recipe doesn't change the conclusions. All three controls
   were re-run with a different FP8 scaling recipe. The CDIs for the tRNA swap,
   window rotation and full scramble moved by at most 0.02, and at
   intermediate shuffle radii by up to 0.06, well inside their intervals.

### Gene identity and the threshold question

Both analyses were written down in
[`docs/plan-v0.2-threshold-and-gene.md`](docs/plan-v0.2-threshold-and-gene.md)
and committed before they were run, and the code applies the decision rules
from that plan (`scripts/05_gene_confound.py`, `scripts/04_threshold_artefact.py`).

For gene identity the confidence intervals come from resampling whole genes:

| Score | AUROC on tRNA variants |
|---|---|
| Evo 2 1B (−ΔL) | 0.824 [0.724, 0.918] |
| Gene prior, ignoring the variant (leave-one-out, biased low) | 0.627 [0.306, 0.802] |
| Gene prior, ignoring the variant (in-sample, biased high) | 0.901 [0.770, 0.969] |
| Evo 2 1B, within-gene pairs only (45 pairs, descriptive) | 0.867 [0.704, 0.938] |

By the plan's rules, gene identity is a substantial part of the AUROC, and this
dataset can't separate gene identity from the effect of the variant.

For the threshold question: shrinking every score toward zero leaves AUROC
exactly where it was, but drags scores across any fixed cut-off. At the
original cut-off (ΔL ≤ −0.0030):

| | Sensitivity | Specificity | AUROC |
|---|---|---|---|
| Native | 0.818 | 0.826 | 0.824 |
| Pure compression (native scores × 0.77) | 0.636 | 0.870 | 0.824 |
| Observed tRNA swap | 0.557 | 0.828 | 0.750 |

Pure compression accounts for 0.697 [0.23, 1.15] of the sensitivity drop. That
is just under the 0.70 the plan set in advance, so the verdict is
indeterminate. Compression would also make specificity rise, and under the swap
it didn't. A simulation of the paper's operating point (synthetic scores, not
the paper's data) shows that compression on its own is enough to turn 65.8%
sensitivity into 5.1% without moving the AUROC.

### Within genes (pre-declared, v0.3)

A second plan,
[`docs/plan-v0.3-within-gene.md`](docs/plan-v0.3-within-gene.md), also
committed before running, looked for variant sets with enough within-gene
pairs (`scripts/06_within_gene.py`; intervals from resampling whole genes).

| Variant set | Within-gene pairs | Gene prior (leave-one-out to in-sample) | Model, within genes |
|---|---|---|---|
| Protein-coding, strict labels | 893 | 0.542 to 0.697 | 0.934 [0.847, 0.978] |
| tRNA, benign widened to ClinVar 1+ star | 544 | 0.693 to 0.834 | 0.789 [0.698, 0.865] |

On protein-coding variants the plan's verdict is that the model separates
pathogenic from benign variants within the same gene, and a gene-only score
gets nowhere near it. The widened tRNA labels are noisier and are a
sensitivity analysis, not the headline. With them, the within-gene
context-dependence index under the tRNA swap is 0.13 [−0.11, 0.43] and under
the full scramble 0.11 [−0.18, 0.53]; both verdicts are indeterminate. So the
within-gene tRNA signal exists, but how much of it depends on context can't be
pinned down at this size.

Everything else, including what was checked along the way and what went wrong,
is in [`docs/findings.md`](docs/findings.md).

## Limitations

- **One small checkpoint.** Only `evo2_1b_base` was tested, on an Ada laptop
  GPU (RTX 4060, 8 GB) running FP8, which the Evo 2 authors document for Hopper
  GPUs. The scores pass sanity checks but haven't been compared with Hopper
  output. The 7B and 40B models from the Evo 2 paper may behave differently.
- **Small samples.** 67 tRNA variants, only 23 of them benign. Most CDI
  intervals are wide and include zero, and a null result at this size is not
  evidence that there is no effect.
- **Labels cluster by gene.** MT-TL1 alone has 13 pathogenic variants and no
  benign ones. A gene prior that ignores the variant scores AUROC 0.627–0.901
  against the model's 0.824, and with only 45 within-gene pairs the strict
  tRNA set can't say how much of the model's signal comes from within genes.
  Every strict-label control result carries this caveat. Wider labels and the
  protein-coding variants show within-gene signal, but the wider tRNA labels
  are single-submitter calls and noisier.
- **Most benign variants are common.** ClinVar-benign mtDNA variants are mostly
  population polymorphisms, so the model could partly be scoring how familiar
  an allele looks. A rough check finds no clear sign of this: among benign
  variants, ΔL correlates with allele frequency at Spearman 0.109
  [−0.022, 0.240].
- **The controls aren't perfect either.** The tRNA swap also changes strand
  context (genes keep their reference orientation) and which tRNAs neighbour
  each other. Window rotation adds an artificial junction. Flank shuffling
  keeps dinucleotides but breaks up longer motifs, so it removes more than
  "context" in the loose sense.
- **Per-variant scores are noisy.** Two FP8 recipes agree only to Spearman
  0.94 (all variants) or 0.95 (tRNA) per variant and flip the sign of 10% of
  variants, even though their AUROCs differ by just 0.007. Sequences are scored
  one at a time because batching changes the scores (`docs/model-choice.md`).
- **One window size.** The headline uses 1,025 bp windows. At 4,097 bp the
  model discriminates no better, rounding noise is about four times larger
  relative to variant effects, and shuffling context more than ~500 bp away has
  no detectable effect (`docs/findings.md`).
- **One scoring run per configuration.** The shuffles use 10 seeds and the
  bootstrap covers variant sampling, but nothing else is repeated.
- **No clinical claims.** This repo measures how a model behaves on a
  benchmark. It doesn't predict anything about patients.

## Status

v0.2, in progress. Done: data loaders, the Evo 2 adapter, the baseline, the
three controls, the trust card and the command-line tool (v0.1, milestones
M0–M7 in `CLAUDE.md`), plus the pre-declared gene-confound, threshold and
within-gene analyses. Not yet implemented: the synonymous-variant control, the
genetic-code check, a second model, and nuclear (ClinVar) variants. The most
useful next step is a larger tRNA set with high-confidence benign labels, so
the within-gene context question can be answered without leaning on
single-submitter calls.

## Citations

- Mathur & Sachidanandam (2026). *Benchmarking DNA Foundation Models:
  Biological Blind Spots in Evo2 Variant-Effect Prediction.* bioRxiv
  [10.64898/2026.03.10.710786](https://www.biorxiv.org/content/10.64898/2026.03.10.710786v1).
  The tRNA permutation design and the 65.8% to 5.1% result this repo builds on.
- Brixi et al. (2026). *Genome modelling and design across all domains of life
  with Evo 2.* Nature
  ([s41586-026-10176-5](https://www.nature.com/articles/s41586-026-10176-5)).
  The model.
- Altschul & Erikson (1985). *Significance of nucleotide sequence alignments.*
  Mol Biol Evol 2:526. The dinucleotide-preserving shuffle.
- Kandel et al. (1996). *Shuffling biological sequences.* Discrete Appl Math
  71:171. The Eulerian-path construction used here.
- MITOMAP (mitomap.org), ClinVar (ncbi.nlm.nih.gov/clinvar) and Ensembl
  (GRCh38) for variants and sequence; versions and checksums are in
  `data/MANIFEST.md`.

## Licence

MIT. See `LICENSE`.
