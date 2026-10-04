# seqcontrol

Control experiments for genomic language models: how much of a DNA model's
variant-effect score comes from the gene it is supposed to be reading, and how
much from the sequence around it?

**Headline.** On 67 mitochondrial tRNA variants, Evo 2 1B separates
pathogenic from benign with AUROC **0.824**, but a substantial part of that is
*which tRNA* a variant sits in. A score that ignores the variant entirely and
uses only its gene's pathogenic fraction reaches **0.627** (leave-one-out,
biased low) to **0.901** (in-sample, biased high), and only 45 pathogenic–benign
pairs share a gene, too few to separate gene identity from variant effect.

With that caveat: scrambling everything around the gene (keeping the gene and
the context's base composition intact) lowers AUROC to **0.710**, a
context-dependence index of **0.35 [0.01, 0.63]**. Individual variant scores
move much more than the ranking does (Spearman 0.47 against native).

All numbers are for **Evo 2 `evo2_1b_base` (1B parameters)**, the smallest
Evo 2 checkpoint, run on one laptop GPU. Larger Evo 2 models were not tested.

![Flank-shuffle sweep](results/flank_sweep.png)

*Left: how closely per-variant scores track the native scores as more real
flanking sequence is kept around each tRNA (the rest is dinucleotide-shuffled).
Right: pathogenic-vs-benign AUROC at the same radii. Bands are 95% bootstrap
intervals over variants; faint dots are the 10 shuffle seeds.*

## Quickstart

Reproduce the headline numbers from the committed per-variant scores, with no
GPU and no downloads:

```bash
pip install -e ".[dev]"
seqcontrol run sweep --from-csv    # flank shuffle: r = 0 row gives AUROC 0.710, CDI 0.351
```

`seqcontrol run swap --from-csv` and `seqcontrol run rotation --from-csv` do the
same for the other two controls, and `seqcontrol card` rebuilds the one-page
[trust card](results/trust_card_tRNA.md).

To re-score from scratch you need a Linux GPU environment with Evo 2
(`scripts/setup_evo2_wsl.sh`; see `docs/model-choice.md`), then:

```bash
seqcontrol run data       # ~200 MB of MITOMAP, ClinVar and Ensembl data; see data/MANIFEST.md
seqcontrol run baseline   # ~2 min on an RTX 4060 Laptop GPU
seqcontrol run swap       # ~10 min
seqcontrol run rotation   # ~7 min
seqcontrol run sweep      # ~40 min
```

## What is measured

**Variants.** 94 pathogenic (MITOMAP, confirmed P/LP) and 228 benign (ClinVar,
2+ review stars) single-base variants in human mtDNA (GRCh38). The controls
target the 67 in tRNA genes (44 pathogenic, 23 benign).

**Score.** ΔL = mean log-likelihood of a 1,025 bp window carrying the
alternate base minus that of the reference window, centred on the variant.
More negative means more damaging.

**Context-dependence index.** The share of above-chance discrimination a
control removes:

> **CDI = (AUROC_native − AUROC_control) / (AUROC_native − 0.5)**

CDI near 0 means the model is reading the gene; near 1 means essentially all
signal came from context. Every estimate carries a 95% interval from 2,000
label-stratified bootstrap resamples of variants, paired between native and
control.

## The controls

Every control keeps the tRNA's own bases byte-identical; this is asserted on
every perturbed sequence and tested.

- **tRNA swap.** Each tRNA is moved into another tRNA's position on the
  chromosome, carrying its own sequence (all 19 cyclic shifts). This follows
  the design of the benchmark that motivated this repo (Mathur &
  Sachidanandam 2026).
- **Window rotation.** The scoring window is rotated so the context on each
  side is rearranged and an artificial junction is introduced (13 offsets,
  never cutting the gene). It changes arrangement, not content.
- **Flank shuffle.** Everything more than *r* bp from the tRNA is replaced by a
  dinucleotide-preserving (Altschul–Erikson) shuffle, for r = 0 to 400 bp and
  10 seeds. Composition stays fixed, so only the order of the context is
  tested.

## Results

tRNA variants, `evo2_1b_base`, Evo 2's shipped FP8 recipe, 1,025 bp windows.
Native AUROC 0.824 [0.713, 0.915].

| Control | AUROC under control | CDI | Spearman ΔL vs native |
|---|---|---|---|
| Flank shuffle, r = 0 (all context scrambled) | 0.710 [0.615, 0.803] | **0.35 [0.01, 0.63]** | 0.47 |
| Flank shuffle, r = 100 bp | 0.765 [0.652, 0.861] | 0.18 [−0.12, 0.44] | 0.77 |
| Flank shuffle, r = 400 bp | 0.797 [0.688, 0.888] | 0.09 [−0.06, 0.22] | 0.93 |
| tRNA swap | 0.750 [0.661, 0.826] | 0.23 [−0.12, 0.45] | 0.56 |
| Window rotation | 0.793 [0.692, 0.877] | 0.09 [−0.05, 0.23] | 0.90 |

Spearman should be read against **0.95**, the agreement between two FP8
rounding recipes with no control applied. Native AUROC on all 322 variants is
0.856 [0.805, 0.903].

What the results say:

1. **Scores are context-sensitive; rankings much less so.** The tRNA swap and
   the stronger flank shuffles move per-variant scores far past the
   rounding-noise floor (only r = 400 bp comes close to it), yet AUROC stays
   well above chance throughout. Only the full scramble (r = 0) has a CDI interval
   that excludes zero, and only just.
2. **The context that matters is local.** Score stability recovers steadily as
   more real flank is kept, reaching 0.93 at 400 bp.
3. **Threshold metrics exaggerate the effect.** Under the tRNA swap,
   sensitivity at a cut-off fixed on native scores falls from 0.82 to 0.56
   while AUROC falls only from 0.824 to 0.750. The source benchmark reports its
   collapse (65.8% to 5.1%) as sensitivity at a fixed threshold; this repo
   reproduces the direction of that drop but not its size. A pre-declared
   decomposition (below) finds that score compression is a large part of the
   mechanism but not all of it.
4. **Robust to the FP8 rounding recipe.** All three controls were re-run with a
   different FP8 scaling recipe. The CDIs for the tRNA swap, window rotation and
   full scramble change by at most 0.02; at intermediate shuffle radii by up to
   0.06, well inside their intervals.

### Gene identity and the threshold artefact (pre-declared analyses)

Both were specified in
[`docs/plan-v0.2-threshold-and-gene.md`](docs/plan-v0.2-threshold-and-gene.md)
and committed before being run; the decision rules are applied in code
(`scripts/05_gene_confound.py`, `scripts/04_threshold_artefact.py`).

**Gene identity.** With confidence intervals from resampling whole genes:

| Score | AUROC on tRNA variants |
|---|---|
| Evo 2 1B (−ΔL) | 0.824 [0.724, 0.918] |
| Gene prior, ignoring the variant (leave-one-out, biased low) | 0.627 [0.306, 0.802] |
| Gene prior, ignoring the variant (in-sample, biased high) | 0.901 [0.770, 0.969] |
| Evo 2 1B, within-gene pairs only (45 pairs, descriptive) | 0.867 [0.704, 0.938] |

Pre-declared verdicts: gene identity is a **substantial** part of the native
AUROC, and this dataset **cannot separate** gene identity from variant effect.

**Threshold artefact.** Shrinking every score toward zero leaves AUROC exactly
unchanged but drags scores across a fixed cut-off. At the native cut-off
(ΔL ≤ −0.0030):

| | Sensitivity | Specificity | AUROC |
|---|---|---|---|
| Native | 0.818 | 0.826 | 0.824 |
| Pure compression (native scores × 0.77) | 0.636 | 0.870 | 0.824 |
| Observed tRNA swap | 0.557 | 0.828 | 0.750 |

Pure compression reproduces **0.697 [0.23, 1.15]** of the sensitivity drop,
just under the pre-declared 0.70 bar, so the verdict is **indeterminate**.
Compression also predicts specificity should rise, and under the swap it did
not. A simulation of the paper's operating point (synthetic scores, not the
paper's data) shows that compression alone can turn 65.8% sensitivity into
5.1% with AUROC unchanged.

The full log, including what was checked and what went wrong, is in
[`docs/findings.md`](docs/findings.md).

## Limitations

- **One small checkpoint.** Only `evo2_1b_base` was tested, on an Ada laptop
  GPU (RTX 4060, 8 GB) using FP8, which the Evo 2 authors document for Hopper.
  Scores were sanity-checked but not compared against Hopper output. The 7B and
  40B models used in the Evo 2 paper may behave differently.
- **Small samples.** 67 tRNA variants, only 23 of them benign. Most CDI
  intervals are wide and include zero; a null result here is not evidence of no
  effect.
- **Labels cluster by gene.** MT-TL1 alone has 13 pathogenic and no benign
  variants. A variant-blind gene prior scores AUROC 0.627–0.901 against the
  model's 0.824, and with only 45 within-gene pairs this dataset cannot say how
  much of the model's signal is within-gene. Every control result above
  inherits this limitation.
- **Benign variants are mostly common.** ClinVar-benign mtDNA variants are
  largely population polymorphisms, so the model could partly be scoring
  familiarity. A weak check finds no clear sign of it: among benign variants,
  ΔL correlates with allele frequency at Spearman 0.109 [−0.022, 0.240].
- **The controls can mislead too.** The tRNA swap also changes strand context
  (genes keep reference orientation) and which tRNAs sit next to each other.
  Window rotation adds an artificial junction. Flank shuffling keeps
  dinucleotides but destroys longer motifs, so it removes more than "context".
- **Per-variant scores carry rounding noise.** Two FP8 recipes agree only to
  Spearman 0.94 (all variants) or 0.95 (tRNA) per variant and flip 10% of variant signs, though AUROC differs
  by 0.007. Scoring is done one sequence at a time because batching changes
  scores (`docs/model-choice.md`).
- **One window size for the headline.** Results are for 1,025 bp windows. At
  4,097 bp the model discriminates no better, rounding noise is about four times
  larger relative to variant effects, and shuffling context beyond ~500 bp has
  no detectable effect (`docs/findings.md`).
- **Single scoring run per configuration.** Shuffle seeds are repeated (10);
  the bootstrap covers variant sampling; nothing else is replicated.
- **No clinical claims.** This repo measures model behaviour on a benchmark.
  It does not predict anything about patients.

## Status

**v0.2, in progress.** Done: data loaders, Evo 2 adapter, baseline, the three
controls, the trust card and the CLI (v0.1, milestones M0–M7 in `CLAUDE.md`),
plus the pre-declared gene-confound and threshold-artefact analyses. Not yet
implemented: the synonymous-variant control, the genetic-code check, a second
model, and nuclear (ClinVar) variants. A dataset with more within-gene
pathogenic–benign pairs is the most important next step.

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
  Mol Biol Evol 2:526. Dinucleotide-preserving shuffle.
- Kandel et al. (1996). *Shuffling biological sequences.* Discrete Appl Math
  71:171. The Eulerian-path construction used here.
- MITOMAP (mitomap.org), ClinVar (ncbi.nlm.nih.gov/clinvar) and Ensembl
  (GRCh38) for variants and sequence; versions and checksums in
  `data/MANIFEST.md`.

## Licence

MIT. See `LICENSE`.
