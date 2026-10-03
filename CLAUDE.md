# CLAUDE.md — `seqcontrol`

> Control experiments for genomic language models.
> **Mission: measure how much of a DNA model's variant-effect score comes from the sequence it claims to be reading, versus the context around it.**

---

## 0. Context for you, the agent

This repo is being built by Aryan Sinha (3rd-year B.Tech Biotechnology, IIT Kharagpur). It will be read by people who evaluate research for a living — reviewers from frontier AI labs. That sets the bar:

- **Honest beats impressive.** A clean negative result, clearly reported, is worth more than an inflated positive one.
- **Every number in the README must be reproducible** by running a command in this repo.
- **Never claim a result the code has not produced.** If a milestone isn't done, the README says "not yet implemented."
- Small and finished beats large and half-built.

Work in small, reviewable commits with real messages. Never backdate commits.

---

## 1. The problem, stated precisely

DNA is a string over {A,C,G,T}. A point mutation changes one letter. Most are harmless; a few cause disease. Distinguishing them is the "variants of uncertain significance" problem — millions sit unresolved in clinical records.

Genomic language models (Evo 2, Nucleotide Transformer, HyenaDNA) are trained autoregressively on genomes. Their variant-effect score is typically a **log-likelihood ratio**: how much less probable the sequence becomes under the mutation. The premise is that evolutionary constraint shows up as model confidence.

**The documented failure.** A March 2026 benchmark stress-tested Evo 2 on mitochondrial DNA (a trustworthy answer key — mtDNA pathogenic variants are well catalogued in MITOMAP). They **cyclically permuted** tRNA genes: rotate the start position so the gene's content is preserved but its genomic context changes. Variant-effect performance **collapsed from 65.8% to 5.1%**.

The model was reading the *address*, not the gene. Related findings from the same and adjacent work:
- Evo 2 selected the human-preferred codon only **24.4%** of the time (≈ chance).
- It called normal mitochondrial start codons pathogenic **100%** of the time — mitochondria use a variant genetic code the model never internalised.
- A Jan 2026 paper found generated sequences drift from real genome structure at long range (a classifier separates real from synthetic with **AUROC 0.93** by 175 kb).

**The gap.** The field benchmarks on accuracy. It almost never runs **controls** — the discipline every wet-lab student learns in week one. The permutation trick above was one experiment, in one paper, on one gene family. Nobody packaged it as a reusable harness.

**What this repo is.** That harness.

---

## 2. Scope

### v0.1 — must ship
1. Score variants with at least one genomic LM and reproduce a sane baseline AUROC on a labelled set.
2. **Context-swap control** (cyclic permutation) with the performance delta reported.
3. **Flank-shuffle control**: randomise context at increasing radii; plot score stability vs radius.
4. A **trust card**: one-page per-region summary with a context-dependence index.
5. CLI + README with reproducible commands and real numbers.

### v0.2 — if time allows
6. **Synonymous control**: variants that provably don't change the protein should mostly score low. Report violation rate.
7. **Genetic-code awareness check**: flag regions using non-standard codes (mitochondria) and quantify systematic error there.
8. A second model, for comparison across architectures.

### Explicitly out of scope
- Training or fine-tuning any genomic model.
- Clinical interpretation or any claim about patient care.
- Beating state-of-the-art variant-effect prediction. **This repo measures models; it does not compete with them.**

---

## 3. Compute reality check — do this FIRST

Before writing analysis code, confirm what actually runs on the available hardware.

1. Try the **smallest** checkpoint first (e.g. an Evo 2 1B-parameter base model). Do not start with 7B or 40B.
2. If local inference is infeasible, evaluate hosted inference APIs, or fall back to a smaller genomic LM (HyenaDNA, Nucleotide Transformer).
3. **Whatever you end up using, name it explicitly everywhere** — README, trust cards, plots. Never let a reader assume a bigger model was tested than actually was.
4. Record the decision and the reasoning in `docs/model-choice.md`, including timings and memory.

Write `scripts/00_smoke_test.py` that loads the model, scores one 1 kb sequence, prints wall-clock time and peak memory, and exits. Commit it before anything else.

---

## 4. Data

Fetch programmatically into `data/raw/` (gitignored). Every loader must cache and be re-runnable.

| Need | Source | Notes |
|---|---|---|
| mtDNA pathogenic/benign variants | MITOMAP | The cleanest answer key; start here |
| Nuclear variant labels | ClinVar (pathogenic vs benign, high review status) | Filter by review stars; document the filter |
| tRNA gene coordinates | GtRNAdb | For permutation targets |
| Reference sequence | Ensembl REST / UCSC | Pin the genome build (GRCh38) in config |

Rules:
- Pin every download with a date and a checksum in `data/MANIFEST.md`.
- Never commit raw genomic data. Commit the loader and the manifest.
- Keep a tiny committed fixture (a few hundred bp + 20 variants) in `tests/fixtures/` so tests run without network.

---

## 5. The controls — exact specifications

Each control lives in `seqcontrol/controls/` and implements a common interface:

```python
class Control(Protocol):
    name: str
    def apply(self, region: Region, rng: np.random.Generator) -> list[Region]: ...
```

### 5.1 Context swap (cyclic permutation)
Rotate the sequence window so the gene's internal sequence and ordering are preserved while its flanking context changes. Sweep rotation offsets. **Critical:** assert the gene's own bases are byte-identical before and after — a rotation that corrupts the gene invalidates the control. Write that assertion as a test.

### 5.2 Flank shuffle
Hold the gene fixed; dinucleotide-shuffle the flanks beyond radius *r*, for *r* in a sweep (e.g. 0, 50, 100, 250, 500, 1000, 2500 bp). Use **dinucleotide-preserving** shuffling (Altschul–Erikson), not naive shuffling — naive shuffling destroys base composition and confounds the result. Repeat with ≥10 random seeds and report variance.

### 5.3 Synonymous control (v0.2)
Within coding regions, generate codon-synonymous substitutions. Expectation: mostly low scores. Report the fraction scored above the pathogenic threshold. Use the correct codon table per region.

### 5.4 Genetic-code check (v0.2)
For regions using a non-standard genetic code (mitochondrial code 2 for human mtDNA), check whether the model's scores reflect the correct code. Report start/stop codon handling specifically.

---

## 6. Metrics

Define once in `seqcontrol/metrics.py` and reuse:

- **Discrimination**: AUROC and AUPRC for pathogenic vs benign. Report both; AUPRC matters because classes are imbalanced.
- **Context-dependence index (CDI)**: the relative drop in discrimination under a control.
  `CDI = (AUROC_native − AUROC_control) / (AUROC_native − 0.5)`
  CDI near 1 means essentially all signal was context. Near 0 means the model is reading the sequence. **Define this formula in the README and keep it fixed.**
- **Score stability**: Spearman correlation between native and perturbed per-variant scores.
- Confidence intervals by bootstrap over variants (≥1000 resamples). **Never report a bare point estimate.**

---

## 7. Repo layout

```
seqcontrol/
  __init__.py
  models/          # thin adapters: load(), score_variants() -> np.ndarray
    base.py        # Protocol all adapters satisfy
    evo2.py
    hyenadna.py
  data/            # loaders: mitomap.py, clinvar.py, gtrnadb.py, reference.py
  controls/        # permute.py, flank_shuffle.py, synonymous.py, genetic_code.py
  metrics.py
  trustcard.py     # renders the one-page summary (markdown + matplotlib figure)
  cli.py           # typer-based: seqcontrol run / card / smoke
scripts/
  00_smoke_test.py
  01_baseline.py
  02_permutation.py
  03_flank_sweep.py
tests/
  fixtures/
docs/
  model-choice.md
  findings.md      # the evolving results log — update as runs complete
results/           # committed: small JSON/CSV summaries + figures (NOT raw data)
README.md
pyproject.toml
```

**Model adapters matter.** Keep `models/base.py` narrow so a new model is ~40 lines. The repo's value is being model-agnostic.

---

## 8. Build order

Each milestone ends with a commit and a green test suite.

- **M0** — `pyproject.toml`, package skeleton, ruff + pytest configured, CI workflow running tests on push, MIT licence, `.gitignore` (exclude `data/raw/`, model weights, `__pycache__`). README stub saying "v0.1, in development."
- **M1** — Model adapter + `00_smoke_test.py` passing. Commit `docs/model-choice.md`.
- **M2** — Data loaders + manifest. Test on the committed fixture.
- **M3** — Variant scoring end to end; produce a **native baseline AUROC** on the mtDNA set. This is the number everything else is measured against. Record it in `docs/findings.md`.
- **M4** — Cyclic permutation control. **Acceptance: the direction and rough magnitude of the published collapse is reproduced on tRNA variants.** If it is not reproduced, that is itself a finding — investigate, document honestly in `docs/findings.md`, and do not quietly tune until it matches.
- **M5** — Flank-shuffle sweep + the radius plot. This plot is the repo's signature figure.
- **M6** — Trust card renderer + CLI.
- **M7** — README with real numbers, figures, limitations. Tag `v0.1.0`.
- **M8+** (v0.2) — Synonymous control, genetic-code check, second model.

---

## 9. README requirements

Written for someone who has 60 seconds. Required sections:

1. **One-sentence what and why**, with the headline finding as a number.
2. **The signature figure** (score stability vs flank radius) near the top.
3. **Quickstart**: install, then one command that reproduces a headline number.
4. **What the controls are**, each in two sentences.
5. **Results table**: model × region × native AUROC × control AUROC × CDI, with CIs.
6. **Limitations** — a real section, not an afterthought. Include: which checkpoint was tested, what wasn't tested, where the controls themselves could mislead, sample sizes.
7. **Status**: "v0.1 — active development," with what's done and what isn't.
8. Citations to the papers in §1.

Avoid hype words: "revolutionary", "breakthrough", "state-of-the-art". State what was measured.

---

## 10. Honesty checklist — run before every README edit

- [ ] Every number traceable to a command in this repo.
- [ ] Model checkpoint and size named wherever results appear.
- [ ] Uncertainty reported alongside point estimates.
- [ ] Negative/null results included, not quietly dropped.
- [ ] No claim about clinical use.
- [ ] Prior work credited; this repo is positioned as *measuring* models, not beating them.
- [ ] Any result from a single run is labelled as such.

---

## 11. Definition of done for v0.1

A reviewer clones the repo, runs two commands, and sees a reproduced, uncertainty-quantified measurement of how much a genomic model's variant-effect signal depends on context rather than sequence — with the limitations stated plainly.
