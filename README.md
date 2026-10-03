# seqcontrol

Control experiments for genomic language models: measuring how much of a DNA
model's variant-effect score comes from the sequence it claims to be reading,
versus the context around it.

**Status: v0.1, in development.** There are no results yet. Every milestone
below is not yet implemented unless marked done.

| Milestone | Status |
|---|---|
| M0: package skeleton, tests, CI | done |
| M1: model adapter + smoke test (Evo 2 `evo2_1b_base` on an 8 GB laptop GPU; see `docs/model-choice.md`) | done |
| M2: data loaders + manifest (94 pathogenic, 228 benign mtDNA SNVs; see `data/MANIFEST.md`) | done |
| M3: native baseline (AUROC 0.856 [0.805, 0.903] on 322 mtDNA variants, `evo2_1b_base`; see `docs/findings.md`) | done |
| M4: context-swap controls (tRNA swap: AUROC 0.824 to 0.750, CDI 0.23 [−0.12, 0.46], n = 67; see `docs/findings.md`) | done |
| M5: flank-shuffle sweep | not yet implemented |
| M6: trust card + CLI | not yet implemented |

## Development

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate    macOS/Linux: source .venv/bin/activate
pip install -e ".[dev]"
ruff check .
pytest
python scripts/fetch_data.py   # downloads ~200 MB into data/raw/, writes data/MANIFEST.md
```

Scoring needs the Linux GPU environment (`scripts/setup_evo2_wsl.sh`), then:

```bash
python scripts/01_baseline.py  # ~2 min on an RTX 4060 Laptop GPU
```

## Model

All results will use Evo 2 `evo2_1b_base` (1B parameters), the smallest Evo 2
checkpoint, not the 7B or 40B models from the Evo 2 paper. Why, and how it was
set up, is in `docs/model-choice.md`.

## Licence

MIT. See `LICENSE`.
