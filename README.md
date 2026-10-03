# seqcontrol

Control experiments for genomic language models: measuring how much of a DNA
model's variant-effect score comes from the sequence it claims to be reading,
versus the context around it.

**Status: v0.1, in development.** There are no results yet. Every milestone
below is not yet implemented unless marked done.

| Milestone | Status |
|---|---|
| M0: package skeleton, tests, CI | done |
| M1: model adapter + smoke test | not yet implemented |
| M2: data loaders + manifest | not yet implemented |
| M3: native baseline AUROC (mtDNA) | not yet implemented |
| M4: cyclic-permutation (context-swap) control | not yet implemented |
| M5: flank-shuffle sweep | not yet implemented |
| M6: trust card + CLI | not yet implemented |

## Development

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate    macOS/Linux: source .venv/bin/activate
pip install -e ".[dev]"
ruff check .
pytest
```

## Licence

MIT. See `LICENSE`.
