"""Does a sequence's Evo 2 score depend on what was scored before it?

Usage (inside the WSL evo2 environment):
    python scripts/check_precision.py

For each precision mode (see seqcontrol/models/evo2.py), loads evo2_1b_base fresh,
then scores one fixed 1,025 bp target window from tests/fixtures/chrM_1_1000.fa
straight after 16 "preceding" sequences of four kinds. With no dependence on history
the target's score is identical every time. The spread is printed next to the median
single-variant effect, for scale. Results are recorded in docs/model-choice.md.
"""

from __future__ import annotations

import gc
from pathlib import Path

import numpy as np
import torch

from seqcontrol.models.base import score_variants
from seqcontrol.models.evo2 import PRECISIONS, Evo2Adapter
from seqcontrol.variants import Variant, variant_windows

FIXTURE = Path(__file__).resolve().parent.parent / "tests" / "fixtures" / "chrM_1_1000.fa"
WINDOW = 501  # the 1 kb fixture cannot hold a 1,025 bp window without wrapping


def main() -> None:
    lines = FIXTURE.read_text().splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    rng = np.random.default_rng(0)
    target = seq[250 : 250 + WINDOW]
    contexts = {
        "none": [],
        "random DNA": ["".join(rng.choice(list("ACGT"), WINDOW)) for _ in range(16)],
        "poly-A": ["A" * WINDOW] * 16,
        "GC-rich": ["".join(rng.choice(list("GGGCCCAT"), WINDOW)) for _ in range(16)],
    }
    variants = []
    for pos in rng.choice(np.arange(260, 740), size=8, replace=False):
        ref = seq[pos - 1]
        variants.append(Variant("MT", int(pos), ref, next(b for b in "ACGT" if b != ref)))
    windows = [variant_windows(seq, v, WINDOW, circular=True) for v in variants]

    for precision in PRECISIONS:
        model = Evo2Adapter("evo2_1b_base", precision=precision)
        model.load()
        target_scores = {}
        for name, before in contexts.items():
            if before:
                model.score_sequences(before)
            target_scores[name] = float(model.score_sequences([target])[0])
        effects = score_variants(model, [r for r, _ in windows], [a for _, a in windows])

        spread = max(target_scores.values()) - min(target_scores.values())
        print(f"\n== {model.label}")
        for name, s in target_scores.items():
            print(f"  target after {name:11s} {s:.6f}")
        print(f"  spread of target score         {spread:.6f}")
        print(f"  median |variant effect|        {np.median(np.abs(effects)):.6f}")
        print(f"  variant effects  {' '.join(f'{x:+.5f}' for x in effects)}")

        del model
        gc.collect()
        torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
