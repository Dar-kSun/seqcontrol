"""Evo 2 adapter. Needs the Linux environment from scripts/setup_evo2_wsl.sh."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

# Context length per checkpoint, from the Evo 2 model cards. Only checkpoints that have
# been run on this repo's hardware are listed (see docs/model-choice.md).
CONTEXT = {"evo2_1b_base": 8_192}


class Evo2Adapter:
    """Scores one sequence per forward pass.

    Batching is deliberately not offered. On evo2_1b_base (FP8), scoring the same
    sequences with batch size 8 or 16 instead of 1 moved scores by up to 3.3e-3, more
    than a typical single-variant effect (median ~9e-4), and flipped some variants'
    signs. Batch size 1 is exactly reproducible. See docs/model-choice.md.
    """

    def __init__(self, name: str = "evo2_1b_base") -> None:
        if name not in CONTEXT:
            raise ValueError(f"unknown or untested Evo 2 checkpoint: {name}")
        self.name = name
        self.max_context = CONTEXT[name]
        self._model = None

    def load(self) -> None:
        from evo2 import Evo2  # imported here so the package works without evo2 installed

        self._model = Evo2(self.name)

    def score_sequences(self, seqs: Sequence[str]) -> np.ndarray:
        if self._model is None:
            raise RuntimeError("call load() before scoring")
        scores = self._model.score_sequences(list(seqs), batch_size=1)
        return np.asarray(scores, dtype=float)
