"""The interface every model adapter satisfies, and the variant score built on it.

An adapter only has to say how to score whole sequences. Variant scoring is the same
for every model, so it lives here, not in the adapters.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Protocol

import numpy as np


class ModelAdapter(Protocol):
    name: str  # exact checkpoint name, e.g. "evo2_1b_base"; printed wherever results appear
    max_context: int  # longest sequence (bp) the model can score

    def load(self) -> None:
        """Load weights. Called once before scoring; may be slow."""
        ...

    def score_sequences(self, seqs: Sequence[str]) -> np.ndarray:
        """Mean log-likelihood per base (natural log) for each sequence, shape (len(seqs),)."""
        ...


def score_variants(
    model: ModelAdapter, ref_windows: Sequence[str], alt_windows: Sequence[str]
) -> np.ndarray:
    """Variant-effect score: mean log-likelihood of the alt window minus that of the ref.

    More negative means the mutation makes the sequence less probable under the model,
    i.e. the model predicts a more damaging variant.
    """
    if len(ref_windows) != len(alt_windows):
        raise ValueError("ref_windows and alt_windows must be the same length")
    for ref, alt in zip(ref_windows, alt_windows, strict=True):
        if len(ref) != len(alt):
            raise ValueError("each ref window must be the same length as its alt window")
        if len(ref) > model.max_context:
            raise ValueError(f"window of {len(ref)} bp exceeds {model.name} context")
    scores = model.score_sequences([*ref_windows, *alt_windows])
    n = len(ref_windows)
    return np.asarray(scores[n:], dtype=float) - np.asarray(scores[:n], dtype=float)
