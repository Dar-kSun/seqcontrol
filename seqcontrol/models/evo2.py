"""Evo 2 adapter. Needs the Linux environment from scripts/setup_evo2_wsl.sh."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

# Context length per checkpoint, from the Evo 2 model cards. Only checkpoints that have
# been run on this repo's hardware are listed (see docs/model-choice.md).
CONTEXT = {"evo2_1b_base": 8_192}

# How the input projections run. Evo 2 ships "fp8-delayed": Transformer Engine's
# DelayedScaling recipe, which sets each layer's FP8 scale from the largest values
# seen over the previous 16 forward passes, so a score can depend on what was scored
# before it. "fp8-current" computes the scale from the current sequence alone.
# "bf16" skips FP8. See docs/model-choice.md for how much each choice matters.
PRECISIONS = ("fp8-delayed", "fp8-current", "bf16")


class Evo2Adapter:
    """Scores one sequence per forward pass.

    Batching is deliberately not offered. On evo2_1b_base (FP8), scoring the same
    sequences with batch size 8 or 16 instead of 1 moved scores by up to 3.3e-3, more
    than a typical single-variant effect (median ~9e-4). See docs/model-choice.md.
    """

    def __init__(self, name: str = "evo2_1b_base", precision: str = "fp8-delayed") -> None:
        if name not in CONTEXT:
            raise ValueError(f"unknown or untested Evo 2 checkpoint: {name}")
        if precision not in PRECISIONS:
            raise ValueError(f"precision must be one of {PRECISIONS}, not {precision!r}")
        self.name = name
        self.precision = precision
        self.max_context = CONTEXT[name]
        self._model = None

    @property
    def label(self) -> str:
        """Checkpoint plus precision, for result files and plots."""
        return f"{self.name} ({self.precision})"

    def load(self) -> None:
        from evo2 import Evo2  # imported here so the package works without evo2 installed

        self._model = Evo2(self.name)
        self.set_precision(self.precision)

    def set_precision(self, precision: str) -> None:
        """Switch precision on the loaded model, without reloading weights."""
        if precision not in PRECISIONS:
            raise ValueError(f"precision must be one of {PRECISIONS}, not {precision!r}")
        if self._model is None:
            raise RuntimeError("call load() before set_precision()")
        from transformer_engine.common.recipe import (
            DelayedScaling,
            Float8CurrentScaling,
            Format,
        )

        net = self._model.model
        use_fp8 = precision != "bf16"
        # The model pads sequences to a multiple of 16 only when this flag is set.
        net.config["use_fp8_input_projections"] = use_fp8
        n = 0
        for module in net.modules():
            if hasattr(module, "use_fp8_input_projections") and hasattr(module, "fp8_recipe"):
                module.use_fp8_input_projections = use_fp8
                if precision == "fp8-current":
                    module.fp8_recipe = Float8CurrentScaling(fp8_format=Format.HYBRID)
                else:  # Evo 2's own recipe (vortex.model.layers.set_format_recipe)
                    module.fp8_recipe = DelayedScaling(
                        fp8_format=Format.HYBRID, amax_history_len=16, amax_compute_algo="max"
                    )
                n += 1
        if n == 0:
            raise RuntimeError("found no FP8 projection layers; has the evo2 package changed?")
        self.precision = precision

    def score_sequences(self, seqs: Sequence[str]) -> np.ndarray:
        if self._model is None:
            raise RuntimeError("call load() before scoring")
        scores = self._model.score_sequences(list(seqs), batch_size=1)
        return np.asarray(scores, dtype=float)
