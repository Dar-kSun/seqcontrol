import numpy as np
import pytest

from seqcontrol.models.base import ModelAdapter, score_variants
from seqcontrol.models.evo2 import Evo2Adapter


class GCModel:
    """Toy adapter for tests: 'log-likelihood' is the GC fraction. No GPU needed."""

    name = "toy-gc"
    max_context = 100

    def load(self) -> None:
        pass

    def score_sequences(self, seqs):
        return np.array([sum(b in "GC" for b in s) / len(s) for s in seqs])


def test_toy_model_satisfies_the_protocol():
    model: ModelAdapter = GCModel()
    assert model.score_sequences(["GGCC"])[0] == 1.0


def test_score_is_alt_minus_ref():
    scores = score_variants(GCModel(), ["AAAA", "GGGG"], ["AAAG", "GGGA"])
    np.testing.assert_allclose(scores, [0.25, -0.25])


def test_window_longer_than_context_is_rejected():
    with pytest.raises(ValueError, match="exceeds"):
        score_variants(GCModel(), ["A" * 101], ["C" + "A" * 100])


def test_mismatched_lengths_are_rejected():
    with pytest.raises(ValueError):
        score_variants(GCModel(), ["AAAA"], [])
    with pytest.raises(ValueError):
        score_variants(GCModel(), ["AAAA"], ["AAA"])


def test_evo2_adapter_constructs_without_evo2_installed():
    adapter = Evo2Adapter("evo2_1b_base")
    assert adapter.max_context == 8_192
    with pytest.raises(RuntimeError, match="load"):
        adapter.score_sequences(["ACGT"])
    with pytest.raises(ValueError):
        Evo2Adapter("evo2_40b")
