from collections import Counter

import numpy as np
import pytest

from seqcontrol.controls.base import Interval, Region
from seqcontrol.controls.flank_shuffle import FlankShuffle, dinucleotide_shuffle


def dinucs(s):
    return Counter(s[i : i + 2] for i in range(len(s) - 1))


def random_dna(n, seed=0, alphabet="ACGT"):
    return "".join(np.random.default_rng(seed).choice(list(alphabet), n))


@pytest.mark.parametrize("seed", range(20))
def test_dinucleotide_shuffle_preserves_counts_and_ends(seed):
    rng = np.random.default_rng(seed)
    # Mix of random, skewed and repetitive inputs.
    s = random_dna(200, seed) + "CGCGCG" + random_dna(50, seed + 1, "AAT")
    t = dinucleotide_shuffle(s, rng)
    assert len(t) == len(s)
    assert dinucs(t) == dinucs(s)
    assert Counter(t) == Counter(s)
    assert (t[0], t[-1]) == (s[0], s[-1])


def test_dinucleotide_shuffle_actually_shuffles():
    s = random_dna(300, seed=5)
    outs = {dinucleotide_shuffle(s, np.random.default_rng(k)) for k in range(10)}
    assert s not in outs
    assert len(outs) == 10


def test_dinucleotide_shuffle_handles_trivial_inputs():
    rng = np.random.default_rng(0)
    for s in ["", "A", "AC", "AAAA", "ACACACAC"]:
        t = dinucleotide_shuffle(s, rng)
        assert dinucs(t) == dinucs(s) and t[:1] == s[:1]


def test_dinucleotide_shuffle_is_reproducible_for_a_seed():
    s = random_dna(150, seed=2)
    a = dinucleotide_shuffle(s, np.random.default_rng(42))
    b = dinucleotide_shuffle(s, np.random.default_rng(42))
    assert a == b


@pytest.fixture
def window():
    return Region(random_dna(1025, seed=9), genes=(Interval("tRNA", 480, 550),))


def test_flank_shuffle_keeps_gene_and_radius_untouched(window):
    r = 100
    new = FlankShuffle.shuffle(window, r, np.random.default_rng(0))
    s, t = window.sequence, new.sequence
    assert t[380:650] == s[380:650]  # gene plus r on each side
    assert t[:380] != s[:380] and t[650:] != s[650:]
    assert dinucs(t[:380]) == dinucs(s[:380])
    assert dinucs(t[650:]) == dinucs(s[650:])


def test_flank_shuffle_radius_beyond_window_is_identity(window):
    new = FlankShuffle.shuffle(window, 600, np.random.default_rng(0))
    assert new.sequence == window.sequence


def test_flank_shuffle_apply_makes_one_copy_per_radius_and_seed(window):
    out = FlankShuffle([0, 50, 200], n_seeds=4).apply(window, np.random.default_rng(1))
    assert [(r.meta["radius"], r.meta["seed"]) for r in out] == [
        (r, k) for r in (0, 50, 200) for k in range(4)
    ]
    assert len({r.sequence for r in out}) == 12
