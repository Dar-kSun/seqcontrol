"""The context-swap controls must move context without touching the genes."""

import numpy as np
import pytest

from seqcontrol.controls.base import Interval, Region, check_genes_preserved
from seqcontrol.controls.permute import TrnaSwap, WindowRotation, merge_overlapping


def random_dna(n, seed=0):
    return "".join(np.random.default_rng(seed).choice(list("ACGT"), n))


@pytest.fixture
def chromosome():
    # Five "tRNAs" of different lengths, two of which overlap, on a 600 bp circle.
    genes = (
        Interval("tA", 20, 80),
        Interval("tB", 150, 215),
        Interval("tC", 212, 270),  # overlaps tB by 3 bp
        Interval("tD", 340, 395),
        Interval("tE", 500, 571),
    )
    return Region(random_dna(600), genes=genes, circular=True)


def assert_origin_exact(original: Region, new: Region):
    """Every original base appears, unchanged, at its mapped index."""
    assert sorted(new.origin) == list(range(len(new.sequence)))
    for i, base in enumerate(original.sequence):
        assert new.sequence[new.origin[i]] == base


def test_overlapping_genes_move_as_one_unit():
    units = merge_overlapping([Interval("a", 0, 10), Interval("b", 8, 20), Interval("c", 25, 30)])
    assert units == [Interval("a+b", 0, 20), Interval("c", 25, 30)]


def test_trna_swap_shift_zero_is_identity(chromosome):
    units = merge_overlapping(chromosome.genes)
    assert TrnaSwap.shift(chromosome, units, 0).sequence == chromosome.sequence


def test_trna_swap_preserves_every_gene_and_every_base(chromosome):
    regions = TrnaSwap().apply(chromosome)
    assert [r.meta["shift"] for r in regions] == [1, 2, 3]  # 4 units after merging tB+tC
    for new in regions:
        assert len(new.sequence) == len(chromosome.sequence)
        assert new.sequence != chromosome.sequence
        assert_origin_exact(chromosome, new)
        check_genes_preserved(chromosome, new)


def test_trna_swap_moves_each_unit_into_the_next_slot(chromosome):
    units = merge_overlapping(chromosome.genes)
    new = TrnaSwap.shift(chromosome, units, 1)
    # Unit 0 (tA) now starts where it would after the gap before slot 1, i.e. it
    # occupies slot 1; the gap between slots is untouched.
    ta = chromosome.sequence[20:80]
    assert new.sequence.find(ta) == new.map_index(20)
    assert new.sequence[:20] == chromosome.sequence[:20]  # gap before slot 0 is fixed
    # Slot 0 now holds the last unit (tE).
    te = chromosome.sequence[500:571]
    assert new.sequence[20 : 20 + len(te)] == te


def test_window_rotation_preserves_genes_and_bases():
    window = Region(random_dna(1025, seed=1), genes=(Interval("tRNA", 480, 550),))
    rot = WindowRotation(step=64)
    offsets = rot.valid_offsets(window)
    assert 512 not in offsets  # inside the gene
    assert 448 in offsets and 576 in offsets
    for new in rot.apply(window):
        assert_origin_exact(window, new)
        check_genes_preserved(window, new)


def test_window_rotation_refuses_to_cut_a_gene():
    window = Region(random_dna(100), genes=(Interval("g", 40, 60),))
    with pytest.raises(ValueError, match="cut a gene"):
        WindowRotation.rotate(window, 50)


def test_check_genes_preserved_catches_a_corrupted_gene(chromosome):
    units = merge_overlapping(chromosome.genes)
    good = TrnaSwap.shift(chromosome, units, 1)
    i = good.map_index(30)  # inside tA
    flipped = "A" if good.sequence[i] != "A" else "C"
    bad_seq = good.sequence[:i] + flipped + good.sequence[i + 1 :]
    bad = Region(bad_seq, chromosome.genes, origin=good.origin, meta=good.meta)
    with pytest.raises(AssertionError, match="tA bases changed"):
        check_genes_preserved(chromosome, bad)
