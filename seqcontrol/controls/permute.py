"""Context-swap controls: keep a gene's own sequence, change what surrounds it.

Two versions, because "cyclic permutation" has been used for two different things:

TrnaSwap (Mathur & Sachidanandam 2026): on the whole chromosome, every tRNA gene
moves into the slot of the tRNA k places further along (cyclically), carrying its
own sequence with it. Everything else stays put. A variant inside a tRNA is then
scored in its usual-size window at the gene's new location, so the gene is the same
but its neighbourhood is another tRNA's. tRNAs that overlap are moved as one unit,
so no base is duplicated or lost. Genes keep their reference-strand orientation.

WindowRotation (CLAUDE.md 5.1): within one scoring window, rotate the sequence by r
bp (bases shifted off one end re-enter at the other). The gene stays contiguous,
but its upstream and downstream context are rearranged and an artificial junction
joins the window's two ends. Offsets that would cut through a gene are skipped.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from seqcontrol.controls.base import Interval, Region, check_genes_preserved


def merge_overlapping(genes: Sequence[Interval]) -> list[Interval]:
    """Join overlapping intervals (sorted by start) into single units named A+B."""
    units: list[Interval] = []
    for g in sorted(genes, key=lambda g: g.start):
        if units and g.start < units[-1].end:
            last = units.pop()
            g = Interval(f"{last.name}+{g.name}", last.start, max(last.end, g.end))
        units.append(g)
    return units


class TrnaSwap:
    name = "trna-swap"

    def __init__(self, shifts: Sequence[int] | None = None) -> None:
        self.shifts = shifts  # None: every shift 1 .. n_units - 1

    def apply(self, region: Region, rng: np.random.Generator | None = None) -> list[Region]:
        units = merge_overlapping(region.genes)
        shifts = self.shifts if self.shifts is not None else range(1, len(units))
        return [self.shift(region, units, k) for k in shifts]

    @staticmethod
    def shift(region: Region, units: list[Interval], k: int) -> Region:
        """Unit i's sequence goes to slot (i + k) mod n. k = 0 returns the original."""
        n, seq = len(units), region.sequence
        # Build the new sequence slot by slot; the gaps between slots never move.
        pieces, origin = [], np.full(len(seq), -1, dtype=np.int64)
        cursor, out_len = 0, 0
        for slot in range(n):
            gap = range(cursor, units[slot].start)
            origin[gap.start : gap.stop] = np.arange(out_len, out_len + len(gap))
            pieces.append(seq[gap.start : gap.stop])
            out_len += len(gap)
            src = units[(slot - k) % n]  # the unit that lands in this slot
            origin[src.start : src.end] = np.arange(out_len, out_len + len(src))
            pieces.append(seq[src.start : src.end])
            out_len += len(src)
            cursor = units[slot].end
        origin[cursor:] = np.arange(out_len, out_len + len(seq) - cursor)
        pieces.append(seq[cursor:])
        new = Region(
            "".join(pieces),
            genes=region.genes,  # names only; positions are found through origin
            circular=region.circular,
            origin=origin,
            meta={"control": "trna-swap", "shift": k},
        )
        check_genes_preserved(region, new)
        return new


class WindowRotation:
    name = "window-rotation"

    def __init__(self, offsets: Sequence[int] | None = None, step: int = 64) -> None:
        self.offsets = offsets
        self.step = step

    def valid_offsets(self, region: Region) -> list[int]:
        """Rotations that leave every gene in one piece: the cut must not fall inside one."""
        n = len(region.sequence)
        candidates = self.offsets if self.offsets is not None else range(self.step, n, self.step)
        return [r for r in candidates if not any(g.start < r < g.end for g in region.genes)]

    def apply(self, region: Region, rng: np.random.Generator | None = None) -> list[Region]:
        return [self.rotate(region, r) for r in self.valid_offsets(region)]

    @staticmethod
    def rotate(region: Region, r: int) -> Region:
        """New sequence = seq[r:] + seq[:r]; original index i moves to (i - r) mod n."""
        n = len(region.sequence)
        if any(g.start < r % n < g.end for g in region.genes):
            raise ValueError(f"rotation by {r} would cut a gene")
        new = Region(
            region.sequence[r:] + region.sequence[:r],
            genes=region.genes,
            circular=False,
            origin=(np.arange(n) - r) % n,
            meta={"control": "window-rotation", "offset": r},
        )
        check_genes_preserved(region, new)
        return new
