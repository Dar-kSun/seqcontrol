"""Flank-shuffle control: hold the gene, scramble its surroundings beyond radius r.

For each radius r, every base more than r bp outside the gene (on either side) is
replaced by a dinucleotide-preserving shuffle of that flank segment. Small r destroys
almost all context; large r leaves the gene's near neighbourhood intact. Shuffling
preserves each segment's exact dinucleotide counts and its first and last base, so
base composition, CpG content and the junction with the untouched part all stay the
same: only the order of the context changes. Naive (mononucleotide) shuffling would
also change composition and confound the result.

Dinucleotide shuffle: Altschul & Erikson (1985), "Significance of nucleotide sequence
alignments", Mol Biol Evol 2:526, via the random Eulerian path construction of
Kandel et al. (1996).
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence

import numpy as np

from seqcontrol.controls.base import Region, check_genes_preserved


def dinucleotide_shuffle(seq: str, rng: np.random.Generator) -> str:
    """A uniformly random sequence with exactly the same dinucleotide counts as `seq`.

    It also keeps the first and last base. Sequences shorter than 3 are returned as is.
    """
    if len(seq) < 3:
        return seq
    # Each adjacent pair is an edge u -> v; any Eulerian path from seq[0] using every
    # edge once spells a sequence with the same dinucleotide counts.
    edges: dict[str, list[str]] = {}
    for u, v in zip(seq, seq[1:], strict=False):
        edges.setdefault(u, []).append(v)
    end = seq[-1]

    # Choose each vertex's last exit at random so the last exits form a tree
    # pointing at `end` (no cycles); this guarantees the walk uses every edge.
    while True:
        last = {u: outs[rng.integers(len(outs))] for u, outs in edges.items() if u != end}
        if all(_reaches(u, end, last) for u in last):
            break

    order: dict[str, list[str]] = {}
    for u, outs in edges.items():
        rest = list(outs)
        if u in last:
            rest.remove(last[u])
        rng.shuffle(rest)
        order[u] = rest + ([last[u]] if u in last else [])

    out, cur = [seq[0]], seq[0]
    pointers = Counter()
    for _ in range(len(seq) - 1):
        nxt = order[cur][pointers[cur]]
        pointers[cur] += 1
        out.append(nxt)
        cur = nxt
    return "".join(out)


def _reaches(u: str, end: str, last: dict[str, str]) -> bool:
    seen = set()
    while u != end:
        if u in seen:
            return False
        seen.add(u)
        u = last[u]
    return True


class FlankShuffle:
    name = "flank-shuffle"

    def __init__(self, radii: Sequence[int], n_seeds: int = 10) -> None:
        self.radii = list(radii)
        self.n_seeds = n_seeds

    def apply(self, region: Region, rng: np.random.Generator) -> list[Region]:
        """One perturbed copy per (radius, seed); each copy's meta records both."""
        out = []
        for r in self.radii:
            for seed in range(self.n_seeds):
                # Seed per (radius, seed) from the caller's rng, so runs are reproducible.
                child = np.random.default_rng(rng.integers(2**63))
                out.append(self.shuffle(region, r, child, meta={"radius": r, "seed": seed}))
        return out

    @staticmethod
    def shuffle(
        region: Region, r: int, rng: np.random.Generator, meta: dict | None = None
    ) -> Region:
        """Shuffle everything more than r bp outside the region's single gene."""
        if len(region.genes) != 1:
            raise ValueError("flank shuffle needs exactly one gene to hold fixed")
        (gene,) = region.genes
        seq = region.sequence
        left_end = max(gene.start - r, 0)  # left flank to shuffle: [0, left_end)
        right_start = min(gene.end + r, len(seq))  # right flank: [right_start, len)
        new_seq = (
            dinucleotide_shuffle(seq[:left_end], rng)
            + seq[left_end:right_start]
            + dinucleotide_shuffle(seq[right_start:], rng)
        )
        new = Region(
            new_seq,
            genes=region.genes,
            circular=region.circular,
            origin=np.arange(len(seq)),
            meta={"control": "flank-shuffle", **(meta or {"radius": r})},
        )
        check_genes_preserved(region, new)
        return new
