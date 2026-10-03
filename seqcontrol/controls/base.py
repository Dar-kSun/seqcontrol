"""The interface every control satisfies.

A control takes an unmutated Region and returns perturbed copies of it. Each copy
carries `origin`, which maps every index of the original sequence to its index in
the copy (-1 if the base was removed). Variants are introduced afterwards, at their
mapped position, so ref and alt are always perturbed identically.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol

import numpy as np


@dataclass(frozen=True)
class Interval:
    """A named span of a sequence, 0-based and half-open."""

    name: str
    start: int
    end: int

    def __len__(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class Region:
    sequence: str
    genes: tuple[Interval, ...] = ()
    circular: bool = False
    origin: np.ndarray | None = field(default=None, compare=False)
    meta: dict = field(default_factory=dict, compare=False)

    def map_index(self, i: int) -> int:
        """Where original index `i` sits in this sequence (identity if unperturbed)."""
        if self.origin is None:
            return i
        j = int(self.origin[i])
        if j < 0:
            raise ValueError(f"original index {i} was removed by {self.meta}")
        return j


class Control(Protocol):
    name: str

    def apply(self, region: Region, rng: np.random.Generator) -> list[Region]: ...


def check_genes_preserved(original: Region, perturbed: Region) -> None:
    """Raise unless every gene's bases are byte-identical and contiguous after the control.

    A control that corrupts the gene it claims to hold fixed is invalid.
    """
    for g in original.genes:
        new = [perturbed.map_index(i) for i in range(g.start, g.end)]
        if new != list(range(new[0], new[0] + len(g))):
            raise AssertionError(f"{g.name} is no longer contiguous after {perturbed.meta}")
        if perturbed.sequence[new[0] : new[0] + len(g)] != original.sequence[g.start : g.end]:
            raise AssertionError(f"{g.name} bases changed after {perturbed.meta}")
