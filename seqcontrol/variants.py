"""Variants and the reference/alternate sequence windows a model scores them on."""

from __future__ import annotations

from dataclasses import dataclass

BASES = frozenset("ACGT")


@dataclass(frozen=True)
class Variant:
    """A single-nucleotide variant. `pos` is 1-based, as in VCF, MITOMAP and ClinVar."""

    chrom: str
    pos: int
    ref: str
    alt: str
    label: int | None = None  # 1 = pathogenic, 0 = benign, None = unlabelled

    def __post_init__(self) -> None:
        if self.ref not in BASES or self.alt not in BASES:
            raise ValueError(f"only single-base A/C/G/T variants are supported: {self}")
        if self.ref == self.alt:
            raise ValueError(f"ref and alt are the same base: {self}")
        if self.pos < 1:
            raise ValueError(f"pos is 1-based and must be >= 1: {self}")


def variant_windows(
    sequence: str, variant: Variant, size: int, circular: bool = False
) -> tuple[str, str]:
    """Return (ref_window, alt_window) of length `size` centred on the variant.

    `sequence` is the whole chromosome (position 1 is sequence[0]). The variant sits at
    index size // 2 of each window, and the two windows differ only there. For a
    circular chromosome (mtDNA) the window wraps around the origin; for a linear one
    it is shifted inwards at the ends, so the variant is then off-centre.
    """
    n = len(sequence)
    if not 1 <= size <= n:
        raise ValueError(f"window size {size} must be between 1 and the sequence length {n}")
    if variant.pos > n:
        raise ValueError(f"{variant} lies beyond the sequence end ({n} bp)")
    i = variant.pos - 1
    if sequence[i].upper() != variant.ref:
        raise ValueError(
            f"reference mismatch at {variant.chrom}:{variant.pos}: "
            f"sequence has {sequence[i]}, variant says {variant.ref}"
        )

    start = i - size // 2
    if circular:
        idx = [(start + k) % n for k in range(size)]
    else:
        start = min(max(start, 0), n - size)
        idx = list(range(start, start + size))
    ref = "".join(sequence[j] for j in idx).upper()
    at = idx.index(i)
    alt = ref[:at] + variant.alt + ref[at + 1 :]
    return ref, alt
