"""Minimal VCF reader: the first eight columns, with INFO parsed into a dict."""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class VcfRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    info: dict[str, str] = field(default_factory=dict)


def read_vcf(path: Path, chrom: str | None = None, encoding: str = "utf-8") -> Iterator[VcfRecord]:
    """Yield records, optionally only those on `chrom`. Reads .vcf or .vcf.gz."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding=encoding) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if chrom is not None and cols[0] != chrom:
                continue
            info = {}
            for item in cols[7].split(";"):
                key, _, value = item.partition("=")
                info[key] = value
            yield VcfRecord(cols[0], int(cols[1]), cols[3], cols[4], info)
