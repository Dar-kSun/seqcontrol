"""The GRCh38 mitochondrial reference sequence and its gene annotation, from Ensembl.

tRNA coordinates come from Ensembl's MT annotation, not GtRNAdb: GtRNAdb covers
nuclear tRNAs, and using the same source as the sequence guarantees they agree.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from seqcontrol import config
from seqcontrol.data.cache import fetch


@dataclass(frozen=True)
class Gene:
    """An annotated gene. `start` and `end` are 1-based and inclusive; strand is 1 or -1."""

    name: str
    biotype: str
    start: int
    end: int
    strand: int

    def contains(self, pos: int) -> bool:
        return self.start <= pos <= self.end


def read_fasta(path: Path) -> str:
    lines = path.read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">")).upper()


def chrm_sequence(raw: Path = config.DATA_RAW) -> str:
    url = f"{config.ENSEMBL_REST}/sequence/region/human/MT:1..{config.CHRM_LENGTH}:1"
    path = fetch(url + "?content-type=text/x-fasta", raw / "chrM_GRCh38.fa")
    seq = read_fasta(path)
    if len(seq) != config.CHRM_LENGTH or set(seq) - set("ACGTN"):
        raise ValueError(f"{path} is not a {config.CHRM_LENGTH} bp DNA sequence; delete it")
    return seq


def parse_genes(records: list[dict]) -> list[Gene]:
    genes = [
        Gene(r.get("external_name") or r["id"], r["biotype"], r["start"], r["end"], r["strand"])
        for r in records
    ]
    return sorted(genes, key=lambda g: g.start)


def chrm_genes(raw: Path = config.DATA_RAW) -> list[Gene]:
    url = (
        f"{config.ENSEMBL_REST}/overlap/region/human/MT:1-{config.CHRM_LENGTH}"
        "?feature=gene;content-type=application/json"
    )
    path = fetch(url, raw / "chrM_GRCh38_genes.json")
    return parse_genes(json.loads(path.read_text()))


def trna_genes(raw: Path = config.DATA_RAW) -> list[Gene]:
    return [g for g in chrm_genes(raw) if g.biotype == "Mt_tRNA"]
