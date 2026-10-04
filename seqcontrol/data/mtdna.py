"""The labelled mtDNA variant set: MITOMAP pathogenic plus ClinVar benign, on GRCh38 chrM."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from seqcontrol import config
from seqcontrol.data import clinvar, mitomap, reference
from seqcontrol.data.reference import Gene
from seqcontrol.variants import Variant


@dataclass
class MtDataset:
    sequence: str
    genes: list[Gene]
    variants: list[Variant]
    report: dict[str, int] = field(default_factory=dict)

    def gene_of(self, v: Variant) -> Gene | None:
        return next((g for g in self.genes if g.contains(v.pos)), None)

    def biotype_counts(self) -> Counter:
        """(biotype, label) -> count, with 'non-coding' for variants outside any gene."""
        counts = Counter()
        for v in self.variants:
            g = self.gene_of(v)
            counts[(g.biotype if g else "non-coding", v.label)] += 1
        return counts


def combine(
    sequence: str, pathogenic: list[Variant], benign: list[Variant]
) -> tuple[list[Variant], dict[str, int]]:
    """Merge the two label sources, dropping duplicates, conflicts and reference mismatches."""
    key = lambda v: (v.pos, v.ref, v.alt)  # noqa: E731
    path_keys = {key(v) for v in pathogenic}
    ben_keys = {key(v) for v in benign}
    conflicting = path_keys & ben_keys

    report = {
        "pathogenic_in": len(pathogenic),
        "benign_in": len(benign),
        "duplicates_dropped": len(pathogenic) + len(benign) - len(path_keys) - len(ben_keys),
        "label_conflicts_dropped": len(conflicting),
        "ref_mismatch_dropped": 0,
    }
    seen, out = set(), []
    for v in sorted(pathogenic + benign, key=lambda v: (v.pos, v.alt)):
        k = key(v)
        if k in conflicting or k in seen:
            continue
        seen.add(k)
        if sequence[v.pos - 1] != v.ref:
            report["ref_mismatch_dropped"] += 1
            continue
        out.append(v)
    report["pathogenic_out"] = sum(v.label == 1 for v in out)
    report["benign_out"] = sum(v.label == 0 for v in out)
    return out, report


LABEL_SETS = {"strict": 2, "expanded": 1}  # minimum ClinVar review stars for benign


def load(raw: Path = config.DATA_RAW, labels: str = "strict") -> MtDataset:
    """Fetch (once) and assemble the labelled mtDNA set.

    labels="strict" (the headline set) takes ClinVar benign at 2+ stars; "expanded"
    also takes 1 star, as a sensitivity analysis (docs/plan-v0.3-within-gene.md).
    """
    if labels not in LABEL_SETS:
        raise ValueError(f"labels must be one of {list(LABEL_SETS)}")
    sequence = reference.chrm_sequence(raw)
    genes = reference.chrm_genes(raw)
    pathogenic = mitomap.pathogenic_variants(mitomap.fetch_disease_vcf(raw))
    benign = clinvar.benign_mt_variants(clinvar.fetch_clinvar_vcf(raw), LABEL_SETS[labels])
    variants, report = combine(sequence, pathogenic, benign)
    return MtDataset(sequence, genes, variants, report)
