"""Benign mtDNA variants from ClinVar.

Label rule: a single-base chrM variant is benign if its ClinVar germline
classification is Benign, Likely_benign or Benign/Likely_benign, with a review status
of at least two stars: criteria provided by multiple submitters with no conflicts,
reviewed by expert panel, or practice guideline.
"""

from __future__ import annotations

from pathlib import Path

from seqcontrol import config
from seqcontrol.data.cache import fetch
from seqcontrol.data.vcf import read_vcf
from seqcontrol.variants import Variant

BENIGN = {"Benign", "Likely_benign", "Benign/Likely_benign"}
MIN_TWO_STARS = {
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "reviewed_by_expert_panel",
    "practice_guideline",
}


def benign_mt_variants(path: Path) -> list[Variant]:
    out = []
    for rec in read_vcf(path, chrom="MT"):
        if rec.info.get("CLNSIG") not in BENIGN:
            continue
        if rec.info.get("CLNREVSTAT") not in MIN_TWO_STARS:
            continue
        if len(rec.ref) == 1 and len(rec.alt) == 1 and rec.alt in "ACGT":
            out.append(Variant("MT", rec.pos, rec.ref, rec.alt, label=0))
    return out


def fetch_clinvar_vcf(raw: Path = config.DATA_RAW) -> Path:
    return fetch(config.CLINVAR_VCF_URLS, raw / f"clinvar_{config.CLINVAR_RELEASE}.vcf.gz")
