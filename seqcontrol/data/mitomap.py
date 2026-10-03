"""Pathogenic mtDNA variants from MITOMAP's disease table.

Label rule: a single-base variant is pathogenic if MITOMAP marks that allele
"Cfrm-[P]" or "Cfrm-[LP]" (confirmed pathogenic or likely pathogenic). Every other
status (Reported, VUS, Cfrm-[VUS*], Conflicting, ...) is left out, not called benign.
On rows with several ALT alleles, DiseaseStatus has one entry per allele.
"""

from __future__ import annotations

from pathlib import Path

from seqcontrol import config
from seqcontrol.data.cache import fetch
from seqcontrol.data.vcf import read_vcf
from seqcontrol.variants import Variant

PATHOGENIC_STATUSES = {"Cfrm-[P]", "Cfrm-[LP]"}


def pathogenic_variants(path: Path) -> list[Variant]:
    out = []
    for rec in read_vcf(path, encoding="latin-1"):
        alts = rec.alt.split(",")
        statuses = rec.info.get("DiseaseStatus", "").split(",")
        if len(statuses) != len(alts):
            continue  # status cannot be matched to an allele; see skipped_rows()
        for alt, status in zip(alts, statuses, strict=True):
            if status in PATHOGENIC_STATUSES and len(rec.ref) == 1 and len(alt) == 1:
                out.append(Variant("MT", rec.pos, rec.ref, alt, label=1))
    return out


def skipped_rows(path: Path) -> list[str]:
    """Rows whose status count does not match their ALT count, so cannot be labelled."""
    return [
        f"MT:{rec.pos} {rec.ref}>{rec.alt} {rec.info.get('DiseaseStatus', '')}"
        for rec in read_vcf(path, encoding="latin-1")
        if len(rec.info.get("DiseaseStatus", "").split(",")) != len(rec.alt.split(","))
    ]


def fetch_disease_vcf(raw: Path = config.DATA_RAW) -> Path:
    dest = raw / "mitomap_disease.vcf"
    try:
        return fetch(config.MITOMAP_DISEASE_VCF, dest)
    except RuntimeError as e:
        # MITOMAP sits behind Cloudflare, which sometimes challenges scripted downloads.
        raise RuntimeError(
            f"{e}\n\nMITOMAP blocked the automatic download. Open "
            f"{config.MITOMAP_DISEASE_VCF} in a browser, save the file as {dest}, "
            "and re-run."
        ) from e
