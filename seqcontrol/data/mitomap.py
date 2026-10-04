"""Pathogenic mtDNA variants from MITOMAP's disease table, and allele frequencies.

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


def allele_frequencies(path: Path) -> dict[tuple[int, str, str], float]:
    """(pos, ref, alt) -> allele frequency as a fraction, from the polymorphism table.

    MITOMAP's AF field is a percentage of ~66,800 full-length GenBank sequences
    (e.g. AC=1169 gives AF=1.7499), so it is divided by 100 here.
    """
    out = {}
    for rec in read_vcf(path, encoding="latin-1"):
        for alt, af in zip(rec.alt.split(","), rec.info.get("AF", "").split(","), strict=False):
            if len(rec.ref) == 1 and len(alt) == 1 and af:
                out[(rec.pos, rec.ref, alt)] = float(af) / 100
    return out


def _fetch_mitomap(url: str, dest: Path) -> Path:
    try:
        return fetch(url, dest)
    except RuntimeError as e:
        # MITOMAP sits behind Cloudflare, which sometimes challenges scripted downloads.
        raise RuntimeError(
            f"{e}\n\nMITOMAP blocked the automatic download. Open {url} in a browser, "
            f"save the file as {dest}, and re-run."
        ) from e


def fetch_disease_vcf(raw: Path = config.DATA_RAW) -> Path:
    return _fetch_mitomap(config.MITOMAP_DISEASE_VCF, raw / "mitomap_disease.vcf")


def fetch_polymorphisms_vcf(raw: Path = config.DATA_RAW) -> Path:
    return _fetch_mitomap(config.MITOMAP_POLYMORPHISMS_VCF, raw / "mitomap_polymorphisms.vcf")
