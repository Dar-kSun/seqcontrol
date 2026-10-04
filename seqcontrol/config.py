"""Pinned data sources. Change a value here and re-run scripts/fetch_data.py."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = ROOT / "data" / "raw"

GENOME_BUILD = "GRCh38"  # chrM in GRCh38 is the revised Cambridge Reference Sequence (rCRS)
CHRM_LENGTH = 16_569
ENSEMBL_REST = "https://rest.ensembl.org"

# MITOMAP serves its current tables with no version number, so each download is
# pinned by date and checksum in data/MANIFEST.md instead.
MITOMAP_DISEASE_VCF = "https://www.mitomap.org/cgi-bin/disease.cgi?format=vcf"
MITOMAP_POLYMORPHISMS_VCF = "https://www.mitomap.org/cgi-bin/polymorphisms.cgi?format=vcf"

CLINVAR_RELEASE = "20260928"
CLINVAR_VCF_URLS = [  # a dated release moves to the archive folder after a few weeks
    f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_{CLINVAR_RELEASE}.vcf.gz",
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/"
    f"{CLINVAR_RELEASE[:4]}/clinvar_{CLINVAR_RELEASE}.vcf.gz",
]
