"""Download every data source into data/raw/ and rewrite data/MANIFEST.md.

Usage:
    python scripts/fetch_data.py

Files already in data/raw/ are reused, not re-downloaded. To refresh a source,
delete its file and re-run. The manifest records each file's URL, download date,
size and SHA-256, plus how many labelled variants the files yield.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

from seqcontrol import config
from seqcontrol.data import mtdna
from seqcontrol.data.cache import sha256

MANIFEST = config.ROOT / "data" / "MANIFEST.md"

SOURCES = [
    ("chrM_GRCh38.fa", "Ensembl REST: GRCh38 chrM sequence (rCRS)", config.ENSEMBL_REST),
    ("chrM_GRCh38_genes.json", "Ensembl REST: chrM gene annotation", config.ENSEMBL_REST),
    ("mitomap_disease.vcf", "MITOMAP disease table (unversioned)", config.MITOMAP_DISEASE_VCF),
    (
        f"clinvar_{config.CLINVAR_RELEASE}.vcf.gz",
        f"ClinVar GRCh38 VCF, release {config.CLINVAR_RELEASE}",
        config.CLINVAR_VCF_URLS[0],
    ),
]


def download_date(path: Path) -> str:
    """The VCF's own ##fileDate if it has one (MITOMAP stamps the serve date), else mtime."""
    if path.suffix == ".vcf":
        with open(path, encoding="latin-1") as f:
            for line in f:
                if not line.startswith("##"):
                    break
                if line.startswith("##fileDate="):
                    d = line.strip().split("=", 1)[1]
                    return f"{d[:4]}-{d[4:6]}-{d[6:8]}"
    return dt.date.fromtimestamp(path.stat().st_mtime).isoformat()


def main() -> None:
    data = mtdna.load()
    r = data.report
    counts = data.biotype_counts()

    lines = [
        "# Data manifest",
        "",
        "Written by `python scripts/fetch_data.py`. Raw files live in `data/raw/`,",
        "which git ignores. Re-running reuses files already there; a different",
        "checksum below means the source changed and results may not match.",
        "",
        f"Genome build: {config.GENOME_BUILD}.",
        "",
        "| File | Source | Downloaded | Size (bytes) | SHA-256 |",
        "|---|---|---|---|---|",
    ]
    for name, desc, url in SOURCES:
        p = config.DATA_RAW / name
        row = [f"`{name}`", f"{desc} ({url})", download_date(p), p.stat().st_size, f"`{sha256(p)}`"]
        lines.append("| " + " | ".join(map(str, row)) + " |")
    lines += [
        "",
        "## Labelled mtDNA variants",
        "",
        "Single-base variants only. Pathogenic: MITOMAP `Cfrm-[P]` or `Cfrm-[LP]`.",
        "Benign: ClinVar Benign / Likely benign with 2+ review stars. Rules are in",
        "`seqcontrol/data/mitomap.py` and `seqcontrol/data/clinvar.py`.",
        "",
        "| Step | Count |",
        "|---|---|",
        f"| Pathogenic from MITOMAP | {r['pathogenic_in']} |",
        f"| Benign from ClinVar | {r['benign_in']} |",
        f"| Duplicates dropped | {r['duplicates_dropped']} |",
        f"| Labelled both ways, dropped | {r['label_conflicts_dropped']} |",
        f"| Reference base mismatch, dropped | {r['ref_mismatch_dropped']} |",
        f"| **Pathogenic in final set** | **{r['pathogenic_out']}** |",
        f"| **Benign in final set** | **{r['benign_out']}** |",
        "",
        "By region:",
        "",
        "| Region | Pathogenic | Benign |",
        "|---|---|---|",
    ]
    for biotype in sorted({b for b, _ in counts}):
        lines.append(f"| {biotype} | {counts[(biotype, 1)]} | {counts[(biotype, 0)]} |")

    MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
