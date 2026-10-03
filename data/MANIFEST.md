# Data manifest

Written by `python scripts/fetch_data.py`. Raw files live in `data/raw/`,
which git ignores. Re-running reuses files already there; a different
checksum below means the source changed and results may not match.

Genome build: GRCh38.

| File | Source | Downloaded | Size (bytes) | SHA-256 |
|---|---|---|---|---|
| `chrM_GRCh38.fa` | Ensembl REST: GRCh38 chrM sequence (rCRS) (https://rest.ensembl.org) | 2026-10-04 | 16878 | `7ee6821c2117e356d9870177997933208d4d5f2e9370443d9989c46964e4575d` |
| `chrM_GRCh38_genes.json` | Ensembl REST: chrM gene annotation (https://rest.ensembl.org) | 2026-10-04 | 15620 | `18d4b981733f95f5bfe9c0ad9361c1fe433459411b26ebd1d9793470b0757aaa` |
| `mitomap_disease.vcf` | MITOMAP disease table (unversioned) (https://www.mitomap.org/cgi-bin/disease.cgi?format=vcf) | 2026-10-03 | 225085 | `23fe09e3707a8438726c126a0d250965ff678e00ba725eeddd472d6ebc368ea2` |
| `clinvar_20260928.vcf.gz` | ClinVar GRCh38 VCF, release 20260928 (https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260928.vcf.gz) | 2026-10-04 | 198125411 | `e0ac8e4d2d8f08d3faf845c5b5baa46b2051a8254cb0be4509a586fb15b81cb0` |

## Labelled mtDNA variants

Single-base variants only. Pathogenic: MITOMAP `Cfrm-[P]` or `Cfrm-[LP]`.
Benign: ClinVar Benign / Likely benign with 2+ review stars. Rules are in
`seqcontrol/data/mitomap.py` and `seqcontrol/data/clinvar.py`.

| Step | Count |
|---|---|
| Pathogenic from MITOMAP | 94 |
| Benign from ClinVar | 228 |
| Duplicates dropped | 0 |
| Labelled both ways, dropped | 0 |
| Reference base mismatch, dropped | 0 |
| **Pathogenic in final set** | **94** |
| **Benign in final set** | **228** |

By region:

| Region | Pathogenic | Benign |
|---|---|---|
| Mt_rRNA | 2 | 12 |
| Mt_tRNA | 44 | 23 |
| protein_coding | 48 | 193 |
