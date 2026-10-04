"""Loader tests on small excerpts of the real MITOMAP and ClinVar files (no network)."""

from pathlib import Path

import pytest

from seqcontrol.data import clinvar, mitomap
from seqcontrol.data.mtdna import combine
from seqcontrol.data.reference import Gene, parse_genes, read_fasta
from seqcontrol.data.vcf import read_vcf
from seqcontrol.variants import Variant

FIXTURES = Path(__file__).parent / "fixtures"
MITOMAP = FIXTURES / "mitomap_disease_sample.vcf"
CLINVAR = FIXTURES / "clinvar_mt_sample.vcf"


def keys(variants):
    return {(v.pos, v.ref, v.alt, v.label) for v in variants}


def test_vcf_reader_parses_info_with_commas_in_values():
    rec = next(r for r in read_vcf(CLINVAR) if r.pos == 1018)
    assert rec.info["CLNREVSTAT"] == "criteria_provided,_multiple_submitters,_no_conflicts"


def test_mitomap_keeps_only_confirmed_single_base_alleles():
    # 72 is only Reported; 3272 is a deletion with Cfrm-[VUS*]: both excluded.
    assert keys(mitomap.pathogenic_variants(MITOMAP)) == {
        (616, "T", "C", 1),  # Cfrm-[LP]; its other allele, T>G, is only Reported-[VUS]
        (1555, "A", "G", 1),
        (3243, "A", "G", 1),  # Cfrm-[P]
        (3243, "A", "T", 1),  # Cfrm-[LP]
        (8344, "A", "G", 1),
    }


def test_mitomap_pairs_each_status_with_its_own_allele():
    alleles_616 = {v.alt for v in mitomap.pathogenic_variants(MITOMAP) if v.pos == 616}
    assert alleles_616 == {"C"}
    assert mitomap.skipped_rows(MITOMAP) == []


def test_clinvar_keeps_two_star_benign_single_base_variants():
    # Excluded: 15 (pathogenic), 235 (one star), 247 (indel).
    assert keys(clinvar.benign_mt_variants(CLINVAR)) == {
        (1018, "G", "A", 0),
        (1420, "T", "C", 0),
        (5293, "G", "A", 0),  # Likely benign, expert panel
    }


def test_combine_drops_conflicts_duplicates_and_reference_mismatches():
    seq = "ACGTACGTAC"
    pathogenic = [
        Variant("MT", 1, "A", "G", 1),
        Variant("MT", 1, "A", "G", 1),  # duplicate
        Variant("MT", 2, "C", "T", 1),  # also benign below: conflict
        Variant("MT", 3, "T", "A", 1),  # sequence has G here: mismatch
    ]
    benign = [Variant("MT", 2, "C", "T", 0), Variant("MT", 4, "T", "C", 0)]
    out, report = combine(seq, pathogenic, benign)
    assert keys(out) == {(1, "A", "G", 1), (4, "T", "C", 0)}
    assert report["duplicates_dropped"] == 1
    assert report["label_conflicts_dropped"] == 1
    assert report["ref_mismatch_dropped"] == 1
    assert (report["pathogenic_out"], report["benign_out"]) == (1, 1)


def test_parse_genes_sorts_and_names():
    genes = parse_genes(
        [
            {"id": "ENSG2", "external_name": "MT-TV", "biotype": "Mt_tRNA",
             "start": 1602, "end": 1670, "strand": 1},
            {"id": "ENSG1", "external_name": "MT-TF", "biotype": "Mt_tRNA",
             "start": 577, "end": 647, "strand": 1},
        ]
    )  # fmt: skip
    assert [g.name for g in genes] == ["MT-TF", "MT-TV"]
    assert Gene("MT-TF", "Mt_tRNA", 577, 647, 1).contains(647)


@pytest.mark.parametrize("pos,ref", [(1, "G"), (16, "A"), (1000, "T")])
def test_reference_fixture_matches_rcrs(pos, ref):
    assert read_fasta(FIXTURES / "chrM_1_1000.fa")[pos - 1] == ref


def test_clinvar_one_star_adds_single_submitter_benign():
    strict = keys(clinvar.benign_mt_variants(CLINVAR))
    widened = keys(clinvar.benign_mt_variants(CLINVAR, min_stars=1))
    assert strict < widened
    assert widened - strict == {(235, "A", "G", 0)}  # single-submitter benign
    with pytest.raises(ValueError):
        clinvar.benign_mt_variants(CLINVAR, min_stars=0)
