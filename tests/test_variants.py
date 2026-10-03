import pytest

from seqcontrol.variants import Variant, variant_windows


def test_fixture_is_rcrs_start(chrm_1kb):
    assert len(chrm_1kb) == 1000
    assert chrm_1kb.startswith("GATCACAGGTCTATCACCC")


def test_windows_differ_only_at_the_centre(chrm_1kb):
    v = Variant("MT", 500, chrm_1kb[499], "A" if chrm_1kb[499] != "A" else "C")
    ref, alt = variant_windows(chrm_1kb, v, size=101)
    assert len(ref) == len(alt) == 101
    assert ref == chrm_1kb[449:550]
    assert [k for k in range(101) if ref[k] != alt[k]] == [50]
    assert alt[50] == v.alt


def test_circular_window_wraps_the_origin(chrm_1kb):
    v = Variant("MT", 1, "G", "A")
    ref, _ = variant_windows(chrm_1kb, v, size=11, circular=True)
    assert ref == chrm_1kb[-5:] + chrm_1kb[:6]
    assert ref[5] == "G"


def test_linear_window_shifts_inwards_at_the_end(chrm_1kb):
    v = Variant("MT", 1, "G", "A")
    ref, alt = variant_windows(chrm_1kb, v, size=11, circular=False)
    assert ref == chrm_1kb[:11]
    assert alt[0] == "A" and alt[1:] == ref[1:]


def test_reference_mismatch_is_rejected(chrm_1kb):
    wrong = "A" if chrm_1kb[0] != "A" else "C"
    with pytest.raises(ValueError, match="reference mismatch"):
        variant_windows(chrm_1kb, Variant("MT", 1, wrong, "T"), size=11)


@pytest.mark.parametrize(
    "kwargs",
    [dict(ref="N", alt="A"), dict(ref="A", alt="A"), dict(ref="AC", alt="A"), dict(pos=0)],
)
def test_invalid_variants_are_rejected(kwargs):
    args = dict(chrom="MT", pos=10, ref="A", alt="G") | kwargs
    with pytest.raises(ValueError):
        Variant(**args)
