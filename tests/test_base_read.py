"""Tests for Read model behavior."""

from scanitd.base.basic import MappingMode
from scanitd.base.basic_read import Read, reverse_complement


def test_read_new_populates_derived_fields() -> None:
    read = Read.new(
        query_name="q1",
        chrom="chr1",
        ref_start=100,
        strand="+",
        cigar_str="5S10M2I5M10N10M15S",
        mapq=60,
        nm=1,
        query_seq="A" * 47,
        query_qualities=None,
    )

    assert read.ref_end == 135
    assert read.sms == (5, 27, 15)
    assert read.simple_mode == MappingMode.MS


def test_read_simple_mode_sm() -> None:
    read = Read.new(
        query_name="q2",
        chrom="chr2",
        ref_start=42,
        strand="-",
        cigar_str="12S20M4S",
        mapq=50,
        nm=0,
        query_seq="C" * 36,
        query_qualities=None,
    )

    assert read.simple_mode == MappingMode.SM


def test_reverse_complement() -> None:
    assert reverse_complement("ACTG") == "CAGT"
