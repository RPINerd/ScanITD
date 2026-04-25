"""Tests for CIGAR parsing behavior."""

from scanitd.base.cigar import parse_cigar


def test_parse_cigar_with_softclip_and_indels() -> None:
    result = parse_cigar("5S10M2I5M10N10M15S")

    assert result.lt_soft_len == 5
    assert result.rt_soft_len == 15
    assert result.read_match == 27
    assert result.ref_match == 35
    assert result.indel_len == 8
    assert result.query_len == 47
    assert result.cigartuples_without_soft == [(0, 10), (1, 2), (0, 5), (3, 10), (0, 10)]


def test_parse_cigar_without_softclip() -> None:
    result = parse_cigar("20M3D7M")

    assert result.lt_soft_len == 0
    assert result.rt_soft_len == 0
    assert result.read_match == 27
    assert result.ref_match == 30
    assert result.indel_len == 3
    assert result.query_len == 27
