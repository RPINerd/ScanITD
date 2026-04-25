"""Tests for core data structures in basic.py."""

import pytest

from scanitd.base.basic import Event, Interval, Intervals, MappingMode, MicroRegion, Strand


def test_microregion_parsing_and_equality() -> None:
    insertion = MicroRegion("+AC")
    homology = MicroRegion("-GG")
    blunt = MicroRegion("")

    assert insertion.micro_type == "microinsertion"
    assert insertion.sequence == "AC"
    assert insertion.length == 2
    assert homology.micro_type == "microhomology"
    assert blunt.micro_type == "blunt_end"
    assert blunt.sequence == ""

    assert insertion == MicroRegion("+AC")
    assert insertion != homology


def test_event_new_sets_end_and_af() -> None:
    event = Event.new(
        event_type="TDUP",
        event_id=("chr2", 100, 25, "ACTG", MicroRegion("+A")),
        oao=9,
        ao=4,
        dp=10,
        ref_allele="A",
        alt_allele="<TDUP>",
    )

    assert event.end == 125
    assert event.af == 0.4
    assert event.event_type == "TDUP"


def test_mapping_mode_helpers() -> None:
    assert MappingMode.from_int(1).is_ms()
    assert MappingMode.from_int(2).is_sm()
    assert MappingMode.MS.reversed() == MappingMode.SM

    with pytest.raises(ValueError):
        MappingMode.from_int(42)


def test_strand_helpers() -> None:
    assert Strand.from_str("+").is_forward()
    assert Strand.from_str("-").is_reverse()
    assert str(Strand.Forward) == "+"

    with pytest.raises(ValueError):
        Strand.from_str("x")


def test_interval_operations_and_containment() -> None:
    left = Interval(10, 20)
    right = Interval(15, 30)

    overlap, union = left.join(right)
    assert overlap == Interval(15, 20)
    assert union == Interval(10, 30)

    assert left + 3 == Interval(13, 23)
    assert right - 5 == Interval(10, 25)
    assert Interval(12, 18) in left
    assert left.contain(Interval(10, 20), same_left=True, same_right=True)


def test_intervals_add_sub_and_introns() -> None:
    exons = Intervals.from_list([(10, 20), (30, 40), (50, 60)])

    shifted = exons + 5
    assert shifted.first == Interval(15, 25)

    unshifted = shifted - 5
    assert unshifted == exons

    introns = exons.introns()
    assert introns is not None
    assert introns.exon_list == [Interval(20, 30), Interval(40, 50)]
