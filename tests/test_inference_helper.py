"""Tests for isolated helper functions in the inference layer."""

from __future__ import annotations

from pathlib import Path

import pytest

from scanitd.base import MappingMode, MicroRegion, Read
from scanitd.base.basic import Event
from scanitd.inference.helper import (
    format_sa_tag,
    get_insertion_reference_pos,
    obtain_bp_region_seq,
    obtain_depth_given_genomic_position,
    obtain_sa_query_seq_from_ra,
    parse_target_genomic_coordinates,
    same_chrom_same_strand_handler,
    self_loop_checker,
    softclipped_length_and_event_size_checker,
    write_events_to_vcf,
)


class _FakeGenomeSlice:
    def __init__(self, seq: str) -> None:
        self.seq = seq


class _FakeGenomeChrom:
    def __init__(self, sequence: str) -> None:
        self._sequence = sequence

    def __getitem__(self, item: slice) -> _FakeGenomeSlice:
        return _FakeGenomeSlice(self._sequence[item])


class _FakeGenome:
    def __init__(self, sequence_by_chrom: dict[str, str]) -> None:
        self._sequence_by_chrom = sequence_by_chrom

    def __getitem__(self, chrom: str) -> _FakeGenomeChrom:
        return _FakeGenomeChrom(self._sequence_by_chrom[chrom])


class _FakeBam:
    def __init__(self) -> None:
        self.calls: list[tuple[str, int, int]] = []

    def count(self, contig: str, start: int, end: int) -> int:
        self.calls.append((contig, start, end))
        return 7


class _FakeLogger:
    def __init__(self) -> None:
        self.messages: list[str] = []

    def info(self, message: str) -> None:
        self.messages.append(message)

    def trace(self, message: str) -> None:
        self.messages.append(message)


def _build_read(
    *,
    ref_start: int,
    cigar_str: str,
    query_seq: str,
    chrom: str = "chr1",
    strand: str = "+",
) -> Read:
    return Read.new(
        query_name=f"{chrom}:{ref_start}",
        chrom=chrom,
        ref_start=ref_start,
        strand=strand,
        cigar_str=cigar_str,
        mapq=60,
        nm=0,
        query_seq=query_seq,
        query_qualities=None,
    )


def test_format_sa_tag_parses_types() -> None:
    assert format_sa_tag("chr7,42,+,10M1I5M,60,2") == ("chr7", 41, "+", "10M1I5M", 60, 2)


def test_obtain_sa_query_seq_from_ra_reverse_complements() -> None:
    assert obtain_sa_query_seq_from_ra("ACTG", "+", "+") == "ACTG"
    assert obtain_sa_query_seq_from_ra("ACTG", "+", "-") == "CAGT"


def test_parse_target_genomic_coordinates_supports_string_and_bed_file(tmp_path: Path) -> None:
    bed_path = tmp_path / "targets.bed"
    bed_path.write_text("chr1\t10\t20\nchr2\t30\t50\n", encoding="utf-8")

    assert parse_target_genomic_coordinates("chr3:5-15") == ["chr3:5-15"]
    assert parse_target_genomic_coordinates(str(bed_path)) == ["chr1:10-20", "chr2:30-50"]


def test_parse_target_genomic_coordinates_rejects_invalid_input() -> None:
    with pytest.raises(ValueError):
        parse_target_genomic_coordinates("chr1:20-10")

    with pytest.raises(ValueError):
        parse_target_genomic_coordinates("nonsense")


def test_self_loop_checker_detects_shifted_duplication() -> None:
    is_loop, offset, combo = self_loop_checker("ATGC", "TTGC", "AACC", allowed_mismatched=0)

    assert is_loop is True
    assert offset >= 0
    assert combo == "TGCA"


def test_self_loop_checker_returns_false_for_novel_sequence() -> None:
    assert self_loop_checker("AAAA", "TTTT", "CCCC", allowed_mismatched=0) == (False, 0, "")


def test_get_insertion_reference_pos_finds_match_and_missing_pattern() -> None:
    assert get_insertion_reference_pos("10M3I5M", 100, 3) == 109
    assert get_insertion_reference_pos("10M2D5M", 100, 3) is None


def test_obtain_depth_given_genomic_position_adjusts_ms_mode() -> None:
    bam = _FakeBam()

    assert obtain_depth_given_genomic_position(bam, "chr1", 50, MappingMode.SM) == 7
    assert bam.calls[-1] == ("chr1", 50, 51)

    assert obtain_depth_given_genomic_position(bam, "chr1", 50, MappingMode.MS) == 7
    assert bam.calls[-1] == ("chr1", 49, 50)


def test_softclipped_length_and_event_size_checker() -> None:
    read = _build_read(ref_start=100, cigar_str="12S20M4S", query_seq="A" * 36)

    assert softclipped_length_and_event_size_checker(read, MappingMode.SM, 10, 0) is False
    assert softclipped_length_and_event_size_checker(read, MappingMode.MS, 3, 0) is False


def test_obtain_bp_region_seq_for_insertion_and_homology() -> None:
    genome = _FakeGenome({"chr1": "A" * 40})
    read = _build_read(ref_start=4, cigar_str="6S10M4S", query_seq="TTTTTTAACCGGTTAACC")

    assert obtain_bp_region_seq(read, MappingMode.SM, 3, genome) == "+TTT"
    assert obtain_bp_region_seq(read, MappingMode.MS, -2, genome) == "-AA"


def test_same_chrom_same_strand_handler_returns_tdup_event() -> None:
    genome = _FakeGenome({"chr1": "A" * 200})
    logger = _FakeLogger()
    left = _build_read(ref_start=10, cigar_str="4S10M", query_seq="TTTTAAAAAAAAAA")
    right = _build_read(ref_start=12, cigar_str="10M4S", query_seq="AAAAAAAAAAGGGG")

    result = same_chrom_same_strand_handler(
        left,
        right,
        left.simple_mode,
        right.simple_mode,
        genome,
        logger,
        microinsertion_cutoff=20,
    )

    assert result is not None
    event_type, positions, read1_info, read2_info, insertion_info, strands = result
    assert event_type == "TDUP"
    assert positions[:2] == ("chr1:10", "chr1:16")
    assert read1_info == (10, 20)
    assert read2_info == (12, 22)
    assert insertion_info == ("-AAAAAA", "-AAAAAA")
    assert strands == (left.strand, right.strand)


def test_write_events_to_vcf_filters_and_writes_file(tmp_path: Path) -> None:
    output_path = tmp_path / "events.vcf"
    logger = _FakeLogger()
    header = {"SQ": [{"SN": "chr1", "LN": 1000}], "PG": [{"ID": "bwa", "CL": "bwa mem ref.fa input.bam"}]}
    events = [
        Event.new("TDUP", ("chr1", 99, 10, "ACTGACTGAA", MicroRegion("+AA")), 8, 5, 20, ref_allele="A", alt_allele="<TDUP>"),
        Event.new("TDUP", ("chr1", 199, 10, "ACTGACTGAA", MicroRegion("")), 4, 1, 20, ref_allele="A", alt_allele="<TDUP>"),
    ]

    write_events_to_vcf(output_path, header, events, logger, min_ao=2, min_depth=10, min_vaf=0.2)

    written = output_path.read_text(encoding="utf-8")
    data_lines = [line for line in written.splitlines() if line and not line.startswith("#")]
    assert "#CHROM\tPOS\tID\tREF\tALT" in written
    assert len(data_lines) == 1
    assert any("Events passing filters" in message for message in logger.messages)