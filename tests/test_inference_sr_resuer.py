"""Tests for sequence rescue helpers."""

from __future__ import annotations

import pytest

from scanitd.base import MappingMode, MicroRegion
from scanitd.inference.sr_resuer import alignment_operation, calculate_variants, update_tdup_ao


class _FakeAlignment:
    def __init__(
        self,
        *,
        reference_start: int,
        reference_end: int,
        read_start: int,
        read_end: int,
        cigar_pair_list: list[tuple[int, str]],
    ) -> None:
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.read_start = read_start
        self.read_end = read_end
        self.cigar_pair_list = cigar_pair_list


class _FakeAlignmentMgr:
    def __init__(self, alignment: _FakeAlignment) -> None:
        self.alignment = alignment
        self.read: str | None = None
        self.reference: str | None = None

    def set_read(self, query_seq: str) -> None:
        self.read = query_seq

    def set_reference(self, reference_seq: str) -> None:
        self.reference = reference_seq

    def align(self, gap_open: int, gap_extension: int) -> _FakeAlignment:
        assert gap_open == 3
        assert gap_extension == 1
        return self.alignment


class _FakeGenomeSlice:
    def __init__(self, seq: str) -> None:
        self.seq = seq


class _FakeGenomeChrom:
    def __init__(self, sequence: str) -> None:
        self.sequence = sequence

    def __getitem__(self, item: slice) -> _FakeGenomeSlice:
        return _FakeGenomeSlice(self.sequence[item])


class _FakeGenome:
    def __init__(self, sequence_by_chrom: dict[str, str]) -> None:
        self.sequence_by_chrom = sequence_by_chrom

    def __getitem__(self, chrom: str) -> _FakeGenomeChrom:
        return _FakeGenomeChrom(self.sequence_by_chrom[chrom])


def test_calculate_variants_counts_mismatches_insertions_and_deletions() -> None:
    variants = calculate_variants(
        [(3, "M"), (2, "I"), (2, "M"), (1, "D")],
        reference_seq="AATCCG",
        query_seq="AACGTCC",
        reference_start=0,
        reference_end=6,
        query_start=0,
        query_end=7,
    )

    assert variants == {"insertions": 1, "deletions": 1, "snvs": 1, "mismatches": 4}


def test_alignment_operation_accepts_matching_sm_alignment() -> None:
    manager = _FakeAlignmentMgr(
        _FakeAlignment(
            reference_start=1,
            reference_end=4,
            read_start=0,
            read_end=3,
            cigar_pair_list=[(4, "M")],
        ),
    )

    assert alignment_operation(manager, "ACGT", "AACGT", MappingMode.SM, mismatches_cutoff=0) is True


def test_alignment_operation_rejects_high_mismatch_alignment() -> None:
    manager = _FakeAlignmentMgr(
        _FakeAlignment(
            reference_start=0,
            reference_end=3,
            read_start=0,
            read_end=3,
            cigar_pair_list=[(4, "M")],
        ),
    )

    assert alignment_operation(manager, "TTTT", "AAAA", MappingMode.MS, mismatches_cutoff=1) is False


def test_update_tdup_ao_counts_rescued_sequences(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[tuple[str, str, MappingMode, int]] = []

    def fake_alignment_operation(align_mgr, query_seq: str, reference_seq: str, read_mode: MappingMode, mismatches_cutoff: int = 5) -> bool:
        calls.append((query_seq, reference_seq, read_mode, mismatches_cutoff))
        return query_seq.endswith("ok")

    monkeypatch.setattr("scanitd.inference.sr_resuer.alignment_operation", fake_alignment_operation)

    genome = _FakeGenome({"chr1": "AACCGGTTAACCGGTTAACCGGTTAACCGGTT"})
    tdup_id = ("chr1", 10, 4, "ACTG", MicroRegion("+GG"))
    to_be_rescued_sequences = {
        ("chr1", 10, MappingMode.SM): ["left-ok", "left-no"],
        ("chr1", 14, MappingMode.MS): ["right-ok"],
    }

    updated = update_tdup_ao(tdup_id, 3, genome, to_be_rescued_sequences, mismatches_cutoff=2)

    assert updated == 5
    assert len(calls) == 3
    assert all(call[3] == 2 for call in calls)