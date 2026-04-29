"""Tests for inference.main using mocked pysam and Fasta objects."""

from __future__ import annotations

from pathlib import Path

import pytest
from pyfaidx import FastaNotFoundError

from scanitd.base import MicroRegion
from scanitd.inference.main import BamScanner, scan_itd


class _FakeHeader:
    def __init__(self, header_dict: dict) -> None:
        self._header_dict = header_dict

    def as_dict(self) -> dict:
        return self._header_dict


class _FakeAlignmentFile:
    def __init__(self, reads: list, header_dict: dict) -> None:
        self._reads = reads
        self.header = _FakeHeader(header_dict)
        self.closed = False

    def fetch(self, region=None):
        return self._reads

    def pileup(self, region=None, stepper=None, truncate=None):
        return []

    def close(self) -> None:
        self.closed = True


class _FakeLogger:
    def __init__(self) -> None:
        self.info_messages: list[str] = []
        self.trace_messages: list[str] = []

    def info(self, message: str) -> None:
        self.info_messages.append(message)

    def trace(self, message: str) -> None:
        self.trace_messages.append(message)

    def warning(self, message: str) -> None:
        self.trace_messages.append(message)


class _FakeRead:
    def __init__(
        self,
        *,
        query_name: str,
        reference_name: str,
        reference_start: int,
        is_reverse: bool,
        cigarstring: str,
        mapping_quality: int,
        query_sequence: str,
        query_qualities: list[int] | None,
        sa_tag: str,
        nm_tag: int = 0,
        has_xa: bool = False,
        is_supplementary: bool = False,
        is_secondary: bool = False,
    ) -> None:
        self.query_name = query_name
        self.reference_name = reference_name
        self.reference_start = reference_start
        self.is_reverse = is_reverse
        self.cigarstring = cigarstring
        self.mapping_quality = mapping_quality
        self.query_sequence = query_sequence
        self.query_qualities = query_qualities
        self._sa_tag = sa_tag
        self._nm_tag = nm_tag
        self._has_xa = has_xa
        self.is_supplementary = is_supplementary
        self.is_secondary = is_secondary

    def has_tag(self, tag: str) -> bool:
        if tag == "SA":
            return True
        if tag == "XA":
            return self._has_xa
        if tag == "NM":
            return True
        return False

    def get_tag(self, tag: str):
        if tag == "SA":
            return self._sa_tag
        if tag == "NM":
            return self._nm_tag
        msg = f"unsupported tag: {tag}"
        raise KeyError(msg)


def _sorted_header() -> dict:
    return {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}


def test_bam_scanner_initializes_with_mocked_pysam_and_fasta(monkeypatch: pytest.MonkeyPatch) -> None:
    fake_alignment_file = _FakeAlignmentFile([], _sorted_header())

    monkeypatch.setattr("scanitd.inference.main.pysam.AlignmentFile", lambda *_args, **_kwargs: fake_alignment_file)
    monkeypatch.setattr("scanitd.inference.main.Fasta", lambda *_args, **_kwargs: object())

    scanner = BamScanner(
        input_bam=Path("input.bam"),
        mapq_cutoff=15,
        ref_genome=Path("~/ref.fa"),
        microinsertion_cutoff=10,
        regions=[None],
        logger=_FakeLogger(),
    )

    assert scanner.header["HD"]["SO"] == "coordinate"
    assert str(scanner.ref_genome).startswith(str(Path.home()))
    assert scanner.in_bam_object is fake_alignment_file


def test_bam_scanner_unsorted_header_returns_false(monkeypatch: pytest.MonkeyPatch) -> None:
    header = {"HD": {"SO": "queryname"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    fake_alignment_file = _FakeAlignmentFile([], header)

    monkeypatch.setattr("scanitd.inference.main.pysam.AlignmentFile", lambda *_args, **_kwargs: fake_alignment_file)
    monkeypatch.setattr("scanitd.inference.main.Fasta", lambda *_args, **_kwargs: object())

    scanner = BamScanner(
        input_bam=Path("input.bam"),
        mapq_cutoff=15,
        ref_genome=Path("ref.fa"),
        microinsertion_cutoff=10,
        regions=[None],
        logger=_FakeLogger(),
    )

    assert scanner._check_bam_sort(header) is False


def test_bam_scanner_raises_when_fasta_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    fake_alignment_file = _FakeAlignmentFile([], _sorted_header())

    monkeypatch.setattr("scanitd.inference.main.pysam.AlignmentFile", lambda *_args, **_kwargs: fake_alignment_file)

    def _raise_fasta_not_found(*_args, **_kwargs):
        raise FastaNotFoundError("missing")

    monkeypatch.setattr("scanitd.inference.main.Fasta", _raise_fasta_not_found)

    with pytest.raises(SystemExit, match="Reference File"):
        BamScanner(
            input_bam=Path("input.bam"),
            mapq_cutoff=15,
            ref_genome=Path("ref.fa"),
            microinsertion_cutoff=10,
            regions=[None],
            logger=_FakeLogger(),
        )


def test_iter_bam_collects_tdup_anchor_with_mocked_read(monkeypatch: pytest.MonkeyPatch) -> None:
    read = _FakeRead(
        query_name="q1",
        reference_name="chr1",
        reference_start=10,
        is_reverse=False,
        cigarstring="4S10M",
        mapping_quality=60,
        query_sequence="TTTTAAAAAAAAAA",
        query_qualities=[30] * 14,
        sa_tag="chr1,13,+,10M4S,60,0;",
    )
    fake_alignment_file = _FakeAlignmentFile([read], _sorted_header())

    monkeypatch.setattr("scanitd.inference.main.pysam.AlignmentFile", lambda *_args, **_kwargs: fake_alignment_file)
    monkeypatch.setattr("scanitd.inference.main.Fasta", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(
        "scanitd.inference.main.same_chrom_same_strand_handler",
        lambda *_args, **_kwargs: (
            "TDUP",
            ("chr1:10", "chr1:16", 2, 1),
            (10, 20),
            (12, 22),
            ("+AA", "+AA"),
            ("+", "+"),
        ),
    )

    logger = _FakeLogger()
    scanner = BamScanner(
        input_bam=Path("input.bam"),
        mapq_cutoff=15,
        ref_genome=Path("ref.fa"),
        microinsertion_cutoff=10,
        regions=[None],
        logger=logger,
    )

    anchors = scanner.iter_bam()

    assert "q1" in anchors
    chrom, start, end, strand, break_point_region = anchors["q1"]
    assert (chrom, start, end, strand) == ("chr1", 10, 16, "+")
    assert isinstance(break_point_region, MicroRegion)
    assert break_point_region.micro_type == "microinsertion"
    assert break_point_region.sequence == "AA"
    assert logger.info_messages


def test_scan_itd_returns_empty_list_and_closes_bam(monkeypatch: pytest.MonkeyPatch) -> None:
    fake_alignment_file = _FakeAlignmentFile([], _sorted_header())

    monkeypatch.setattr("scanitd.inference.main.parse_target_genomic_coordinates", lambda _target: [])
    monkeypatch.setattr("scanitd.inference.main.pysam.AlignmentFile", lambda *_args, **_kwargs: fake_alignment_file)
    monkeypatch.setattr("scanitd.inference.main.Fasta", lambda *_args, **_kwargs: {"chr1": "A" * 200})

    events, header = scan_itd(
        in_bam_path="input.bam",
        mapq_cutoff=15,
        ref_genome="ref.fa",
        target_file="",
        itd_length_cutoff=10,
        allowed_mismatches_for_sr_rescue=1,
        allowed_mismatches_for_insertion=2,
        logger=_FakeLogger(),
    )

    assert events == []
    assert header["HD"]["SO"] == "coordinate"
    assert fake_alignment_file.closed is True
