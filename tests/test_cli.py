"""Tests for CLI entrypoint and validators."""

from pathlib import Path

import pytest
import typer
from typer.testing import CliRunner

from scanitd.cli.cli import app, itd_len_type


runner = CliRunner()


def test_itd_len_type_rejects_negative() -> None:
    with pytest.raises(typer.BadParameter):
        itd_len_type(-1)


def test_itd_len_type_accepts_positive() -> None:
    assert itd_len_type(10) == 10


def test_cli_main_success(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    bam_file = tmp_path / "sample.bam"
    fasta_file = tmp_path / "ref.fa"
    bam_file.write_text("bam", encoding="utf-8")
    fasta_file.write_text(">chr1\nACGT\n", encoding="utf-8")

    def fake_scan_itd(**_: object):
        return [], {"SQ": [{"SN": "chr1", "LN": 1000}]}

    calls: dict[str, object] = {}

    def fake_write_events_to_vcf(output: str, bam_header: dict[str, object], event_list: list[object], logger: object, **kwargs: object) -> None:
        calls["output"] = output
        calls["bam_header"] = bam_header
        calls["event_list"] = event_list
        calls["kwargs"] = kwargs

    monkeypatch.setattr("scanitd.cli.cli.scan_itd", fake_scan_itd)
    monkeypatch.setattr("scanitd.cli.cli.write_events_to_vcf", fake_write_events_to_vcf)

    result = runner.invoke(
        app,
        [
            "--input",
            str(bam_file),
            "--ref",
            str(fasta_file),
            "--output",
            str(tmp_path / "out.vcf"),
        ],
    )

    assert result.exit_code == 0
    assert calls["event_list"] == []
