"""Tests for CLI entrypoint and validators."""

from pathlib import Path

import pytest
import typer
from typer.testing import CliRunner

from scanitd import __version__
from scanitd.cli.cli import LogLevel, app, itd_len_type, version_callback


runner = CliRunner()


def test_itd_len_type_rejects_negative() -> None:
    with pytest.raises(typer.BadParameter):
        itd_len_type(-1)


def test_itd_len_type_accepts_positive() -> None:
    assert itd_len_type(10) == 10


def test_version_callback_exits_and_prints(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(typer.Exit):
        version_callback(True)

    captured = capsys.readouterr()
    assert f"ScanITD version: {__version__}" in captured.out


def test_version_callback_noop_when_false() -> None:
    assert version_callback(False) is None


def test_log_level_values_are_expected() -> None:
    assert LogLevel.INFO.value == "info"
    assert LogLevel.WARNING.value == "warning"
    assert LogLevel.ERROR.value == "error"
    assert LogLevel.DEBUG.value == "debug"
    assert LogLevel.TRACE.value == "trace"


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
