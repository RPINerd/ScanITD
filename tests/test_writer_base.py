"""Tests for the abstract writer base class."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from scanitd.writer.writer import Writer


class _DummyWriter(Writer):
    def write_data(self, data_object: Any, object_id: str):
        return (data_object, object_id)

    def write_line(self, line: str):
        return line

    def open(self, mode: str = "w"):
        return mode

    def close(self):
        return None


def test_writer_initializes_with_path(tmp_path: Path) -> None:
    out = tmp_path / "result.vcf"

    writer = _DummyWriter(str(out))

    assert writer.file_path == out
    assert writer.io is None


def test_writer_warns_when_target_exists(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    out = tmp_path / "result.vcf"
    out.write_text("existing", encoding="utf-8")
    calls: list[str] = []

    monkeypatch.setattr("scanitd.writer.writer.logger.warning", lambda message: calls.append(message))

    _DummyWriter(str(out))

    assert calls
    assert str(out) in calls[0]
