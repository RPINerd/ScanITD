"""Tests for logger type protocol module imports and method contract shape."""

from __future__ import annotations

from typing import Any

from scanitd.mtype import LoggerType


class _LoggerImpl:
    def trace(self, msg: str) -> None:
        self._last = msg

    def debug(self, msg: str) -> None:
        self._last = msg

    def info(self, msg: str) -> None:
        self._last = msg

    def warning(self, msg: str) -> None:
        self._last = msg

    def error(self, msg: str) -> None:
        self._last = msg

    def critical(self, msr: str) -> None:
        self._last = msr

    def success(self, msg: str) -> None:
        self._last = msg

    def complete(self) -> Any:
        return {"status": "done"}


def _consume_logger(logger: LoggerType) -> Any:
    logger.trace("t")
    logger.debug("d")
    logger.info("i")
    logger.warning("w")
    logger.error("e")
    logger.critical("c")
    logger.success("s")
    return logger.complete()


def test_logger_type_usage_contract() -> None:
    logger = _LoggerImpl()

    result = _consume_logger(logger)

    assert result == {"status": "done"}
