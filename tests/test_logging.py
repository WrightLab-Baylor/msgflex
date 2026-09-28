"""Tests for the rewritten core.logging module (modes, verbosity, per-sample)."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from msgflex.core.logging import (
    setup_logging,
    _level_for,
)

@pytest.fixture(autouse=True)
def reset_msgflex_logger():
    """Reset the msgflex logger after every test so they don't pollute each other."""
    yield
    root = logging.getLogger("msgflex")
    root.handlers.clear()
    root.setLevel(logging.NOTSET)
    rx = logging.getLogger("rxflow")
    rx.handlers.clear()
    rx.setLevel(logging.NOTSET)


class TestVerbosityLevels:
    """The three named verbosities map to standard logging levels."""

    def test_quiet_is_warning(self):
        assert _level_for("quiet") == logging.WARNING

    def test_normal_is_info(self):
        assert _level_for("normal") == logging.INFO

    def test_verbose_is_debug(self):
        assert _level_for("verbose") == logging.DEBUG


class TestModes:
    """CLI mode adds a stream handler; GUI mode does not."""

    def test_cli_mode_adds_stream_handler(self):
        setup_logging(mode="cli", verbosity="normal")
        root = logging.getLogger("msgflex")
        stream_handlers = [
            h for h in root.handlers if isinstance(h, logging.StreamHandler)
            and not isinstance(h, logging.FileHandler)
        ]
        assert len(stream_handlers) == 1

    def test_gui_mode_no_stream_handler(self):
        setup_logging(mode="gui", verbosity="normal")
        root = logging.getLogger("msgflex")
        stream_handlers = [
            h for h in root.handlers if isinstance(h, logging.StreamHandler)
            and not isinstance(h, logging.FileHandler)
        ]
        assert len(stream_handlers) == 0

    def test_logfile_added_in_both_modes(self, tmp_path):
        log = tmp_path / "test.log"
        setup_logging(mode="gui", verbosity="normal", logfile=log)
        root = logging.getLogger("msgflex")
        file_handlers = [
            h for h in root.handlers if isinstance(h, logging.FileHandler)
        ]
        assert len(file_handlers) == 1

class TestBackwardsCompat:
    """Old kwargs (`level`, `quiet`) keep working."""

    def test_level_debug_maps_to_verbose(self):
        setup_logging(level="DEBUG")
        # In CLI mode (default), the stream handler should be at DEBUG level
        root = logging.getLogger("msgflex")
        stream = next(
            h for h in root.handlers
            if isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
        )
        assert stream.level == logging.DEBUG

    def test_quiet_true_lifts_threshold(self):
        setup_logging(quiet=True)
        root = logging.getLogger("msgflex")
        stream = next(
            h for h in root.handlers
            if isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler)
        )
        assert stream.level == logging.WARNING
