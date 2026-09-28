"""
Logging configuration for MSGFLEX.

Single setup function used by both the CLI and GUI. All stages log
through the ``msgflex`` logger tree; XPECTRA also writes to its own
``rxflow`` logger which we adopt as a child of ours.

Three verbosity levels for the terminal:

  - ``"quiet"``    — WARNING and above only
  - ``"normal"``   — INFO and above (default; same content as the GUI live log)
  - ``"verbose"``  — DEBUG and above

GUI mode silences the terminal handler entirely. The GUI's live log
pane attaches its own handler via the runner machinery, so terminal
output would be redundant.

Per-sample log files are written by the per-sample helpers in this
module (``sample_logger`` context manager) — used by the sparx and
quantix pipelines.

Usage at the CLI entry point::

    from msgflex.core.logging import setup_logging
    setup_logging(mode="cli", verbosity="normal",
                  logfile=paths.base / "msgflex.log")

Usage in a library module::
    import logging
    log = logging.getLogger(__name__)
    log.info("...")
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Iterator, Literal, Optional

_FMT_CONSOLE = "%(asctime)s  %(levelname)-7s  %(name)s  %(message)s"
_FMT_FILE    = "%(asctime)s  %(levelname)-7s  %(name)s  %(message)s"
_DATEFMT     = "%Y-%m-%d %H:%M:%S"

Mode = Literal["cli", "gui"]
Verbosity = Literal["quiet", "normal", "verbose"]

def _level_for(verbosity: Verbosity) -> int:
    return {
        "quiet":   logging.WARNING,
        "normal":  logging.INFO,
        "verbose": logging.DEBUG,
    }[verbosity]

def setup_logging(
    *,
    mode: Mode = "cli",
    verbosity: Verbosity = "normal",
    logfile: Optional[Path] = None,
    level: Optional[str] = None,
    quiet: Optional[bool] = None,
) -> logging.Logger:
    """Configure the root msgflex logger.

    Idempotent — safe to call multiple times; old handlers are cleared.

    Args:
        mode:      "cli" enables a terminal stream handler; "gui"
                   suppresses it (the GUI bridge handles its own output).
        verbosity: "quiet" → WARNING+, "normal" → INFO+ (default),
                   "verbose" → DEBUG+.
        logfile:   Optional persistent log file path (always at level INFO+
                   regardless of terminal verbosity, so a full run record
                   is preserved on disk).
        level:     Backwards-compat — pre-rename code passed a level string.
                   If provided, it overrides ``verbosity``.
        quiet:     Backwards-compat — pre-rename ``quiet=True`` flag.
                   If True, sets verbosity to "quiet" unless ``verbosity``
                   is also passed.

    Returns:
        The configured ``msgflex`` logger.
    """
    # Honour back-compat kwargs without breaking existing callers
    if level is not None:
        lv = level.upper()
        if lv == "DEBUG":
            verbosity = "verbose"
        elif lv == "WARNING" or lv == "ERROR":
            verbosity = "quiet"
        else:
            verbosity = "normal"
    if quiet:
        verbosity = "quiet"

    root = logging.getLogger("msgflex")
    root.handlers.clear()
    root.setLevel(logging.DEBUG)
    root.propagate = False

    term_level = _level_for(verbosity)

    # Terminal handler — only when in CLI mode
    if mode == "cli":
        ch = logging.StreamHandler(sys.stderr)
        ch.setLevel(term_level)
        ch.setFormatter(logging.Formatter(_FMT_CONSOLE, _DATEFMT))
        root.addHandler(ch)

    # Persistent log file always at INFO so on-disk record is complete
    if logfile is not None:
        logfile = Path(logfile).expanduser().resolve()
        logfile.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(logfile)
        fh.setLevel(min(logging.INFO, term_level))
        fh.setFormatter(logging.Formatter(_FMT_FILE, _DATEFMT))
        root.addHandler(fh)

    rx = logging.getLogger("rxflow")
    rx.handlers = list(root.handlers)
    rx.setLevel(logging.DEBUG)
    rx.propagate = False

    return root

def get_logger(name: str) -> logging.Logger:
    """Short-hand for ``logging.getLogger(f"msgflex.{name}")``."""
    if not name.startswith("msgflex"):
        name = f"msgflex.{name}"
    return logging.getLogger(name)
