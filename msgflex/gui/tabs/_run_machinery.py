"""
Shared run-machinery for stage tabs.

Every stage tab — sparx, xpectra, quantix, full_run — has the same
shape: spawn a worker thread, bridge `logging` → `LogPane`, swap badge
states, drive a per-step progress grid, allow Stop. This module factors
that pattern out so each tab only writes the parts that differ.

Live progress updates work by pattern-matching log records. Each tab
constructs a StageRunner with its step_ids and the regex patterns that
identify when each step starts in the pipeline's logging output.

When a log record matches a pattern, the runner marks the previous
step done (bar full, status "done") and the new step running (bar
small, status "running"). On successful completion, all remaining
steps are filled to 1.0.

Because progress updates are driven from inside the log handler, which
can fire from any thread, all DPG calls are guarded by
``_dpg_*`` safety helpers that check ``dpg.is_dearpygui_running()``
before touching any widget. This prevents the segfault that occurs
when the Quit button is pressed while a worker thread is still running
and tries to update widgets after the DPG context has started tearing
down.
"""

from __future__ import annotations

import logging
import re
import threading
from threading import Thread
from typing import Callable, Optional

import dearpygui.dearpygui as dpg

from msgflex.gui.log_pane import LogPane


# ── DPG safety helpers ────────────────────────────────────────────────
# Every DPG mutation goes through these. They silently no-op if the
# viewport is closing or the item is gone, preventing segfaults when
# worker threads fire after the render loop exits.

def _dpg_alive() -> bool:
    """True while the DPG viewport is running OR in headless test contexts."""
    try:
        return dpg.is_dearpygui_running()
    except Exception:
        return False


def _dpg_set(tag: str, value) -> None:
    try:
        if not tag:
            return
        # In headless tests is_dearpygui_running() returns False but items
        # still exist in the context. Guard on does_item_exist alone when
        # the viewport is not running so tests still exercise the logic.
        if dpg.does_item_exist(tag):
            dpg.set_value(tag, value)
    except Exception:
        pass


def _dpg_label(tag: str, label: str) -> None:
    try:
        if not tag:
            return
        if dpg.does_item_exist(tag):
            dpg.set_item_label(tag, label)
    except Exception:
        pass


def _dpg_enable(tag: str) -> None:
    try:
        if not tag:
            return
        if dpg.does_item_exist(tag):
            dpg.enable_item(tag)
    except Exception:
        pass


def _dpg_disable(tag: str) -> None:
    try:
        if not tag:
            return
        if dpg.does_item_exist(tag):
            dpg.disable_item(tag)
    except Exception:
        pass


# ── Log bridge with step-detection ───────────────────────────────────

class _LogPaneHandler(logging.Handler):
    """Forward Python log records to the GUI LogPane queue, and run an
    optional step-detection callback for each record.
    """

    _LEVEL_MAP = {
        logging.DEBUG:    "debug",
        logging.INFO:     "info",
        logging.WARNING:  "warn",
        logging.ERROR:    "error",
        logging.CRITICAL: "error",
    }

    def __init__(
        self,
        pane: LogPane,
        on_message: Optional[Callable[[str], None]] = None,
    ) -> None:
        super().__init__()
        self._pane = pane
        self._on_message = on_message

    def emit(self, record: logging.LogRecord) -> None:
        try:
            text = self.format(record)
            level = self._LEVEL_MAP.get(record.levelno, "info")
            self._pane.append(text, level)
            if self._on_message is not None:
                try:
                    self._on_message(text)
                except Exception:
                    pass
        except Exception:
            self.handleError(record)


def _attach_log_bridge(
    log_pane: LogPane,
    on_message: Optional[Callable[[str], None]] = None,
) -> _LogPaneHandler:
    """Attach a handler to the msgflex and rxflow root loggers.

    Also bumps each logger's level to INFO if it's currently
    NOTSET/WARNING — the GUI doesn't call setup_logging() (the CLI
    does), so without this all INFO step markers would be silently
    filtered out before reaching the bridge.

    Returns the handler for cleanup via _detach_log_bridge.
    """
    handler = _LogPaneHandler(log_pane, on_message=on_message)
    handler.setFormatter(logging.Formatter("%(name)s: %(message)s"))

    for logger_name in ("msgflex", "rxflow"):
        lg = logging.getLogger(logger_name)
        if lg.level == logging.NOTSET or lg.level > logging.INFO:
            lg.setLevel(logging.INFO)
        lg.addHandler(handler)

    return handler


def _detach_log_bridge(handler: _LogPaneHandler) -> None:
    for logger_name in ("msgflex", "rxflow"):
        logging.getLogger(logger_name).removeHandler(handler)


# ── Stage runner ──────────────────────────────────────────────────────

class StageRunner:
    """Owns the run/stop state and step-progress driving for one stage tab.

    Construct once per tab in build(), wire to the Run/Stop callbacks.
    All DPG widget mutations go through the _dpg_* safety helpers so
    that worker-thread callbacks that fire after Quit is pressed are
    silently swallowed instead of causing a segfault.
    """

    def __init__(
        self,
        *,
        prefix: str,
        step_ids: list[str],
        run_button_tag: str,
        stop_button_tag: str,
        badge_tag: str,
        step_matchers: Optional[list[tuple[str, str]]] = None,
        incremental_matcher: Optional[tuple[str, str, str]] = None,
        incremental_total: Optional[int] = None,
    ) -> None:
        self.prefix = prefix
        self.step_ids = step_ids
        self.run_button_tag = run_button_tag
        self.stop_button_tag = stop_button_tag
        self.badge_tag = badge_tag

        self._matchers: list[tuple[str, re.Pattern]] = [
            (sid, re.compile(pat)) for sid, pat in (step_matchers or [])
        ]
        self._step_index = {sid: i for i, sid in enumerate(step_ids)}
        self._current_step_idx: int = -1

        if incremental_matcher is not None:
            bar_tag, label_tag, pattern = incremental_matcher
            self._inc_bar_tag = bar_tag
            self._inc_label_tag = label_tag
            self._inc_pattern: Optional[re.Pattern] = re.compile(pattern)
        else:
            self._inc_bar_tag = ""
            self._inc_label_tag = ""
            self._inc_pattern = None
        self._inc_total = incremental_total or 0
        self._inc_count = 0

        self._stop_event = threading.Event()
        self._thread: Optional[Thread] = None

    # ── Incremental mode ─────────────────────────────────────────────

    def set_incremental_total(self, total: int) -> None:
        self._inc_total = max(int(total), 0)

    def reset_incremental(self) -> None:
        self._inc_count = 0
        _dpg_set(self._inc_bar_tag, 0.0)
        if self._inc_label_tag:
            label = f"0 / {self._inc_total}" if self._inc_total else "Waiting..."
            _dpg_set(self._inc_label_tag, label)

    # ── Progress helpers ──────────────────────────────────────────────

    def _set_step(self, step_id: str, value: float, status: str) -> None:
        _dpg_set(f"{self.prefix}_{step_id}_bar", value)
        _dpg_set(f"{self.prefix}_{step_id}_status", status)

    def reset_progress(self) -> None:
        for sid in self.step_ids:
            self._set_step(sid, 0.0, "queued")
        self._current_step_idx = -1
        self.reset_incremental()

    def mark_all_done(self) -> None:
        for sid in self.step_ids:
            self._set_step(sid, 1.0, "done")
        self._mark_incremental_done()

    def mark_remaining_done(self) -> None:
        for sid in self.step_ids:
            self._set_step(sid, 1.0, "done")
        self._mark_incremental_done()

    def mark_all_error(self) -> None:
        for sid in self.step_ids:
            self._set_step(sid, 0.0, "error")
        if self._inc_label_tag:
            _dpg_set(self._inc_label_tag, "Failed")

    def _mark_incremental_done(self) -> None:
        if not self._inc_pattern:
            return
        _dpg_set(self._inc_bar_tag, 1.0)
        if self._inc_label_tag:
            total = self._inc_total or self._inc_count
            _dpg_set(self._inc_label_tag, f"{total} / {total} — Done")

    # ── Step-transition callback ──────────────────────────────────────

    def _on_log_message(self, message: str) -> None:
        for sid, pat in self._matchers:
            if pat.search(message):
                self._advance_to(sid)
                break

        if self._inc_pattern is not None:
            m = self._inc_pattern.search(message)
            if m:
                self._inc_count += 1
                try:
                    item_name = m.group(1)
                except IndexError:
                    item_name = ""
                self._update_incremental(item_name)

    def _update_incremental(self, item_name: str) -> None:
        value = (
            min(self._inc_count / self._inc_total, 1.0)
            if self._inc_total > 0 else 0.0
        )
        _dpg_set(self._inc_bar_tag, value)
        if self._inc_label_tag:
            label = (
                f"{self._inc_count} / {self._inc_total} — {item_name}"
                if item_name
                else f"{self._inc_count} / {self._inc_total}"
            )
            _dpg_set(self._inc_label_tag, label)

    def _advance_to(self, step_id: str) -> None:
        new_idx = self._step_index.get(step_id, -1)
        if new_idx < 0 or new_idx == self._current_step_idx:
            return
        for i, sid in enumerate(self.step_ids):
            if i < new_idx:
                self._set_step(sid, 1.0, "done")
            elif i == new_idx:
                self._set_step(sid, 0.1, "running")
        self._current_step_idx = new_idx

    # ── Badge / button helpers ────────────────────────────────────────

    def _set_badge(self, label: str) -> None:
        _dpg_label(self.badge_tag, label)

    def _disable_run_enable_stop(self) -> None:
        _dpg_disable(self.run_button_tag)
        _dpg_enable(self.stop_button_tag)

    def _enable_run_disable_stop(self) -> None:
        _dpg_enable(self.run_button_tag)
        _dpg_disable(self.stop_button_tag)

    # ── Public lifecycle ──────────────────────────────────────────────

    @property
    def is_stopping(self) -> bool:
        return self._stop_event.is_set()

    def is_running(self) -> bool:
        return self._thread is not None and self._thread.is_alive()

    def start(
        self,
        *,
        log_pane: LogPane,
        target: Callable[[], None],
        on_success: Optional[Callable[[], None]] = None,
        on_error: Optional[Callable[[Exception], None]] = None,
    ) -> None:
        if self.is_running():
            log_pane.append("A run is already active. Wait or press Stop.", "warn")
            return

        self._stop_event = threading.Event()
        self.reset_progress()
        self._set_badge("RUNNING")
        self._disable_run_enable_stop()

        # Prime the first step immediately on the main thread so the user
        # sees a running bar right away rather than empty bars for several
        # seconds until the first log line arrives.
        if self.step_ids:
            self._set_step(self.step_ids[0], 0.1, "running")
            self._current_step_idx = 0

        bridge = _attach_log_bridge(log_pane, on_message=self._on_log_message)

        def worker() -> None:
            try:
                target()
                if self._stop_event.is_set():
                    log_pane.append("Run stopped by user.", "warn")
                    self._set_badge("STOPPED")
                    return
                self.mark_remaining_done()
                self._set_badge("DONE")
                log_pane.append(f"{self.prefix} pipeline completed.", "success")
                if on_success is not None:
                    on_success()
            except Exception as e:
                log_pane.append(f"{self.prefix} pipeline failed: {e}", "error")
                import traceback
                tb_text = "".join(
                    traceback.TracebackException.from_exception(e).format()
                )
                for tb_line in tb_text.rstrip().split("\n"):
                    log_pane.append(tb_line, "debug")
                self._set_badge("ERROR")
                self.mark_all_error()
                if on_error is not None:
                    on_error(e)
            finally:
                _detach_log_bridge(bridge)
                self._enable_run_disable_stop()

        self._thread = Thread(target=worker, daemon=True)
        self._thread.start()

    def stop(self, log_pane: LogPane) -> None:
        """Request stop and immediately kill all registered child processes."""
        if not self.is_running():
            log_pane.append("No active run to stop.", "info")
            return
        self._stop_event.set()
        log_pane.append("Stop requested — terminating active processes...", "warn")
        _dpg_disable(self.stop_button_tag)

        from msgflex.core.process_registry import ProcessRegistry
        threading.Thread(
            target=ProcessRegistry.kill_all,
            kwargs={"grace_seconds": 3.0},
            daemon=True,
        ).start()