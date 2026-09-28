"""
Live log pane — the 280px-wide right column shown on every stage tab.

Uses a DPG child_window with auto-scroll. Each log line is a colored
text widget tagged by its level for styling.

Thread safety: append() may be called from any thread (e.g. pipeline
worker threads). Messages are put into a queue.Queue and flushed by
_drain(), which app.py calls once per frame inside the manual render
loop — always on the main DPG thread.
"""

from __future__ import annotations

import queue
from dataclasses import dataclass, field
from typing import Literal

import dearpygui.dearpygui as dpg

from msgflex.gui.theme import TOKENS, get_stage_color

LogLevel = Literal["info", "warn", "error", "sample", "debug", "success"]


@dataclass
class LogPane:
    """
    Scrollable log display bound to a parent container.
    """

    parent: int | str
    state: object = None          # AppState reference for project path resolution
    _window_id: int = 0
    _content_id: int = 0
    _max_lines: int = 2000
    _lines: list[int] = field(default_factory=list)
    _queue: queue.Queue = field(default_factory=queue.Queue)

    def build(self) -> None:
        with dpg.group(parent=self.parent):
            with dpg.group(horizontal=True):
                dpg.add_text("LIVE LOG", color=TOKENS["text_muted"])
                dpg.add_spacer(width=-1)

            self._window_id = dpg.add_child_window(
                autosize_x=True,
                height=-60,
                border=True,
                horizontal_scrollbar=True,
            )
            self._content_id = self._window_id

            with dpg.group(horizontal=True):
                dpg.add_button(
                    label="Clear log",
                    width=-110,
                    callback=lambda: self.clear(),
                )
                dpg.add_button(
                    label="Save log",
                    width=-1,
                    callback=lambda: self._save_to_file(),
                )

    # ── Public API ────────────────────────────────────────────────────

    def append(self, text: str, level: LogLevel = "info") -> None:
        """Queue one log line — safe to call from any thread."""
        self._queue.put((text, level))

    def _drain(self) -> None:
        if self._content_id == 0:
            return
        try:
            if not dpg.is_dearpygui_running():
                return
        except Exception:
            return
        try:
            while True:
                text, level = self._queue.get_nowait()
                self._write(text, level)
        except queue.Empty:
            pass
        except Exception:
            pass

    def clear(self) -> None:
        while not self._queue.empty():
            try:
                self._queue.get_nowait()
            except queue.Empty:
                break
        if self._content_id and dpg.does_item_exist(self._content_id):
            for line_id in self._lines:
                if dpg.does_item_exist(line_id):
                    dpg.delete_item(line_id)
            self._lines.clear()

    # ── Internal helpers ──────────────────────────────────────────────

    def _write(self, text: str, level: LogLevel) -> None:
        """Create the DPG text widget — main thread only."""
        try:
            if not dpg.is_dearpygui_running():
                return
        except Exception:
            return
        try:
            color  = self._color_for_level(level)
            prefix = self._prefix_for_level(level)

            line_id = dpg.add_text(
                f"{prefix}{text}",
                color=color,
                parent=self._content_id,
                wrap=0,
            )
            self._lines.append(line_id)

            while len(self._lines) > self._max_lines:
                old = self._lines.pop(0)
                if dpg.does_item_exist(old):
                    dpg.delete_item(old)

            dpg.set_y_scroll(self._content_id, -1.0)
        except Exception:
            pass

    def _save_to_file(self) -> None:
        if not self._lines:
            return
        import datetime, pathlib

        # Save to project base dir if a project is loaded, otherwise home dir
        if self.state is not None and getattr(self.state, "paths", None) is not None:
            log_path = pathlib.Path(self.state.paths.base) / "msgflex.log"
        else:
            log_path = pathlib.Path.home() / "msgflex.log"

        texts = []
        for line_id in self._lines:
            if dpg.does_item_exist(line_id):
                texts.append(dpg.get_value(line_id))
        try:
            with log_path.open("w", encoding="utf-8") as fh:
                fh.write(f"# msgflex log — {datetime.datetime.now():%Y-%m-%d %H:%M:%S}\n")
                fh.write("\n".join(texts))
                fh.write("\n")
            self.append(f"Log saved → {log_path}", "success")
        except OSError as exc:
            self.append(f"Save failed: {exc}", "error")

    @staticmethod
    def _color_for_level(level: LogLevel) -> tuple[int, int, int, int]:
        if level == "info":
            return TOKENS["info_text"]
        if level == "warn":
            return TOKENS["warning_text"]
        if level == "error":
            return TOKENS["error_text"]
        if level == "success":
            return TOKENS["success_text"]
        if level == "sample":
            return get_stage_color("sparx", "600")
        return TOKENS["text_muted"]

    @staticmethod
    def _prefix_for_level(level: LogLevel) -> str:
        return {
            "info": "[INFO]  ",
            "warn": "[WARN]  ",
            "error": "[ERROR] ",
            "success": "[OK]    ",
            "sample": "        ",
            "debug": "[DEBUG] ",
        }.get(level, "        ")