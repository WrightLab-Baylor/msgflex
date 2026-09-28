"""Check tab — preflight validator wrapped in GUI."""

from __future__ import annotations

import dearpygui.dearpygui as dpg

from msgflex.core.preflight import Status, run_all_checks
from msgflex.gui.state import AppState
from msgflex.gui.log_pane import LogPane
from msgflex.gui.theme import TOKENS


_STATUS_COLOR = {
    Status.OK:   TOKENS["success_text"],
    Status.WARN: TOKENS["warning_text"],
    Status.FAIL: TOKENS["error_text"],
    Status.SKIP: TOKENS["text_faint"],
}

_STATUS_GLYPH = {
    Status.OK:   "✓",
    Status.WARN: "!",
    Status.FAIL: "✗",
    Status.SKIP: "-",
}


def _run_check(state: AppState, log: LogPane, stage: str):
    def handler():
        log.append(f"Running preflight (stage={stage})...", "info")

        # Clear previous results
        if dpg.does_item_exist("check_results_container"):
            dpg.delete_item("check_results_container", children_only=True)

        paths = state.paths  # may be None for env+tools only
        try:
            reports = run_all_checks(paths=paths, stage=stage)
        except Exception as e:
            log.append(f"Preflight error: {e}", "error")
            return

        total_ok = total_warn = total_fail = 0
        container = "check_results_container"

        for report in reports:
            # Category header
            dpg.add_spacer(height=6, parent=container)
            dpg.add_text(
                report.category,
                color=TOKENS["text"],
                parent=container,
            )
            dpg.add_separator(parent=container)

            for result in report.results:
                color = _STATUS_COLOR[result.status]
                glyph = _STATUS_GLYPH[result.status]
                with dpg.group(horizontal=True, parent=container):
                    dpg.add_text(f" {glyph} ", color=color)
                    dpg.add_text(result.name, color=TOKENS["text"])
                    dpg.add_text(f"  {result.message}", color=TOKENS["text_muted"])

                for detail in result.details:
                    dpg.add_text(
                        f"      {detail}",
                        color=TOKENS["text_faint"],
                        parent=container,
                    )

                if result.status == Status.OK:
                    total_ok += 1
                elif result.status == Status.WARN:
                    total_warn += 1
                elif result.status == Status.FAIL:
                    total_fail += 1

        # Summary
        dpg.add_spacer(height=10, parent=container)
        with dpg.group(horizontal=True, parent=container):
            dpg.add_text(f"✓ {total_ok} passed", color=TOKENS["success_text"])
            dpg.add_text("  ·  ")
            dpg.add_text(f"! {total_warn} warnings", color=TOKENS["warning_text"])
            dpg.add_text("  ·  ")
            dpg.add_text(f"✗ {total_fail} failed", color=TOKENS["error_text"])

        if total_fail:
            log.append(f"Preflight finished: {total_fail} failures", "error")
        elif total_warn:
            log.append(f"Preflight finished: {total_warn} warnings", "warn")
        else:
            log.append("Preflight finished: all checks passed", "success")

    return handler


def build(
    parent: int | str,
    state: AppState,
    log_pane: LogPane,
    themes: dict,
) -> None:
    with dpg.group(parent=parent):
        dpg.add_text("Preflight check", color=TOKENS["text"])
        dpg.add_text(
            "Validates environment, tools, and project layout without running the pipeline.",
            color=TOKENS["text_muted"], wrap=520,
        )
        dpg.add_spacer(height=12)

        dpg.add_text("Scope", color=TOKENS["text_muted"])
        dpg.add_radio_button(
            items=["all", "sparx", "xpectra", "quantix"],
            default_value="all",
            horizontal=True,
            tag="check_stage",
        )
        dpg.add_spacer(height=10)

        with dpg.group(horizontal=True):
            btn = dpg.add_button(
                label="▶  Run preflight",
                width=160, height=34,
                callback=lambda: _run_check(
                    state, log_pane,
                    dpg.get_value("check_stage"),
                )(),
            )
            dpg.bind_item_theme(btn, themes["neutral"])

        dpg.add_spacer(height=14)
        dpg.add_separator()
        dpg.add_spacer(height=6)

        dpg.add_child_window(
            tag="check_results_container",
            autosize_x=True,
            height=-1,
            border=False,
        )
