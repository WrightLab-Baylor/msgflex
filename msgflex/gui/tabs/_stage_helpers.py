"""
Shared UI pattern for stage tabs (SPARX, XPECTRA, QUANTIX).
This module provides build_stage_header(), build_progress_grid(), and
build_action_bar() that individual stage tabs compose.
"""

from __future__ import annotations

from typing import Callable

import dearpygui.dearpygui as dpg

from msgflex.gui.state import AppState
from msgflex.gui.theme import TOKENS, get_stage_color


def build_stage_header(
    title: str,
    description: str,
    stage: str,
    badge_tag: str,
    themes: dict,
) -> None:
    """Render stage title + inline status badge + description."""
    with dpg.group(horizontal=True):
        dpg.add_text(title, color=get_stage_color(stage, "800"))
        # Status badge — hidden until a run starts
        status_btn = dpg.add_button(
            label="IDLE",
            tag=badge_tag,
            small=True,
            enabled=False,
        )
        dpg.bind_item_theme(status_btn, themes[f"badge_{stage}"])
    dpg.add_text(description, color=TOKENS["text_muted"], wrap=520)
    dpg.add_spacer(height=14)


def build_progress_grid(
    steps: list[tuple[str, str]],
    progress_prefix: str,
    stage: str,
    themes: dict,
) -> None:
    """
    Render a 3-column progress grid (step name, bar, status text).
    """
    dpg.add_text("Pipeline progress", color=TOKENS["text"])
    dpg.add_spacer(height=6)

    with dpg.table(
        header_row=False,
        borders_innerH=False, borders_outerH=False,
        borders_innerV=False, borders_outerV=False,
        policy=dpg.mvTable_SizingFixedFit,
        row_background=False,
    ):
        dpg.add_table_column(width_fixed=True, init_width_or_weight=110)
        dpg.add_table_column(width_stretch=True)
        dpg.add_table_column(width_fixed=True, init_width_or_weight=80)

        for step_id, step_name in steps:
            with dpg.table_row():
                dpg.add_text(step_name, color=TOKENS["text_muted"])
                bar = dpg.add_progress_bar(
                    default_value=0.0,
                    width=-1,
                    tag=f"{progress_prefix}_{step_id}_bar",
                )
                dpg.bind_item_theme(bar, themes[stage])
                dpg.add_text(
                    "queued",
                    tag=f"{progress_prefix}_{step_id}_status",
                    color=TOKENS["text_faint"],
                )


def build_action_bar(
    stage: str,
    run_callback: Callable,
    dry_run_callback: Callable,
    stop_callback: Callable,
    themes: dict,
    run_button_tag: str,
    stop_button_tag: str,
) -> None:
    """Render the Run / Dry run / Stop button row."""
    dpg.add_spacer(height=12)
    dpg.add_separator()
    dpg.add_spacer(height=10)

    with dpg.group(horizontal=True):
        run_btn = dpg.add_button(
            label=f"▶  Run {stage}",
            callback=run_callback,
            tag=run_button_tag,
            width=150,
            height=34,
        )
        dpg.bind_item_theme(run_btn, themes[stage])

        dpg.add_button(
            label="Dry run",
            callback=dry_run_callback,
            width=100,
            height=34,
        )

        dpg.add_spacer(width=8)

        stop_btn = dpg.add_button(
            label="■  Stop",
            callback=stop_callback,
            tag=stop_button_tag,
            width=100,
            height=34,
            enabled=False,
        )
        dpg.bind_item_theme(stop_btn, themes["danger"])


def require_project(state: AppState, log_pane) -> bool:
    """Guard helper — call before starting any stage. Returns True if OK."""
    if not state.has_project:
        log_pane.append(
            "No project selected. Set a base directory on the Project tab first.",
            "error",
        )
        return False
    return True
