"""SPARX tab — MS-GF+ conventional/binning search."""

from __future__ import annotations

import dearpygui.dearpygui as dpg

from msgflex.gui.log_pane import LogPane
from msgflex.gui.state import AppState
from msgflex.gui.theme import TOKENS
from msgflex.gui.tabs._run_machinery import StageRunner
from msgflex.gui.tabs._stage_helpers import (
    build_action_bar,
    build_progress_grid,
    build_stage_header,
    require_project,
)

SPARX_STEPS = [
    ("step0", "Step 0 (conv.)"),
    ("msgf",  "MS-GF+"),
    ("tsv",   "MzID → TSV"),
    ("phrp",  "PHRP"),
    ("masic", "MASIC"),
    ("merge", "Merger"),
]

_runner = StageRunner(
    prefix="sparx",
    step_ids=[s[0] for s in SPARX_STEPS],
    run_button_tag="sparx_run_btn",
    stop_button_tag="sparx_stop_btn",
    badge_tag="sparx_status_badge",
    # Match log lines emitted by msgflex.pipelines.conventional / binning.
    # First pattern that matches a record wins; order matters.
    step_matchers=[
        ("step0", r"=== Step 0:|Step-0 complete|Step 0:"),
        ("msgf",  r"=== Step 1: Running MS-GF\+|=== Step 1: MS-GF\+|Step-1 complete"
                  r"|Running MS-GF\+:.*\[part\s*\d+"),
        ("tsv",   r"=== Step 2:|Step-1\.5:|Step-1\.5 complete\.|Step-1\.6 complete\."),
        ("phrp",  r"=== Step 3:|Step-3 already|Starting PHRP"),
        ("masic", r"=== Step 4:|Starting MASIC|Step-4 already"),
        ("merge", r"=== Step 5:|Step-5 already|All intensity merging"),
    ],
)

def _on_run(state: AppState, log: LogPane):
    def handler():
        if not require_project(state, log):
            return
        if not state.has_params:
            log.append(
                "No MSGFPlus params file selected. Set one on the Project tab.",
                "error",
            )
            return

        mode = dpg.get_value("sparx_mode")
        paths = state.paths
        conf = state.params_file
        java_mem = state.java_mem

        log.append(f" base:   {paths.base}",  "debug")
        log.append(f" params: {conf}", "debug")
        log.append(f" mode:   {mode}", "debug")
        log.append(f" heap:   {java_mem}", "debug")

        try:
            if mode == "conventional":
                from msgflex.pipelines.conventional import run as run_sparx
            elif mode == "binning":
                from msgflex.pipelines.binning import run as run_sparx
            else:
                log.append(f"Unknown saprx mode: {mode}", "error")
                return
        except ImportError as e:
            log.append(f"Failed to load sparx pipeline: {e}", "error")
            return

        log.append("Step 0: file conversion + QC plots...", "info")
        _runner.start(
            log_pane=log,
            target=lambda: run_sparx(
                paths=paths, 
                conf_loc=conf, 
                java_mem=java_mem,
                stop_event=_runner._stop_event,
            ),
        )

    return handler

def _on_dry_run(state: AppState, log: LogPane):
    def handler():
        if not require_project(state, log):
            return
        log.append("Dry run: would execute `msgflex sparx` with:", "info")
        mode   = dpg.get_value("sparx_mode") if dpg.does_item_exist("sparx_mode") else "conventional"
        params = state.params_file or "(none)"
        log.append(
            f"  msgflex sparx -b {state.paths.base} -c {params} "
            f"--mode {mode} --java-mem {state.java_mem}",
            "debug",
        )
    return handler


def _on_stop(log: LogPane):
    def handler():
        _runner.stop(log)
    return handler


def build(
    parent: int | str,
    state: AppState,
    log_pane: LogPane,
    themes: dict,
) -> None:
    with dpg.group(parent=parent):
        build_stage_header(
            title="SPARX — MS-GF+ search",
            description=(
                "Converts RAW → mzML, runs MS-GF+, PHRP, MASIC, then merges "
                "per-sample outputs into SICdir."
            ),
            stage="sparx",
            badge_tag="sparx_status_badge",
            themes=themes,
        )

        with dpg.group(horizontal=True):
            with dpg.group():
                dpg.add_text("Mode", color=TOKENS["text_muted"])
                dpg.add_radio_button(
                    items=["conventional", "binning"],
                    default_value="conventional",
                    horizontal=True,
                    tag="sparx_mode",
                )
            dpg.add_spacer(width=24)
            with dpg.group():
                dpg.add_text("Java heap (from Project tab)", color=TOKENS["text_muted"])
                dpg.add_text(
                    state.java_mem,
                    tag="sparx_heap_display",
                    color=TOKENS["text"],
                )

        dpg.add_spacer(height=12)
        dpg.add_separator()
        dpg.add_spacer(height=10)

        build_progress_grid(
            SPARX_STEPS,
            progress_prefix="sparx",
            stage="sparx",
            themes=themes,
        )

        build_action_bar(
            stage="sparx",
            run_callback=_on_run(state, log_pane),
            dry_run_callback=_on_dry_run(state, log_pane),
            stop_callback=_on_stop(log_pane),
            themes=themes,
            run_button_tag="sparx_run_btn",
            stop_button_tag="sparx_stop_btn",
        )

    state.subscribe(
        lambda s: dpg.set_value("sparx_heap_display", s.java_mem)
        if dpg.does_item_exist("sparx_heap_display")
        else None
    )
