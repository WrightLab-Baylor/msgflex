"""Full run tab — chains SPARX → [XPECTRA] → QUANTIX via the orchestrator."""

from __future__ import annotations

from pathlib import Path

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

FULL_RUN_STEPS = [
    ("sparx", "SPARX"),
    ("xpectra",  "XPECTRA"),
    ("quantix",  "QUANTIX"),
]

_runner = StageRunner(
    prefix="run",
    step_ids=[s[0] for s in FULL_RUN_STEPS],
    run_button_tag="run_run_btn",
    stop_button_tag="run_stop_btn",
    badge_tag="run_status_badge",
    # Match orchestrator stage transitions
    step_matchers=[
        ("sparx", r"Stage 1: SPARX|Stage 1: SKIPPED"),
        ("xpectra", r"Stage 2: XPECTRA"),
        ("quantix", r"Stage 3: QUANTIX"),
    ],
)

def _on_run(state: AppState, log: LogPane):
    def handler():
        if not require_project(state, log):
            return
        if not state.has_params:
            log.append(
                "No MSGFPlus params file selected. Set one on Project tab.",
                "error",
            )
            return

        rescore       = dpg.get_value("run_rescore")
        skip_sparx = dpg.get_value("run_skip_sparx")
        mode          = dpg.get_value("sparx_mode") if dpg.does_item_exist("sparx_mode") else "conventional"

        # Inherit QUANTIX settings from its tab so the user only configures once
        score_field = (
            dpg.get_value("quantix_score_field")
            if dpg.does_item_exist("quantix_score_field")
            else ("xpec_q" if rescore else "MSMSScore")
        )
        threshold = (
            dpg.get_value("quantix_threshold")
            if dpg.does_item_exist("quantix_threshold")
            else (0.01 if rescore else 10.0)
        )
        num_pep = (
            dpg.get_value("quantix_num_pep")
            if dpg.does_item_exist("quantix_num_pep")
            else 2
        )
        output_path = Path(
            dpg.get_value("quantix_output_tsv")
            if dpg.does_item_exist("quantix_output_tsv")
            and dpg.get_value("quantix_output_tsv")
            else state.paths.base / "final.tsv"
        )

        log.append(f" base: {state.paths.base}", "debug")
        log.append(f" mode: {mode}", "debug")
        log.append(f" java_mem: {state.java_mem}", "debug")
        log.append(f" rescore: {rescore}", "debug")
        log.append(f" skip_sparx:{skip_sparx}", "debug")
        log.append(f" output: {output_path}", "debug")

        try:
            from msgflex.orchestrator import workflow
        except ImportError as e:
            log.append(f"Failed to load orchestrator: {e}", "error")
            return

        log.append("Starting full pipeline...", "info")
        _runner.start(
            log_pane=log,
            target=lambda: workflow.run_full(
                paths=state.paths,
                conf_loc=str(state.params_file),
                mode=mode,
                java_mem=state.java_mem,
                rescore=rescore,
                skip_sparx=skip_sparx,
                final_output=output_path,
                score_field=score_field,
                threshold=float(threshold),
                num_pep=int(num_pep),
            ),
        )

    return handler

def _on_dry_run(state: AppState, log: LogPane):
    def handler():
        if not require_project(state, log):
            return
        rescore = " --rescore" if dpg.get_value("run_rescore") else ""
        skip    = " --skip-sparx" if dpg.get_value("run_skip_sparx") else ""
        params  = state.params_file or "(none)"
        mode    = dpg.get_value("saprx_mode") if dpg.does_item_exist("saprx_mode") else "conventional"
        log.append("Dry run: would execute `msgflex run` with:", "info")
        log.append(
            f"  msgflex run -b {state.paths.base} -c {params} "
            f"--mode {mode} --java-mem {state.java_mem}{rescore}{skip}",
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
            title="Full run — SPARX → [XPECTRA] → QUANTIX",
            description=(
                "Chains all three stages. Optionally insert XPECTRA rescoring, "
                "or skip sparx to reuse an existing SICdir. Per-stage settings "
                "are read from their tabs."
            ),
            stage="neutral",
            badge_tag="run_status_badge",
            themes=themes,
        )

        with dpg.group(horizontal=True):
            dpg.add_checkbox(
                label="Insert XPECTRA rescoring",
                tag="run_rescore",
                default_value=False,
            )
            dpg.add_spacer(width=24)
            dpg.add_checkbox(
                label="Skip Sparx (reuse SICdir)",
                tag="run_skip_sparx",
                default_value=False,
            )

        dpg.add_spacer(height=12)
        dpg.add_separator()
        dpg.add_spacer(height=10)

        build_progress_grid(
            FULL_RUN_STEPS,
            progress_prefix="run",
            stage="neutral",
            themes=themes,
        )

        build_action_bar(
            stage="neutral",
            run_callback=_on_run(state, log_pane),
            dry_run_callback=_on_dry_run(state, log_pane),
            stop_callback=_on_stop(log_pane),
            themes=themes,
            run_button_tag="run_run_btn",
            stop_button_tag="run_stop_btn",
        )
