"""XPECTRA tab — ML-based PSM rescoring.

Uses **per-sample progress** instead of a 4-step grid because rxflow
processes each sample sequentially through its own
features→rescore→stats→write cycle. The 4-step grid was misleading: it
sat at "features" for the entire run because every sample fired the
"Processing sample" log marker which we mapped to that step.
"""

from __future__ import annotations

from pathlib import Path

import dearpygui.dearpygui as dpg

from msgflex.gui.log_pane import LogPane
from msgflex.gui.state import AppState
from msgflex.gui.theme import TOKENS, get_stage_color
from msgflex.gui.widgets import build_browse_row
from msgflex.gui.tabs._run_machinery import StageRunner
from msgflex.gui.tabs._stage_helpers import (
    build_action_bar,
    build_stage_header,
)
from msgflex.xpectra.rxflow import rxflow

# Empty step grid — XPECTRA uses incremental mode only.
_runner = StageRunner(
    prefix="xpectra",
    step_ids=[],  # incremental only, no fixed phases
    run_button_tag="xpectra_run_btn",
    stop_button_tag="xpectra_stop_btn",
    badge_tag="xpectra_status_badge",
    incremental_matcher=(
        "xpectra_progress_bar",
        "xpectra_progress_label",
        # Capture sample name from "Processing sample: <name>"
        r"Processing sample:\s+(\S+)",
    ),
)

def _count_mzml_files(mzml_dir: str) -> int:
    """Count *.mzML files. Used to set the incremental total before run."""
    try:
        p = Path(mzml_dir)
        if not p.is_dir():
            return 0
        return len(list(p.glob("*.mzML")))
    except Exception:
        return 0

def _on_run(state: AppState, log: LogPane):
    def handler():
        mzml_dir = dpg.get_value("xpectra_mzml_dir")
        psm_dir = dpg.get_value("xpectra_psm_dir")
        out_dir = dpg.get_value("xpectra_rescored_dir") 
        ensemble = dpg.get_value("xpectra_ensemble")
        folds = dpg.get_value("xpectra_folds")
        log_level = dpg.get_value("xpectra_log_level")

        if not mzml_dir or not Path(mzml_dir).is_dir():
            log.append(f"mzML directory missing or invalid: {mzml_dir}", "error")
            return
        if not psm_dir or not Path(psm_dir).is_dir():
            log.append(f"PSM directory missing or invalid: {psm_dir}", "error")
            return
        if not out_dir:
            log.append("Output directory not set.", "error")
            return

        n_samples = _count_mzml_files(mzml_dir)
        if n_samples == 0:
            log.append(f"No *.mzML files in {mzml_dir}", "error")
            return

        Path(out_dir).mkdir(parents=True, exist_ok=True)  # ← create it if needed

        _runner.set_incremental_total(n_samples)
        log.append(f"Starting XPECTRA rescoring on {n_samples} sample(s)...", "info")
        _runner.start(
            log_pane=log,
            target=lambda: rxflow(
                mzml_dir=mzml_dir,
                msgf_dir=psm_dir,
                out_dir=out_dir,            
                ensemble=ensemble,
                folds=folds,
                log_level=log_level,
            ),
        )

    return handler


def _on_dry_run(state: AppState, log: LogPane):
    def handler():
        ens   = dpg.get_value("xpectra_ensemble")
        folds = dpg.get_value("xpectra_folds")
        mzml  = dpg.get_value("xpectra_mzml_dir")
        psm   = dpg.get_value("xpectra_psm_dir")
        resc  = dpg.get_value("xpectra_rescored_dir")
        n_samples = _count_mzml_files(mzml) if mzml else 0
        log.append("Dry run: would execute XPECTRA with:", "info")
        log.append(
            f"  rxflow(mzml_dir={mzml!r}, msgf_dir={psm!r}, "
            f"out_dir={str(Path(resc).parent) if resc else '?'!r}, "
            f"ensemble={ens}, folds={folds})",
            "debug",
        )
        log.append(f"  Samples to process: {n_samples}", "debug")
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
            title="XPECTRA — ML rescoring",
            description=(
                "Rescores PSMs from upstream using spectral features and an "
                "XGBoost ensemble. Optional — skip to feed SICdir directly "
                "into QUANTIX."
            ),
            stage="xpectra",
            badge_tag="xpectra_status_badge",
            themes=themes,
        )

        # ── Directory pickers ──────────────────────────────────────
        build_browse_row(
            label="mzML directory",
            tag="xpectra_mzml_dir",
            kind="directory",
            help_text="Directory containing *.mzML spectra (defaults to base/data)",
        )
        dpg.add_spacer(height=8)

        build_browse_row(
            label="PSM directory (SICdir)",
            tag="xpectra_psm_dir",
            kind="directory",
            help_text="Directory with *_fht_PlusSICStats.txt from sparx",
        )
        dpg.add_spacer(height=8)

        build_browse_row(
            label="Output directory",
            tag="xpectra_rescored_dir",
            kind="directory",
            help_text="Base output folder. Subdirectories rescored/, features/, "
                    "stats/, and logs/ will be created inside it.",
        )
        dpg.add_spacer(height=12)
        dpg.add_separator()
        dpg.add_spacer(height=10)

        # ── Model parameters ───────────────────────────────────────
        dpg.add_text("Model parameters", color=TOKENS["text"])
        dpg.add_spacer(height=6)
        with dpg.group(horizontal=True):
            with dpg.group():
                dpg.add_text("Ensemble size", color=TOKENS["text_muted"])
                dpg.add_input_int(
                    default_value=5, min_value=1, max_value=20,
                    width=90, tag="xpectra_ensemble",
                    min_clamped=True, max_clamped=True,
                )
            dpg.add_spacer(width=20)
            with dpg.group():
                dpg.add_text("CV folds", color=TOKENS["text_muted"])
                dpg.add_input_int(
                    default_value=5, min_value=2, max_value=10,
                    width=90, tag="xpectra_folds",
                    min_clamped=True, max_clamped=True,
                )
            dpg.add_spacer(width=20)
            with dpg.group():
                dpg.add_text("Log level", color=TOKENS["text_muted"])
                dpg.add_combo(
                    items=["INFO", "DEBUG", "WARNING"],
                    default_value="INFO",
                    width=120, tag="xpectra_log_level",
                )

        dpg.add_spacer(height=14)
        dpg.add_separator()
        dpg.add_spacer(height=10)

        # ── Per-sample progress ────────────────────────────────────
        dpg.add_text("Sample progress", color=TOKENS["text"])
        dpg.add_spacer(height=6)
        dpg.add_progress_bar(
            tag="xpectra_progress_bar",
            default_value=0.0,
            width=-1, height=22,
            overlay="",
        )
        dpg.add_spacer(height=4)
        dpg.add_text(
            "Waiting...",
            tag="xpectra_progress_label",
            color=get_stage_color("xpectra", "600"),
        )
        dpg.add_spacer(height=14)
        dpg.add_separator()
        dpg.add_spacer(height=10)

        build_action_bar(
            stage="xpectra",
            run_callback=_on_run(state, log_pane),
            dry_run_callback=_on_dry_run(state, log_pane),
            stop_callback=_on_stop(log_pane),
            themes=themes,
            run_button_tag="xpectra_run_btn",
            stop_button_tag="xpectra_stop_btn",
        )
