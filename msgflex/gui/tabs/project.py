"""Project tab — sets the base directory, params file, inputfile.tsv, and JVM heap."""

from __future__ import annotations

from pathlib import Path

import dearpygui.dearpygui as dpg

from msgflex.gui.log_pane import LogPane
from msgflex.gui.state import AppState
from msgflex.gui.theme import TOKENS
from msgflex.gui.widgets import build_browse_row, set_browse_value_if_empty


def _on_base_change(state: AppState, log: LogPane):
    def handler(path: str):
        try:
            state.set_base(path)
            log.append(f"Project base set: {path}", "success")
            # Ensure the input widget shows the value (covers typed + browsed paths)
            if dpg.does_item_exist("project_base"):
                dpg.set_value("project_base", path)
            # Auto-fill default paths for stage tabs if they're empty
            base = Path(path)
            set_browse_value_if_empty("project_inputfile", base / "inputfile.tsv")
            set_browse_value_if_empty("xpectra_mzml_dir", base / "data")
            set_browse_value_if_empty("xpectra_psm_dir", base / "SICdir")
            set_browse_value_if_empty(
                "xpectra_rescored_dir", base / "xpectra" # removed rescored
            )
            set_browse_value_if_empty("quantix_input_dir", base / "SICdir")
            set_browse_value_if_empty("quantix_database_dir", base / "database")
            set_browse_value_if_empty("quantix_output_tsv", base / "final.tsv")
        except Exception as e:
            log.append(f"Error setting base: {e}", "error")
    return handler


def _on_params_change(state: AppState, log: LogPane):
    def handler(path: str):
        state.set_params(path)
        log.append(f"Params file set: {path}", "success")
    return handler


def _refresh_summary(state: AppState) -> None:
    if dpg.does_item_exist("project_summary"):
        if state.paths:
            data = state.paths.data
            db = state.paths.database
            raws = len(list(data.glob("*.raw"))) if data.is_dir() else 0
            mzmls = len(list(data.glob("*.mzML"))) if data.is_dir() else 0
            fastas = (
                len(list(db.glob("*.fasta")) + list(db.glob("*.faa")))
                if db.is_dir() else 0
            )
            dpg.set_value(
                "project_summary",
                f"{raws} RAW · {mzmls} mzML · {fastas} FASTA",
            )
        else:
            dpg.set_value("project_summary", "")


def build(
    parent: int | str,
    state: AppState,
    log_pane: LogPane,
    themes: dict,
) -> None:
    with dpg.group(parent=parent):
        dpg.add_text("Project setup", color=TOKENS["text"])
        dpg.add_text(
            "Choose the project directory containing data/, database/, and inputfile.tsv.",
            color=TOKENS["text_muted"],
        )
        dpg.add_spacer(height=14)

        # ── Base directory ─────────────────────────────────────────
        build_browse_row(
            label="Base directory",
            tag="project_base",
            kind="directory",
            on_change=_on_base_change(state, log_pane),
            help_text="Must contain data/, database/, and inputfile.tsv",
        )
        dpg.add_spacer(height=10)

        # ── MSGFPlus params ────────────────────────────────────────
        build_browse_row(
            label="MSGFPlus params file",
            tag="project_params",
            kind="file_open",
            extensions=[".txt", ".*"],
            on_change=_on_params_change(state, log_pane),
            help_text="e.g. MSGFPlus_Params.txt — search tolerances, enzymes, modifications",
        )
        dpg.add_spacer(height=10)

        # ── inputfile.tsv ──────────────────────────────────────────
        build_browse_row(
            label="inputfile.tsv (decoder)",
            tag="project_inputfile",
            kind="file_open",
            extensions=[".tsv", ".txt", ".*"],
            help_text="Auto-detected from base/inputfile.tsv — override only if non-standard",
        )
        dpg.add_spacer(height=14)

        # ── Java heap ──────────────────────────────────────────────
        dpg.add_text("Java heap (MS-GF+ memory)", color=TOKENS["text_muted"])
        dpg.add_combo(
            items=["4G","8G", "16G", "32G", "64G"],
            default_value=state.java_mem,
            width=140,
            callback=lambda s, v: state.set_java_mem(v),
        )
        dpg.add_spacer(height=16)

        dpg.add_separator()
        dpg.add_spacer(height=12)

        # ── Summary readout ────────────────────────────────────────
        dpg.add_text("Project summary", color=TOKENS["text"])
        dpg.add_text("", tag="project_summary", color=TOKENS["text_muted"])

    state.subscribe(lambda s: _refresh_summary(s))
