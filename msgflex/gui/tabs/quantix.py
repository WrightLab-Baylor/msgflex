"""QUANTIX tab — downstream peptide/protein LFQ."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import dearpygui.dearpygui as dpg

from msgflex.gui.log_pane import LogPane
from msgflex.gui.state import AppState
from msgflex.gui.theme import TOKENS
from msgflex.gui.widgets import build_browse_row
from msgflex.gui.tabs._run_machinery import StageRunner
from msgflex.gui.tabs._stage_helpers import (
    build_action_bar,
    build_progress_grid,
    build_stage_header,
)

QUANTIX_STEPS = [
    ("sics",     "Prepare SICs"),
    ("merge",    "Merge shared"),
    ("filter",   "Filter PSMs"),
    ("crosstab", "Peptide crosstab"),
    ("annot",    "Annotate"),
    ("coverage", "Coverage map"),
    ("rollup",   "Protein rollup"),
]

_Q_LIKE_FIELDS = {"xpec_q", "xpec_pepq", "QValue", "PepQValue"}


_runner = StageRunner(
    prefix="quantix",
    step_ids=[s[0] for s in QUANTIX_STEPS],
    run_button_tag="quantix_run_btn",
    stop_button_tag="quantix_stop_btn",
    badge_tag="quantix_status_badge",
    # Match log lines emitted by msgflex.quantix.pipeline
    step_matchers=[
        ("sics",     r"Step 1: Preparing SICs"),
        ("merge",    r"Step 2: Merging shared"),
        ("filter",   r"Step 3: Filtering PSMs"),
        ("crosstab", r"Step 4: Generating peptide crosstab"),
        ("annot",    r"Step 5: Annotating"),
        ("coverage", r"Step 6: Building peptide-protein"),
        ("rollup",   r"Step 7: Protein rollup"),
    ],
)

def _on_from_xpectra_toggle(state: AppState):
    """Switch defaults + input dir hint when --from-xpectra is flipped."""
    def handler(sender, value):
        if value:
            dpg.set_value("quantix_score_field", "xpec_q")
            dpg.set_value("quantix_threshold", 0.01)
            if state.paths:
                dpg.set_value("quantix_input_dir", str(state.paths.xpectra_rescored))
            if dpg.does_item_exist("quantix_input_dir_hint"):
                dpg.set_value(
                    "quantix_input_dir_hint",
                    "XPECTRA rescored *_rescored.tsv files",
                )
        else:
            dpg.set_value("quantix_score_field", "MSMSScore")
            dpg.set_value("quantix_threshold", 10.0)
            if state.paths:
                dpg.set_value("quantix_input_dir", str(state.paths.sicdir))
            if dpg.does_item_exist("quantix_input_dir_hint"):
                dpg.set_value(
                    "quantix_input_dir_hint",
                    "SICdir *_fht_PlusSICStats.txt files from sparx",
                )
    return handler

def _on_score_field_change(sender, value):
    """Switch threshold default when user picks a q-like vs score-like field."""
    if value in _Q_LIKE_FIELDS:
        dpg.set_value("quantix_threshold", 0.01)
    else:
        dpg.set_value("quantix_threshold", 10.0)

def _on_run(state: AppState, log: LogPane):
    def handler():
        # QUANTIX does NOT require a project base if the user has filled
        # in the input/database/output boxes themselves. We construct an
        # ad-hoc ProjectPaths if needed so the pipeline still gets one.
        from_xpectra = dpg.get_value("quantix_from_xpectra")
        score_field  = dpg.get_value("quantix_score_field")
        threshold    = dpg.get_value("quantix_threshold")
        score_field2_raw = dpg.get_value("quantix_score_field2")
        threshold2_raw   = dpg.get_value("quantix_threshold2")
        mode  = dpg.get_value("quantix_mode")
        num_pep = dpg.get_value("quantix_num_pep")
        rollup  = dpg.get_value("quantix_rollup")
        outlier_alpha_raw = dpg.get_value("quantix_outlier_alpha")
        input_dir = dpg.get_value("quantix_input_dir")
        database_dir = dpg.get_value("quantix_database_dir")
        output_path  = dpg.get_value("quantix_output_tsv")

        # Normalise optional fields — '(none)' / 0.0 → None
        score_field2: Optional[str] = (
            None if not score_field2_raw or score_field2_raw == "(none)"
            else score_field2_raw
        )
        threshold2: Optional[float] = (
            None if score_field2 is None else float(threshold2_raw)
        )
        outlier_alpha: Optional[float] = (
            float(outlier_alpha_raw) if outlier_alpha_raw and outlier_alpha_raw > 0
            else None
        )

        # Validate
        if not input_dir or not Path(input_dir).is_dir():
            log.append(f"Input directory missing or invalid: {input_dir}", "error")
            return
        if not database_dir or not Path(database_dir).is_dir():
            log.append(f"Database directory missing or invalid: {database_dir}", "error")
            return

        # Resolve paths: prefer real project, else synthesise one from input dir
        if state.paths is not None:
            paths = state.paths
        else:
            from msgflex.core import ProjectPaths
            # Use parent of input dir as ad-hoc base — workdirs (SICs/fdr_estd
            # /fdr_filt) will land there, alongside the user's actual inputs.
            ad_hoc_base = Path(input_dir).parent
            paths = ProjectPaths.from_base(ad_hoc_base)
            # Override sicdir/database to the user's selections (in case
            # they don't follow the default <base>/SICdir, <base>/database
            # convention)
            paths = type(paths)(
                **{**paths.__dict__,
                   "sicdir": Path(input_dir) if not from_xpectra else paths.sicdir,
                   "xpectra_rescored": Path(input_dir) if from_xpectra else paths.xpectra_rescored,
                   "database": Path(database_dir),
                   }
            )
            log.append(
                f"No project base set — using ad-hoc base: {ad_hoc_base}",
                "info",
            )

        if not output_path:
            output_path = str(paths.base / "final.tsv")
            log.append(f"No output TSV set — defaulting to: {output_path}", "info")

        log.append(f" input_dir:    {input_dir}",     "debug")
        log.append(f" database:     {database_dir}",  "debug")
        log.append(f" output_tsv:   {output_path}",   "debug")
        log.append(f" from_xpectra: {from_xpectra}",  "debug")
        log.append(f" score_field:  {score_field}",   "debug")
        log.append(f" threshold:    {threshold}",     "debug")
        log.append(f" score_field2: {score_field2}",  "debug")
        log.append(f" threshold2:   {threshold2}",    "debug")
        log.append(f" mode:         {mode}",          "debug")
        log.append(f" num_pep:      {num_pep}",       "debug")
        log.append(f" rollup:       {rollup}",        "debug")
        log.append(f" outlier_alpha:{outlier_alpha}", "debug")

        try:
            from msgflex import quantix as quantix_stage
        except ImportError as e:
            log.append(f"Failed to load QUANTIX: {e}", "error")
            return

        log.append("Starting QUANTIX...", "info")
        _runner.start(
            log_pane=log,
            target=lambda: quantix_stage.run(
                paths=paths,
                final_output=Path(output_path),
                from_xpectra=from_xpectra,
                score_field=score_field,
                threshold=float(threshold),
                score_field2=score_field2,
                threshold2=threshold2,
                mode=mode,
                num_pep=int(num_pep),
                rollup=rollup,
                outlier_alpha=outlier_alpha,
            ),
        )

    return handler

def _on_dry_run(state: AppState, log: LogPane):
    def handler():
        sf = dpg.get_value("quantix_score_field")
        th = dpg.get_value("quantix_threshold")
        np_ = dpg.get_value("quantix_num_pep")
        rollup = dpg.get_value("quantix_rollup")
        out = dpg.get_value("quantix_output_tsv") or "(default)"
        fx = "--from-xpectra " if dpg.get_value("quantix_from_xpectra") else ""
        base = state.paths.base if state.paths else "(no project)"
        log.append("Dry run: would execute `msgflex quantix` with:", "info")
        log.append(
            f"  msgflex quantix -b {base} {fx}"
            f"--score-field {sf} --threshold {th} --num-pep {np_} "
            f"--rollup {rollup} -o {out}",
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
            title="QUANTIX — LFQ quantification",
            description=(
                "Filters PSMs, builds peptide/protein crosstabs, and rolls up "
                "to final LFQ intensities. Auto-switches input and defaults "
                "based on source."
            ),
            stage="quantix",
            badge_tag="quantix_status_badge",
            themes=themes,
        )

        # ── Input source checkbox ──────────────────────────────────
        dpg.add_checkbox(
            label="Use XPECTRA rescored PSMs as input",
            tag="quantix_from_xpectra",
            default_value=False,
            callback=_on_from_xpectra_toggle(state),
        )
        dpg.add_spacer(height=10)

        # ── Directory / file pickers ───────────────────────────────
        build_browse_row(
            label="Input directory",
            tag="quantix_input_dir",
            kind="directory",
            help_text="SICdir *_fht_PlusSICStats.txt files from sparx",
        )
        dpg.add_text(
            "",
            tag="quantix_input_dir_hint",
            color=TOKENS["text_faint"],
            show=False,
        )
        dpg.add_spacer(height=8)

        build_browse_row(
            label="Database directory",
            tag="quantix_database_dir",
            kind="directory",
            help_text="Contains *.fasta / *.faa files for peptide→protein mapping",
        )
        dpg.add_spacer(height=8)

        build_browse_row(
            label="Final output TSV",
            tag="quantix_output_tsv",
            kind="file_save",
            extensions=[".tsv", ".txt", ".*"],
            help_text="Type a new filename in the dialog's File Name box, "
                      "or pick an existing file to overwrite. "
                      "Final LFQ table with peptide intensities + metadata.",
        )
        dpg.add_spacer(height=14)
        dpg.add_separator()
        dpg.add_spacer(height=10)

        # ── Filter settings ────────────────────────────────────────
        dpg.add_text("Filter settings", color=TOKENS["text"])
        dpg.add_spacer(height=6)
        with dpg.group(horizontal=True):
            with dpg.group():
                dpg.add_text("Score field", color=TOKENS["text_muted"])
                dpg.add_combo(
                    items=["MSMSScore", "MSGFScore", "QValue", "PepQValue",
                           "xpec_score", "xpec_q", "xpec_pepq", "rescore_prob"],
                    default_value="MSMSScore",
                    width=160, tag="quantix_score_field",
                    callback=_on_score_field_change,
                )
            dpg.add_spacer(width=16)
            with dpg.group():
                dpg.add_text("Threshold", color=TOKENS["text_muted"])
                dpg.add_input_float(
                    default_value=10.0, width=110,
                    tag="quantix_threshold", step=0,
                )

        dpg.add_spacer(height=10)

        # ── Optional second filter ────────────────────────────────
        with dpg.group(horizontal=True):
            with dpg.group():
                dpg.add_text("Score field 2 (optional)", color=TOKENS["text_muted"])
                dpg.add_combo(
                    items=["(none)", "MSMSScore", "MSGFScore", "QValue", "PepQValue",
                           "xpec_score", "xpec_q", "xpec_pepq", "rescore_prob"],
                    default_value="(none)",
                    width=160, tag="quantix_score_field2",
                )
            dpg.add_spacer(width=16)
            with dpg.group():
                dpg.add_text("Threshold 2 (optional)", color=TOKENS["text_muted"])
                dpg.add_input_float(
                    default_value=0.0, width=110,
                    tag="quantix_threshold2", step=0,
                )
        dpg.add_text(
            "Set 'Score field 2' to apply an additional filter (e.g. PepQValue ≤ 0.01 alongside MSMSScore ≥ 10).",
            color=TOKENS["text_faint"],
        )

        dpg.add_spacer(height=10)

        # ── Matching mode + min peptides + rollup + outlier alpha ─
        with dpg.group(horizontal=True):
            with dpg.group():
                dpg.add_text("Min peptides / protein", color=TOKENS["text_muted"])
                dpg.add_input_int(
                    default_value=2, min_value=1, max_value=20,
                    width=90, tag="quantix_num_pep",
                    min_clamped=True, max_clamped=True,
                )
            dpg.add_spacer(width=16)
            with dpg.group():
                dpg.add_text("Rollup method", color=TOKENS["text_muted"])
                dpg.add_combo(
                    items=["sum", "rrollup", "zrollup", "qrollup", "dwrollup"],
                    default_value="sum",
                    width=140, tag="quantix_rollup",
                )
            dpg.add_spacer(width=16)
            with dpg.group():
                dpg.add_text("Peptide-protein mode", color=TOKENS["text_muted"])
                dpg.add_combo(
                    items=["all_matches", "unique_only", "requires_unique"],
                    default_value="all_matches",
                    width=160, tag="quantix_mode",
                )
            dpg.add_spacer(width=16)
            with dpg.group():
                dpg.add_text("Outlier α (rrollup only)", color=TOKENS["text_muted"])
                dpg.add_input_float(
                    default_value=0.0, width=100,
                    tag="quantix_outlier_alpha", step=0,
                    min_value=0.0, max_value=1.0,
                )

        dpg.add_spacer(height=14)
        dpg.add_separator()
        dpg.add_spacer(height=10)

        build_progress_grid(
            QUANTIX_STEPS,
            progress_prefix="quantix",
            stage="quantix",
            themes=themes,
        )

        build_action_bar(
            stage="quantix",
            run_callback=_on_run(state, log_pane),
            dry_run_callback=_on_dry_run(state, log_pane),
            stop_callback=_on_stop(log_pane),
            themes=themes,
            run_button_tag="quantix_run_btn",
            stop_button_tag="quantix_stop_btn",
        )
