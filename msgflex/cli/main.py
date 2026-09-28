"""
MSGFLEX command-line interface.
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional

import typer

from msgflex import __version__
from msgflex.core import (
    MsgflexError,
    ProjectPaths,
    get_logger,
    setup_logging,
)

app = typer.Typer(
    name="msgflex",
    help="Unified MSGFLEX proteomics pipeline (SPARX + XPECTRA + QUANTIX).",
    no_args_is_help=True,
    add_completion=True,
    context_settings={"help_option_names": ["-h", "--help"]},
    pretty_exceptions_show_locals=False,
)

#  Shared options
def _base_opt() -> typer.Option:
    return typer.Option(
        ...,
        "-b", "--base",
        help="Project base directory (contains data/, database/, inputfile.tsv).",
        exists=True, file_okay=False, dir_okay=True, resolve_path=True,
    )


def _optional_base_opt() -> typer.Option:
    """Like _base_opt but optional — for stages that can run standalone.

    XPECTRA and QUANTIX can operate with explicit directory arguments
    instead of a project base. If neither is given, the command errors
    with a helpful message.
    """
    return typer.Option(
        None,
        "-b", "--base",
        help="Project base directory (optional). If omitted, provide "
             "explicit directories via the stage-specific flags.",
        exists=True, file_okay=False, dir_okay=True, resolve_path=True,
    )

def _level_opt() -> typer.Option:
    return typer.Option(
        "INFO", "--log-level",
        help="Logging verbosity (advanced; prefer --quiet / --verbose).",
        case_sensitive=False,
    )

def _quiet_opt() -> typer.Option:
    return typer.Option(
        False, "-q", "--quiet",
        help="Quiet mode — only WARNING and ERROR on the terminal.",
    )


def _verbose_opt() -> typer.Option:
    return typer.Option(
        False, "-v", "--verbose",
        help="Verbose mode — include DEBUG output on the terminal.",
    )

def _resolve_verbosity(quiet: bool, verbose: bool, log_level: str) -> str:
    """Reconcile the three verbosity flags.

    Mutually-exclusive priority: --quiet beats --verbose beats --log-level.
    Returns one of "quiet", "normal", "verbose".
    """
    if quiet and verbose:
        raise typer.BadParameter("Cannot pass both --quiet and --verbose.")
    if quiet:
        return "quiet"
    if verbose:
        return "verbose"
    lv = log_level.upper()
    if lv == "DEBUG":
        return "verbose"
    if lv in ("WARNING", "ERROR"):
        return "quiet"
    return "normal"

# ── Top-level callback: version + global flags ────────────────────────
def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"msgflex {__version__}")
        raise typer.Exit()


@app.callback()
def _main(
    version: bool = typer.Option(
        None, "--version", "-V",
        callback=_version_callback, is_eager=True,
        help="Show version and exit.",
    ),
) -> None:
    """MSGFLEX — unified proteomics pipeline."""
    pass

# ── sparx ──────────────────────────────────────────────────────────
@app.command()
def sparx(
    base: Path = _base_opt(),
    config: Path = typer.Option(
        ..., "-c", "--config",
        help="MS-GF+ parameters file (MSGFPlus_Params.txt).",
        exists=True, dir_okay=False, resolve_path=True,
    ),
    mode: str = typer.Option(
        "conventional", "--mode",
        help="Search mode: `conventional` or `binning`.",
        case_sensitive=False,
    ),
    java_mem: str = typer.Option(
        "4G", "--java-mem",
        help="Java heap per MS-GF+ process (e.g. 4G, 8G, 16G).",
    ),
    quiet: bool = _quiet_opt(),
    verbose: bool = _verbose_opt(),
    log_level: str = _level_opt(),
) -> None:
    
    """
    Run the MS-GF+ upstream pipeline (search → MASIC → PHRP → merger).
    """
    setup_logging(
        mode="cli",
        verbosity=_resolve_verbosity(quiet, verbose, log_level),
        logfile=base / "msgflex.log",
    )
    log = get_logger("cli.sparx")

    paths = ProjectPaths.from_base(base)
    log.info("Base:  %s", paths.base)
    log.info("Mode:  %s", mode.lower())
    log.info("Heap:  %s", java_mem)
    log.info("Params: %s", config)

    mode_norm = mode.strip().lower()
    if mode_norm == "conventional":
        from msgflex.pipelines import conventional as pipeline
    elif mode_norm == "binning":
        from msgflex.pipelines import binning as pipeline
    else:
        raise typer.BadParameter(
            f"Unknown mode: {mode!r}. Expected 'conventional' or 'binning'."
        )

    try:
        pipeline.run(paths=paths, conf_loc=str(config), java_mem=java_mem)
    except MsgflexError as e:
        log.exception("SPARX failed")
        raise typer.Exit(code=2)

# ── xpectra ───────────────────────────────────────────────────────────
@app.command()
def xpectra(
    base: Optional[Path] = _optional_base_opt(),
    mzml_dir: Optional[Path] = typer.Option(
        None, "--mzml-dir",
        help="Directory containing *.mzML files. "
             "Defaults to <base>/data when --base is given.",
        exists=True, file_okay=False, dir_okay=True, resolve_path=True,
    ),
    psm_dir: Optional[Path] = typer.Option(
        None, "--psm-dir",
        help="Directory containing *_fht_PlusSICStats.txt files. "
             "Defaults to <base>/SICdir when --base is given.",
        exists=True, file_okay=False, dir_okay=True, resolve_path=True,
    ),
    out_dir: Optional[Path] = typer.Option(
        None, "--out-dir",
        help="XPECTRA output directory (will contain features/, rescored/, "
             "stats/, logs/). Defaults to <base>/xpectra when --base is given.",
        resolve_path=True,
    ),
    ensemble: int = typer.Option(
        5, "--ensemble",
        help="Number of models in the XPECTRA ensemble.",
        min=1,
    ),
    folds: int = typer.Option(
        5, "--folds",
        help="Cross-validation folds for rescoring.",
        min=2,
    ),
    log_level: str = _level_opt(),
) -> None:
    """Run XPECTRA ML-based PSM rescoring (optional stage).

    Either pass --base to use the standard project layout, or pass
    --mzml-dir + --psm-dir + --out-dir for a standalone run.
    """
    # Validate input combination
    if base is None and (mzml_dir is None or psm_dir is None or out_dir is None):
        raise typer.BadParameter(
            "Provide either --base for a standard project, OR all three "
            "of --mzml-dir, --psm-dir, --out-dir for a standalone run."
        )

    log_dir = base if base is not None else (out_dir or Path.cwd())
    setup_logging(level=log_level, logfile=log_dir / "msgflex.log")
    log = get_logger("cli.xpectra")

    if base is not None:
        # Standard project mode
        paths = ProjectPaths.from_base(base)
        paths.require_xpectra_inputs()
        paths.ensure_xpectra_dirs()

        from msgflex import xpectra as xpectra_stage
        try:
            xpectra_stage.run(
                paths=paths,
                ensemble=ensemble,
                folds=folds,
                log_level=log_level,
            )
        except MsgflexError as e:
            log.exception("XPECTRA rescoring failed")
            raise typer.Exit(code=2)
    else:
        # Standalone mode — call rxflow directly
        log.info("Standalone mode: mzml=%s, psm=%s, out=%s",
                 mzml_dir, psm_dir, out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            from msgflex.xpectra.rxflow import rxflow
            rxflow(
                mzml_dir=str(mzml_dir),
                msgf_dir=str(psm_dir),
                out_dir=str(out_dir),
                ensemble=ensemble,
                folds=folds,
                log_level=log_level,
            )
        except Exception as e:
            log.exception("XPECTRA rescoring failed")
            raise typer.Exit(code=2)

# ── quantix ───────────────────────────────────────────────────────────
@app.command()
def quantix(
    base: Optional[Path] = _optional_base_opt(),
    input_dir: Optional[Path] = typer.Option(
        None, "--input-dir",
        help="Directory containing PSM files. Defaults to <base>/SICdir, "
             "or <base>/xpectra when --from-xpectra is set.",
        exists=True, file_okay=False, dir_okay=True, resolve_path=True,
    ),
    database_dir: Optional[Path] = typer.Option(
        None, "--database-dir",
        help="Directory with FASTA files. Defaults to <base>/database.",
        exists=True, file_okay=False, dir_okay=True, resolve_path=True,
    ),
    output: Path = typer.Option(
        None, "-o", "--output",
        help="Final LFQ TSV path. Default: <base>/final.tsv (or "
             "<input_dir>/../final.tsv in standalone mode).",
        resolve_path=True,
    ),
    from_xpectra: bool = typer.Option(
        False, "--from-xpectra",
        help="Use XPECTRA's rescored PSMs instead of SICdir. "
             "When set, defaults --score-field=xpec_q and --threshold=0.01.",
    ),
    score_field: Optional[str] = typer.Option(
        None, "--score-field",
        help="Primary PSM score field for filtering. "
             "Default: MSMSScore (standard) or xpec_q (when --from-xpectra). "
             "Accepted: MSMSScore, QValue, xpec_q, xpec_score, rescore_prob.",
    ),
    threshold: Optional[float] = typer.Option(
        None, "--threshold",
        help="Threshold for the primary score field. "
             "Default: 10.0 (MSMSScore/xpec_score) or 0.01 (xpec_q/QValue).",
    ),
    score_field2: Optional[str] = typer.Option(
        None, "--score-field2",
        help="Optional secondary score field for filtering.",
    ),
    threshold2: Optional[float] = typer.Option(
        None, "--threshold2",
        help="Optional secondary threshold.",
    ),
    num_pep: int = typer.Option(
        2, "--num-pep",
        help="Minimum peptides per protein for rollup.",
        min=1,
    ),
    rollup_mode: str = typer.Option(
        "all_matches", "--mode",
        help="Peptide-protein matching mode.",
    ),
    rollup: str = typer.Option(
        "sum", "--rollup",
        help="Protein rollup method.",
    ),
    outlier_alpha: Optional[float] = typer.Option(
        None, "--outlier-alpha",
        help="Grubbs alpha for outlier detection (rrollup only).",
    ),
    coverage_tsv: Optional[Path] = typer.Option(
        None, "--coverage-tsv",
        help="Optional grouped_coverage.tsv.",
        resolve_path=True,
    ),
    quiet: bool = _quiet_opt(),
    verbose: bool = _verbose_opt(),
    log_level: str = _level_opt(),
) -> None:
    """Run QUANTIX downstream peptide/protein LFQ quantification.

    Either pass --base for the standard project layout, or pass
    --input-dir + --database-dir for a standalone run.
    """
    if base is None and (input_dir is None or database_dir is None):
        raise typer.BadParameter(
            "Provide either --base for a standard project, OR both "
            "--input-dir and --database-dir for a standalone run."
        )

    log_dir = base if base is not None else input_dir.parent
    setup_logging(
        mode="cli",
        verbosity=_resolve_verbosity(quiet, verbose, log_level),
        logfile=log_dir / "msgflex.log",
    )
    log = get_logger("cli.quantix")

    # Auto-select sensible defaults based on input mode
    if score_field is None:
        score_field = "xpec_q" if from_xpectra else "MSMSScore"
        log.info("Defaulting --score-field to %s (from_xpectra=%s)", score_field, from_xpectra)

    if threshold is None:
        _q_like = {"xpec_q", "xpec_pepq", "QValue", "PepQValue"}
        threshold = 0.01 if score_field in _q_like else 10.0
        log.info("Defaulting --threshold to %s", threshold)

    # Build ProjectPaths — real or synthesized
    if base is not None:
        paths = ProjectPaths.from_base(base)
        paths.require_quantix_inputs(from_xpectra=from_xpectra)
    else:
        # Standalone mode — synthesize paths from input dir's parent
        log.info("Standalone mode: input=%s, database=%s", input_dir, database_dir)
        ad_hoc_base = input_dir.parent
        paths = ProjectPaths.from_base(ad_hoc_base)
        # Override key dirs with the explicit inputs
        overrides = {"database": database_dir}
        if from_xpectra:
            overrides["xpectra_rescored"] = input_dir
        else:
            overrides["sicdir"] = input_dir
        paths = type(paths)(**{**paths.__dict__, **overrides})

    final_output = output if output is not None else (paths.base / "final.tsv")

    from msgflex import quantix as quantix_stage

    try:
        quantix_stage.run(
            paths=paths,
            final_output=final_output,
            from_xpectra=from_xpectra,
            score_field=score_field,
            threshold=threshold,
            score_field2=score_field2,
            threshold2=threshold2,
            num_pep=num_pep,
            mode=rollup_mode,
            rollup=rollup,
            outlier_alpha=outlier_alpha,
            coverage_tsv=coverage_tsv,
        )
    except MsgflexError as e:
        log.exception("QUANTIX quantification failed")
        raise typer.Exit(code=2)

# ── run (the whole chain) ─────────────────────────────────────────────
@app.command()
def run(
    base: Path = _base_opt(),
    config: Path = typer.Option(
        ..., "-c", "--config",
        help="MS-GF+ parameters file.",
        exists=True, dir_okay=False, resolve_path=True,
    ),
    mode: str = typer.Option(
        "conventional", "--mode",
        help="Upstream search mode.",
    ),
    java_mem: str = typer.Option(
        "4G", "--java-mem",
        help="Java heap per MS-GF+ process (e.g. 4G, 8G, 16G).",
    ),
    rescore: bool = typer.Option(
        False, "--rescore",
        help="Insert XPECTRA rescoring between SPARX and QUANTIX.",
    ),
    skip_sparx: bool = typer.Option(
        False, "--skip-sparx",
        help="Reuse existing SICdir, don't re-run sparx. Useful for re-runs.",
    ),
    output: Path = typer.Option(
        None, "-o", "--output",
        help="Final LFQ TSV path. Default: <base>/final.tsv",
        resolve_path=True,
    ),
    # Pass-through args for downstream stages (most common knobs only;
    # for the full knob set use the per-stage subcommands).
    ensemble: int = typer.Option(5, "--ensemble", min=1),
    folds: int = typer.Option(5, "--folds", min=2),
    score_field: Optional[str] = typer.Option(
        None, "--score-field",
        help="PSM score field for QUANTIX filtering. "
             "Default: MSMSScore (standard) or xpec_q (with --rescore).",
    ),
    threshold: Optional[float] = typer.Option(
        None, "--threshold",
        help="Threshold for the score field. "
             "Default: 10.0 (MSMSScore) or 0.01 (xpec_q).",
    ),
    num_pep: int = typer.Option(2, "--num-pep", min=1),
    quiet: bool = _quiet_opt(),
    verbose: bool = _verbose_opt(),
    log_level: str = _level_opt(),
) -> None:
    """Run the full pipeline: sparx → [xpectra] → quantix."""
    setup_logging(
        mode="cli",
        verbosity=_resolve_verbosity(quiet, verbose, log_level),
        logfile=base / "msgflex.log",
    )
    log = get_logger("cli.run")

    # Auto-select sensible defaults for QUANTIX based on rescore flag
    if score_field is None:
        score_field = "xpec_q" if rescore else "MSMSScore"
        log.info("Defaulting --score-field to %s (rescore=%s)", score_field, rescore)
    if threshold is None:
        _q_like = {"xpec_q", "xpec_pepq", "QValue", "PepQValue"}
        threshold = 0.01 if score_field in _q_like else 10.0
        log.info("Defaulting --threshold to %s", threshold)

    paths = ProjectPaths.from_base(base)
    final_output = output if output is not None else (paths.base / "final.tsv")

    from msgflex.orchestrator import workflow

    try:
        workflow.run_full(
            paths=paths,
            conf_loc=str(config),
            mode=mode,
            java_mem=java_mem,
            rescore=rescore,
            skip_sparx=skip_sparx,
            final_output=final_output,
            ensemble=ensemble,
            folds=folds,
            score_field=score_field,
            threshold=threshold,
            num_pep=num_pep,
            log_level=log_level,
        )
    except MsgflexError as e:
        log.error("Pipeline failed: %s", e)
        raise typer.Exit(code=2)
    except KeyboardInterrupt:
        log.warning("Interrupted by user.")
        raise typer.Exit(code=130)


# ── check (preflight validator) ───────────────────────────────────────

# Status → icon/color helpers (lazy — only used inside check)
_STATUS_GLYPH = {
    "OK":   "✓",
    "WARN": "!",
    "FAIL": "✗",
    "SKIP": "-",
}
_STATUS_COLOR = {
    "OK":   typer.colors.GREEN,
    "WARN": typer.colors.YELLOW,
    "FAIL": typer.colors.RED,
    "SKIP": typer.colors.BRIGHT_BLACK,
}

def _render_report(report) -> tuple[int, int, int]:
    """Print one PreflightReport. Returns (n_ok, n_warn, n_fail)."""
    from msgflex.core.preflight import Status

    typer.echo()
    typer.secho(f"── {report.category} ──", bold=True)

    n_ok = n_warn = n_fail = 0
    for r in report.results:
        glyph = _STATUS_GLYPH[r.status.value]
        color = _STATUS_COLOR[r.status.value]
        typer.secho(f"  {glyph} ", fg=color, nl=False, bold=True)
        typer.secho(f"{r.name:<40} ", nl=False)
        typer.echo(r.message)
        for detail in r.details:
            typer.secho(f"      {detail}", fg=typer.colors.BRIGHT_BLACK)

        if r.status == Status.OK:
            n_ok += 1
        elif r.status == Status.WARN:
            n_warn += 1
        elif r.status == Status.FAIL:
            n_fail += 1

    return n_ok, n_warn, n_fail

@app.command()
def check(
    base: Optional[Path] = typer.Option(
        None, "-b", "--base",
        help="Project base directory. Omit to check environment + tools only.",
        exists=True, file_okay=False, dir_okay=True, resolve_path=True,
    ),
    stage: str = typer.Option(
        "all", "--stage",
        help="Which stage(s) to validate: all | sparx | xpectra | quantix.",
        case_sensitive=False,
    ),
    from_xpectra: bool = typer.Option(
        False, "--from-xpectra",
        help="For quantix checks: assume input is XPECTRA's rescored output.",
    ),
) -> None:
    """
    Run preflight checks — validates environment, tools, and project layout.
    """
    from msgflex.core.preflight import run_all_checks

    stage = stage.lower().strip()
    if stage not in ("all", "sparx", "xpectra", "quantix"):
        raise typer.BadParameter(
            f"Unknown stage: {stage!r}. Expected all|sparx|xpectra|quantix."
        )

    paths = ProjectPaths.from_base(base) if base is not None else None

    typer.secho(
        f"\nMSGFLEX preflight (msgflex {__version__})",
        bold=True, fg=typer.colors.CYAN,
    )
    if paths:
        typer.echo(f"Project: {paths.base}")
        typer.echo(f"Stage:   {stage}")
    else:
        typer.echo("Mode:    environment + tools only (no project supplied)")

    reports = run_all_checks(paths=paths, stage=stage, from_xpectra=from_xpectra)

    total_ok = total_warn = total_fail = 0
    for report in reports:
        n_ok, n_warn, n_fail = _render_report(report)
        total_ok += n_ok
        total_warn += n_warn
        total_fail += n_fail

    # Final summary
    typer.echo()
    typer.secho("── Summary ──", bold=True)
    typer.secho(f"  ✓ {total_ok} passed", fg=typer.colors.GREEN)
    if total_warn:
        typer.secho(f"  ! {total_warn} warnings", fg=typer.colors.YELLOW)
    if total_fail:
        typer.secho(f"  ✗ {total_fail} failed", fg=typer.colors.RED)
        typer.echo()
        typer.secho(
            "Fix the failures above before running the pipeline.",
            fg=typer.colors.RED, bold=True,
        )
        raise typer.Exit(code=1)
    else:
        typer.echo()
        typer.secho("All preflight checks passed.", fg=typer.colors.GREEN, bold=True)

# ── gui ───────────────────────────────────────────────────────────────
@app.command()
def gui() -> None:
    """Launch the MSGFLEX graphical user interface (Dear PyGui)."""
    try:
        from msgflex.gui import run_gui
    except ImportError as e:
        typer.secho(
            "Dear PyGui is not installed. Install it via:\n"
            "    pip install dearpygui\n"
            "or add it to your conda environment.",
            fg=typer.colors.RED,
        )
        typer.secho(f"Import error: {e}", fg=typer.colors.BRIGHT_BLACK)
        raise typer.Exit(code=1)

    run_gui()

# ── Entry point ───────────────────────────────────────────────────────
def main() -> None:
    """Module entry point (used by pyproject's [project.scripts])."""
    app()

def _check_shim() -> None:
    """Entry point for the `msgflex-check` command (alias for `msgflex check`)."""
    import sys as _sys
    _sys.argv = [_sys.argv[0], "check", *_sys.argv[1:]]
    app()

if __name__ == "__main__":
    main()
