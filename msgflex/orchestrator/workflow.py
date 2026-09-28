#!/usr/bin/env python3
"""
Orchestrator: runs the full pipeline.

Chains upstream ==> [XPECTRA] ==> QUANTIX. The top-level `msgflex run`
subcommand delegates here; the per-stage subcommands (`upstream`,
`xpectra`, `quantix`) bypass this module entirely.
"""

from __future__ import annotations

from msgflex.core import (
    MsgflexError,
    ProjectPaths,
    SparxError,
    get_logger,
)

log = get_logger("orchestrator")

def run_full(
    paths: ProjectPaths,
    *,
    conf_loc,
    mode = "conventional",
    java_mem = "4G",
    rescore = False,
    skip_sparx = False,
    final_output = None,
    # xpectra knobs
    ensemble = 5,
    folds = 5,
    # quantix knobs
    score_field = "MSMSScore",
    threshold = 10.0,
    num_pep = 2,
    log_level = "INFO",
):
    """
    Run SPARX ==> [XPECTRA] ==> QUANTIX.
    """
    log.info("=" * 60)
    step = "XPECTRA" if rescore else "(skip)"
    log.info(f"Pipeline: SPARX ({mode}) ==> {step} ==> quantix")
    log.info("=" * 66)

    #  Stage 1: SPARX
    if skip_sparx:
        log.info(f"Stage 1: SKIPPED (--skip-sparx). Reusing {paths.sicdir}")
        paths.require_xpectra_inputs() 
    else:
        log.info(f"Stage 1: sparx {mode}  heap={java_mem}")
        mode_norm = mode.strip().lower()
        if mode_norm == "conventional":
            from msgflex.pipelines import conventional as pipeline
        elif mode_norm == "binning":
            from msgflex.pipelines import binning as pipeline
        else:
            raise SparxError(
                f"Unknown sparx mode: {mode!r}. Expected 'conventional' or 'binning'."
            )
        pipeline.run(paths=paths, conf_loc=conf_loc, java_mem=java_mem)

    # Stage 2: XPECTRA
    if rescore:
        log.info(f"Stage 2: XPECTRA rescoring")
        from msgflex import xpectra as xpectra_stage
        xpectra_stage.run(
            paths=paths,
            ensemble=ensemble,
            folds=folds,
            log_level=log_level,
        )
    else:
        log.info(f"Stage 2: SKIPPED (no --rescore)")

    # Stage 3: QUANTIX
    log.info("Stage 3: QUANTIX quantification")
    from msgflex import quantix as quantix_stage

    out = final_output if final_output is not None else (paths.base / "final.tsv")
    quantix_stage.run(
        paths=paths,
        final_output=out,
        from_xpectra=rescore,
        score_field=score_field,
        threshold=threshold,
        num_pep=num_pep,
    )

    log.info("=" * 60)
    log.info(f"Pipeline done -> final output: {out}")
    log.info("=" * 60)
    return out
