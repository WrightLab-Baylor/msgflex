#!/usr/bin/env python3
"""
XPECTRA model engine initiation
"""

from __future__ import annotations
from pathlib import Path

from msgflex.core import ProjectPaths, XpectraError, get_logger
from msgflex.xpectra.rxflow import rxflow
from msgflex.xpectra.features import (
    extract_spectral_features,
    merge_features_with_psms,
    extract_spectral_features_mzml,
)
from msgflex.xpectra.model import cross_validated_rescoring
from msgflex.xpectra.metrics import (
    compare_gains,
    estimate_fdr,
    estimate_pep_fdr,
)
from msgflex.xpectra.utils import _infer_is_decoy

__all__ = [
    "run",
    "rxflow",
    "extract_spectral_features",
    "merge_features_with_psms",
    "extract_spectral_features_mzml",
    "cross_validated_rescoring",
    "compare_gains",
    "estimate_fdr",
    "estimate_pep_fdr",
]

log = get_logger("xpectra")

def run(
    paths: ProjectPaths,
    *,
    ensemble: int = 5,
    folds: int = 5,
    log_level: str = "INFO",
) -> Path:
    """
    Run XPECTRA rescoring on sparx outputs.
    """
    paths.require_xpectra_inputs()
    paths.ensure_xpectra_dirs()

    log.info(f"XPECTRA: rescoring PSMs in {paths.sicdir}")
    log.info(f"XPECTRA: mzML directory = {paths.data}")
    log.info(f"XPECTRA: output = {paths.xpectra_out}")
    log.info(f"XPECTRA: ensemble={ensemble} folds={folds}")

    try:
        rxflow(
            mzml_dir=str(paths.data),
            msgf_dir=str(paths.sicdir),
            out_dir=str(paths.xpectra_out),
            ensemble=ensemble,
            folds=folds,
            log_level=log_level,
            quiet=False,
        )
    except SystemExit as e:
        code = e.code if isinstance(e.code, int) else 1
        if code != 0:
            raise XpectraError(
                f"rxflow exited with code {code} — see logs for details."
            ) from e
    except Exception as e:
        raise XpectraError(f"rxflow raised: {type(e).__name__}: {e}") from e

    log.info(f"XPECTRA: done -- {paths.xpectra_rescored}")
    return paths.xpectra_rescored
