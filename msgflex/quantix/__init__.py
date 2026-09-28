"""
QUANTIX -- downstream peptide/protein LFQ quantification.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from msgflex.core import ProjectPaths, QuantixError, get_logger

log = get_logger("quantix")


def run(
    paths: ProjectPaths,
    *,
    final_output: Path,
    from_xpectra: bool = False,
    score_field: str = "MSMSScore",
    threshold: float = 10.0,
    score_field2: Optional[str] = None,
    threshold2: Optional[float] = None,
    num_pep: int = 2,
    mode: str = "all_matches",
    rollup: str = "sum",
    outlier_alpha: Optional[float] = None,
    coverage_tsv: Optional[Path] = None,
) -> Path:
    """
    Run the QUANTIX downstream quantification pipeline.
    """
    paths.require_quantix_inputs(from_xpectra=from_xpectra)

    log.info(
        "QUANTIX: from_xpectra=%s  score_field=%s  threshold=%s  output=%s",
        from_xpectra, score_field, threshold, final_output,
    )

    from msgflex.quantix.pipeline import QuantixPipeline

    try:
        pipeline = QuantixPipeline(paths)
        result = pipeline.run(
            final_output=final_output,
            from_xpectra=from_xpectra,
            score_field=score_field,
            threshold=threshold,
            score_field2=score_field2,
            threshold2=threshold2,
            num_pep=num_pep,
            mode=mode,
            rollup=rollup,
            outlier_alpha=outlier_alpha,
            coverage_tsv=coverage_tsv,
        )
    except QuantixError:
        raise
    except Exception as e:
        raise QuantixError(f"QUANTIX pipeline failed: {type(e).__name__}: {e}") from e

    log.info("QUANTIX: done — %s", result)
    return result
