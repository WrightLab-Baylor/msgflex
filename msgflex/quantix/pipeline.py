"""
QUANTIX -- downstream peptide/protein LFQ pipeline.
Orchestrating Pipeline
"""
#!/usr/bin/env python3
from __future__ import annotations

import shutil
import time
from pathlib import Path
from typing import Optional

import pandas as pd

from msgflex.core import ProjectPaths, QuantixError, get_logger

log = get_logger("quantix.pipeline")

class QuantixPipeline:
    """
    Orchestrates the QUANTIX downstream quantification steps.
    """

    def __init__(self, paths: ProjectPaths) -> None:
        self.paths = paths

    # Helpers 
    def _prepare_dir(self, path: Path) -> None:
        if path.exists() and path.is_dir():
            log.info("Cleaning existing directory: %s", path)
            shutil.rmtree(path)
        path.mkdir(parents=True, exist_ok=True)

    # Step 1: Prepare SICs 
    def prepare_sics(self, sic_input_dir: Path, sic_output_dir: Path) -> None:
        log.info("Step 1: Preparing SICs")
        self._prepare_dir(sic_output_dir)
        from msgflex.quantix.prepare_sics import prepare_SICs
        prepare_SICs(input_dir=str(sic_input_dir), output_dir=str(sic_output_dir))

    # Step 2: Merge shared peptides from PHRPOut 
    def merge_shared(
        self, sic_output_dir: Path, phrp_dir: Path
    ) -> None:
        log.info("Step 2: Merging shared peptides from PHRPOut")
        from msgflex.quantix.merge_shared import append_syn_to_fdr

        suffix_map = {
            "_PlusSICStats.tsv": sorted(sic_output_dir.glob("*_PlusSICStats.tsv")),
            "_rescored.tsv": sorted(sic_output_dir.glob("*_rescored.tsv")),
        }
        suffix, sic_files = max(suffix_map.items(), key=lambda kv: len(kv[1]))

        if not sic_files and list(sic_output_dir.glob("*_withsyn.tsv")):
            log.info("Already merged — skipping")
            return

        total_added = 0
        for sic_file in sic_files:
            base = sic_file.name.replace(suffix, "")
            syn_file = phrp_dir / f"{base}_syn.txt"
            if syn_file.exists():
                added, unique_pairs = append_syn_to_fdr(sic_file, syn_file)
                log.info("%s: +%d rows (unique pairs=%d)", sic_file.name, added, unique_pairs)
                total_added += added
            else:
                log.warning("No matching SYN for %s", sic_file.name)

        log.info("Total new rows appended: %d", total_added)

    # Step 3: Filter PSMs
    def filter_psms(
        self,
        input_dir: Path,
        output_dir: Path,
        *,
        score_field: str = "MSMSScore",
        threshold: float = 10.0,
        score_field2: Optional[str] = None,
        threshold2: Optional[float] = None,
    ) -> None:
        log.info("Step 3: Filtering PSMs (score=%s, threshold=%s)", score_field, threshold)
        self._prepare_dir(output_dir)
        
        from msgflex.quantix.filter_psms import process_fdrdir

        process_fdrdir(
            input_directory=str(input_dir),
            output_directory=str(output_dir),
            metric1=score_field,
            threshold1=threshold,
            metric2=score_field2,
            threshold2=threshold2,
        )

    #  Step 4: Peptide crosstab 
    def peptide_crosstab(
        self, input_dir: Path, output_file: Path
    ) -> None:
        log.info("Step 4: Generating peptide crosstab")
        from msgflex.quantix.peptide_crosstab import merge_peptide_ctab

        tsvs = sorted(input_dir.glob("*_filtered.tsv"))
        if not tsvs:
            raise QuantixError(f"No *_filtered.tsv files in {input_dir}")

        ctab = merge_peptide_ctab(tsvs)
        rename_map = {
            c: c[: -len("_rescored")]
            for c in ctab.columns
            if c not in ("Peptide", "PeptideFlanked", "Protein") and c.endswith("_rescored")
        }
        if rename_map:
            log.info("Stripping '_rescored' suffix from %d sample column(s)", len(rename_map))
            ctab = ctab.rename(columns=rename_map)

        output_file.parent.mkdir(parents=True, exist_ok=True)
        ctab.to_csv(output_file, sep="\t", index=False)
        log.info(
            "Peptide crosstab: %d peptides × %d samples → %s",
            len(ctab), len(tsvs), output_file,
            )

    #  Step 5: Annotate 
    def annotate(
        self,
        crosstab_file: Path,
        database_path: Path,
        output_file: Path,
    ) -> None:
        log.info("Step 5: Annotating peptide crosstab")
        from msgflex.quantix.annotator import build_protein_info_map, get_protein_function

        df = pd.read_csv(crosstab_file, sep="\t", dtype=str)
        for c in ("Peptide", "PeptideFlanked", "Protein"):
            if c not in df.columns:
                raise QuantixError(f"Missing required column in crosstab: {c}")
            df[c] = df[c].astype(str).str.strip()

        prot2info = build_protein_info_map(database_path)

        info_series = df["Protein"].apply(
            lambda pid: get_protein_function(str(pid), prot2info)
        )
        df["Function"] = info_series.apply(lambda d: d.get("function", ""))
        df["Gene"] = info_series.apply(lambda d: d.get("gene", ""))

        # Uniqueness flags
        df["_nprot_"] = df.groupby("PeptideFlanked")["Protein"].transform("nunique")
        df["Unique"] = df["_nprot_"] == 1
        df.drop(columns=["_nprot_"], inplace=True)

        df["Shared"] = ~df["Unique"]

        # Coerce sample columns to numeric
        reserved = {
            "Peptide", "PeptideFlanked", "Protein", "Unique", "Shared", "Gene", "Function",
        }
        sample_cols = [c for c in df.columns if c not in reserved]
        for c in sample_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

        # Reorder
        front = ["Peptide", "PeptideFlanked", "Protein"]
        tail = [ "Unique", "Shared", "Gene", "Function"]
        middle = [c for c in df.columns if c not in set(front + tail)]
        df = df[front + middle + tail]

        output_file.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_file, sep="\t", index=False)
        log.info("Annotated crosstab: %d rows → %s", len(df), output_file)

    # Step 6: Peptide-protein coverage map 
    def coverage_map(
        self,
        annotated_file: Path,
        database_dir: Path,
        *,
        num_pep: int = 2,
        mode: str = "all_matches",
        sic_dir: Optional[Path] = None,
    ) -> None:
        log.info("Step 6: Building peptide-protein maps")
        from msgflex.quantix.coverage import (
            ensure_required_columns,
            detect_sample_columns,
            load_fasta_for_samples,
            prepare_subtable_for_sample,
            compute_coverage_for_sample,
        )

        outdir = annotated_file.parent / "map_files"
        outdir.mkdir(parents=True, exist_ok=True)

        df = pd.read_csv(annotated_file, sep="\t", low_memory=False, dtype=str)
        ensure_required_columns(df)
        sample_cols = detect_sample_columns(df)

        if not sample_cols:
            raise QuantixError("No sample columns detected in annotated crosstab.")

        # Locate inputfile.tsv for sample→FASTA mapping
        inputfile_path = None
        if sic_dir is not None:
            candidate = Path(sic_dir).parent / "inputfile.tsv"
            if candidate.exists():
                inputfile_path = str(candidate)
        if inputfile_path is None and self.paths.inputfile.exists():
            inputfile_path = str(self.paths.inputfile)

        fasta_dict, valid_sample_cols = load_fasta_for_samples(
            database_dir, sample_cols, inputfile_path=inputfile_path,
        )

        per_sample_tables = []
        for sample_col in valid_sample_cols:
            sub = prepare_subtable_for_sample(df, sample_col, mode)
            cov_tbl = compute_coverage_for_sample(
                sub, fasta_dict[sample_col], num_pep, sample_col,
            )
            if not cov_tbl.empty:
                map_tbl = cov_tbl.rename(columns={sample_col: "Coverage"})
                map_tbl.to_csv(outdir / f"{sample_col}_mapfile.tsv", sep="\t", index=False)
                per_sample_tables.append(cov_tbl[["Protein", sample_col]].copy())

        if per_sample_tables:
            merged = per_sample_tables[0]
            for tbl in per_sample_tables[1:]:
                merged = merged.merge(tbl, on="Protein", how="outer")
            merged.fillna(0, inplace=True)
            merged.to_csv(outdir / "grouped_coverage.tsv", sep="\t", index=False)
            log.info("Coverage map: %d proteins → %s", len(merged), outdir)

            #sample coverage heatmap generation
            from msgflex.quantix.density_plots import generate_coverage_heatmap
            heatmap_pdf = generate_coverage_heatmap(merged, outdir)
            if heatmap_pdf:
                log.info("Coverage heatmap written to: %s", heatmap_pdf)
        else:
            log.warning("No proteins passed coverage filters.")

    # ── Step 7: Protein rollup ────────────────────────────────────────
    def protein_rollup(
        self,
        annotated_file: Path,
        output_file: Path,
        *,
        rollup: str = "sum",
        mode: str = "all_matches",
        outlier_alpha: Optional[float] = None,
        coverage_tsv: Optional[Path] = None,
    ) -> None:
        log.info("Step 7: Protein rollup (%s)", rollup)
        from msgflex.quantix.rollup import run_rollup

        cov = str(coverage_tsv) if coverage_tsv else None
        run_rollup(
            input_tsv=str(annotated_file),
            output_protein_tsv=str(output_file),
            rollup=rollup,
            mode=mode,
            outlier_alpha=outlier_alpha,
            coverage_tsv=cov,
        )
        log.info("Protein rollup → %s", output_file)

    # ── Full pipeline ─────────────────────────────────────────────────
    def run(
        self,
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
        Run the full QUANTIX pipeline.
        """
        start = time.time()
        paths = self.paths

        if from_xpectra:
            log.info("Mode: XPECTRA rescored input")
            # Rescored files go directly to filtering
            self.merge_shared(paths.xpectra_rescored, paths.phrp_out)
            filter_input = paths.xpectra_rescored
        else:
            log.info("Mode: standard sparx input")

            # Step 1: Prepare SICs
            sic_output_dir = paths.quantix_sics
            self.prepare_sics(paths.sicdir, sic_output_dir)

            # Step 2: Merge shared peptides
            self.merge_shared(sic_output_dir, paths.phrp_out)

            filter_input = sic_output_dir

        # Step 3: Filter PSMs
        filter_output = paths.quantix_fdr_filt
        self.filter_psms(
            filter_input, filter_output,
            score_field=score_field,
            threshold=threshold,
            score_field2=score_field2,
            threshold2=threshold2,
        )

        # Step 4: Peptide crosstab
        self.peptide_crosstab(filter_output, paths.quantix_peptide_ctab)

        # Step 5: Annotate
        self.annotate(
            paths.quantix_peptide_ctab,
            paths.database,
            paths.quantix_annotated,
        )

        # Step 6: Coverage map
        self.coverage_map(
            paths.quantix_annotated,
            paths.database,
            num_pep=num_pep,
            mode=mode,
            sic_dir=paths.sicdir if not from_xpectra else None,
        )

        # Step 7: Protein rollup
        self.protein_rollup(
            paths.quantix_annotated,
            final_output,
            rollup=rollup,
            mode=mode,
            outlier_alpha=outlier_alpha,
            coverage_tsv=coverage_tsv,
        )

        elapsed = time.time() - start
        log.info("QUANTIX pipeline completed in %.1f seconds", elapsed)
        log.info("Final output: %s", final_output)
        return final_output