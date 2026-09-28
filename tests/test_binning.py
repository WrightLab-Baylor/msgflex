"""
Tests for binning v0.5 helpers (E-value correction, FDR recomputation)
and the --java-mem threading from CLI → orchestrator → pipeline.run.
"""
from __future__ import annotations

import inspect

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from msgflex.cli.main import app

runner = CliRunner()

# Skip the whole module if user-provided upstream modules are missing
# (masic_merger, dbcurator, MASIC_wrapper). They're optional in CI/sandbox.
binning = pytest.importorskip(
    "msgflex.pipelines.binning",
    reason="msgflex.modules.{masic_merger,dbcurator,MASIC_wrapper} not present",
    exc_type=ImportError,
)
conventional = pytest.importorskip(
    "msgflex.pipelines.conventional",
    reason="msgflex.modules.{masic_merger,dbcurator,...} not present",
    exc_type=ImportError,
)

# =============================================================================
# Binning helper unit tests
# =============================================================================

class TestBinningHelpers:

    def test_count_amino_acids(self, tmp_path):
        count_amino_acids = binning.count_amino_acids
        fasta = tmp_path / "tiny.fasta"
        fasta.write_text(
            ">P1\n"
            "ACDEFG\n"
            ">P2\n"
            "MKVLWA\n"
            "ALLVT\n"
        )
        # 6 + 6 + 5 = 17 (header lines ignored)
        assert count_amino_acids(str(fasta)) == 17

    def test_count_amino_acids_skips_blank_lines(self, tmp_path):
        count_amino_acids = binning.count_amino_acids
        fasta = tmp_path / "blanks.fasta"
        fasta.write_text(">P1\n\nACDEF\n\n>P2\nGHIK\n")
        assert count_amino_acids(str(fasta)) == 9

    def test_get_db_aa_counts_split(self, tmp_path):
        get_db_aa_counts = binning.get_db_aa_counts
        # Three split parts
        for i, seq in enumerate(["ACDE", "FGHIK", "LMNPQR"], 1):
            (tmp_path / f"db_part{i}.fasta").write_text(f">P{i}\n{seq}\n")
        split_dbs = sorted((tmp_path).glob("db_part*.fasta"))
        full, parts = get_db_aa_counts("db.fasta", [str(p) for p in split_dbs], str(tmp_path))
        assert parts == [4, 5, 6]
        assert full == 15

    def test_get_db_aa_counts_no_split(self, tmp_path):
        get_db_aa_counts = binning.get_db_aa_counts
        (tmp_path / "single.fasta").write_text(">P1\nACDEFGHIK\n")
        full, parts = get_db_aa_counts("single.fasta", [], str(tmp_path))
        assert full == 9
        assert parts == [9]


class TestSpecEValueCorrection:

    def test_factor_is_n_full_over_n_part(self):
        correct_specevalue = binning.correct_specevalue
        df = pd.DataFrame({
            "MSGFDB_SpecEValue": [1e-5, 1e-3],
            "EValue":            [1e-2, 1e-1],
        })
        out = correct_specevalue(df, n_full_aa=1000, n_part_aa=250)
        # factor = 4
        assert np.allclose(out["MSGFDB_SpecEValue"], [4e-5, 4e-3])
        assert np.allclose(out["EValue"],            [4e-2, 4e-1])
        assert (out["CorrectionFactor"] == 4.0).all()

    def test_no_evalue_column_is_ok(self):
        correct_specevalue = binning.correct_specevalue
        df = pd.DataFrame({"MSGFDB_SpecEValue": [1e-5]})
        out = correct_specevalue(df, n_full_aa=10, n_part_aa=5)
        assert "EValue" not in out.columns
        assert np.isclose(out["MSGFDB_SpecEValue"].iloc[0], 2e-5)

    def test_correction_is_conservative(self):
        """Corrected value should always be >= split value (factor >= 1)."""
        correct_specevalue = binning.correct_specevalue
        df = pd.DataFrame({"MSGFDB_SpecEValue": [1e-5]})
        out = correct_specevalue(df, n_full_aa=100, n_part_aa=10)
        assert out["MSGFDB_SpecEValue"].iloc[0] > df["MSGFDB_SpecEValue"].iloc[0]


class TestIsDecoyProtein:

    def test_pure_target(self):
        is_decoy_protein = binning.is_decoy_protein
        assert not is_decoy_protein("sp|P12345|GENE")
        assert not is_decoy_protein("P1;P2;P3")

    def test_pure_decoy(self):
        is_decoy_protein = binning.is_decoy_protein
        assert is_decoy_protein("XXX_P12345")
        assert is_decoy_protein("XXX_P1;XXX_P2")

    def test_mixed_is_target(self):
        """If even one protein in the group is a target, PSM is target."""
        is_decoy_protein = binning.is_decoy_protein
        assert not is_decoy_protein("XXX_P1;P2")


class TestQValueComputation:

    def test_monotone_non_decreasing(self):
        """QValue must never decrease as score worsens (after sort)."""
        compute_qvalue = binning.compute_qvalue
        df = pd.DataFrame({
            "MSGFDB_SpecEValue": [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4],
            "_IsDecoy":          [False, False, True,  False, True,  False, True],
        })
        q = compute_qvalue(df)
        # Reorder by score ascending and check monotone
        order = df["MSGFDB_SpecEValue"].argsort()
        q_sorted = q[order]
        assert all(q_sorted[i] <= q_sorted[i + 1] for i in range(len(q_sorted) - 1))

    def test_pure_target_qvalue_zero(self):
        compute_qvalue = binning.compute_qvalue
        df = pd.DataFrame({
            "MSGFDB_SpecEValue": [1e-10, 1e-9, 1e-8],
            "_IsDecoy":          [False, False, False],
        })
        q = compute_qvalue(df)
        assert all(q == 0.0)


class TestPepQValue:

    def test_collapses_to_best_psm_per_peptide(self):
        compute_pep_qvalue = binning.compute_pep_qvalue
        df = pd.DataFrame({
            "Peptide":           ["K.PEPTIDER.K", "K.PEPTIDER.K", "K.OTHER.R"],
            "MSGFDB_SpecEValue": [1e-9, 1e-3, 1e-7],
            "_IsDecoy":          [False, False, False],
        })
        pq = compute_pep_qvalue(df)
        # Same peptide → same PepQValue
        assert pq[0] == pq[1]

# =============================================================================
# CLI threading: --java-mem reaches the pipeline
# =============================================================================
class TestJavaMemThreading:

    def test_sparx_help_advertises_java_mem(self):
        result = runner.invoke(app, ["sparx", "--help"])
        assert result.exit_code == 0
        assert "--java-mem" in result.stdout

    def test_run_help_advertises_java_mem(self):
        result = runner.invoke(app, ["run", "--help"])
        assert result.exit_code == 0
        assert "--java-mem" in result.stdout

    def test_conventional_run_signature_takes_java_mem(self):
        conv_run = conventional.run
        sig = inspect.signature(conv_run)
        assert "java_mem" in sig.parameters
        assert sig.parameters["java_mem"].default == "4G"

    def test_binning_run_signature_takes_java_mem(self):
        bin_run = binning.run
        sig = inspect.signature(bin_run)
        assert "java_mem" in sig.parameters
        assert sig.parameters["java_mem"].default == "4G"

    def test_workflow_run_full_signature_takes_java_mem(self):
        from msgflex.orchestrator.workflow import run_full
        sig = inspect.signature(run_full)
        assert "java_mem" in sig.parameters

    def test_run_msgfplus_signature_takes_java_mem(self):
        run_msgfplus = conventional.run_msgfplus
        sig = inspect.signature(run_msgfplus)
        assert "java_mem" in sig.parameters

# =============================================================================
# Imports
# =============================================================================

class TestBinningImports:

    def test_run_callable(self):
        assert callable(binning.run)

    def test_helpers_exposed(self):
        assert callable(binning.count_amino_acids)
        assert callable(binning.correct_specevalue)
        assert callable(binning.compute_qvalue)
        assert callable(binning.compute_pep_qvalue)
        assert callable(binning.is_decoy_protein)
