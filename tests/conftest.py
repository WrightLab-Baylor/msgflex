"""
Shared pytest fixtures for MSGFLEX tests.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def empty_project(tmp_path) -> Path:
    """A base directory that exists but has no subdirs."""
    return tmp_path


@pytest.fixture
def project_with_data(tmp_path) -> Path:
    """
    Project skeleton with valid data/ and database/ but no upstream outputs yet.

    Layout:
        tmp_path/
          ├── data/
          │   ├── sample_A.raw   (zero-byte)
          │   └── sample_A.mzML  (zero-byte)
          ├── database/
          │   └── ref.fasta      (one entry)
          └── inputfile.tsv
    """
    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "sample_A.raw").touch()
    (tmp_path / "data" / "sample_A.mzML").touch()

    (tmp_path / "database").mkdir()
    fasta = tmp_path / "database" / "ref.fasta"
    fasta.write_text(">sp|P12345|TEST_PROT Test protein OS=Example\nMKVLWAALLVTFLAGCQA\n")

    df = pd.DataFrame({
        "msfilename": ["sample_A"],
        "database": ["ref.fasta"],
        "QCplot": ["sample_A_tic"],
    })
    df.to_csv(tmp_path / "inputfile.tsv", sep="\t", index=False)

    return tmp_path


@pytest.fixture
def project_after_sparx(project_with_data) -> Path:
    """Project state after upstream completed — SICdir and PHRPOut populated."""
    base = project_with_data
    (base / "SICdir").mkdir()
    (base / "SICdir" / "sample_A_fht_PlusSICStats.txt").write_text(
        "Scan\tPeptide\tProtein\tMSGFDB_SpecEValue\tQValue\n"
        "100\tK.PEPTIDE.R\tsp|P12345|TEST_PROT\t1e-10\t0.001\n"
    )
    (base / "results").mkdir()
    (base / "results" / "PHRPOut").mkdir()
    (base / "results" / "PHRPOut" / "sample_A_syn.txt").write_text(
        "Scan\tPeptide\tProtein\n100\tK.PEPTIDE.R\tsp|P12345|TEST_PROT\n"
    )
    (base / "results" / "MasicOut").mkdir()

    return base


@pytest.fixture
def project_after_xpectra(project_after_sparx) -> Path:
    """Project state after XPECTRA rescoring completed."""
    base = project_after_sparx
    (base / "xpectra").mkdir()
    (base / "xpectra" / "features").mkdir()
    (base / "xpectra" / "rescored").mkdir()
    (base / "xpectra" / "stats").mkdir()
    (base / "xpectra" / "logs").mkdir()
    (base / "xpectra" / "rescored" / "sample_A_rescored.tsv").write_text(
        "ScanNum\tPeptide\tProtein\tis_decoy\txpec_q\txpec_score\n"
        "100\tK.PEPTIDE.R\tsp|P12345|TEST_PROT\t0\t0.001\t0.95\n"
    )
    return base


@pytest.fixture
def mock_tools_dir(tmp_path, monkeypatch) -> Path:
    """
    Build a fake tools/ directory with empty files for every expected tool,
    and point MSGFLEX_TOOLS_DIR at it.
    """
    tools = tmp_path / "tools"
    tools.mkdir()

    # Files
    (tools / "MSGFPlus.jar").touch()
    (tools / "MASICParameters.xml").write_text("<MASICParameters/>")

    # Subdirectories
    (tools / "MASIC").mkdir()
    (tools / "MASIC" / "MASIC_Console.exe").touch()

    (tools / "PHRP").mkdir()
    (tools / "PHRP" / "PeptideHitResultsProcRunner.exe").touch()
    (tools / "PHRP" /"MSGFDB_Mods.txt").touch()
    (tools / "PHRP" / "Mass_Correction_Tags.txt").touch()

    (tools / "ThermoRawFileParser").mkdir()
    (tools / "ThermoRawFileParser" / "ThermoRawFileParser.exe").touch()

    (tools / "MzidMerger").mkdir()
    (tools / "MzidMerger" / "net8.0").mkdir()
    (tools / "MzidMerger" / "net8.0" / "MzidMerger.exe").touch()

    monkeypatch.setenv("MSGFLEX_TOOLS_DIR", str(tools))
    return tools
