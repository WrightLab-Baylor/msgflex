"""
Tests for ProjectPaths — the single source of truth for directory layout.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from msgflex.core import ConfigurationError, ProjectPaths


class TestProjectPathsResolution:
    """All 17+ path attributes resolve correctly from a base dir."""

    def test_basic_paths(self, empty_project):
        p = ProjectPaths.from_base(empty_project)
        assert p.base == empty_project.resolve()
        assert p.data == p.base / "data"
        assert p.database == p.base / "database"
        assert p.inputfile == p.base / "inputfile.tsv"

    def test_sparx_output_paths(self, empty_project):
        p = ProjectPaths.from_base(empty_project)
        assert p.qcdir == p.base / "QCdir"
        assert p.sicdir == p.base / "SICdir"
        assert p.results == p.base / "results"
        assert p.phrp_out == p.results / "PHRPOut"
        assert p.masic_out == p.results / "MasicOut"

    def test_xpectra_paths(self, empty_project):
        p = ProjectPaths.from_base(empty_project)
        assert p.xpectra_out == p.base / "xpectra"
        assert p.xpectra_features == p.xpectra_out / "features"
        assert p.xpectra_rescored == p.xpectra_out / "rescored"
        assert p.xpectra_stats == p.xpectra_out / "stats"
        assert p.xpectra_logs == p.xpectra_out / "logs"

    def test_quantix_paths(self, empty_project):
        p = ProjectPaths.from_base(empty_project)
        assert p.quantix_sics == p.base / "SICs"
        assert p.quantix_fdr_filt == p.base / "fdr_filt"
        assert p.quantix_peptide_ctab == p.base / "peptide_crosstab.tsv"

    def test_accepts_path_and_str(self, empty_project):
        p1 = ProjectPaths.from_base(str(empty_project))
        p2 = ProjectPaths.from_base(empty_project)
        assert p1.base == p2.base

    def test_resolves_tilde(self):
        p = ProjectPaths.from_base("~/nonexistent_test_dir")
        assert "~" not in str(p.base)
        assert str(p.base).startswith("/")


class TestProjectPathsImmutability:
    """ProjectPaths is frozen — can't accidentally mutate."""

    def test_frozen(self, empty_project):
        p = ProjectPaths.from_base(empty_project)
        with pytest.raises((AttributeError, Exception)):
            p.base = Path("/other")  # type: ignore[misc]


class TestSparxValidation:

    def test_missing_data_dir_fails(self, empty_project):
        p = ProjectPaths.from_base(empty_project)
        with pytest.raises(ConfigurationError, match="Data directory"):
            p.require_sparx_inputs()

    def test_missing_database_dir_fails(self, empty_project):
        (empty_project / "data").mkdir()
        p = ProjectPaths.from_base(empty_project)
        with pytest.raises(ConfigurationError, match="Database directory"):
            p.require_sparx_inputs()

    def test_valid_sparx_passes(self, project_with_data):
        p = ProjectPaths.from_base(project_with_data)
        p.require_sparx_inputs()  # should not raise


class TestXpectraValidation:

    def test_missing_sicdir_fails(self, project_with_data):
        p = ProjectPaths.from_base(project_with_data)
        with pytest.raises(ConfigurationError, match="SICdir"):
            p.require_xpectra_inputs()

    def test_empty_sicdir_fails(self, project_with_data):
        (project_with_data / "SICdir").mkdir()
        p = ProjectPaths.from_base(project_with_data)
        with pytest.raises(ConfigurationError, match="_fht_PlusSICStats"):
            p.require_xpectra_inputs()

    def test_valid_xpectra_passes(self, project_after_sparx):
        p = ProjectPaths.from_base(project_after_sparx)
        p.require_xpectra_inputs()


class TestQuantixValidation:

    def test_standard_requires_sicdir(self, project_with_data):
        p = ProjectPaths.from_base(project_with_data)
        (project_with_data / "results" / "PHRPOut").mkdir(parents=True)
        with pytest.raises(ConfigurationError, match="SICdir"):
            p.require_quantix_inputs(from_xpectra=False)

    def test_from_xpectra_requires_rescored(self, project_after_sparx):
        p = ProjectPaths.from_base(project_after_sparx)
        with pytest.raises(ConfigurationError, match="rescored"):
            p.require_quantix_inputs(from_xpectra=True)

    def test_valid_standard_passes(self, project_after_sparx):
        p = ProjectPaths.from_base(project_after_sparx)
        p.require_quantix_inputs(from_xpectra=False)

    def test_valid_from_xpectra_passes(self, project_after_xpectra):
        p = ProjectPaths.from_base(project_after_xpectra)
        p.require_quantix_inputs(from_xpectra=True)


class TestEnsureDirs:

    def test_creates_sparx_dirs(self, empty_project):
        p = ProjectPaths.from_base(empty_project)
        p.ensure_dirs()
        assert p.qcdir.is_dir()
        assert p.sicdir.is_dir()
        assert p.phrp_out.is_dir()
        assert p.masic_out.is_dir()

    def test_creates_xpectra_dirs(self, empty_project):
        p = ProjectPaths.from_base(empty_project)
        p.ensure_xpectra_dirs()
        assert p.xpectra_out.is_dir()
        assert p.xpectra_features.is_dir()
        assert p.xpectra_rescored.is_dir()
        assert p.xpectra_stats.is_dir()
        assert p.xpectra_logs.is_dir()

    def test_idempotent(self, empty_project):
        p = ProjectPaths.from_base(empty_project)
        p.ensure_dirs()
        p.ensure_dirs()  # second call must not raise
