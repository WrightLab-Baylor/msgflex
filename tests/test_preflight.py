"""
Tests for the preflight validator module.
"""

from __future__ import annotations

from msgflex.core import ProjectPaths
from msgflex.core.preflight import (
    Status,
    check_environment,
    check_project_quantix,
    check_project_sparx,
    check_project_xpectra,
    check_tools,
    run_all_checks,
)

class TestEnvironmentCheck:

    def test_returns_report(self):
        r = check_environment()
        assert r.category == "Environment"
        assert len(r.results) > 0

    def test_core_python_deps_present(self):
        r = check_environment()
        # These must be importable — we depend on them directly
        deps = {res.name for res in r.results if res.status == Status.OK}
        assert "python:typer" in deps
        assert "python:pandas" in deps
        assert "python:numpy" in deps

class TestToolsCheck:

    def test_missing_env_var(self, monkeypatch):
        monkeypatch.delenv("MSGFLEX_TOOLS_DIR", raising=False)
        # With msgflex.tools possibly still resolving, this may pass or fail,
        # but either way it must not crash
        r = check_tools()
        assert r.category == "Tools"

    def test_valid_tools_dir(self, mock_tools_dir):
        r = check_tools()
        # MSGFPlus.jar and all REQUIRED_TOOLS should show OK
        ok_names = {res.name for res in r.results if res.status == Status.OK}
        assert "tools:MSGFPlus.jar" in ok_names
        assert "tools:MASIC/MASIC_Console.exe" in ok_names
        assert "tools:MASICParameters.xml" in ok_names

    def test_missing_tools_flagged(self, tmp_path, monkeypatch):
        # Create empty tools dir — every required tool should FAIL
        tools = tmp_path / "empty_tools"
        tools.mkdir()
        monkeypatch.setenv("MSGFLEX_TOOLS_DIR", str(tools))
        r = check_tools()
        failures = {res.name for res in r.results if res.status == Status.FAIL}
        assert "tools:MSGFPlus.jar" in failures

class TestSparxProjectCheck:

    def test_empty_project_fails(self, empty_project):
        p = ProjectPaths.from_base(empty_project)
        r = check_project_sparx(p)
        assert r.has_failures

    def test_valid_project_passes(self, project_with_data):
        p = ProjectPaths.from_base(project_with_data)
        r = check_project_sparx(p)
        assert not r.has_failures

    def test_inputfile_cross_check_catches_missing_db(self, project_with_data):
        # Rewrite inputfile.tsv to reference a non-existent database
        import pandas as pd
        df = pd.DataFrame({
            "msfilename": ["sample_A"],
            "database": ["nonexistent.fasta"],
            "QCplot": ["sample_A_tic"],
        })
        df.to_csv(project_with_data / "inputfile.tsv", sep="\t", index=False)

        p = ProjectPaths.from_base(project_with_data)
        r = check_project_sparx(p)

        inputfile_result = next(
            res for res in r.results if res.name == "project:inputfile.tsv"
        )
        assert inputfile_result.status in (Status.WARN, Status.FAIL)

class TestXpectraProjectCheck:

    def test_requires_sicdir(self, project_with_data):
        p = ProjectPaths.from_base(project_with_data)
        r = check_project_xpectra(p)
        assert r.has_failures

    def test_passes_after_sparx(self, project_after_sparx):
        p = ProjectPaths.from_base(project_after_sparx)
        r = check_project_xpectra(p)
        assert not r.has_failures

class TestQuantixProjectCheck:

    def test_standard_mode_after_sparx(self, project_after_sparx):
        p = ProjectPaths.from_base(project_after_sparx)
        r = check_project_quantix(p, from_xpectra=False)
        assert not r.has_failures

    def test_xpectra_mode_before_rescore_fails(self, project_after_sparx):
        p = ProjectPaths.from_base(project_after_sparx)
        r = check_project_quantix(p, from_xpectra=True)
        assert r.has_failures

    def test_xpectra_mode_after_rescore_passes(self, project_after_xpectra):
        p = ProjectPaths.from_base(project_after_xpectra)
        r = check_project_quantix(p, from_xpectra=True)
        assert not r.has_failures

class TestRunAllChecks:

    def test_no_paths_skips_project(self):
        """Default stage='all' → environment + tools, no project checks."""
        reports = run_all_checks()
        categories = [r.category for r in reports]
        assert "Environment" in categories
        assert "Tools" in categories  # stage='all' includes upstream, so tools checked
        assert not any("Project" in c for c in categories)

    def test_xpectra_only_skips_tools(self, project_after_sparx):
        """XPECTRA is pure Python — no Java/Mono needed."""
        from msgflex.core import ProjectPaths
        p = ProjectPaths.from_base(project_after_sparx)
        reports = run_all_checks(paths=p, stage="xpectra")
        categories = [r.category for r in reports]
        assert "Environment" in categories
        assert "Tools" not in categories  # skipped

    def test_quantix_only_skips_tools(self, project_after_sparx):
        """QUANTIX is pure Python — no Java/Mono needed."""
        from msgflex.core import ProjectPaths
        p = ProjectPaths.from_base(project_after_sparx)
        reports = run_all_checks(paths=p, stage="quantix")
        categories = [r.category for r in reports]
        assert "Tools" not in categories

    def test_all_stages_with_paths(self, mock_tools_dir, project_after_xpectra):
        p = ProjectPaths.from_base(project_after_xpectra)
        reports = run_all_checks(paths=p, stage="all", from_xpectra=True)
        categories = [r.category for r in reports]
        assert any("sparx" in c.lower() for c in categories)
        assert any("XPECTRA" in c for c in categories)
        assert any("QUANTIX" in c for c in categories)

    def test_single_stage_narrows_scope(self, mock_tools_dir, project_after_sparx):
        p = ProjectPaths.from_base(project_after_sparx)
        reports = run_all_checks(paths=p, stage="sparx")
        categories = [r.category for r in reports]
        # Only upstream project check — no xpectra or quantix
        project_reports = [c for c in categories if "Project" in c]
        assert len(project_reports) == 1
        assert "Sparx" in project_reports[0]
