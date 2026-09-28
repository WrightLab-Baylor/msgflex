"""Tests for the Typer CLI — argument parsing and error handling."""

from __future__ import annotations

from typer.testing import CliRunner

from msgflex import __version__
from msgflex.cli.main import app


runner = CliRunner()


class TestTopLevel:

    def test_version(self):
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert __version__ in result.stdout

    def test_help(self):
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        # All subcommands listed
        for cmd in ("sparx", "xpectra", "quantix", "run", "check"):
            assert cmd in result.stdout

    def test_no_args_shows_help(self):
        result = runner.invoke(app, [])
        assert "Usage" in result.stdout

    def test_unknown_command_fails(self):
        result = runner.invoke(app, ["nonexistent"])
        assert result.exit_code != 0


class TestSparxCommand:

    def test_help(self):
        result = runner.invoke(app, ["sparx", "--help"])
        assert result.exit_code == 0
        assert "--base" in result.stdout
        assert "--config" in result.stdout
        assert "--mode" in result.stdout

    def test_missing_required_args(self):
        result = runner.invoke(app, ["sparx"])
        assert result.exit_code != 0

    def test_invalid_base_fails(self):
        result = runner.invoke(
            app, ["sparx", "-b", "/nonexistent/path/xyz", "-c", "params.txt"],
        )
        assert result.exit_code != 0


class TestXpectraCommand:

    def test_help(self):
        result = runner.invoke(app, ["xpectra", "--help"])
        assert result.exit_code == 0
        assert "--ensemble" in result.stdout
        assert "--folds" in result.stdout

    def test_invalid_ensemble_value(self, project_after_sparx):
        result = runner.invoke(
            app, ["xpectra", "-b", str(project_after_sparx), "--ensemble", "0"],
        )
        assert result.exit_code != 0

    def test_invalid_folds_value(self, project_after_sparx):
        result = runner.invoke(
            app, ["xpectra", "-b", str(project_after_sparx), "--folds", "1"],
        )
        assert result.exit_code != 0


class TestQuantixCommand:

    def test_help(self):
        result = runner.invoke(app, ["quantix", "--help"])
        assert result.exit_code == 0
        assert "--from-xpectra" in result.stdout
        assert "--score-field" in result.stdout
        assert "--threshold" in result.stdout


class TestRunCommand:

    def test_help(self):
        result = runner.invoke(app, ["run", "--help"])
        assert result.exit_code == 0
        assert "--rescore" in result.stdout
        assert "--skip-sparx" in result.stdout


class TestCheckCommand:
    """The new preflight command."""

    def test_help(self):
        result = runner.invoke(app, ["check", "--help"])
        assert result.exit_code == 0
        assert "--stage" in result.stdout
        assert "--from-xpectra" in result.stdout

    def test_no_base_runs_env_tools_only(self, mock_tools_dir):
        # Without --base, only environment + tools get checked.
        # This should either pass or fail based on env, but must not crash.
        result = runner.invoke(app, ["check"])
        assert "MSGFLEX preflight" in result.stdout
        assert "Environment" in result.stdout
        assert "Tools" in result.stdout

    def test_invalid_stage(self, project_with_data):
        result = runner.invoke(
            app, ["check", "-b", str(project_with_data), "--stage", "nonsense"],
        )
        assert result.exit_code != 0

    def test_with_base_runs_project_checks(self, mock_tools_dir, project_with_data):
        result = runner.invoke(
            app, ["check", "-b", str(project_with_data), "--stage", "sparx"],
        )
        assert "Project (Sparx)" in result.stdout

    def test_preflight_fails_when_sicdir_missing_for_xpectra(
        self, mock_tools_dir, project_with_data
    ):
        # project_with_data has no SICdir yet → xpectra check must fail
        result = runner.invoke(
            app, ["check", "-b", str(project_with_data), "--stage", "xpectra"],
        )
        assert result.exit_code == 1
        assert "FAIL" in result.stdout or "✗" in result.stdout
