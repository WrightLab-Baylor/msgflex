"""Tests for the GUI theme system and app state.
"""

from __future__ import annotations

import pytest


# Skip the whole module if Dear PyGui isn't available
dpg = pytest.importorskip("dearpygui.dearpygui")


@pytest.fixture
def dpg_context():
    """Create and tear down a DPG context for each test."""
    dpg.create_context()
    yield
    dpg.destroy_context()


class TestTheme:

    def test_tokens_defined(self):
        from msgflex.gui.theme import TOKENS
        # Base chrome
        for key in ("surface", "panel", "border", "accent", "text"):
            assert key in TOKENS, f"Missing token: {key}"
        # Each must be a 4-tuple (RGBA)
        for key, value in TOKENS.items():
            assert len(value) == 4, f"{key} is not RGBA"
            for component in value:
                assert 0 <= component <= 255, f"{key} component out of range"

    def test_stage_accents_present(self):
        from msgflex.gui.theme import STAGE_ACCENTS
        for stage in ("sparx", "xpectra", "quantix", "neutral"):
            assert stage in STAGE_ACCENTS
            shades = STAGE_ACCENTS[stage]
            for shade in ("50", "200", "400", "600", "800"):
                assert shade in shades, f"Missing shade {shade} for {stage}"

    def test_get_stage_color(self):
        from msgflex.gui.theme import get_stage_color
        assert get_stage_color("sparx", "600") == (83, 74, 183, 255)
        assert get_stage_color("xpectra", "600") == (29, 158, 117, 255)
        assert get_stage_color("quantix", "600") == (153, 60, 29, 255)

    def test_get_stage_color_unknown_falls_back(self):
        from msgflex.gui.theme import get_stage_color
        result = get_stage_color("bogus", "600")  # type: ignore[arg-type]
        # Falls back to neutral
        assert result == (83, 74, 183, 255)

    def test_apply_global_theme(self, dpg_context):
        from msgflex.gui.theme import apply_global_theme
        themes = apply_global_theme()

        required_keys = {
            "base", "sparx", "xpectra", "quantix", "neutral",
            "danger", "muted",
            "badge_sparx", "badge_xpectra", "badge_quantix", "badge_neutral",
        }
        assert required_keys.issubset(themes.keys())

        # Every theme ID should exist as a DPG item
        for key, theme_id in themes.items():
            assert dpg.does_item_exist(theme_id), f"Theme {key} invalid"


class TestState:

    def test_empty_state(self):
        from msgflex.gui.state import AppState
        state = AppState()
        assert not state.has_project
        assert not state.has_params
        assert state.java_mem == "4G"

    def test_set_base_populates_paths(self, tmp_path):
        from msgflex.gui.state import AppState
        state = AppState()
        state.set_base(tmp_path)
        assert state.has_project
        assert state.paths.base == tmp_path.resolve()

    def test_subscribers_called(self, tmp_path):
        from msgflex.gui.state import AppState
        state = AppState()
        calls = []
        state.subscribe(lambda s: calls.append(s.paths))
        state.set_base(tmp_path)
        assert len(calls) == 1
        assert calls[0] is not None

    def test_broken_subscriber_doesnt_crash(self, tmp_path):
        from msgflex.gui.state import AppState
        state = AppState()
        state.subscribe(lambda s: (_ for _ in ()).throw(RuntimeError("boom")))
        # Should not raise — broken subscribers are logged, not propagated
        state.set_base(tmp_path)

    def test_java_mem_update(self):
        from msgflex.gui.state import AppState
        state = AppState()
        state.set_java_mem("16G")
        assert state.java_mem == "16G"


class TestFullBuild:
    """Ensure the full GUI tree can be assembled."""

    def test_build_app(self, dpg_context):
        from msgflex.gui.app import build_app
        from msgflex.gui.log_pane import LogPane
        from msgflex.gui.state import AppState
        from msgflex.gui.theme import apply_global_theme
        from msgflex.gui.fonts import register_fonts
        import msgflex.gui.app as app_module

        # create_viewport is required before building windows
        dpg.create_viewport(title="test", width=1240, height=720)

        fonts = register_fonts(base_size=16, mono_size=14)
        state = AppState()
        themes = apply_global_theme()
        app_module._log_pane = LogPane(parent=0)

        build_app(state, themes, fonts=fonts)

        # All six tabs exist
        assert dpg.does_item_exist("main_tabs")
        tabs = dpg.get_item_children("main_tabs", 1)
        assert len(tabs) == 6

        # Header + status bar items exist
        assert dpg.does_item_exist("header_version")
        assert dpg.does_item_exist("header_project_path")
        assert dpg.does_item_exist("statusbar_text")

        # New browse inputs exist
        assert dpg.does_item_exist("project_base")
        assert dpg.does_item_exist("project_params")
        assert dpg.does_item_exist("project_inputfile")
        assert dpg.does_item_exist("xpectra_mzml_dir")
        assert dpg.does_item_exist("xpectra_psm_dir")
        assert dpg.does_item_exist("xpectra_rescored_dir")
        assert dpg.does_item_exist("quantix_input_dir")
        assert dpg.does_item_exist("quantix_database_dir")
        assert dpg.does_item_exist("quantix_output_tsv")

        # Quit button lives in the header, always visible
        assert dpg.does_item_exist("header_quit_btn")

    def test_state_updates_reflected_in_header(self, dpg_context, tmp_path):
        from msgflex.gui.app import build_app
        from msgflex.gui.log_pane import LogPane
        from msgflex.gui.state import AppState
        from msgflex.gui.theme import apply_global_theme
        from msgflex.gui.fonts import register_fonts
        import msgflex.gui.app as app_module

        dpg.create_viewport(title="test", width=1240, height=720)
        fonts = register_fonts(base_size=16, mono_size=14)
        state = AppState()
        themes = apply_global_theme()
        app_module._log_pane = LogPane(parent=0)
        build_app(state, themes, fonts=fonts)

        assert dpg.get_value("header_project_path") == "(no project)"

        state.set_base(tmp_path)
        assert dpg.get_value("header_project_path") == str(tmp_path.resolve())

    def test_project_base_autofills_stage_browse_boxes(self, dpg_context, tmp_path):
        """Selecting a project base should populate empty XPECTRA/QUANTIX paths."""
        from msgflex.gui.app import build_app
        from msgflex.gui.log_pane import LogPane
        from msgflex.gui.state import AppState
        from msgflex.gui.theme import apply_global_theme
        from msgflex.gui.fonts import register_fonts
        from msgflex.gui.tabs.project import _on_base_change
        import msgflex.gui.app as app_module

        dpg.create_viewport(title="test", width=1240, height=720)
        fonts = register_fonts(base_size=16, mono_size=14)
        state = AppState()
        themes = apply_global_theme()
        log = LogPane(parent=0)
        app_module._log_pane = log
        build_app(state, themes, fonts=fonts)

        # Simulate the base directory being picked
        handler = _on_base_change(state, log)
        handler(str(tmp_path))

        assert str(tmp_path) in dpg.get_value("project_inputfile")
        assert "data" in dpg.get_value("xpectra_mzml_dir")
        assert "SICdir" in dpg.get_value("xpectra_psm_dir")
        assert "rescored" in dpg.get_value("xpectra_rescored_dir")
        assert "SICdir" in dpg.get_value("quantix_input_dir")
        assert "database" in dpg.get_value("quantix_database_dir")
        assert "final.tsv" in dpg.get_value("quantix_output_tsv")


class TestFonts:
    """Verify font registration doesn't crash and yields a dict."""

    def test_register_fonts(self, dpg_context):
        from msgflex.gui.fonts import register_fonts
        fonts = register_fonts(base_size=16, mono_size=14)
        assert isinstance(fonts, dict)
        assert set(fonts.keys()) == {"sans", "sans_bold", "mono"}

    def test_register_fonts_hidpi_scale(self, dpg_context):
        from msgflex.gui.fonts import register_fonts
        fonts = register_fonts(base_size=16, mono_size=14, global_scale=1.5)
        # Doesn't crash; both fonts still registered
        assert isinstance(fonts, dict)


class TestCliIntegration:

    def test_gui_subcommand_exists(self):
        from typer.testing import CliRunner
        from msgflex.cli.main import app

        runner = CliRunner()
        result = runner.invoke(app, ["gui", "--help"])
        assert result.exit_code == 0
        assert "graphical user interface" in result.stdout.lower()
