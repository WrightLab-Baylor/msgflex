"""
Tests for StageRunner step-pattern matching → progress updates.
"""

from __future__ import annotations

import pytest

dpg = pytest.importorskip("dearpygui.dearpygui")


@pytest.fixture
def dpg_context():
    dpg.create_context()
    dpg.create_viewport(title="t", width=800, height=600)
    yield
    dpg.destroy_context()


def _make_runner(prefix, step_ids, matchers):
    """Build a StageRunner with dummy DPG widgets so set_value works."""
    from msgflex.gui.tabs._run_machinery import StageRunner

    with dpg.window(label="t"):
        for sid in step_ids:
            dpg.add_progress_bar(default_value=0.0, tag=f"{prefix}_{sid}_bar")
            dpg.add_text("queued", tag=f"{prefix}_{sid}_status")
        dpg.add_button(label="Run", tag=f"{prefix}_run_btn")
        dpg.add_button(label="Stop", tag=f"{prefix}_stop_btn", enabled=False)
        dpg.add_button(label="IDLE", tag=f"{prefix}_status_badge")

    return StageRunner(
        prefix=prefix,
        step_ids=step_ids,
        run_button_tag=f"{prefix}_run_btn",
        stop_button_tag=f"{prefix}_stop_btn",
        badge_tag=f"{prefix}_status_badge",
        step_matchers=matchers,
    )


class TestStepAdvance:

    def test_first_match_advances_to_first_step(self, dpg_context):
        runner = _make_runner(
            "test", ["a", "b", "c"],
            matchers=[("a", r"start a"), ("b", r"start b"), ("c", r"start c")],
        )
        runner._on_log_message("start a")
        assert runner._current_step_idx == 0
        # First step is running
        assert dpg.get_value("test_a_status") == "running"
        assert dpg.get_value("test_b_status") == "queued"

    def test_advances_through_sequential_matches(self, dpg_context):
        runner = _make_runner(
            "u", ["a", "b", "c"],
            matchers=[("a", r"start a"), ("b", r"start b"), ("c", r"start c")],
        )
        runner._on_log_message("start a")
        runner._on_log_message("start b")
        # Previous step → done, new step → running
        assert dpg.get_value("u_a_status") == "done"
        assert dpg.get_value("u_a_bar") == 1.0
        assert dpg.get_value("u_b_status") == "running"

    def test_unmatched_message_does_not_advance(self, dpg_context):
        runner = _make_runner(
            "v", ["a", "b"],
            matchers=[("a", r"start a"), ("b", r"start b")],
        )
        runner._on_log_message("start a")
        runner._on_log_message("nothing matches this")
        assert runner._current_step_idx == 0  # still on a
        assert dpg.get_value("v_a_status") == "running"

    def test_forward_jump_marks_skipped_steps_done(self, dpg_context):
        """If pipeline checkpoints skip steps, intermediate steps must
        still show as done (no orphaned 'queued' bars before the active
        step)."""
        runner = _make_runner(
            "w", ["a", "b", "c", "d"],
            matchers=[
                ("a", r"do a"), ("b", r"do b"),
                ("c", r"do c"), ("d", r"do d"),
            ],
        )
        runner._on_log_message("do a")
        runner._on_log_message("do d")  # jumps over b and c
        assert dpg.get_value("w_a_status") == "done"
        assert dpg.get_value("w_b_status") == "done"
        assert dpg.get_value("w_c_status") == "done"
        assert dpg.get_value("w_d_status") == "running"

    def test_repeated_match_for_same_step_is_idempotent(self, dpg_context):
        runner = _make_runner(
            "x", ["a", "b"],
            matchers=[("a", r"do a"), ("b", r"do b")],
        )
        runner._on_log_message("do a")
        runner._on_log_message("do a")  # same step matched twice
        assert runner._current_step_idx == 0
        assert dpg.get_value("x_a_status") == "running"


class TestLogBridgeLevelFix:
    """Regression test for the silent-INFO-drop bug.

    Without setup_logging() (which the GUI doesn't call), the msgflex
    logger has level NOTSET, falling back to root's WARNING. INFO records
    from child loggers like msgflex.pipelines.conventional would be
    silently filtered out before reaching our bridge handler.

    _attach_log_bridge must lower the level so step markers actually
    arrive at the GUI.
    """

    def test_attach_bridge_lowers_level_when_unset(self, dpg_context):
        import logging
        from msgflex.gui.log_pane import LogPane
        from msgflex.gui.tabs._run_machinery import (
            _attach_log_bridge,
            _detach_log_bridge,
        )

        # Simulate fresh GUI start — no setup_logging() called yet
        msgflex_root = logging.getLogger("msgflex")
        original_level = msgflex_root.level
        msgflex_root.setLevel(logging.NOTSET)
        try:
            pane = LogPane(parent=0)
            handler = _attach_log_bridge(pane)

            # Bridge should have raised level to INFO
            assert msgflex_root.level == logging.INFO

            # And a child logger emitting INFO should now be captured
            captured = []
            class _Probe(logging.Handler):
                def emit(self, rec):
                    captured.append(rec.getMessage())
            probe = _Probe()
            msgflex_root.addHandler(probe)

            child = logging.getLogger("msgflex.pipelines.conventional")
            child.info("=== Step 0: ===")
            assert "=== Step 0: ===" in captured

            msgflex_root.removeHandler(probe)
            _detach_log_bridge(handler)
        finally:
            msgflex_root.setLevel(original_level)

    def test_attach_bridge_keeps_explicit_debug(self, dpg_context):
        """If user explicitly enabled DEBUG, attach should not raise to INFO."""
        import logging
        from msgflex.gui.log_pane import LogPane
        from msgflex.gui.tabs._run_machinery import (
            _attach_log_bridge,
            _detach_log_bridge,
        )

        msgflex_root = logging.getLogger("msgflex")
        original_level = msgflex_root.level
        msgflex_root.setLevel(logging.DEBUG)
        try:
            pane = LogPane(parent=0)
            handler = _attach_log_bridge(pane)
            # DEBUG was already more verbose than INFO; don't raise it
            assert msgflex_root.level == logging.DEBUG
            _detach_log_bridge(handler)
        finally:
            msgflex_root.setLevel(original_level)


class TestRealPipelinePatterns:
    """Ensure each tab's matchers fire on the actual log strings."""

    def test_upstream_conventional_pattern(self, dpg_context):
        from msgflex.gui.tabs.sparx import _runner
        _runner.reset_progress()
        for line, expected_idx in [
            ("=== Step 0: File conversion + QC plots ===", 0),
            ("=== Step 1: Running MS-GF+ ===",             1),
            ("=== Step 2: mzid → TSV conversion ===",      2),
            ("=== Step 3: PHRP processing ===",            3),
            ("=== Step 4: MASIC SIC generation ===",       4),
            ("=== Step 5: MASIC merger ===",               5),
        ]:
            _runner._on_log_message(line)
            assert _runner._current_step_idx == expected_idx, line

    def test_quantix_pattern(self, dpg_context):
        from msgflex.gui.tabs.quantix import _runner
        _runner.reset_progress()
        for line, expected_idx in [
            ("Step 1: Preparing SICs",                  0),
            ("Step 2: Merging shared peptides",         1),
            ("Step 3: Filtering PSMs",                  2),
            ("Step 4: Generating peptide crosstab",     3),
            ("Step 5: Annotating peptide crosstab",     4),
            ("Step 6: Building peptide-protein maps",   5),
            ("Step 7: Protein rollup (sum)",            6),
        ]:
            _runner._on_log_message(line)
            assert _runner._current_step_idx == expected_idx, line

    def test_xpectra_incremental_pattern(self, dpg_context):
        """XPECTRA uses per-sample incremental mode, not a fixed 4-step grid.

        Each "Processing sample: <name>" log line should advance the
        incremental counter and update the bar/label.
        """
        from msgflex.gui.tabs.xpectra import _runner
        # Need to call set_incremental_total + reset before testing
        _runner.set_incremental_total(3)
        _runner.reset_progress()
        assert _runner._inc_count == 0

        # First sample
        _runner._on_log_message("Processing sample: sample_A")
        assert _runner._inc_count == 1

        # Second sample
        _runner._on_log_message("Processing sample: sample_B")
        assert _runner._inc_count == 2

        # An unrelated log line shouldn't advance
        _runner._on_log_message("Loaded 12345 PSMs")
        assert _runner._inc_count == 2

        # Third sample
        _runner._on_log_message("Processing sample: sample_C")
        assert _runner._inc_count == 3

    def test_xpectra_progress_widgets_update(self, dpg_context):
        """Verify the progress bar and label widgets get updated values."""
        # Build the app so the xpectra widgets exist
        import msgflex.gui.app as app_module
        from msgflex.gui.log_pane import LogPane
        from msgflex.gui.state import AppState
        from msgflex.gui.theme import apply_global_theme
        from msgflex.gui.fonts import register_fonts

        register_fonts(base_size=16, mono_size=14)
        state = AppState()
        themes = apply_global_theme()
        app_module._log_pane = LogPane(parent=0)
        app_module.build_app(state, themes, fonts={})

        from msgflex.gui.tabs.xpectra import _runner
        _runner.set_incremental_total(2)
        _runner.reset_progress()

        # Initial state
        assert dpg.get_value("xpectra_progress_bar") == 0.0
        assert "0 / 2" in dpg.get_value("xpectra_progress_label")

        _runner._on_log_message("Processing sample: alpha")
        assert dpg.get_value("xpectra_progress_bar") == 0.5
        label = dpg.get_value("xpectra_progress_label")
        assert "1 / 2" in label and "alpha" in label

        _runner._on_log_message("Processing sample: beta")
        assert dpg.get_value("xpectra_progress_bar") == 1.0
        label = dpg.get_value("xpectra_progress_label")
        assert "2 / 2" in label and "beta" in label

    def test_full_run_pattern(self, dpg_context):
        from msgflex.gui.tabs.full_run import _runner
        _runner.reset_progress()
        for line, expected_idx in [
            ("Stage 1: SPARX (conventional)  heap=4G", 0),
            ("Stage 2: XPECTRA rescoring", 1),
            ("Stage 3: QUANTIX quantification", 2),
        ]:
            _runner._on_log_message(line)
            assert _runner._current_step_idx == expected_idx, line
