"""
MSGFLEX GUI theme — slate-blue-purple with stage-aware accents.

Exposes:
    TOKENS          : raw color tokens (RGB tuples for Dear PyGui)
    STAGE_ACCENTS   : per-stage color families (SPARX/XPECTRA/QUANTIX)
    build_base_theme() -> int
    build_stage_theme(stage) -> int
    apply_global_theme()
    get_stage_color(stage, shade) -> tuple

Dear PyGui color format: (R, G, B, A) as 0-255 ints.
"""

from __future__ import annotations

from typing import Literal

import dearpygui.dearpygui as dpg


# ── Token definitions (RGB as 0-255 tuples) ──────────────────────────

def _hex(h: str, alpha: int = 255) -> tuple[int, int, int, int]:
    """'#F7F6FA' -> (247, 246, 250, 255)"""
    h = h.lstrip("#")
    return (int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16), alpha)


TOKENS = {
    # Base chrome — slate blue purple
    "surface":        _hex("#F7F6FA"),
    "surface_alt":    _hex("#FFFFFF"),
    "panel":          _hex("#ECEAF2"),
    "border":         _hex("#D9D6E3"),
    "border_strong":  _hex("#BDB8CE"),
    "accent":         _hex("#534AB7"),
    "accent_dark":    _hex("#3C3489"),
    "accent_light":   _hex("#EEEDFE"),
    "text":           _hex("#2A2340"),
    "text_muted":     _hex("#6B6680"),
    "text_faint":     _hex("#9D99AE"),

    # Semantic
    "success":        _hex("#639922"),
    "success_light":  _hex("#EAF3DE"),
    "success_text":   _hex("#3B6D11"),
    "warning":        _hex("#BA7517"),
    "warning_light":  _hex("#FAEEDA"),
    "warning_text":   _hex("#854F0B"),
    "error":          _hex("#A32D2D"),
    "error_light":    _hex("#FCEBEB"),
    "error_text":     _hex("#A32D2D"),
    "info":           _hex("#378ADD"),
    "info_light":     _hex("#E6F1FB"),
    "info_text":      _hex("#0C447C"),
}


# Stage color families — the "surprise" accent system
STAGE_ACCENTS = {
    "sparx": {
        "50":  _hex("#EEEDFE"),
        "200": _hex("#AFA9EC"),
        "400": _hex("#7F77DD"),
        "600": _hex("#534AB7"),
        "800": _hex("#3C3489"),
    },
    "xpectra": {
        "50":  _hex("#E1F5EE"),
        "200": _hex("#9FE1CB"),
        "400": _hex("#5DCAA5"),
        "600": _hex("#1D9E75"),
        "800": _hex("#0F6E56"),
    },
    "quantix": {
        "50":  _hex("#FAECE7"),
        "200": _hex("#F5C4B3"),
        "400": _hex("#D85A30"),
        "600": _hex("#993C1D"),
        "800": _hex("#712B13"),
    },
    # Neutral — Project, Check, Run, Log tabs
    "neutral": {
        "50":  _hex("#ECEAF2"),
        "200": _hex("#BDB8CE"),
        "400": _hex("#6B6680"),
        "600": _hex("#534AB7"),   # falls back to base accent
        "800": _hex("#3C3489"),
    },
}

StageName = Literal["sparx", "xpectra", "quantix", "neutral"]


def get_stage_color(stage: StageName, shade: str = "600") -> tuple[int, int, int, int]:
    return STAGE_ACCENTS.get(stage, STAGE_ACCENTS["neutral"])[shade]


# ── Base theme (applied globally) ─────────────────────────────────────

def build_base_theme() -> int:
    """Build the global theme — slate chrome, applied to every widget.

    Must be called after dpg.create_context().

    Returns:
        Theme ID. Bind with dpg.bind_theme(theme_id).
    """
    with dpg.theme() as theme_id:
        # ── Window / child / popup backgrounds ─────────────────────────
        with dpg.theme_component(dpg.mvAll):
            dpg.add_theme_color(dpg.mvThemeCol_WindowBg,       TOKENS["surface"])
            dpg.add_theme_color(dpg.mvThemeCol_ChildBg,        TOKENS["surface"])
            dpg.add_theme_color(dpg.mvThemeCol_PopupBg,        TOKENS["surface_alt"])
            dpg.add_theme_color(dpg.mvThemeCol_Border,         TOKENS["border"])
            dpg.add_theme_color(dpg.mvThemeCol_BorderShadow,   (0, 0, 0, 0))

            # Text
            dpg.add_theme_color(dpg.mvThemeCol_Text,           TOKENS["text"])
            dpg.add_theme_color(dpg.mvThemeCol_TextDisabled,   TOKENS["text_faint"])

            # Frame bg (inputs, combos, checkboxes, plots)
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg,        TOKENS["surface_alt"])
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, TOKENS["panel"])
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive,  TOKENS["accent_light"])

            # Title bar (for popups/dialogs)
            dpg.add_theme_color(dpg.mvThemeCol_TitleBg,         TOKENS["panel"])
            dpg.add_theme_color(dpg.mvThemeCol_TitleBgActive,   TOKENS["panel"])
            dpg.add_theme_color(dpg.mvThemeCol_TitleBgCollapsed, TOKENS["panel"])

            # Menu bar
            dpg.add_theme_color(dpg.mvThemeCol_MenuBarBg,      TOKENS["panel"])

            # Scrollbars
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarBg,      TOKENS["surface"])
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarGrab,    TOKENS["border_strong"])
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarGrabHovered, TOKENS["text_faint"])
            dpg.add_theme_color(dpg.mvThemeCol_ScrollbarGrabActive,  TOKENS["accent"])

            # Checkmark / slider / progress
            dpg.add_theme_color(dpg.mvThemeCol_CheckMark,         TOKENS["accent"])
            dpg.add_theme_color(dpg.mvThemeCol_SliderGrab,        TOKENS["accent"])
            dpg.add_theme_color(dpg.mvThemeCol_SliderGrabActive,  TOKENS["accent_dark"])

            # Buttons (default style — stage tabs override per-tab)
            dpg.add_theme_color(dpg.mvThemeCol_Button,         TOKENS["surface_alt"])
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered,  TOKENS["accent_light"])
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive,   TOKENS["accent"])

            # Header (collapsing headers, selectables)
            dpg.add_theme_color(dpg.mvThemeCol_Header,         TOKENS["panel"])
            dpg.add_theme_color(dpg.mvThemeCol_HeaderHovered,  TOKENS["accent_light"])
            dpg.add_theme_color(dpg.mvThemeCol_HeaderActive,   TOKENS["accent_light"])

            # Separator
            dpg.add_theme_color(dpg.mvThemeCol_Separator,        TOKENS["border"])
            dpg.add_theme_color(dpg.mvThemeCol_SeparatorHovered, TOKENS["border_strong"])
            dpg.add_theme_color(dpg.mvThemeCol_SeparatorActive,  TOKENS["accent"])

            # Tab strip (our main nav)
            dpg.add_theme_color(dpg.mvThemeCol_Tab,              TOKENS["panel"])
            dpg.add_theme_color(dpg.mvThemeCol_TabHovered,       TOKENS["accent_light"])
            dpg.add_theme_color(dpg.mvThemeCol_TabActive,        TOKENS["surface"])
            dpg.add_theme_color(dpg.mvThemeCol_TabUnfocused,     TOKENS["panel"])
            dpg.add_theme_color(dpg.mvThemeCol_TabUnfocusedActive, TOKENS["surface"])

            # Tables (used on Project tab)
            dpg.add_theme_color(dpg.mvThemeCol_TableHeaderBg,     TOKENS["panel"])
            dpg.add_theme_color(dpg.mvThemeCol_TableBorderStrong, TOKENS["border"])
            dpg.add_theme_color(dpg.mvThemeCol_TableBorderLight,  TOKENS["border"])
            dpg.add_theme_color(dpg.mvThemeCol_TableRowBg,        TOKENS["surface"])
            dpg.add_theme_color(dpg.mvThemeCol_TableRowBgAlt,     TOKENS["panel"])

            # Plot colors (for future per-sample progress plots)
            dpg.add_theme_color(dpg.mvThemeCol_PlotHistogram,        TOKENS["accent"])
            dpg.add_theme_color(dpg.mvThemeCol_PlotHistogramHovered, TOKENS["accent_dark"])

            # ── Rounding, padding, spacing ─────────────────────────────
            dpg.add_theme_style(dpg.mvStyleVar_WindowRounding,     6)
            dpg.add_theme_style(dpg.mvStyleVar_ChildRounding,      6)
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding,      4)
            dpg.add_theme_style(dpg.mvStyleVar_PopupRounding,      6)
            dpg.add_theme_style(dpg.mvStyleVar_GrabRounding,       3)
            dpg.add_theme_style(dpg.mvStyleVar_TabRounding,        4)
            dpg.add_theme_style(dpg.mvStyleVar_ScrollbarRounding,  6)
            dpg.add_theme_style(dpg.mvStyleVar_FrameBorderSize,    1)
            dpg.add_theme_style(dpg.mvStyleVar_WindowBorderSize,   0)
            dpg.add_theme_style(dpg.mvStyleVar_ChildBorderSize,    1)
            dpg.add_theme_style(dpg.mvStyleVar_WindowPadding,      12, 12)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding,       10, 6)
            dpg.add_theme_style(dpg.mvStyleVar_ItemSpacing,        8, 6)
            dpg.add_theme_style(dpg.mvStyleVar_ItemInnerSpacing,   6, 4)
            dpg.add_theme_style(dpg.mvStyleVar_CellPadding,        8, 4)
            dpg.add_theme_style(dpg.mvStyleVar_ScrollbarSize,      12)
            dpg.add_theme_style(dpg.mvStyleVar_GrabMinSize,        12)
            dpg.add_theme_style(dpg.mvStyleVar_TabBorderSize,      0)

    return theme_id


# ── Stage-accent themes ───────────────────────────────────────────────

def build_stage_theme(stage: StageName) -> int:
    """Build an accent theme for a specific stage.

    Bind this to widgets you want to stage-tint — stage title text,
    progress bars, 'selected' buttons, running-status dots.

    Returns:
        Theme ID.
    """
    accent = get_stage_color(stage, "600")
    accent_dark = get_stage_color(stage, "800")
    accent_light = get_stage_color(stage, "50")

    with dpg.theme() as theme_id:
        with dpg.theme_component(dpg.mvAll):
            # Button in this stage's accent
            dpg.add_theme_color(dpg.mvThemeCol_Button,         accent_light)
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered,  accent)
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive,   accent_dark)
            dpg.add_theme_color(dpg.mvThemeCol_Text,           accent_dark)

            # Progress bar fill
            dpg.add_theme_color(dpg.mvThemeCol_PlotHistogram,  accent)
            dpg.add_theme_color(dpg.mvThemeCol_CheckMark,      accent)
            dpg.add_theme_color(dpg.mvThemeCol_SliderGrab,     accent)

            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding,  4)

    return theme_id


# ── Reusable mini-themes ──────────────────────────────────────────────

def build_danger_button_theme() -> int:
    """Red 'Stop' button style."""
    with dpg.theme() as theme_id:
        with dpg.theme_component(dpg.mvButton):
            dpg.add_theme_color(dpg.mvThemeCol_Button,        TOKENS["error"])
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, TOKENS["error_text"])
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive,  TOKENS["error_text"])
            dpg.add_theme_color(dpg.mvThemeCol_Text,          TOKENS["surface_alt"])
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 4)
    return theme_id


def build_muted_text_theme() -> int:
    """Secondary/muted text color."""
    with dpg.theme() as theme_id:
        with dpg.theme_component(dpg.mvAll):
            dpg.add_theme_color(dpg.mvThemeCol_Text, TOKENS["text_muted"])
    return theme_id


def build_badge_theme(stage: StageName) -> int:
    """Pill-style status badge (e.g., 'RUNNING' label next to stage title)."""
    accent_50 = get_stage_color(stage, "50")
    accent_800 = get_stage_color(stage, "800")
    with dpg.theme() as theme_id:
        with dpg.theme_component(dpg.mvAll):
            dpg.add_theme_color(dpg.mvThemeCol_Button,        accent_50)
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, accent_50)
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive,  accent_50)
            dpg.add_theme_color(dpg.mvThemeCol_Text,          accent_800)
            dpg.add_theme_style(dpg.mvStyleVar_FrameRounding, 10)
            dpg.add_theme_style(dpg.mvStyleVar_FramePadding,  8, 2)
    return theme_id


# ── Global application entry ──────────────────────────────────────────

def apply_global_theme() -> dict:
    """Build all themes and bind the base theme globally.

    Returns:
        Dict of all theme IDs: {'base', 'saprx', 'xpectra', 'quantix',
        'neutral', 'danger', 'muted', 'badge_sparx', ...}
    """
    themes = {
        "base":      build_base_theme(),
        "sparx":  build_stage_theme("sparx"),
        "xpectra":   build_stage_theme("xpectra"),
        "quantix":   build_stage_theme("quantix"),
        "neutral":   build_stage_theme("neutral"),
        "danger":    build_danger_button_theme(),
        "muted":     build_muted_text_theme(),
        "badge_sparx": build_badge_theme("saprx"),
        "badge_xpectra":  build_badge_theme("xpectra"),
        "badge_quantix":  build_badge_theme("quantix"),
        "badge_neutral":  build_badge_theme("neutral"),
    }
    dpg.bind_theme(themes["base"])
    return themes
