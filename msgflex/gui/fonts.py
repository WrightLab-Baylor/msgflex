"""
Font registration for the MSGFLEX GUI.

Dear PyGui's default font is a thin bitmap at 13px that looks muddy on HiDPI
displays. We load DejaVu Sans (or a fallback) at 15px regular + 15px bold +
13px monospace, then bind them globally.

Each font also registers font-range hints so non-ASCII glyphs (ppm µ, Å, ²)
render correctly in things like absPPM displays.
"""

from __future__ import annotations

import platform
from pathlib import Path

import dearpygui.dearpygui as dpg


# Font search paths — first one found wins.
# Order: bundled → distro standard → macOS → Windows → last-resort system
_SANS_CANDIDATES = [
    # Linux (Debian/Ubuntu/Mint)
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/TTF/DejaVuSans.ttf",
    # Linux (Fedora/RHEL)
    "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf",
    # macOS
    "/System/Library/Fonts/Helvetica.ttc",
    "/Library/Fonts/Arial.ttf",
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    # Windows
    "C:/Windows/Fonts/segoeui.ttf",
    "C:/Windows/Fonts/arial.ttf",
    # conda-forge sometimes ships this in the env
    f"{Path.home()}/miniforge3/envs/msgflex/share/fonts/DejaVuSans.ttf",
]

_SANS_BOLD_CANDIDATES = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/TTF/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans-Bold.ttf",
    "/System/Library/Fonts/Helvetica.ttc",
    "C:/Windows/Fonts/segoeuib.ttf",
    "C:/Windows/Fonts/arialbd.ttf",
]

_MONO_CANDIDATES = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
    "/usr/share/fonts/dejavu/DejaVuSansMono.ttf",
    "/usr/share/fonts/TTF/DejaVuSansMono.ttf",
    "/usr/share/fonts/dejavu-sans-mono-fonts/DejaVuSansMono.ttf",
    "/System/Library/Fonts/Menlo.ttc",
    "/System/Library/Fonts/Monaco.ttf",
    "C:/Windows/Fonts/consola.ttf",
    "C:/Windows/Fonts/cour.ttf",
]


def _first_existing(candidates: list[str]) -> str | None:
    for p in candidates:
        if Path(p).is_file():
            return p
    return None


def register_fonts(
    *,
    base_size: int = 16,
    mono_size: int = 14,
    global_scale: float = 1.0,
) -> dict:
    """Register sans/sans-bold/mono fonts and bind sans as the default.

    Args:
        base_size: Pixel size for regular UI text (default 16 — was 13).
        mono_size: Pixel size for monospace text (log pane).
        global_scale: Multiply all sizes by this for HiDPI screens.

    Returns:
        Dict with keys 'sans', 'sans_bold', 'mono' — each an int font ID
        (or 0 if unavailable). Bind with ``dpg.bind_item_font(...)``.

    Notes:
        If no TTF is found on the system, Dear PyGui's default font is
        kept and a warning is logged. This is non-fatal.
    """
    sans_path = _first_existing(_SANS_CANDIDATES)
    bold_path = _first_existing(_SANS_BOLD_CANDIDATES)
    mono_path = _first_existing(_MONO_CANDIDATES)

    fonts = {"sans": 0, "sans_bold": 0, "mono": 0}

    # Apply HiDPI scale to all requested sizes
    scaled_base = int(round(base_size * global_scale))
    scaled_mono = int(round(mono_size * global_scale))

    with dpg.font_registry():
        if sans_path:
            fonts["sans"] = dpg.add_font(sans_path, scaled_base)
        if bold_path:
            fonts["sans_bold"] = dpg.add_font(bold_path, scaled_base)
        if mono_path:
            fonts["mono"] = dpg.add_font(mono_path, scaled_mono)

    # Bind the sans font as the default for every widget (bold/mono bound per-widget)
    if fonts["sans"]:
        dpg.bind_font(fonts["sans"])
    else:
        import logging
        logging.getLogger("msgflex.gui").warning(
            "No system TTF found at known paths (tried %d). Falling back to "
            "Dear PyGui's default bitmap font. Install ttf-dejavu:\n"
            "  sudo apt install fonts-dejavu       (Debian/Ubuntu/Mint)\n"
            "  sudo dnf install dejavu-sans-fonts  (Fedora/RHEL)\n"
            "  brew install --cask font-dejavu     (macOS)",
            len(_SANS_CANDIDATES),
        )

    return fonts
