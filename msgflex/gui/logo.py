"""
MSGFLEX logo — M built from spectrum peaks.

Renders the letter M as four vertical peaks over an m/z baseline.
Pure drawlist primitives so there's no image asset to ship.

Usage:
    # In a toolbar layout:
    draw_logo(parent_id, size=28)

    # Or just give it a dimension and let it inline:
    with dpg.drawlist(width=28, height=28, parent=my_group):
        _draw_logo_on(0, 0, 28)
"""

from __future__ import annotations

import dearpygui.dearpygui as dpg

from msgflex.gui.theme import TOKENS


# Logo geometry in a 120×120 unit box (matches the mockup)
# Each tuple is (x, y, width, height, color_key)
# Bar positions: 28, 44, 60, 76, 92   (8px wide, 16px apart)
# Outer bars (shorter): y=42, height=50
# Inner bars (taller):  y=28, height=64
# Middle bar (dip):     y=50, height=42
_LOGO_BARS = [
    (28, 42, 8, 50, "text"),      # outer-left — dark
    (44, 28, 8, 64, "accent"),    # up-stroke
    (60, 50, 8, 42, "accent"),    # middle dip
    (76, 28, 8, 64, "accent"),    # up-stroke
    (92, 42, 8, 50, "text"),      # outer-right — dark
]
_BASELINE = (22, 96, 104, 96)  # x1, y1, x2, y2


def _draw_logo_on(
    draw_parent: int | str,
    *,
    x: float = 0,
    y: float = 0,
    size: float = 120,
    dark_on_dark: bool = False,
) -> None:
    """Draw the logo into an existing drawlist or draw_node.

    Args:
        draw_parent: ID of a drawlist, draw_node, or viewport drawlist.
        x, y: Top-left offset within the drawlist.
        size: Overall size in pixels (logo is square).
        dark_on_dark: If True, use lighter colors for dark backgrounds.
    """
    scale = size / 120.0

    if dark_on_dark:
        bar_color = (141, 133, 224, 255)   # #8D85E0
        baseline_color = (141, 133, 224, 255)
        accent_color = (175, 169, 236, 255)  # #AFA9EC
    else:
        bar_color = TOKENS["text"]
        baseline_color = TOKENS["text"]
        accent_color = TOKENS["accent"]

    for bx, by, bw, bh, key in _LOGO_BARS:
        color = accent_color if key == "accent" else bar_color
        dpg.draw_rectangle(
            pmin=(x + bx * scale, y + by * scale),
            pmax=(x + (bx + bw) * scale, y + (by + bh) * scale),
            color=color,
            fill=color,
            rounding=1,
            parent=draw_parent,
        )

    # Baseline (m/z axis)
    x1, y1, x2, y2 = _BASELINE
    dpg.draw_line(
        p1=(x + x1 * scale, y + y1 * scale),
        p2=(x + x2 * scale, y + y2 * scale),
        color=baseline_color,
        thickness=max(1, 1.5 * scale),
        parent=draw_parent,
    )


def draw_logo(parent: int | str, *, size: int = 28) -> int:
    """Create a drawlist of the given size and render the logo into it.

    Args:
        parent: DPG parent (group, window, etc.) to attach the drawlist to.
        size:   Pixel size of the logo (square).

    Returns:
        The drawlist ID.
    """
    dl = dpg.add_drawlist(width=size, height=size, parent=parent)
    _draw_logo_on(dl, size=size)
    return dl
