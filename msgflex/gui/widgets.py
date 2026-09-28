"""
Reusable browse-row widget for the MSGFLEX GUI.

Every "browse for a file/directory" form row is built with this helper
so behavior stays consistent across tabs. Includes a default-path
subscriber that keeps each input auto-filled from project paths, which
the user can override at any time.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Literal, Optional

import dearpygui.dearpygui as dpg

from msgflex.gui.theme import TOKENS


BrowseKind = Literal["file_open", "file_save", "directory"]


def build_browse_row(
    *,
    label: str,
    tag: str,
    kind: BrowseKind,
    default_value: str = "",
    extensions: Optional[list[str]] = None,
    on_change: Optional[Callable[[str], None]] = None,
    help_text: Optional[str] = None,
) -> None:
    """Render one browse row into the current DPG container.

    Args:
        label:          Left-side field label.
        tag:            DPG tag for the input text widget. Use ``get_value(tag)``
                        to read the current path elsewhere.
        kind:           'file_open' (pick existing), 'file_save' (new or
                        existing), or 'directory' (folder picker).
        default_value:  Initial path to display (can be edited or replaced
                        via browse).
        extensions:     For file_open/file_save, list of extensions to show
                        (e.g. ['.tsv', '.txt']). Ignored for directory.
        on_change:      Optional callback invoked with the new path string
                        whenever the user picks or types a new path.
        help_text:      Small hint under the input. Omitted if None.
    """
    dpg.add_text(label, color=TOKENS["text_muted"])

    def _dialog_callback(sender, app_data):
        # Directory picker → app_data['file_path_name']
        # File picker (file_open) → app_data['selections'] is {name: fullpath} dict
        # File save → use file_path_name directly so a typed filename works
        #             without requiring the user to first pick an existing file
        if kind == "directory":
            path = app_data.get("file_path_name", "")
        elif kind == "file_save":
            # Always honour the typed filename for save dialogs
            path = app_data.get("file_path_name", "")
        else:  # file_open
            sels = app_data.get("selections") or {}
            path = next(iter(sels.values())) if sels else app_data.get("file_path_name", "")

        if path:
            dpg.set_value(tag, path)
            if on_change is not None:
                try:
                    on_change(path)
                except Exception as e:
                    import logging
                    logging.getLogger("msgflex.gui").warning(
                        "Browse on_change for %s raised: %s", tag, e
                    )

    dialog_tag = f"dlg_{tag}"
    dpg.add_file_dialog(
        directory_selector=(kind == "directory"),
        show=False,
        modal=True,
        tag=dialog_tag,
        callback=_dialog_callback,
        width=700, height=450,
    )

    # Register extensions (file pickers only)
    if kind != "directory":
        # Dear PyGui defaults file-extension text to near-white which
        # disappears against our light slate background. We force dark
        # text on primary extensions and a muted grey on the '.*' wildcard
        # so the primary filter visually dominates.
        primary_color = (42, 35, 64, 255)    # TOKENS["text"]
        wildcard_color = (107, 102, 128, 255)  # TOKENS["text_muted"]
        for ext in (extensions or [".*"]):
            color = wildcard_color if ext == ".*" else primary_color
            dpg.add_file_extension(ext, parent=dialog_tag, color=color)

    with dpg.group(horizontal=True):
        dpg.add_input_text(
            tag=tag,
            default_value=default_value,
            width=-120,
            callback=(lambda s, v: on_change(v)) if on_change else None,
        )
        dpg.add_button(
            label="Browse...",
            width=110,
            callback=lambda: dpg.show_item(dialog_tag),
        )

    if help_text:
        dpg.add_text(help_text, color=TOKENS["text_faint"])


def set_browse_value_if_empty(tag: str, default_path: str | Path) -> None:
    """Update a browse input's value only if it's currently empty.

    Used when the user picks a base project dir — we fill in the XPECTRA
    and QUANTIX input/output defaults without overriding any path the user
    already typed manually.
    """
    if not dpg.does_item_exist(tag):
        return
    current = dpg.get_value(tag)
    if not current or not current.strip():
        dpg.set_value(tag, str(default_path))
