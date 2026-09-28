"""
Shared application state for the MSGFLEX GUI.

Holds:
  - current ProjectPaths (or None until the user picks a base dir)
  - current MSGFPlus params file
  - current JVM heap size
  - subscriber list for change notifications

Tabs register callbacks here and rerender when the project changes
(e.g., Preflight tab re-runs checks when base dir is updated).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

from msgflex.core import ProjectPaths


@dataclass
class AppState:
    paths: Optional[ProjectPaths] = None
    params_file: Optional[Path] = None
    java_mem: str = "4G"
    _subscribers: list[Callable[["AppState"], None]] = field(default_factory=list)

    def set_base(self, base: str | Path) -> None:
        self.paths = ProjectPaths.from_base(base)
        self._notify()

    def set_params(self, path: str | Path) -> None:
        self.params_file = Path(path)
        self._notify()

    def set_java_mem(self, value: str) -> None:
        self.java_mem = value
        self._notify()

    def clear(self) -> None:
        self.paths = None
        self.params_file = None
        self._notify()

    def subscribe(self, callback: Callable[["AppState"], None]) -> None:
        self._subscribers.append(callback)

    def _notify(self) -> None:
        for cb in self._subscribers:
            try:
                cb(self)
            except Exception as e:
                # Never let a broken subscriber take down the app
                import logging
                logging.getLogger("msgflex.gui").warning(
                    "State subscriber error: %s", e
                )

    @property
    def has_project(self) -> bool:
        return self.paths is not None

    @property
    def has_params(self) -> bool:
        return self.params_file is not None and self.params_file.is_file()
