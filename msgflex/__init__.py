"""
MSGFLEX — Unified proteomics pipeline (SPARX + XPECTRA + QUANTIX).
"""

__version__ = "1.1.0"
__author__ = "Tulasi Rao Relangi, PhD"

from msgflex.core import (
    ConfigurationError,
    MsgflexError,
    ProjectPaths,
    QuantixError,
    SparxError,
    XpectraError,
    get_logger,
    setup_logging,
)

__all__ = [
    "__version__",
    "ProjectPaths",
    "setup_logging",
    "get_logger",
    "MsgflexError",
    "ConfigurationError",
    "SparxError",
    "XpectraError",
    "QuantixError",
]
