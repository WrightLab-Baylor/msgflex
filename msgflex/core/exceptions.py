"""
Exception hierarchy for MSGFLEX.

All stages raise a subclass of MsgflexError on failure so the CLI can
catch a single base class and report cleanly. Never use sys.exit() from
library code — let the CLI layer translate exceptions to exit codes.
"""


class MsgflexError(Exception):
    """Base class for all MSGFLEX errors."""


class ConfigurationError(MsgflexError):
    """Invalid config, missing required files, wrong directory layout."""


class SparxError(MsgflexError):
    """Failure in the MS-GF+ conventional or binning pipeline."""


class XpectraError(MsgflexError):
    """Failure in the XPECTRA rescoring stage."""


class QuantixError(MsgflexError):
    """Failure in the QUANTIX downstream quantification stage."""


class ToolsNotFoundError(ConfigurationError):
    """An external binary or tool file could not be located."""
