"""
msgflex/tools.py
~~~~~~~~~~~~~~~
Central resolver for all external tool paths used by MSGFLEX pipelines.
"""

import os
from pathlib import Path

# Tool registry 
TOOLS = {
    "msgfplus_jar": "MSGFPlus.jar",
    "thermorawfileparser": "ThermoRawFileParser/ThermoRawFileParser.exe",
    "mzidmerger": "MzidMerger/net8.0/MzidMerger.exe",
    "phrp": "PHRP/PeptideHitResultsProcRunner.exe",
    "masic": "MASIC/MASIC_Console.exe",
    "masic_params":"MASICParameters.xml",
    "phrp_mods": "PHRP/MSGFDB_Mods.txt",
    "phrp_mass_corr": "PHRP/Mass_Correction_Tags.txt",
}

#  Core resolving
def tools_dir():
    """
    Return the MSGFLEX tools directory as a Path.
    """
    env = os.environ.get("MSGFX_TOOLS_DIR", "").strip()
    if env:
        p = Path(env)
        if p.is_dir():
            return p
        import warnings
        warnings.warn(
            f"MSGFX_TOOLS_DIR='{env}' is set but does not exist.",
            RuntimeWarning,
            stacklevel=2,
        )

    repo_tools = Path(__file__).parent.parent / "tools"
    if repo_tools.is_dir():
        return repo_tools.resolve()

    raise RuntimeError(
        "Cannot locate MSGFLEX tools directory.\n"
        "Set MSGFX_TOOLS_DIR=/path/to/msgflex/tools and re-launch."
    )

def tool(relative_path):
    """
    Resolve a tool path relative to the MSGFLEX tools directory.
    """
    resolved = tools_dir() / relative_path
    if not resolved.exists():
        raise FileNotFoundError(
            f"\n"
            f" Tool not found : {resolved}\n"
            f" Tools dir : {tools_dir()}\n"
            f"\n"
            f" Fix: set MSGFX_TOOLS_DIR to the directory containing\n"
            f" MSGFPlus.jar, MASIC/, PHRP/, etc.\n"
            f"\n"
            f" Example:\n"
            f"  export MSGFX_TOOLS_DIR=/path/to/msgflex/tools\n"
            f"  Or re-run install.sh which sets this automatically.\n"
        )
    return str(resolved)
