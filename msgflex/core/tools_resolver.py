
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path
from typing import Optional

from msgflex.core.exceptions import ToolNotFoundError

# ── Constants ─────────────────────────────────────────────────────────────────
_PLATFORM: str = sys.platform          # 'linux', 'darwin', 'win32'
_IS_WINDOWS: bool = _PLATFORM == "win32"

def _search_bases() -> list[Path]:
    bases: list[Path] = []

    # 0. Explicit tools dir set by install.sh or user (highest priority base)
    tools_dir_env = os.environ.get("MSGFLEX_TOOLS_DIR", "")
    if tools_dir_env:
        bases.append(Path(tools_dir_env))

    # 1. Conda-deployed share directory
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix:
        bases.append(Path(conda_prefix) / "share" / "msgflex" / "tools")

    # 2. Dev bundled tree
    dev_tools = Path(__file__).resolve().parents[2] / "tools"
    if dev_tools.is_dir():
        bases.append(dev_tools)

    return bases

def _find_relative(relative: str | Path) -> Optional[Path]:
    """Return the first hit of `relative` across all search bases, or None."""
    for base in _search_bases():
        candidate = base / relative
        if candidate.exists():
            return candidate.resolve()
    return None

# ── Java ──────────────────────────────────────────────────────────────────────
def resolve_java() -> Path:
    """Return the path to a ``java`` executable (≥ 11), or raise."""
    java_home = os.environ.get("JAVA_HOME")
    if java_home:
        candidate = Path(java_home) / "bin" / ("java.exe" if _IS_WINDOWS else "java")
        if candidate.is_file():
            return candidate.resolve()

    found = shutil.which("java")
    if found:
        return Path(found).resolve()

    raise ToolNotFoundError(
        "Java runtime not found.\n"
        "  Install with:  conda install -c conda-forge openjdk>=11\n"
        "  Or set JAVA_HOME to an existing JDK/JRE installation."
    )

# ── MSGFPlus ──────────────────────────────────────────────────────────────────
def resolve_msgfplus() -> Path:
    """Return the path to ``MSGFPlus.jar``, or raise."""
    # Env override
    env = os.environ.get("MSGFLEX_MSGFPLUS_JAR")
    if env:
        p = Path(env)
        if p.is_file():
            return p.resolve()
        raise ToolNotFoundError(
            f"MSGFLEX_MSGFPLUS_JAR is set but the file was not found:\n  {p}"
        )

    # Search standard locations
    jar = _find_relative("MSGFPlus.jar")
    if jar:
        return jar

    raise ToolNotFoundError(
        "MSGFPlus.jar not found.\n"
        "Remediation options (choose one):\n"
        "  1. Run:  msgflex setup-tools   (downloads the 2024 release automatically)\n"
        "  2. Set:  MSGFLEX_MSGFPLUS_JAR=/path/to/MSGFPlus.jar\n"
        "  3. Place MSGFPlus.jar in:  <conda_prefix>/share/msgflex/tools/"
    )

def resolve_msgfplus_cmd() -> list[str]:
    """Return ``[java, -jar, /path/to/MSGFPlus.jar]`` ready for subprocess."""
    return [str(resolve_java()), "-jar", str(resolve_msgfplus())]

# ── MASIC ─────────────────────────────────────────────────────────────────────
def resolve_masic() -> Path:
    """Return the path to ``MASIC.exe``, or raise with platform-aware hint."""
    env = os.environ.get("MSGFLEX_MASIC_PATH")
    if env:
        p = Path(env)
        if p.is_file():
            return p.resolve()
        raise ToolNotFoundError(
            f"MSGFLEX_MASIC_PATH is set but the file was not found:\n  {p}"
        )

    # Search standard locations
    for candidate_name in ("MASIC/MASIC.exe", "MASIC.exe"):
        hit = _find_relative(candidate_name)
        if hit:
            return hit

    # Platform-specific error message
    if _IS_WINDOWS:
        raise ToolNotFoundError(
            "MASIC.exe not found.\n"
            "  1. Run:  msgflex setup-tools\n"
            "  2. Set:  MSGFLEX_MASIC_PATH=C:\\path\\to\\MASIC.exe\n"
            "  3. Place MASIC/ folder under <conda_prefix>/share/msgflex/tools/"
        )
    else:
        raise ToolNotFoundError(
            "MASIC.exe not found.\n"
            "MASIC is a .NET application.  On Linux/macOS it runs via Mono:\n"
            "  conda install -c conda-forge mono\n"
            "Then either:\n"
            "  1. Run:  msgflex setup-tools\n"
            "  2. Set:  MSGFLEX_MASIC_PATH=/path/to/MASIC.exe"
        )

def resolve_masic_cmd() -> list[str]:
    """
    Return the full command prefix to invoke MASIC.

    Windows  → ``[MASIC.exe, ...]``
    Linux/macOS → ``[mono, MASIC.exe, ...]``
    """
    exe = resolve_masic()
    if _IS_WINDOWS:
        return [str(exe)]
    return [_require_mono(), str(exe)]

# ── PHRP ──────────────────────────────────────────────────────────────────────
def resolve_phrp() -> Path:
    """Return the path to the PHRP executable, or raise with platform-aware hint."""
    env = os.environ.get("MSGFLEX_PHRP_PATH")
    if env:
        p = Path(env)
        if p.is_file():
            return p.resolve()
        raise ToolNotFoundError(
            f"MSGFLEX_PHRP_PATH is set but the file was not found:\n  {p}"
        )

    # PHRP ships under different names depending on release
    for candidate_name in (
        "PHRP/PeptideHitResultsProcessor.exe",
        "PHRP/PHRP.exe",
        "PeptideHitResultsProcessor.exe",
        "PHRP.exe",
    ):
        hit = _find_relative(candidate_name)
        if hit:
            return hit

    if _IS_WINDOWS:
        raise ToolNotFoundError(
            "PHRP (PeptideHitResultsProcessor) not found.\n"
            "  1. Run:  msgflex setup-tools\n"
            "  2. Set:  MSGFLEX_PHRP_PATH=C:\\path\\to\\PHRP.exe"
        )
    else:
        raise ToolNotFoundError(
            "PHRP (PeptideHitResultsProcessor) not found.\n"
            "PHRP is a .NET application.  On Linux/macOS it runs via Mono:\n"
            "  conda install -c conda-forge mono\n"
            "Then either:\n"
            "  1. Run:  msgflex setup-tools\n"
            "  2. Set:  MSGFLEX_PHRP_PATH=/path/to/PHRP.exe"
        )

def resolve_phrp_cmd() -> list[str]:
    """Return the full command prefix to invoke PHRP (Mono-wrapped on non-Windows)."""
    exe = resolve_phrp()
    if _IS_WINDOWS:
        return [str(exe)]
    return [_require_mono(), str(exe)]

# ── Internal helpers ──────────────────────────────────────────────────────────

def _require_mono() -> str:
    """Return the path to the ``mono`` executable or raise helpfully."""
    mono = shutil.which("mono")
    if mono:
        return mono
    raise ToolNotFoundError(
        "Mono runtime not found (required to run .NET tools on Linux/macOS).\n"
        "  Install with:  conda install -c conda-forge mono"
    )

# ── Bulk health-check (used by `msgflex check` and setup wizard) ──────────────
_TOOL_RESOLVERS: dict[str, object] = {
    "java":     resolve_java,
    "msgfplus": resolve_msgfplus,
    "masic":    resolve_masic,
    "phrp":     resolve_phrp,
}

def check_all_tools(*, strict: bool = False) -> dict[str, str]:

    status: dict[str, str] = {}
    for name, resolver in _TOOL_RESOLVERS.items():
        try:
            path = resolver()            # type: ignore[operator]
            status[name] = f"ok:{path}"
        except ToolNotFoundError as exc:
            first_line = str(exc).splitlines()[0]
            status[name] = f"missing:{first_line}"
            if strict:
                raise
    return status

def assert_tools_available(*names: str) -> None:
    missing: list[str] = []
    for name in names:
        if name not in _TOOL_RESOLVERS:
            raise ValueError(f"Unknown tool name: {name!r}. "
                             f"Valid names: {sorted(_TOOL_RESOLVERS)}")
        try:
            _TOOL_RESOLVERS[name]()      # type: ignore[operator]
        except ToolNotFoundError as exc:
            missing.append(f"  [{name}] {str(exc).splitlines()[0]}")

    if missing:
        raise ToolNotFoundError(
            f"{len(missing)} required tool(s) not found:\n" + "\n".join(missing)
        )
