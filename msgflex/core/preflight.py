"""
Preflight validation for MSGFLEX.

Runs exhaustive read-only checks and reports pass/warn/fail for every
precondition needed to run the SPARX, XPECTRA, and QUANTIX stages.
No external binaries are invoked.

Result collection is structured so the CLI layer can format as text
and the future GUI can render as a table.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import List, Optional

import pandas as pd

from msgflex.core import ProjectPaths

class Status(str, Enum):
    OK = "OK"
    WARN = "WARN"
    FAIL = "FAIL"
    SKIP = "SKIP"


@dataclass
class CheckResult:
    """One check's outcome."""
    name: str
    status: Status
    message: str = ""
    details: List[str] = field(default_factory=list)

    @property
    def is_fatal(self) -> bool:
        return self.status == Status.FAIL


@dataclass
class PreflightReport:
    """Aggregated results across all checks in a category."""
    category: str
    results: List[CheckResult] = field(default_factory=list)

    def add(self, result: CheckResult) -> None:
        self.results.append(result)

    @property
    def has_failures(self) -> bool:
        return any(r.is_fatal for r in self.results)

    @property
    def has_warnings(self) -> bool:
        return any(r.status == Status.WARN for r in self.results)

    def count(self, status: Status) -> int:
        return sum(1 for r in self.results if r.status == status)


# ── Environment checks ────────────────────────────────────────────────

def check_environment() -> PreflightReport:
    """Verify Python dependencies and external runtimes."""
    report = PreflightReport("Environment")

    # Python packages
    packages = [
        ("typer", "CLI framework"),
        ("pandas", "data handling"),
        ("numpy", "numerical ops"),
        ("pyteomics", "XPECTRA spectral features"),
        ("xgboost", "XPECTRA classifier"),
        ("sklearn", "XPECTRA ML"),
        ("statsmodels", "XPECTRA RT calibration (LOWESS)"),
        ("matplotlib", "QC plots"),
        ("seaborn", "QUANTIX plots"),
        ("Bio", "Biopython for FASTA handling"),
    ]
    for pkg, purpose in packages:
        try:
            __import__(pkg)
            report.add(CheckResult(
                name=f"python:{pkg}", status=Status.OK,
                message=f"importable ({purpose})",
            ))
        except ImportError as e:
            report.add(CheckResult(
                name=f"python:{pkg}", status=Status.FAIL,
                message=f"not importable ({purpose})",
                details=[str(e)],
            ))

    # External binaries
    java = shutil.which("java")
    if java:
        try:
            out = subprocess.run(
                ["java", "-version"], capture_output=True, text=True, timeout=5,
            )
            version_line = (out.stderr or out.stdout).splitlines()[0] if (out.stderr or out.stdout) else "?"
            report.add(CheckResult(
                name="binary:java", status=Status.OK,
                message=f"{java}", details=[version_line],
            ))
        except (subprocess.TimeoutExpired, OSError) as e:
            report.add(CheckResult(
                name="binary:java", status=Status.WARN,
                message=f"found at {java} but version check failed: {e}",
            ))
    else:
        report.add(CheckResult(
            name="binary:java", status=Status.FAIL,
            message="not on PATH (required for MS-GF+)",
        ))

    mono = shutil.which("mono")
    if mono:
        try:
            out = subprocess.run(
                ["mono", "--version"], capture_output=True, text=True, timeout=5,
            )
            version_line = out.stdout.splitlines()[0] if out.stdout else "?"
            report.add(CheckResult(
                name="binary:mono", status=Status.OK,
                message=f"{mono}", details=[version_line],
            ))
        except (subprocess.TimeoutExpired, OSError) as e:
            report.add(CheckResult(
                name="binary:mono", status=Status.WARN,
                message=f"found at {mono} but version check failed: {e}",
            ))
    else:
        report.add(CheckResult(
            name="binary:mono", status=Status.FAIL,
            message="not on PATH (required for MASIC, PHRP, ThermoRawFileParser)",
        ))

    return report


# ── Tools directory ───────────────────────────────────────────────────

REQUIRED_TOOLS = [
    "MSGFPlus.jar",
    "MASIC/MASIC_Console.exe",
    "PHRP/PeptideHitResultsProcRunner.exe",
    "PHRP/MSGFDB_Mods.txt",
    "PHRP/Mass_Correction_Tags.txt",
    "ThermoRawFileParser/ThermoRawFileParser.exe",
    "MASICParameters.xml",
]

BINNING_TOOLS = [
    "MzidMerger/net8.0/MzidMerger.exe",
]


def check_tools() -> PreflightReport:
    """Verify the tools/ directory is populated."""
    report = PreflightReport("Tools")

    tools_dir = os.environ.get("MSGFLEX_TOOLS_DIR")
    if not tools_dir:
        # Fallback: look next to the installed package
        try:
            from msgflex.tools import tool
            # Try resolving a known tool; any exception → missing
            tool("MSGFPlus.jar")
            # If tool() worked, fish the dir out by re-resolving
            _probe = Path(tool("MSGFPlus.jar")).parent
            tools_dir = str(_probe)
        except Exception:
            report.add(CheckResult(
                name="tools:MSGFLEX_TOOLS_DIR",
                status=Status.FAIL,
                message="environment variable not set and tool() cannot resolve",
                details=[
                    "Set MSGFLEX_TOOLS_DIR to the directory containing MSGFPlus.jar",
                    "Or reinstall via install.sh which sets this automatically",
                ],
            ))
            return report

    tools_path = Path(tools_dir).expanduser().resolve()
    if not tools_path.is_dir():
        report.add(CheckResult(
            name="tools:MSGFLEX_TOOLS_DIR", status=Status.FAIL,
            message=f"directory does not exist: {tools_path}",
        ))
        return report

    report.add(CheckResult(
        name="tools:MSGFLEX_TOOLS_DIR", status=Status.OK,
        message=str(tools_path),
    ))

    for tool_rel in REQUIRED_TOOLS:
        full = tools_path / tool_rel
        if full.is_file():
            size_mb = full.stat().st_size / (1024 * 1024)
            report.add(CheckResult(
                name=f"tools:{tool_rel}", status=Status.OK,
                message=f"{size_mb:.1f} MB",
            ))
        else:
            report.add(CheckResult(
                name=f"tools:{tool_rel}", status=Status.FAIL,
                message=f"missing: {full}",
            ))

    for tool_rel in BINNING_TOOLS:
        full = tools_path / tool_rel
        if full.is_file():
            report.add(CheckResult(
                name=f"tools:{tool_rel}", status=Status.OK,
                message="present (needed for binning mode)",
            ))
        else:
            report.add(CheckResult(
                name=f"tools:{tool_rel}", status=Status.WARN,
                message="missing — binning mode will fail",
            ))

    return report


# ── Project layout ────────────────────────────────────────────────────

def check_project_sparx(paths: ProjectPaths) -> PreflightReport:
    """Verify project directory layout and inputs for sparx."""
    report = PreflightReport(f"Project (Sparx) — {paths.base}")

    # data/
    if paths.data.is_dir():
        raws = list(paths.data.glob("*.raw"))
        mzmls = list(paths.data.glob("*.mzML"))
        if not raws and not mzmls:
            report.add(CheckResult(
                name="project:data", status=Status.FAIL,
                message=f"no *.raw or *.mzML files in {paths.data}",
            ))
        else:
            msg = f"{len(raws)} RAW, {len(mzmls)} mzML"
            report.add(CheckResult(
                name="project:data", status=Status.OK, message=msg,
            ))
    else:
        report.add(CheckResult(
            name="project:data", status=Status.FAIL,
            message=f"missing directory: {paths.data}",
        ))

    # database/
    if paths.database.is_dir():
        fastas = list(paths.database.glob("*.fasta")) + list(paths.database.glob("*.faa"))
        if not fastas:
            report.add(CheckResult(
                name="project:database", status=Status.FAIL,
                message=f"no FASTA files in {paths.database}",
            ))
        else:
            report.add(CheckResult(
                name="project:database", status=Status.OK,
                message=f"{len(fastas)} FASTA",
                details=[f.name for f in fastas[:5]],
            ))
    else:
        report.add(CheckResult(
            name="project:database", status=Status.FAIL,
            message=f"missing directory: {paths.database}",
        ))

    # inputfile.tsv (optional — generated if missing for conventional)
    if paths.inputfile.is_file():
        try:
            df = pd.read_csv(paths.inputfile, sep="\t")
            required_cols = {"msfilename", "database", "QCplot"}
            missing = required_cols - set(df.columns)
            if missing:
                report.add(CheckResult(
                    name="project:inputfile.tsv", status=Status.FAIL,
                    message=f"missing columns: {missing}",
                ))
            else:
                # Cross-check sample→database mapping
                missing_db = []
                if paths.database.is_dir():
                    for _, row in df.iterrows():
                        db = paths.database / str(row["database"])
                        if not db.exists():
                            missing_db.append(str(row["database"]))
                if missing_db:
                    report.add(CheckResult(
                        name="project:inputfile.tsv",
                        status=Status.WARN,
                        message=f"{len(df)} samples, {len(missing_db)} reference databases not found",
                        details=missing_db[:5],
                    ))
                else:
                    # Cross-check sample→raw mapping
                    if paths.data.is_dir():
                        missing_raw = []
                        for _, row in df.iterrows():
                            name = str(row["msfilename"])
                            has_raw = (paths.data / f"{name}.raw").exists()
                            has_mzml = (paths.data / f"{name}.mzML").exists()
                            if not (has_raw or has_mzml):
                                missing_raw.append(name)
                        if missing_raw:
                            report.add(CheckResult(
                                name="project:inputfile.tsv",
                                status=Status.WARN,
                                message=f"{len(df)} samples, {len(missing_raw)} have no RAW/mzML",
                                details=missing_raw[:5],
                            ))
                        else:
                            report.add(CheckResult(
                                name="project:inputfile.tsv",
                                status=Status.OK,
                                message=f"{len(df)} samples, all RAW and DB references valid",
                            ))
                    else:
                        report.add(CheckResult(
                            name="project:inputfile.tsv",
                            status=Status.OK,
                            message=f"{len(df)} samples, DB references valid",
                        ))
        except Exception as e:
            report.add(CheckResult(
                name="project:inputfile.tsv", status=Status.FAIL,
                message=f"cannot parse: {e}",
            ))
    else:
        report.add(CheckResult(
            name="project:inputfile.tsv", status=Status.WARN,
            message="not present — will be auto-generated for conventional mode; required for binning",
        ))

    # Output directories (should be creatable)
    for out_path, label in [
        (paths.qcdir, "QCdir"),
        (paths.sicdir, "SICdir"),
        (paths.results, "results/"),
    ]:
        if out_path.exists():
            report.add(CheckResult(
                name=f"project:{label}", status=Status.OK,
                message="exists (will overwrite)",
            ))
        else:
            parent = out_path.parent
            writable = parent.exists() and os.access(parent, os.W_OK)
            if writable:
                report.add(CheckResult(
                    name=f"project:{label}", status=Status.OK,
                    message=f"will create at {out_path}",
                ))
            else:
                report.add(CheckResult(
                    name=f"project:{label}", status=Status.FAIL,
                    message=f"parent not writable: {parent}",
                ))

    return report

def check_project_xpectra(paths: ProjectPaths) -> PreflightReport:
    """Verify preconditions for XPECTRA stage."""
    report = PreflightReport(f"Project (XPECTRA) — {paths.base}")

    if not paths.data.is_dir():
        report.add(CheckResult(
            name="xpectra:data", status=Status.FAIL,
            message=f"mzML directory missing: {paths.data}",
        ))
    else:
        mzmls = list(paths.data.glob("*.mzML"))
        if not mzmls:
            report.add(CheckResult(
                name="xpectra:data", status=Status.FAIL,
                message=f"no *.mzML files in {paths.data}",
            ))
        else:
            report.add(CheckResult(
                name="xpectra:data", status=Status.OK,
                message=f"{len(mzmls)} mzML files",
            ))

    if not paths.sicdir.is_dir():
        report.add(CheckResult(
            name="xpectra:sicdir", status=Status.FAIL,
            message=f"SICdir missing — run sparx first: {paths.sicdir}",
        ))
    else:
        sic_stats = list(paths.sicdir.glob("*_fht_PlusSICStats.txt"))
        if not sic_stats:
            report.add(CheckResult(
                name="xpectra:sicdir", status=Status.FAIL,
                message=f"no *_fht_PlusSICStats.txt files in {paths.sicdir}",
            ))
        else:
            report.add(CheckResult(
                name="xpectra:sicdir", status=Status.OK,
                message=f"{len(sic_stats)} SIC stats files",
            ))

    return report

def check_project_quantix(paths: ProjectPaths, from_xpectra: bool = False) -> PreflightReport:
    """Verify preconditions for QUANTIX stage."""
    report = PreflightReport(f"Project (QUANTIX) — {paths.base}")

    if not paths.database.is_dir():
        report.add(CheckResult(
            name="quantix:database", status=Status.FAIL,
            message=f"database directory missing: {paths.database}",
        ))
    else:
        fastas = list(paths.database.glob("*.fasta")) + list(paths.database.glob("*.faa"))
        report.add(CheckResult(
            name="quantix:database",
            status=Status.OK if fastas else Status.FAIL,
            message=f"{len(fastas)} FASTA",
        ))

    if not paths.phrp_out.is_dir():
        report.add(CheckResult(
            name="quantix:phrp_out", status=Status.FAIL,
            message=f"PHRPOut missing — run sparx first: {paths.phrp_out}",
        ))
    else:
        syns = list(paths.phrp_out.glob("*_syn.txt"))
        report.add(CheckResult(
            name="quantix:phrp_out",
            status=Status.OK if syns else Status.WARN,
            message=f"{len(syns)} SYN files",
        ))

    if from_xpectra:
        if not paths.xpectra_rescored.is_dir():
            report.add(CheckResult(
                name="quantix:xpectra_rescored", status=Status.FAIL,
                message=f"missing — run `msgflex xpectra` first: {paths.xpectra_rescored}",
            ))
        else:
            rescored = list(paths.xpectra_rescored.glob("*_rescored.tsv"))
            report.add(CheckResult(
                name="quantix:xpectra_rescored",
                status=Status.OK if rescored else Status.FAIL,
                message=f"{len(rescored)} rescored TSVs",
            ))
    else:
        if not paths.sicdir.is_dir():
            report.add(CheckResult(
                name="quantix:sicdir", status=Status.FAIL,
                message=f"SICdir missing: {paths.sicdir}",
            ))
        else:
            sic_stats = list(paths.sicdir.glob("*_fht_PlusSICStats.txt"))
            report.add(CheckResult(
                name="quantix:sicdir",
                status=Status.OK if sic_stats else Status.FAIL,
                message=f"{len(sic_stats)} SIC stats files",
            ))

    return report


# ── Full run ──────────────────────────────────────────────────────────
# NEW:
def run_all_checks(
    paths: Optional[ProjectPaths] = None,
    *,
    stage: str = "all",
    from_xpectra: bool = False,
) -> List[PreflightReport]:
    """Run every preflight check.

    Tools are only checked when the upstream stage is involved —
    XPECTRA and QUANTIX are pure Python and don't need Java/Mono binaries.
    """
    reports = [check_environment()]

    # Tools only needed for upstream (Java for MS-GF+, Mono for MASIC/PHRP)
    if stage in ("all", "sparx"):
        reports.append(check_tools())

    if paths is None:
        return reports

    if stage in ("all", "sparx"):
        reports.append(check_project_sparx(paths))
    if stage in ("all", "xpectra"):
        reports.append(check_project_xpectra(paths))
    if stage in ("all", "quantix"):
        reports.append(check_project_quantix(paths, from_xpectra=from_xpectra))

    return reports