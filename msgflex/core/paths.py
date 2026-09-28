"""
ProjectPaths — single source of truth for directory layout.

Every stage (SPARX, XPECTRA, QUANTIX) reads and writes paths derived
from one base directory. No more .parent walking, no more hardcoded
subdir names scattered across modules.

Standard layout under `base`:

    base/
    ├── data/                       # RAW and mzML
    ├── database/                   # FASTA files
    ├── inputfile.tsv               # decoder (msfilename, database, QCplot)
    ├── QCdir/                      # TIC plots
    ├── SICdir/                     # sparx FINAL output (PHRP + SIC merged)
    ├── results/
    │   ├── PHRPOut/                # PHRP per-sample outputs
    │   └── MasicOut/               # MASIC per-sample outputs
    ├── xpectra/                    # XPECTRA output (optional)
    │   ├── features/               #   → spectral features per sample
    │   ├── rescored/               #   → *_rescored.tsv  (fed into QUANTIX)
    │   ├── stats/                  #   → per-sample stats + plots
    │   └── logs/                   #   → per-sample logs
    └── final.tsv                   # QUANTIX final LFQ table (default name)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from msgflex.core.exceptions import ConfigurationError

@dataclass(frozen=True)
class ProjectPaths:
    """
    Immutable view of a project's directory layout.

    Construct with ``ProjectPaths.from_base(base)``. All attributes are
    absolute ``Path`` objects.
    """

    base: Path

    data: Path
    database: Path
    inputfile: Path

    qcdir: Path
    sicdir: Path
    results: Path
    phrp_out: Path
    masic_out: Path

    xpectra_out: Path
    xpectra_features: Path
    xpectra_rescored: Path
    xpectra_stats: Path
    xpectra_logs: Path

    quantix_sics: Path          
    quantix_fdr_estd: Path      
    quantix_fdr_filt: Path      
    quantix_peptide_ctab: Path  
    quantix_annotated: Path     

    @classmethod
    def from_base(cls, base: str | Path) -> "ProjectPaths":
        """Build a ProjectPaths from a single base directory."""
        base = Path(base).expanduser().resolve()
        results = base / "results"
        xpectra_out = base / "xpectra"

        return cls(
            base=base,
            data=base / "data",
            database=base / "database",
            inputfile=base / "inputfile.tsv",
            qcdir=base / "QCdir",
            sicdir=base / "SICdir",
            results=results,
            phrp_out=results / "PHRPOut",
            masic_out=results / "MasicOut",
            xpectra_out=xpectra_out,
            xpectra_features=xpectra_out / "features",
            xpectra_rescored=xpectra_out / "rescored",
            xpectra_stats=xpectra_out / "stats",
            xpectra_logs=xpectra_out / "logs",
            quantix_sics=base / "SICs",
            quantix_fdr_estd=base / "fdr_estd",
            quantix_fdr_filt=base / "fdr_filt",
            quantix_peptide_ctab=base / "peptide_crosstab.tsv",
            quantix_annotated=base / "peptide_crosstab_annotated.tsv",
        )

    # ── Pre-flight validators ─────────────────────────────────────────
    def require_sparx_inputs(self) -> None:
        """Raise ConfigurationError if sparx is missing required inputs."""
        if not self.data.is_dir():
            raise ConfigurationError(f"Data directory not found: {self.data}")
        if not self.database.is_dir():
            raise ConfigurationError(f"Database directory not found: {self.database}")
        # inputfile.tsv is generated if missing, so don't require it here

    def require_xpectra_inputs(self) -> None:
        """Raise if XPECTRA preconditions aren't met (sparx must have run)."""
        if not self.data.is_dir():
            raise ConfigurationError(f"mzML data directory not found: {self.data}")
        if not self.sicdir.is_dir():
            raise ConfigurationError(
                f"SICdir not found at {self.sicdir}. "
                f"Run `msgflex sparx (upstream)` first."
            )
        has_sic_stats = any(self.sicdir.glob("*_fht_PlusSICStats.txt"))
        if not has_sic_stats:
            raise ConfigurationError(
                f"No *_fht_PlusSICStats.txt files in {self.sicdir}. "
                f"Upstream pipeline may not have completed."
            )

    def require_quantix_inputs(self, *, from_xpectra: bool) -> None:
        """Raise if QUANTIX preconditions aren't met."""
        if not self.database.is_dir():
            raise ConfigurationError(f"Database directory not found: {self.database}")
        if not self.phrp_out.is_dir():
            raise ConfigurationError(
                f"PHRPOut not found at {self.phrp_out}. Run sparx first."
            )
        if from_xpectra:
            if not self.xpectra_rescored.is_dir():
                raise ConfigurationError(
                    f"XPECTRA rescored dir not found at {self.xpectra_rescored}. "
                    f"Run `msgflex xpectra` first or omit --rescore."
                )
            if not any(self.xpectra_rescored.glob("*_rescored.tsv")):
                raise ConfigurationError(
                    f"No *_rescored.tsv files in {self.xpectra_rescored}."
                )
        else:
            if not self.sicdir.is_dir():
                raise ConfigurationError(
                    f"SICdir not found at {self.sicdir}. Run sparx first."
                )

    def ensure_dirs(self) -> None:
        """Create all output directories (idempotent)."""
        for p in (
            self.qcdir,
            self.sicdir,
            self.results,
            self.phrp_out,
            self.masic_out,
        ):
            p.mkdir(parents=True, exist_ok=True)

    def ensure_xpectra_dirs(self) -> None:
        for p in (
            self.xpectra_out,
            self.xpectra_features,
            self.xpectra_rescored,
            self.xpectra_stats,
            self.xpectra_logs,
        ):
            p.mkdir(parents=True, exist_ok=True)
