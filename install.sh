#!/usr/bin/env bash
# =============================================================================
# MSGFLEX Pipeline v1.1.0 — Installer
# =============================================================================
VERSION="$(sed -n 's/^__version__ *= *"\(.*\)"/\1/p' "$SCRIPT_DIR/msgflex/__init__.py")"
echo "  MSGFLEX Pipeline v${VERSION} Installer"

# Usage:
#   bash install.sh # default: installs to ~/msgflex
#   bash install.sh --prefix /opt # installs to /opt/msgflex
#   bash install.sh --help
# =============================================================================

set -euo pipefail

# ── Colours ───────────────────────────────────────────────────────────────────
RED='\033[0;31m';  GREEN='\033[0;32m'; YELLOW='\033[1;33m'
BLUE='\033[0;34m'; BOLD='\033[1m';     RESET='\033[0m'

info()  { echo -e "${BLUE}[INFO]${RESET}  $*"; }
success()   { echo -e "${GREEN}[OK]${RESET}    $*"; }
warn()  { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
error() { echo -e "${RED}[ERROR]${RESET} $*" >&2; }
die()   { error "$*"; exit 1; }
header()    { echo -e "\n${BOLD}── $* ──${RESET}"; }

# ── Defaults ──────────────────────────────────────────────────────────────────
INSTALL_PREFIX="$HOME/msgflex"
ENV_NAME="${MSGFLEX_ENV_NAME:-msgflex}"
TOOLS_VAR="MSGFLEX_TOOLS_DIR"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERSION="$(sed -n 's/^__version__ *= *"\(.*\)"/\1/p' "$SCRIPT_DIR/msgflex/__init__.py")"
MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)   INSTALL_PREFIX="$2"; shift 2 ;;
        --prefix=*) INSTALL_PREFIX="${1#*=}"; shift ;;
        --help|-h)
            echo "Usage: bash install.sh [--prefix /install/path]"
            echo "  Default prefix: $HOME/msgflex"
            exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

TOOLS_DIR="$SCRIPT_DIR/tools"

# =============================================================================
# Banner
# =============================================================================
echo -e "${BOLD}"
VERSION="$(sed -n 's/^__version__ *= *"\(.*\)"/\1/p' "$SCRIPT_DIR/msgflex/__init__.py")"
echo "  MSGFLEX Pipeline v${VERSION} Installer"
echo -e "${RESET}"
echo "  Install prefix : $INSTALL_PREFIX"
echo "  Tools dir : $TOOLS_DIR"
echo "  Conda env : $ENV_NAME"
echo ""

# =============================================================================
# Step 1 — System requirements
# =============================================================================
header "Step 1/6  System requirements"

# Bash version (need 4+ for associative arrays)
if (( BASH_VERSINFO[0] < 4 )); then
    die "Bash 4+ required (found $BASH_VERSION). Update via: sudo apt install bash"
fi
success "Bash $BASH_VERSION"

# OS check — Linux only (Mono is Linux-only in this pipeline)
if [[ "$(uname -s)" != "Linux" ]]; then
    die "MSGFLEX currently supports Linux only (Mono .exe tools require Linux)."
fi
success "Linux detected"

# Java (will be installed via conda if missing)
if command -v java &>/dev/null; then
    JAVA_VER=$(java -version 2>&1 | awk -F '"' '/version/ {print $2}')
    success "Java found: $JAVA_VER"
else
    warn "Java not found in PATH — will be installed via conda (openjdk)."
fi

# Mono (will be installed via conda if missing)
if command -v mono &>/dev/null; then
    MONO_VER=$(mono --version 2>&1 | head -1)
    success "Mono found: $MONO_VER"
else
    warn "Mono not found in PATH — will be installed via conda."
fi

# Downloader
if command -v curl &>/dev/null; then
    success "curl available"
elif command -v wget &>/dev/null; then
    success "wget available"
else
    die "Neither curl nor wget found. Install one: sudo apt install curl"
fi

# Tools directory
TOOLS_URL="https://github.com/thulasis/msgflex/releases/latest/download/tools.tar.gz"
if [[ ! -d "$TOOLS_DIR" ]]; then
    info "Tools directory not found — downloading from GitHub release..."
    TOOLS_ARCHIVE="/tmp/msgflex-tools.tar.gz"
    if command -v curl &>/dev/null; then
        curl -fL "$TOOLS_URL" -o "$TOOLS_ARCHIVE"
    else
        wget -q "$TOOLS_URL" -O "$TOOLS_ARCHIVE"
    fi
    [[ -s "$TOOLS_ARCHIVE" ]] || die "Download failed or returned an empty file. Check $TOOLS_URL"
    info "Extracting tools to $SCRIPT_DIR..."
    tar -xzf "$TOOLS_ARCHIVE" -C "$SCRIPT_DIR"
    rm -f "$TOOLS_ARCHIVE"
    success "Tools downloaded and extracted to $TOOLS_DIR"
fi

# Spot-check critical tools (same list as preflight)
REQUIRED_TOOLS=(
    "MSGFPlus.jar"
    "MASIC/MASIC_Console.exe"
    "PHRP/PeptideHitResultsProcRunner.exe"
    "PHRP/MSGFDB_Mods.txt"
    "PHRP/Mass_Correction_Tags.txt"
    "ThermoRawFileParser/ThermoRawFileParser.exe"
    "MASICParameters.xml"
)
MISSING_TOOLS=()
for t in "${REQUIRED_TOOLS[@]}"; do
    if [[ ! -f "$TOOLS_DIR/$t" ]]; then
        MISSING_TOOLS+=("$t")
    fi
done

if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
    error "The following tools are missing from $TOOLS_DIR :"
    for t in "${MISSING_TOOLS[@]}"; do
        error "  ✗  $t"
    done
    die "Cannot continue. Add the missing tools and re-run install.sh."
fi
success "All required tools found in $TOOLS_DIR"

# =============================================================================
# Step 2 — Locate or install Miniforge (mamba)
# =============================================================================
header "Step 2/6  Conda / Mamba"

CONDA_BASE=""

if command -v mamba &>/dev/null; then
    PKG_MGR="mamba"
    CONDA_BASE="$(dirname "$(dirname "$(readlink -f "$(command -v mamba)")")")"
    success "mamba found: $(which mamba)"
elif command -v conda &>/dev/null; then
    PKG_MGR="conda"
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    warn "conda found but not mamba. Will use conda (slower). Consider installing mamba."
    success "conda found: $(which conda)"
else
    info "No conda/mamba found. Installing Miniforge3..."
    MINIFORGE_INSTALLER="/tmp/Miniforge3-installer.sh"

    if command -v curl &>/dev/null; then
        curl -fsSL "$MINIFORGE_URL" -o "$MINIFORGE_INSTALLER"
    else
        wget -q "$MINIFORGE_URL" -O "$MINIFORGE_INSTALLER"
    fi

    bash "$MINIFORGE_INSTALLER" -b -p "$HOME/miniforge3"
    rm -f "$MINIFORGE_INSTALLER"

    export PATH="$HOME/miniforge3/bin:$PATH"
    PKG_MGR="mamba"
    CONDA_BASE="$HOME/miniforge3"
    success "Miniforge3 installed at $CONDA_BASE"
fi

CONDA_INIT_SCRIPT="$CONDA_BASE/etc/profile.d/conda.sh"
if [[ ! -f "$CONDA_INIT_SCRIPT" ]]; then
    die "Cannot find conda init script at $CONDA_INIT_SCRIPT"
fi
# shellcheck source=/dev/null
set +u; source "$CONDA_INIT_SCRIPT"; set -u

if [[ "$PKG_MGR" == "conda" ]] && command -v mamba &>/dev/null; then
    PKG_MGR="mamba"
    success "mamba now available — switching to mamba"
fi

# =============================================================================
# Step 3 — Create conda environment
# =============================================================================
header "Step 3/6  Conda environment ($ENV_NAME)"

ENV_YML="$SCRIPT_DIR/environment.yml"
if [[ ! -f "$ENV_YML" ]]; then
    die "environment.yml not found at $SCRIPT_DIR"
fi

if $PKG_MGR env list | grep -qE "^${ENV_NAME}\s"; then
    warn "Environment '$ENV_NAME' already exists."
    GUI_CONFIRM="n"; [[ -t 0 ]] && read -rp "  Install GUI support (dearpygui)? [y/N]: "
    if [[ "$CONFIRM" =~ ^[Yy]$ ]]; then
        info "Updating environment '$ENV_NAME'..."
        $PKG_MGR env update -f "$ENV_YML" --prune
        success "Environment updated"
    else
        info "Skipping environment update — using existing '$ENV_NAME'"
    fi
else
    info "Creating environment '$ENV_NAME' (this may take 3-5 minutes)..."
    $PKG_MGR env create -f "$ENV_YML"
    success "Environment '$ENV_NAME' created"
fi

# =============================================================================
# Step 4 — Install MSGFLEX package (CRITICAL: registers CLI entry points)
# =============================================================================
header "Step 4/6  Installing MSGFLEX package"

set +u; conda activate "$ENV_NAME"; set -u

if [[ ! -f "$SCRIPT_DIR/pyproject.toml" ]]; then
    die "pyproject.toml not found at $SCRIPT_DIR — cannot install package."
fi

info "Running: pip install -e $SCRIPT_DIR"
pip install -e "$SCRIPT_DIR" --quiet
success "MSGFLEX package installed (editable)"

if ! command -v msgflex &>/dev/null; then
    die "After pip install, 'msgflex' still not on PATH. Installation failed."
fi
success "CLI registered: $(which msgflex)"

# ── GUI (optional) ────────────────────────────────────────────────────
echo ""
read -rp "  Install GUI support (dearpygui)? [y/N]: " GUI_CONFIRM
GUI_INSTALLED=false
if [[ "$GUI_CONFIRM" =~ ^[Yy]$ ]]; then
    info "Installing dearpygui..."
    pip install -e "$SCRIPT_DIR[gui]" --quiet
    success "GUI installed — launch with: msgflex-gui"
    GUI_INSTALLED=true
else
    info "Skipping GUI — install later with: pip install 'msgflex[gui]'"
fi

# =============================================================================
# Step 5 — Set MSGFLEX_TOOLS_DIR via conda activation hooks
# =============================================================================
header "Step 5/6  Configuring $TOOLS_VAR"

# Conda activation hooks: scripts in env/etc/conda/activate.d/ run
# automatically on every `conda activate msgflex`
ENV_PREFIX="$(python -c 'import sys; print(sys.prefix)')"
ACTIVATE_DIR="$ENV_PREFIX/etc/conda/activate.d"
DEACTIVATE_DIR="$ENV_PREFIX/etc/conda/deactivate.d"
mkdir -p "$ACTIVATE_DIR" "$DEACTIVATE_DIR"

cat > "$ACTIVATE_DIR/msgflex_env_vars.sh" << EOF
#!/bin/sh
# Auto-generated by MSGFLEX install.sh — do not edit manually
export ${TOOLS_VAR}="${TOOLS_DIR}"
EOF

cat > "$DEACTIVATE_DIR/msgflex_env_vars.sh" << EOF
#!/bin/sh
# Auto-generated by MSGFLEX install.sh — do not edit manually
unset ${TOOLS_VAR}
EOF

chmod +x "$ACTIVATE_DIR/msgflex_env_vars.sh"
chmod +x "$DEACTIVATE_DIR/msgflex_env_vars.sh"

# Export for the rest of this script session too
export "${TOOLS_VAR}"="${TOOLS_DIR}"

success "$TOOLS_VAR=${TOOLS_DIR}"
success "Will be set automatically on every: conda activate $ENV_NAME"

# =============================================================================
# Step 6 — Shell launchers (optional convenience for when env isn't active)
# =============================================================================
header "Step 6/6  Shell launchers"

BIN_DIR="$HOME/.local/bin"
mkdir -p "$BIN_DIR"

_write_launcher() {
    local cmd="$1"
    local launcher="$BIN_DIR/$cmd"
    cat > "$launcher" << EOF
#!/usr/bin/env bash
# $cmd launcher — auto-generated by MSGFLEX install.sh
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"
exec ${cmd} "\$@"
EOF
    chmod +x "$launcher"
    success "Shell launcher created: $launcher"
}

_write_launcher "msgflex"
#_write_launcher "msgflex-check"
#_write_launcher "msgflex-conventional"
#_write_launcher "msgflex-binning"

if [[ ":$PATH:" != *":$BIN_DIR:"* ]]; then
    warn "$BIN_DIR is not in your PATH."
    warn "Add this to ~/.bashrc or ~/.zshrc:"
    warn "  export PATH=\"\$HOME/.local/bin:\$PATH\""
fi

# =============================================================================
# Verify installation with preflight
# =============================================================================
header "Running preflight check"

# `msgflex check` runs the comprehensive preflight validator
info "Running: msgflex check"
if msgflex check; then
    VERIFY_STATUS="${GREEN}PASSED${RESET}"
else
    # Exit code 1 from check = some FAILs found (e.g., mono not in conda yet)
    # We still report the install as "completed" but flag it
    VERIFY_STATUS="${YELLOW}PARTIAL (some preflight checks failed)${RESET}"
fi

# =============================================================================
# Summary
# =============================================================================
echo ""
echo -e "${BOLD}══════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Installation complete  —  status: $VERIFY_STATUS${BOLD}${RESET}"
echo -e "${BOLD}══════════════════════════════════════════════════${RESET}"
echo ""
echo "  Conda environment : $ENV_NAME"
echo "  Tools directory   : $TOOLS_DIR"
echo "  $TOOLS_VAR set automatically on activate"
echo ""
echo "  How to use:"
echo ""
echo -e "  ${BOLD}Activate the environment first:${RESET}"
echo "    conda activate $ENV_NAME"
echo ""
echo -e "  ${BOLD}Or use the shell launchers (if ~/.local/bin is in PATH):${RESET}"
echo "    msgflex --help"
echo ""
echo -e "  ${BOLD}Preflight check (validates everything is ready):${RESET}"
echo "    msgflex check                          # env + tools only"
echo "    msgflex check -b /path/to/project       # + project layout"
echo ""
echo -e "  ${BOLD}Subcommands:${RESET}"
echo "    msgflex sparx  -b /project -c params.txt [--mode conventional|binning]"
echo "    msgflex xpectra   -b /project [--ensemble 5] [--folds 5]"
echo "    msgflex quantix   -b /project [--score-field MSMSScore --threshold 10]"
echo "    msgflex run       -b /project -c params.txt [--rescore]"
echo ""
if [[ "$GUI_INSTALLED" == true ]]; then
    echo "  GUI               : installed (msgflex-gui)"
else
    echo "  GUI               : not installed (pip install 'msgflex[gui]')"
fi
VERSION="$(sed -n 's/^__version__ *= *"\(.*\)"/\1/p' "$SCRIPT_DIR/msgflex/__init__.py")"
echo "  MSGFLEX Pipeline v${VERSION} Installer"
echo ""
