#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# post-link.sh  —  Linux / macOS
# File:  recipes/msgflex/post-link.sh
#
# Called by conda after the package is extracted into the environment.
# Pip-installs dearpygui (GUI extra) because it has no conda-forge recipe.
#
# Design constraints:
#   • Must be silent on success (conda captures output to a log)
#   • Must NOT fail hard — a missing GUI is acceptable for headless installs
#   • Uses $PREFIX/bin/pip to stay inside the activated environment
# ══════════════════════════════════════════════════════════════════════════════
set -euo pipefail

PIP="${PREFIX}/bin/pip"

# Allow users to skip GUI install (e.g. HPC headless nodes)
if [[ "${MSGFLEX_SKIP_GUI:-0}" == "1" ]]; then
    echo "[msgflex post-link] MSGFLEX_SKIP_GUI=1 — skipping dearpygui install."
    exit 0
fi

echo "[msgflex post-link] Installing dearpygui (GUI extra) via pip..."
if "${PIP}" install --quiet "dearpygui>=2.0"; then
    echo "[msgflex post-link] dearpygui installed successfully."
else
    cat <<'EOF'
[msgflex post-link] WARNING: dearpygui could not be installed automatically.
  The CLI (msgflex) and all pipeline functions work without it.
  To enable the GUI later, run:
      pip install "dearpygui>=2.0"
  On headless systems, set MSGFLEX_SKIP_GUI=1 to suppress this warning.
EOF
fi
