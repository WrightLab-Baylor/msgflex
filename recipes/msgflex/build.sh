#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# Bioconda build script  —  Linux / macOS
# File:  recipes/msgflex/build.sh
# ══════════════════════════════════════════════════════════════════════════════
set -euo pipefail

TOOLS_DEST="${PREFIX}/share/msgflex/tools"
mkdir -p "${TOOLS_DEST}"

# ── 1. Install Python package ─────────────────────────────────────────────────
pushd src
"${PYTHON}" -m pip install . \
    --no-deps \
    --no-build-isolation \
    --ignore-installed \
    -vv
popd

# ── 2. Bundle MSGFPlus.jar ────────────────────────────────────────────────────
cp tools_src/msgfplus/MSGFPlus.jar "${TOOLS_DEST}/"

# ── 3. Bundle MASIC (.NET) ────────────────────────────────────────────────────
mkdir -p "${TOOLS_DEST}/MASIC"
unzip -q tools_src/masic/*.zip -d "${TOOLS_DEST}/MASIC" || \
    cp -r tools_src/masic/. "${TOOLS_DEST}/MASIC/"

# ── 4. Bundle PHRP (.NET) ────────────────────────────────────────────────────
mkdir -p "${TOOLS_DEST}/PHRP"
unzip -q tools_src/phrp/*.zip -d "${TOOLS_DEST}/PHRP" || \
    cp -r tools_src/phrp/. "${TOOLS_DEST}/PHRP/"

echo "Tool installation complete:"
echo "  ${TOOLS_DEST}/MSGFPlus.jar"
ls -lh "${TOOLS_DEST}/MSGFPlus.jar"
echo "  ${TOOLS_DEST}/MASIC/"
echo "  ${TOOLS_DEST}/PHRP/"
