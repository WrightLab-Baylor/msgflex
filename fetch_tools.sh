#!/usr/bin/env bash
# =============================================================================
# fetch_tools.sh downloads and extracts the MSGFLEX tools/ bundle
# Used by install.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_URL="https://github.com/thulasis/msgflex/releases/latest/download/tools.tar.gz"
TOOLS_ARCHIVE="/tmp/msgflex-tools.tar.gz"

if [[ -d "$SCRIPT_DIR/tools" ]]; then
    echo "[OK]    tools/ already present at $SCRIPT_DIR/tools — skipping download"
    exit 0
fi

echo "[INFO]  Downloading tools bundle from $TOOLS_URL ..."
if command -v curl &>/dev/null; then
    curl -fL "$TOOLS_URL" -o "$TOOLS_ARCHIVE"
elif command -v wget &>/dev/null; then
    wget -q "$TOOLS_URL" -O "$TOOLS_ARCHIVE"
else
    echo "[ERROR] Neither curl nor wget found. Install one and retry." >&2
    exit 1
fi

[[ -s "$TOOLS_ARCHIVE" ]] || { echo "[ERROR] Download failed or returned an empty file. Check $TOOLS_URL" >&2; exit 1; }

echo "[INFO]  Extracting to $SCRIPT_DIR ..."
tar -xzf "$TOOLS_ARCHIVE" -C "$SCRIPT_DIR"
rm -f "$TOOLS_ARCHIVE"
echo "[OK]    tools/ downloaded and extracted"
