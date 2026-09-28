@echo off
:: ══════════════════════════════════════════════════════════════════════════
:: Bioconda build script  —  Windows
:: File:  recipes/msgflex/bld.bat
:: ══════════════════════════════════════════════════════════════════════════

SET TOOLS_DEST=%PREFIX%\share\msgflex\tools
IF NOT EXIST "%TOOLS_DEST%" MKDIR "%TOOLS_DEST%"

:: ── 1. Install Python package ─────────────────────────────────────────────
PUSHD src
"%PYTHON%" -m pip install . ^
    --no-deps ^
    --no-build-isolation ^
    --ignore-installed ^
    -vv
IF ERRORLEVEL 1 EXIT 1
POPD

:: ── 2. Bundle MSGFPlus.jar ────────────────────────────────────────────────
COPY "tools_src\msgfplus\MSGFPlus.jar" "%TOOLS_DEST%\" || EXIT 1

:: ── 3. Bundle MASIC ───────────────────────────────────────────────────────
IF NOT EXIST "%TOOLS_DEST%\MASIC" MKDIR "%TOOLS_DEST%\MASIC"
tar -xf tools_src\masic\*.zip -C "%TOOLS_DEST%\MASIC" 2>nul || ^
    XCOPY /E /I /Q "tools_src\masic" "%TOOLS_DEST%\MASIC\"
IF ERRORLEVEL 1 EXIT 1

:: ── 4. Bundle PHRP ───────────────────────────────────────────────────────
IF NOT EXIST "%TOOLS_DEST%\PHRP" MKDIR "%TOOLS_DEST%\PHRP"
tar -xf tools_src\phrp\*.zip -C "%TOOLS_DEST%\PHRP" 2>nul || ^
    XCOPY /E /I /Q "tools_src\phrp" "%TOOLS_DEST%\PHRP\"
IF ERRORLEVEL 1 EXIT 1

:: ── 5. Install dearpygui via pip (Windows post-link equivalent) ───────────
"%PREFIX%\Scripts\pip.exe" install --quiet "dearpygui>=2.0"
IF ERRORLEVEL 1 (
    ECHO [msgflex] WARNING: dearpygui could not be installed.
    ECHO           Run: pip install "dearpygui>=2.0" to enable the GUI.
)

ECHO Build complete.
EXIT 0
