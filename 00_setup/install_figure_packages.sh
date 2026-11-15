#!/usr/bin/env bash
# ==============================================================================
# Script: install_figure_packages.sh
# Purpose: One-click local installer for all R figure packages and PyMOL check
# ==============================================================================

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "${SCRIPT_DIR}"

echo "======================================================================="
echo "   Setting Up Publication Figure Environment (R & PyMOL)"
echo "======================================================================="

# Locate Rscript
RSCRIPT_BIN=""
if [ -x "/usr/local/bin/Rscript" ]; then
    RSCRIPT_BIN="/usr/local/bin/Rscript"
elif command -v Rscript &> /dev/null; then
    RSCRIPT_BIN="$(command -v Rscript)"
else
    echo "[ERROR] Rscript not found. Please ensure R is installed."
    exit 1
fi

echo "[INFO] Using Rscript at: ${RSCRIPT_BIN}"
"${RSCRIPT_BIN}" --version

echo ""
echo "[INFO] Running R package installation..."
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/install_figure_packages.R"

echo ""
echo "======================================================================="
echo "   Checking PyMOL Environment"
echo "======================================================================="

PYMOL_BIN=""
if command -v pymol &> /dev/null; then
    PYMOL_BIN="$(command -v pymol)"
elif [ -x "$HOME/miniconda3/envs/pymol_env/bin/pymol" ]; then
    PYMOL_BIN="$HOME/miniconda3/envs/pymol_env/bin/pymol"
elif [ -x "$HOME/miniforge3/envs/pymol_env/bin/pymol" ]; then
    PYMOL_BIN="$HOME/miniforge3/envs/pymol_env/bin/pymol"
fi

if [ -n "${PYMOL_BIN}" ]; then
    echo "[OK] PyMOL detected at: ${PYMOL_BIN}"
    "${PYMOL_BIN}" -c -q -d "print('>> PyMOL is operational in headless mode <<')"
else
    echo "[WARNING] PyMOL binary not found in standard paths."
    echo "          PyMOL is required for 3D AlphaFold ray-traced rendering (Figure 4)."
fi

echo ""
echo "======================================================================="
echo "   Environment Setup Complete!"
echo "======================================================================="
