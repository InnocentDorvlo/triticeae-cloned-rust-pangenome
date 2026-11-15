#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Reproducible Pipeline: Regenerate Statistical Audits, Figures 1-5,
# Standardized Tables, and Supplementary Workbook (.xlsx) from results/ datasets.
# ==============================================================================

STAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$STAGE_DIR/.." && pwd)"
cd "$PROJECT_DIR"

RSCRIPT_BIN="${RSCRIPT_BIN:-$(command -v Rscript || true)}"
if [[ -z "$RSCRIPT_BIN" || ! -x "$RSCRIPT_BIN" ]]; then
  for cand in /usr/local/bin/Rscript "$HOME/miniforge3/envs/rust31/bin/Rscript" "$HOME/miniconda3/envs/rust31/bin/Rscript"; do
    [[ -x "$cand" ]] && { RSCRIPT_BIN="$cand"; break; }
  done
fi

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "$PYTHON_BIN" || ! -x "$PYTHON_BIN" ]]; then
  for cand in "$HOME/miniforge3/envs/rust31/bin/python3" "$HOME/miniconda3/envs/rust31/bin/python3" "$(command -v python3 || true)"; do
    [[ -n "$cand" && -x "$cand" ]] && { PYTHON_BIN="$cand"; break; }
  done
fi

echo "Using Rscript: $RSCRIPT_BIN"
echo "Using Python : $PYTHON_BIN"
echo "Working dir  : $PROJECT_DIR"
echo

echo "[1/8] Running R statistical audits (Sections 3.4 & 3.5)..."
"$RSCRIPT_BIN" 04_section3.4_selection_and_3d/stats_and_fig4_3.4.R
"$RSCRIPT_BIN" 05_section3.5_haplotype_diversity/stats_and_fig5_3.5.R

echo "[2/8] Generating Figure 1 (Pangenome inventory & CNV waterfall)..."
"$RSCRIPT_BIN" 01_section3.1_pangenome_inventory/make_unified_master_fig1.R

echo "[3/8] Generating Figure 2 (Chromosomal position, triads, SRI & TE)..."
"$RSCRIPT_BIN" 02_section3.2_chromosomal_synteny/make_fig2_review.R

echo "[4/8] Generating Figure 3 (Tandem kinase origins & phase-0 junctions)..."
"$RSCRIPT_BIN" 03_section3.3_tandem_kinases/make_publication_fig3.R

echo "[5/8] Generating Figure 4 (Selection regimes & 3D AlphaFold mapping)..."
"$RSCRIPT_BIN" 04_section3.4_selection_and_3d/make_publication_fig4.R

echo "[6/8] Generating Figure 5 (Haplotype networks, pi_aa / theta_W, & candidates)..."
"$RSCRIPT_BIN" 05_section3.5_haplotype_diversity/make_publication_fig5.R

echo "[7/8] Syncing Table 1 and Supplementary Tables (results/ -> tables/)..."
"$PYTHON_BIN" 06_tables_and_downstream/make_tables.py

echo "[8/8] Generating Master Supplementary Excel Workbook (Supplementary_Data_Tables.xlsx)..."
"$PYTHON_BIN" 06_tables_and_downstream/build_complete_supplementary_workbook.py


echo
echo "Done! All figures (figures/), tables (tables/), and supplementary workbook (.xlsx) are up to date."
