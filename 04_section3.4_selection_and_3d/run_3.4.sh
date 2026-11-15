#!/usr/bin/env bash
# run_3.4.sh - runs the Section 3.4 & Methods 5.6 pipeline (rust34.py + Fig. 4 + tables).
#
#   bash run_3.4.sh status      # show which Section 3.4 outputs exist
#   bash run_3.4.sh             # run all Section 3.4 steps + Fig. 4
#   bash run_3.4.sh align       # run one step: align | trees | codeml | hyphy | domains | asymmetry | structures | stats | figure
#   bash run_3.4.sh --force all # force re-run
set -eo pipefail
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$STAGE_DIR/.." && pwd)"
cd "$ROOT"
ROOT="$PWD"; W="$ROOT/work"; R="$ROOT/results"; L="$ROOT/logs"
mkdir -p "$L" "$R" "$ROOT/figures" "$ROOT/tables"

FORCE=0
[ "${1:-}" = "--force" ] && { FORCE=1; shift; }
STEP="${1:-all}"
THREADS="${THREADS:-4}"

say()  { printf '%s %s\n' "$(date +'[%H:%M:%S]')" "$*"; }
have() { [ -s "$1" ]; }
mark() { have "$1" && echo "  [done]    $2" || echo "  [missing] $2"; }

status() {
  echo "--- Section 3.4 state in $ROOT"
  mark "$R/Supplementary_Table_S10a_codeml_m8.tsv"           "1 codeml (PAML M1a/M2a, M7/M8 & BEB sites)     (S10a)"
  mark "$R/Supplementary_Table_S10b_hyphy_selection.tsv"     "2 hyphy (FUBAR purifying & MEME episodic)      (S10b)"
  mark "$R/Supplementary_Table_S10c_domain_partitioned_omega.tsv" "3 domains (domain-partitioned omega across classes) (S10c)"
  mark "$R/Supplementary_Table_S10d_tkp_module_asymmetry.tsv" "4 asymmetry (Class III TKP active vs decoy asymmetry) (S10d)"
  mark "$R/Supplementary_Table_S10e_af_surface_mapping.tsv" "5 structures (AlphaFold coordinates & DSSP RSA) (S10e)"
  mark "$R/section3.4_values.tsv"                            "6 stats (all Section 3.4 computed values)      (section3.4_values.tsv)"
  mark "$ROOT/figures/Fig4_selection_structures.pdf"         "8 figure (Fig. 4a-g multi-panel plot)          (figures/)"
  mark "$R/Supplementary_Table_S10f_R_statistical_audit_3.4.tsv" "9 r_audit (R statistical audit & R Fig. 4)     (S10f, Fig4_R)"
}
[ "$STEP" = "status" ] && { status; exit 0; }

# Activate conda env rust31
CONDA="${CONDA:-$(command -v conda 2>/dev/null || true)}"
for c in "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/conda" "$HOME/miniconda3/bin/conda" \
         "$HOME/anaconda3/bin/conda" /opt/miniforge3/bin/conda /opt/miniconda3/bin/conda; do
  [ -z "$CONDA" ] && [ -x "$c" ] && CONDA="$c"
done
[ -z "$CONDA" ] && { echo "ERROR: conda not found"; exit 1; }
eval "$("$CONDA" shell.bash hook)"
conda activate rust31

caffeinate -i -s -w $$ 2>/dev/null &

run_one() {
  local s="$1"
  case "$s" in
    codeml)      [ "$FORCE" = 0 ] && have "$R/Supplementary_Table_S10a_codeml_m8.tsv"           && { say "codeml already done - skipped";      return; } ;;
    hyphy)       [ "$FORCE" = 0 ] && have "$R/Supplementary_Table_S10b_hyphy_selection.tsv"     && { say "hyphy already done - skipped";       return; } ;;
    domains)     [ "$FORCE" = 0 ] && have "$R/Supplementary_Table_S10c_domain_partitioned_omega.tsv" && { say "domains already done - skipped"; return; } ;;
    asymmetry)   [ "$FORCE" = 0 ] && have "$R/Supplementary_Table_S10d_tkp_module_asymmetry.tsv" && { say "asymmetry already done - skipped"; return; } ;;
    structures)  [ "$FORCE" = 0 ] && have "$R/Supplementary_Table_S10e_af_surface_mapping.tsv" && { say "structures already done - skipped"; return; } ;;
    stats)       [ "$FORCE" = 0 ] && have "$R/section3.4_values.tsv"                            && { say "stats already done - skipped";       return; } ;;
    figure)      [ "$FORCE" = 0 ] && have "$ROOT/figures/Fig4_selection_structures.pdf"         && { say "figure already done - skipped";      return; } ;;
  esac
  say "===== rust34 $s ====="
  python "$STAGE_DIR/rust34.py" "$s" --threads "$THREADS" 2>&1 | tee -a "$L/rust34_${s}.log"
}

if [ "$STEP" = "all" ]; then
  for s in align trees codeml hyphy domains asymmetry structures stats figure; do
    run_one "$s"
  done
  say "===== R statistical audit & R Fig. 4 ====="
  Rscript "$STAGE_DIR/stats_and_fig4_3.4.R" 2>&1 | tee -a "$L/rust34_r_audit.log"
  python "$ROOT/06_tables_and_downstream/make_tables.py"
  status
else
  run_one "$STEP"
fi
