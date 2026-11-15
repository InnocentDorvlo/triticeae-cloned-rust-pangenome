#!/usr/bin/env bash
# run_3.2.sh - runs the Section 3.2 & Methods 5.4 pipeline (rust32.py + Fig. 2 + tables).
#
#   bash run_3.2.sh status      # show which Section 3.2 outputs exist
#   bash run_3.2.sh             # run all Section 3.2 steps + update pruning anchor support + Fig. 2
#   bash run_3.2.sh positions   # run one step: positions | triads | sri | te | stats | figures
#   bash run_3.2.sh --force sri # force re-run of a step
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
  echo "--- Section 3.2 state in $ROOT"
  mark "$R/S4b_locusgroup_positions.tsv"         "1 positions (d_norm across 56 subgenomes)   (S4, S4b)"
  mark "$R/S5b_orthogroup_triad_summary.tsv"     "2 triads (hexaploid A:B:D & subgenome CNV)  (S5, S5b, S5d)"
  mark "$R/S6b_sri_pairwise.tsv"                 "3 sri (+/-500 kb Synteny Retention Index)   (S6a, S6b, S6c)"
  mark "$R/S7a_flanking_te_windows.tsv"          "4 te (+/-200 kb TE class profiling)         (S7a, S7b)"
  mark "$R/section3.2_values.tsv"                "5 stats (Section 3.2 summary metrics)       (section3.2_values.tsv)"
  mark "$ROOT/figures/Fig2_chromosomal_distribution.pdf" "7 figures (Fig. 2a-f, Suppl. Fig. 1 & 2)    (figures/)"
  mark "$R/S7c_R_statistical_audit_3.2.tsv"      "8 R audit & figures (stats_and_fig2_3.2.R)  (S7c, Fig2_R)"
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
    positions) [ "$FORCE" = 0 ] && have "$R/S4b_locusgroup_positions.tsv"     && { say "positions already done - skipped"; return; } ;;
    triads)    [ "$FORCE" = 0 ] && have "$R/S5b_orthogroup_triad_summary.tsv" && { say "triads already done - skipped";    return; } ;;
    sri)       [ "$FORCE" = 0 ] && have "$R/S6b_sri_pairwise.tsv"             && { say "sri already done - skipped";       return; } ;;
    te)        [ "$FORCE" = 0 ] && have "$R/S7a_flanking_te_windows.tsv"      && { say "te already done - skipped";        return; } ;;
    stats)     [ "$FORCE" = 0 ] && have "$R/section3.2_values.tsv"            && { say "stats already done - skipped";     return; } ;;
    figures)
      say "===== figures & tables (3.2) ====="
      # Now that S6b_sri_pairwise.tsv exists, re-run prune_paralogues.py to fill flanking-anchor support
      if have "$R/S6b_sri_pairwise.tsv"; then
        python "$ROOT/01_section3.1_pangenome_inventory/prune_paralogues.py" --threads "$THREADS" 2>&1 | tee -a "$L/prune.log"
      fi
      Rscript "$STAGE_DIR/stats_and_fig2_3.2.R" 2>&1 | tee -a "$L/fig2_R.log"
      Rscript "$STAGE_DIR/make_fig2_review.R" 2>&1 | tee -a "$L/fig2.log"
      python "$ROOT/06_tables_and_downstream/make_tables.py" 2>&1 | tee -a "$L/tables.log"
      return
      ;;
  esac
  say "===== rust32 $s ====="
  python "$STAGE_DIR/rust32.py" "$s" --threads "$THREADS" 2>&1 | tee -a "$L/rust32_${s}.log"
}

if [ "$STEP" = "all" ]; then
  for s in positions triads sri te stats figures; do
    run_one "$s"
  done
  status
else
  run_one "$STEP"
fi

