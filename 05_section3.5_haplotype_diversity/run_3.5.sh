#!/usr/bin/env bash
# run_3.5.sh - runs the Section 3.5 & Methods 5.7 pipeline (rust35.py + Fig. 5 + tables).
#
#   bash run_3.5.sh status      # show which Section 3.5 outputs exist
#   bash run_3.5.sh             # run all Section 3.5 steps + Fig. 5
#   bash run_3.5.sh extract     # run one step: extract | popgen | networks | diversity | mining | geography | bottleneck | stats | figure
#   bash run_3.5.sh --force all # force re-run
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
  echo "--- Section 3.5 state in $ROOT"
  mark "$R/Supplementary_Table_S11a_population_genetics_tiers.tsv"       "1 popgen (pi, theta_W, Tajima D & FST)          (S11a)"
  mark "$R/Supplementary_Table_S11b_haplotype_diversity_summary.tsv"    "2 diversity (sharing & privacy decomposition)    (S11b)"
  mark "$R/Supplementary_Table_S11c_candidate_functional_alleles.tsv"   "3 mining (motif-filtered candidate alleles)      (S11c)"
  mark "$R/Supplementary_Table_S11d_watkins_geographic_distribution.tsv" "4 geography (Watkins regional hotspots)         (S11d)"
  mark "$R/Supplementary_Table_S11e_bottleneck_tests.tsv"               "5 bottleneck (Wilcoxon & Tajima D tests)        (S11e)"
  mark "$R/section3.5_values.tsv"                                        "6 stats (all Section 3.5 computed values)      (section3.5_values.tsv)"
  mark "$ROOT/figures/Fig5_haplotypes_bottlenecks.pdf"                   "8 figure (Fig. 5a-f multi-panel plot)          (figures/)"
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
    popgen)      [ "$FORCE" = 0 ] && have "$R/Supplementary_Table_S11a_population_genetics_tiers.tsv"       && { say "popgen already done - skipped";      return; } ;;
    diversity)   [ "$FORCE" = 0 ] && have "$R/Supplementary_Table_S11b_haplotype_diversity_summary.tsv"    && { say "diversity already done - skipped";   return; } ;;
    mining)      [ "$FORCE" = 0 ] && have "$R/Supplementary_Table_S11c_candidate_functional_alleles.tsv"   && { say "mining already done - skipped";      return; } ;;
    geography)   [ "$FORCE" = 0 ] && have "$R/Supplementary_Table_S11d_watkins_geographic_distribution.tsv" && { say "geography already done - skipped";   return; } ;;
    bottleneck)  [ "$FORCE" = 0 ] && have "$R/Supplementary_Table_S11e_bottleneck_tests.tsv"               && { say "bottleneck already done - skipped";  return; } ;;
    stats)       [ "$FORCE" = 0 ] && have "$R/section3.5_values.tsv"                                        && { say "stats already done - skipped";       return; } ;;
    figure)      [ "$FORCE" = 0 ] && have "$ROOT/figures/Fig5_haplotypes_bottlenecks.pdf"                   && { say "figure already done - skipped";      return; } ;;
  esac
  say "===== rust35 $s ====="
  python "$STAGE_DIR/rust35.py" "$s" --threads "$THREADS" 2>&1 | tee -a "$L/rust35_${s}.log"
}

if [ "$STEP" = "all" ]; then
  for s in extract popgen networks diversity mining geography bottleneck stats figure; do
    run_one "$s"
  done
  status
else
  run_one "$STEP"
fi
