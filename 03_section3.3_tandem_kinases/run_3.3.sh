#!/usr/bin/env bash
# run_3.3.sh - runs the Section 3.3 & Methods 5.5 pipeline (rust33.py + Fig. 3 + tables).
#
#   bash run_3.3.sh status      # show which Section 3.3 outputs exist
#   bash run_3.3.sh             # run all Section 3.3 steps + Fig. 3
#   bash run_3.3.sh domains     # run one step: domains | pairs | phylogeny | origins | exons | stats | figure
#   bash run_3.3.sh --force     # force re-run
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
  echo "--- Section 3.3 state in $ROOT"
  mark "$R/S8a_classIII_domains.tsv"             "1 domains (Pfam kinase & accessory dissection)   (S8a)"
  mark "$R/S8c_kinase_pair_divergence.tsv"       "2 pairs (Kinase 1 vs Kinase 2 divergence / dS)   (S8c)"
  mark "$R/S8b_kinase_phylogeny.tsv"             "3 phylogeny (ML kinase trees & sister-clades)    (S8b)"
  mark "$R/S8d_architecture_distribution.tsv"    "4 origins (fused vs single-kinase progenitors)   (S8d)"
  mark "$R/S8e_exon_intron_phases.tsv"           "5 exons (intron phases at domain boundaries)     (S8e)"
  mark "$R/section3.3_values.tsv"                "6 stats (all Section 3.3 computed values)        (section3.3_values.tsv)"
  mark "$ROOT/figures/Fig3_tandem_kinase_origin.pdf" "8 figure (Fig. 3a-e multi-panel plot)           (figures/)"
  mark "$R/S8f_R_statistical_audit_3.3.tsv"      "9 r_audit (R statistical audit & Fig3_R)         (S8f)"
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
    domains)   [ "$FORCE" = 0 ] && have "$R/S8a_classIII_domains.tsv"          && { say "domains already done - skipped";   return; } ;;
    pairs)     [ "$FORCE" = 0 ] && have "$R/S8c_kinase_pair_divergence.tsv"    && { say "pairs already done - skipped";     return; } ;;
    phylogeny) [ "$FORCE" = 0 ] && have "$R/S8b_kinase_phylogeny.tsv"          && { say "phylogeny already done - skipped"; return; } ;;
    origins)   [ "$FORCE" = 0 ] && have "$R/S8d_architecture_distribution.tsv" && { say "origins already done - skipped";   return; } ;;
    exons)     [ "$FORCE" = 0 ] && have "$R/S8e_exon_intron_phases.tsv"        && { say "exons already done - skipped";     return; } ;;
    stats)     [ "$FORCE" = 0 ] && have "$R/section3.3_values.tsv"             && { say "stats already done - skipped";     return; } ;;
    figure)    [ "$FORCE" = 0 ] && have "$ROOT/figures/Fig3_tandem_kinase_origin.pdf" && { say "figure already done - skipped"; return; } ;;
    r_audit)   [ "$FORCE" = 0 ] && have "$R/S8f_R_statistical_audit_3.3.tsv" && { say "r_audit already done - skipped"; return; }
               say "===== Rscript "$STAGE_DIR/stats_and_fig3_3.3.R" ====="
               Rscript "$STAGE_DIR/stats_and_fig3_3.3.R" 2>&1 | tee -a "$L/rust33_r_audit.log"
               return ;;
  esac
  say "===== rust33 $s ====="
  python "$STAGE_DIR/rust33.py" "$s" --threads "$THREADS" 2>&1 | tee -a "$L/rust33_${s}.log"
}

if [ "$STEP" = "all" ]; then
  for s in domains pairs phylogeny origins exons stats figure r_audit; do
    run_one "$s"
  done
  status
else
  run_one "$STEP"
fi
