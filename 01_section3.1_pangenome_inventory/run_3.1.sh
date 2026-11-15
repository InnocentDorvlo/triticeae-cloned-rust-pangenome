#!/usr/bin/env bash
# run_3.1.sh - Section 3.1 pipeline. Every step is resumable: just re-run after an interruption.
#   bash run_3.1.sh              # all steps, then Heaps' law
#   bash run_3.1.sh genome       # one step only: setup | prep | genome | annot | loci | orthofinder | stats | heaps
set -eo pipefail   # no -u: conda activation scripts use unset variables
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$STAGE_DIR/.." && pwd)"
cd "$ROOT"
STEP="${1:-all}"
CONDA="${CONDA:-$(command -v conda 2>/dev/null || true)}"
for c in "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/conda" "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" \
         /opt/miniforge3/bin/conda /opt/miniconda3/bin/conda /usr/local/Caskroom/miniforge/base/bin/conda; do
  [ -z "$CONDA" ] && [ -x "$c" ] && CONDA="$c"
done
[ -z "$CONDA" ] && { echo "conda not found - install Miniforge (README step 1)"; exit 1; }
eval "$("$CONDA" shell.bash hook)"
conda activate rust31 || { echo "env rust31 missing - run: conda env create -f env_3.1.yml"; exit 1; }
[ -x "$CONDA_PREFIX/lib/jvm/bin/java" ] && export JAVA_HOME="$CONDA_PREFIX/lib/jvm"
caffeinate -i -s -w $$ 2>/dev/null &        # keep the Mac awake while this script runs
mkdir -p logs
run_step() {
  if [ "$1" = "heaps" ]; then Rscript "$STAGE_DIR/heaps_3.1.R" results 1000 1000 2>&1 | tee -a logs/heaps.log
  else python "$STAGE_DIR/rust31.py" "$1" 2>&1 | tee -a "logs/$1.log"; fi
}
if [ "$STEP" = "all" ]; then
  for s in setup prep genome annot loci orthofinder stats heaps; do echo "===== $s ====="; run_step "$s"; done
else
  run_step "$STEP"
fi
echo "finished: $STEP"
