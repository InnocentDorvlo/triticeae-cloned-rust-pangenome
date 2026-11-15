#!/usr/bin/env bash
# run_meme_3.4.sh - runs HyPhy MEME (episodic diversifying selection) on every Section 3.4 orthogroup.
# Usage: bash run_meme_3.4.sh [n_parallel]
set -eo pipefail
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$STAGE_DIR/.."
HYPHY="${HYPHY:-$(command -v hyphy || true)}"
if [[ -z "$HYPHY" || ! -x "$HYPHY" ]]; then
  for cand in "$HOME/miniforge3/envs/rust31/bin/hyphy" "$HOME/miniconda3/envs/rust31/bin/hyphy"; do
    [[ -x "$cand" ]] && { HYPHY="$cand"; break; }
  done
fi
NPAR="${1:-4}"
run_one() {
  d="$1"
  [ -s "$d/meme.json" ] && { echo "skip $(basename "$d")"; return; }
  "$HYPHY" meme --alignment "$d/aln.fna" --tree "$d/tree.nwk" --output "$d/meme.json" > "$d/meme.log" 2>&1 \
    && echo "done $(basename "$d")" || echo "FAILED $(basename "$d")"
}
export -f run_one; export HYPHY
ls -d work/sel_3.4/OG*/ | sed 's#/$##' | xargs -P "$NPAR" -I{} bash -c 'run_one {}'
