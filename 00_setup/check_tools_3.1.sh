#!/usr/bin/env bash
# check_tools_3.1.sh - checks that everything needed for the Section 3.1 pipeline is installed.
# Run from the rust31 folder:   bash check_tools_3.1.sh
# Optional:  ENV_NAME=myenv bash check_tools_3.1.sh     (default env name: rust31)
#            CONDA=/path/to/conda bash check_tools_3.1.sh
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$STAGE_DIR/.." && pwd)"
cd "$ROOT"
ENV_NAME="${ENV_NAME:-rust31}"
NFAIL=0; NWARN=0
ok()   { printf "  [ OK ]  %-28s %s\n" "$1" "$2"; }
bad()  { printf "  [MISS]  %-28s %s\n" "$1" "$2"; NFAIL=$((NFAIL+1)); }
warn() { printf "  [WARN]  %-28s %s\n" "$1" "$2"; NWARN=$((NWARN+1)); }

echo "== 1. Computer"
ARCH=$(uname -m); OSV=$(sw_vers -productVersion 2>/dev/null || uname -sr)
ok "macOS / CPU" "$OSV / $ARCH"
CORES=$(sysctl -n hw.ncpu 2>/dev/null || echo "?"); MEM=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1073741824 ))
ok "CPU cores / RAM" "$CORES cores / ${MEM} GB"
[ "$MEM" -eq 0 ] || [ "$MEM" -ge 16 ] || warn "RAM" "16 GB recommended (largest wheat chromosome ~0.9 Gb)"
FREE=$(df -g . | awk 'NR==2{print $4}')
if [ "${FREE:-0}" -ge 25 ]; then ok "free disk here" "${FREE} GB (need >= 25 GB)"; else bad "free disk here" "${FREE} GB - free up space (need >= 25 GB)"; fi
for t in curl gzip tar; do
  if command -v $t >/dev/null; then ok "$t" "$(command -v $t)"; else bad "$t" "install Xcode command line tools: xcode-select --install"; fi
done

echo "== 2. Pipeline files in this folder"
for f in config_3.1.json panel_3.1.tsv catalogue_3.1.tsv queries_3.1.faa env_3.1.yml \
         01_section3.1_pangenome_inventory/rust31.py 01_section3.1_pangenome_inventory/heaps_3.1.R 01_section3.1_pangenome_inventory/run_3.1.sh; do
  if [ -s "$f" ]; then ok "$f" ""; else bad "$f" "missing - re-extract rust31_section3.1_scripts.tar.gz"; fi
done

echo "== 3. conda"
if [ -z "$CONDA" ]; then
  CONDA=$(command -v conda 2>/dev/null)
  for c in "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/conda" "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" \
           /opt/miniforge3/bin/conda /opt/miniconda3/bin/conda /opt/anaconda3/bin/conda \
           /usr/local/Caskroom/miniforge/base/bin/conda /usr/local/miniforge3/bin/conda; do
    [ -z "$CONDA" ] && [ -x "$c" ] && CONDA="$c"
  done
fi
if [ -z "$CONDA" ] || ! "$CONDA" --version >/dev/null 2>&1; then
  bad "conda" "not found/working. Install Miniforge (see README, step 1)"
  [ -e /opt/anaconda3/_conda ] && [ ! -e /opt/anaconda3/bin/conda ] && \
    warn "/opt/anaconda3" "incomplete Anaconda install (no bin/conda) - ignore it, use Miniforge"
  echo; echo "RESULT: $NFAIL missing, $NWARN warnings. Fix conda first, then run this check again."; exit 1
fi
ok "conda" "$CONDA ($("$CONDA" --version 2>&1))"
PREFIX="${ENV_PREFIX:-}"
[ -z "$PREFIX" ] && PREFIX=$("$CONDA" env list 2>/dev/null | awk -v n="$ENV_NAME" '$1==n {print $NF}')
if [ -z "$PREFIX" ]; then   # fallback: <conda root>/envs/<name>
  for d in "$(dirname "$(dirname "$CONDA")")/envs/$ENV_NAME" "$HOME/miniforge3/envs/$ENV_NAME" "$HOME/.conda/envs/$ENV_NAME"; do
    [ -z "$PREFIX" ] && [ -x "$d/bin/python" ] && PREFIX="$d"
  done
fi
if [ -z "$PREFIX" ] || [ ! -d "$PREFIX" ]; then
  bad "env '$ENV_NAME'" "not created yet:  $CONDA env create -f env_3.1.yml"
  echo; echo "RESULT: $NFAIL missing, $NWARN warnings. Create the env, then run this check again."; exit 1
fi
ok "env '$ENV_NAME'" "$PREFIX"
B="$PREFIX/bin"

echo "== 4. Bioinformatics tools (inside env '$ENV_NAME')"
[ -x "$B/python" ]   && ok "python" "$("$B/python" --version 2>&1)" || bad "python" "conda install -n $ENV_NAME python=3.11"
"$B/python" -c "import numpy, scipy; print(numpy.__version__, scipy.__version__)" >/dev/null 2>&1 \
  && ok "numpy + scipy" "$("$B/python" -c 'import numpy, scipy; print("numpy", numpy.__version__, "| scipy", scipy.__version__)')" \
  || bad "numpy + scipy" "conda install -n $ENV_NAME -c conda-forge numpy scipy"
[ -x "$B/miniprot" ]  && ok "miniprot" "$("$B/miniprot" --version 2>&1)"           || bad "miniprot" "conda install -n $ENV_NAME -c bioconda miniprot"
[ -x "$B/diamond" ]   && ok "diamond" "$("$B/diamond" version 2>&1 | head -1)"      || bad "diamond" "conda install -n $ENV_NAME -c bioconda diamond"
[ -x "$B/hmmscan" ] && [ -x "$B/hmmpress" ] && ok "hmmer (hmmscan, hmmpress)" "$("$B/hmmscan" -h 2>&1 | grep -m1 HMMER | sed 's/^# //')" \
  || bad "hmmer" "conda install -n $ENV_NAME -c bioconda hmmer=3.4"
[ -x "$B/bedtools" ]  && ok "bedtools" "$("$B/bedtools" --version 2>&1)"           || bad "bedtools" "conda install -n $ENV_NAME -c bioconda bedtools"
if [ -x "$B/orthofinder" ]; then
  OFV=$(env -u PYTHONSAFEPATH PATH="$B:$PATH" "$B/orthofinder" -h 2>&1 | grep -m1 -i "OrthoFinder version")
  [ -n "$OFV" ] && ok "orthofinder" "$OFV" || bad "orthofinder" "installed but does not start: conda install -n $ENV_NAME -c bioconda orthofinder=2.5.5 --force-reinstall"
  for d in mcl fastme; do [ -x "$B/$d" ] && ok "  $d (OrthoFinder dep.)" "" || bad "  $d" "conda install -n $ENV_NAME -c bioconda $d"; done
else bad "orthofinder" "conda install -n $ENV_NAME -c bioconda orthofinder=2.5.5"; fi

echo "== 5. Java (only for NLR-Annotator; off by default)"
J=""
for j in "$PREFIX/lib/jvm/bin/java" "$PREFIX/bin/java" "${JAVA_HOME:+$JAVA_HOME/bin/java}"; do
  [ -z "$J" ] && [ -n "$j" ] && [ -x "$j" ] && "$j" -version >/dev/null 2>&1 && J="$j"
done
if [ -n "$J" ]; then ok "java" "$("$J" -version 2>&1 | head -1)"
else
  NLR_ON=$(grep -A1 '"nlr_annotator"' config_3.1.json | grep -c '"enabled": true')
  if [ "$NLR_ON" -gt 0 ]; then bad "java" "NLR-Annotator is enabled: conda install -n $ENV_NAME -c conda-forge openjdk=17"
  else warn "java" "not found; fine while NLR-Annotator is disabled in config_3.1.json"; fi
fi

echo "== 6. R (Heaps' law script)"
R_BIN=""; [ -x "$B/Rscript" ] && R_BIN="$B/Rscript"
[ -z "$R_BIN" ] && command -v Rscript >/dev/null && R_BIN=$(command -v Rscript)
if [ -n "$R_BIN" ]; then
  ok "Rscript" "$R_BIN ($("$R_BIN" --version 2>&1 | head -1))"
  for pkg in micropan vegan; do
    "$R_BIN" -e "suppressMessages(library($pkg)); cat(as.character(packageVersion('$pkg')))" >/tmp/_r31 2>&1 \
      && ok "  R package $pkg" "$(tail -1 /tmp/_r31)" \
      || bad "  R package $pkg" "conda install -n $ENV_NAME -c conda-forge r-$pkg"
  done
else bad "Rscript" "conda install -n $ENV_NAME -c conda-forge r-base r-micropan r-vegan"; fi

echo "== 7. Internet access to the data sources"
net() { local code; code=$(curl -s -o /dev/null -w "%{http_code}" --max-time 30 "$2")
        case "$code" in 2*|3*) ok "$1" "HTTP $code";; *) bad "$1" "HTTP ${code:-none} - check internet / VPN / firewall: $2";; esac; }
net "Ensembl Plants FTP (EBI)" "https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-63/"
net "NCBI genomes FTP"         "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/057/428/605/"
net "NCBI Datasets API"        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_018294505.1/dataset_report"
net "Pfam (EBI)"               "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt"
# GitHub is only needed for NLR-Annotator (off by default)
NLR_ON=$(grep -A1 '"nlr_annotator"' config_3.1.json | grep -c '"enabled": true')
GH1=$(curl -s -o /dev/null -w "%{http_code}" --max-time 30 "https://raw.githubusercontent.com/steuernb/NLR-Annotator/NLR-Annotator-2/readme.md")
GH2=$(curl -sL -o /dev/null -w "%{http_code}" --max-time 30 "https://github.com/steuernb/NLR-Annotator/raw/NLR-Annotator-2/readme.md")
if [ -s work/db/NLR-Annotator-v2.1b.jar ] && [ -s work/db/mot.txt ] && [ -s work/db/store.txt ]; then
  ok "NLR-Annotator files" "already in work/db"
elif echo "$GH1 $GH2" | grep -qE "(^| )(2|3)[0-9][0-9]"; then
  ok "GitHub (NLR-Annotator)" "raw=$GH1 github=$GH2"
elif [ "$NLR_ON" -gt 0 ]; then
  bad "GitHub (NLR-Annotator)" "unreachable (raw=$GH1 github=$GH2). Download the 3 files by browser (README section 5)"
else
  warn "GitHub (NLR-Annotator)" "unreachable, but NLR-Annotator is OFF in config_3.1.json - not needed"
fi

echo "== 8. Quick function test (miniprot + DIAMOND + hmmscan on 2 query proteins)"
T="$PWD/_check_tmp"; rm -rf "$T"; mkdir -p "$T"
"$B/python" - "$T" <<'EOF' 2>/dev/null
import sys, random
random.seed(1); T = sys.argv[1]
seqs, name = {}, None
for l in open("queries_3.1.faa"):
    if l.startswith(">"): name = l[1:].split()[0]; seqs[name] = ""
    else: seqs[name] += l.strip()
q = {k: seqs[k] for k in ("Sr35", "Lr67")}
open(f"{T}/q.faa", "w").write("".join(f">{k}\n{v}\n" for k, v in q.items()))
# back-translate Sr35 into a fake chromosome with one intron
cod = {"A":"GCT","R":"CGT","N":"AAT","D":"GAT","C":"TGT","Q":"CAA","E":"GAA","G":"GGT","H":"CAT","I":"ATT","L":"CTT",
       "K":"AAA","M":"ATG","F":"TTT","P":"CCT","S":"TCT","T":"ACT","W":"TGG","Y":"TAT","V":"GTT","X":"GCT"}
cds = "".join(cod[a] for a in q["Sr35"]) + "TAA"
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
g = rnd(5000) + cds[:900] + "GTAAGT" + rnd(1200) + "TTTCAG" + cds[900:] + rnd(5000)
open(f"{T}/g.fa", "w").write(">1A\n" + g + "\n")
EOF
if [ -s "$T/g.fa" ]; then
  MP=$("$B/miniprot" -t 2 -G 30k --gff --trans "$T/g.fa" "$T/q.faa" 2>/dev/null | grep -c $'\tmRNA\t')
  [ "$MP" = "1" ] && ok "miniprot spliced alignment" "Sr35 found across an intron" || bad "miniprot test" "expected 1 alignment, got $MP"
  "$B/diamond" makedb --in "$T/q.faa" -d "$T/q" --quiet >/dev/null 2>&1 && \
  DM=$("$B/diamond" blastp -q "$T/q.faa" -d "$T/q" --ultra-sensitive --quiet -f 6 2>/dev/null | wc -l | tr -d ' ')
  [ "${DM:-0}" -ge 2 ] && ok "DIAMOND ultra-sensitive" "$DM hits" || bad "DIAMOND test" "no hits"
  HMM=$(ls work/db/diag.hmm.h3m 2>/dev/null)
  if [ -n "$HMM" ]; then
    "$B/hmmscan" --domtblout "$T/d.tbl" -o /dev/null work/db/diag.hmm "$T/q.faa" >/dev/null 2>&1
    grep -q "NB-ARC" "$T/d.tbl" && ok "hmmscan (Pfam subset)" "NB-ARC found in Sr35" || bad "hmmscan test" "NB-ARC not found"
  else warn "hmmscan test" "Pfam subset not built yet (made by: python 01_section3.1_pangenome_inventory/rust31.py setup)"; fi
else bad "function test" "could not build test files (python in env broken?)"; fi
rm -rf "$T"

echo
if [ "$NFAIL" -eq 0 ]; then echo "RESULT: all required tools OK ($NWARN warnings). You can run:  bash 01_section3.1_pangenome_inventory/run_3.1.sh"
else echo "RESULT: $NFAIL item(s) missing, $NWARN warning(s). Fix the [MISS] lines above, then run this check again."; fi
exit $NFAIL
