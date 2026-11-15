#!/usr/bin/env python3
"""
rust33.py - Section 3.3 / Methods 5.5: origin of Class III tandem kinases (TKPs) and kinase fusions.

Every number written by this script is computed from local data in work/ and results/.
Nothing is typed in by hand.

Steps (run in order, or `all`):
  domains    split every Class III query and every retained pangenome homologue into
             Kinase 1 / Kinase 2 / accessory domains (Pfam hmmscan domtbl), and score the
             catalytic motifs (b3 VAIK lysine, HRD catalytic loop, DFG activation loop)
             -> results/S8a_classIII_domains.tsv
  pairs      Kinase 1 vs Kinase 2 divergence inside the same polypeptide
               * amino-acid identity for every pangenome TKP homologue (protein level)
               * Nei-Gojobori pS/pN with Jukes-Cantor dS/dN from real nucleotide CDS:
                 reference CDS (queries_3.1.cds.fna, NCBI) + genomic CDS extracted from any
                 genome FASTA still on disk (genomes/<asm>.fa.gz)
             -> results/S8c_kinase_pair_divergence.tsv
  phylogeny  MAFFT L-INS-i -> trimAl -automated1 -> IQ-TREE (-m MFP -B 1000 -alrt 1000) on
             dissected Kinase 1 / Kinase 2 / single-kinase progenitor / outgroup kinase domains
             (Brachypodium, Hordeum from work/tkp_3.3/outgroup_kinases.faa); separate K1-only and
             K2-only trees for the tanglegram; sister-clade test per TKP with Bio.Phylo
             -> results/S8b_kinase_phylogeny.tsv, work/tkp_3.3/*.treefile
  origins    per assembly x Class III gene: fused architecture vs unfused single-kinase progenitor
             vs standalone accessory-domain progenitor
             -> results/S8d_architecture_distribution.tsv
  exons      exon-intron structure and intron phase at the domain junctions (reference CDS from
             NCBI + best full-length Miniprot genomic model), plus flanking TE from S7a
             -> results/S8e_exon_intron_phases.tsv
  stats      all Section 3.3 values -> results/section3.3_values.tsv, tables/Section3.3_Summary_Statistics.md
  figure     figures/Fig3_tandem_kinase_origin.pdf/.png

Usage: python rust33.py all [--threads 4]
"""
import argparse, collections, glob, itertools, math, os, re, subprocess, sys
from pathlib import Path

import numpy as np
import scipy.stats as ss

STAGE_DIR = Path(__file__).resolve().parent
HERE = STAGE_DIR.parent
sys.path.insert(0, str(HERE / "01_section3.1_pangenome_inventory"))
os.environ["PATH"] = f"{Path(sys.executable).parent}:{os.environ.get('PATH', '')}"
from rust31 import log, load_cfg, read_tsv, write_tsv, read_fasta, catalogue, panel_rows, which_or_die

KIN = ("Pkinase", "PK_Tyr_Ser-Thr")
CLASS3 = ["Yr15", "Sr60", "Sr62", "Lr9", "Yr36", "Sr43", "Rpg1", "WTK4"]
DOM_E = 1e-3
MIN_KIN_LEN = 150          # aa; a kinase module shorter than this is not counted as a full domain

CODON = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ----------------------------------------------------------------------------- helpers
def fasta_dict(path):
    return {n: s for n, _, s in read_fasta(path)}


def cds_with_headers(path):
    out = {}
    for n, d, s in read_fasta(path):
        m = re.search(r"\[location=([^\]]+)\]", d)
        out[n] = (s.upper(), m.group(1) if m else "")
    return out


def translate(nt):
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def merge_kinase(ivs):
    """Merge overlapping Pkinase / PK_Tyr_Ser-Thr hits into kinase modules."""
    merged = []
    for s, e in sorted(ivs):
        if merged and min(e, merged[-1][1]) - max(s, merged[-1][0]) > 30:
            merged[-1] = (min(s, merged[-1][0]), max(e, merged[-1][1]))
        else:
            merged.append((s, e))
    return merged


def read_domtbl(path, keep=None):
    hits = collections.defaultdict(list)
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 21 or (keep is not None and f[3] not in keep):
            continue
        if float(f[12]) <= DOM_E:
            hits[f[3]].append((f[0], int(f[17]), int(f[18])))
    return hits


def architecture(dom_hits):
    """-> (list of kinase modules >= MIN_KIN_LEN, accessory domain names, accessory span)"""
    kin = [iv for iv in merge_kinase([(s, e) for d, s, e in dom_hits if d in KIN]) if iv[1] - iv[0] + 1 >= MIN_KIN_LEN]
    acc = [(d, s, e) for d, s, e in dom_hits if d not in KIN]
    names = sorted({d for d, _, _ in acc})
    span = (min(s for _, s, _ in acc), max(e for _, _, e in acc)) if acc else None
    return kin, names, span


def catalytic(seq):
    """b3-strand lysine (VAIK), catalytic-loop HRD, activation-loop DFG."""
    k = bool(re.search(r"[VIL]A[VILM]K", seq))
    hrd = bool(re.search(r"H[RKGN]D[LIVMFA][KR]", seq))
    dfg = bool(re.search(r"DFG", seq))
    return ("active" if (k and hrd and dfg) else "pseudokinase"), f"K={int(k)};HRD={int(hrd)};DFG={int(dfg)}"


def aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score, a.extend_gap_score = -10, -0.5
    a.target_end_gap_score = a.query_end_gap_score = 0
    return a


def align_pair(al, p1, p2):
    aln = al.align(p1, p2)[0]
    return str(aln[0]), str(aln[1])


def aa_identity(a1, a2):
    cols = [(x, y) for x, y in zip(a1, a2) if x != "-" and y != "-"]
    return (sum(x == y for x, y in cols) / len(cols) if cols else float("nan")), len(cols)


def _syn_sites(c):
    aa = CODON.get(c)
    if aa in (None, "*"):
        return None
    s = 0.0
    for i in range(3):
        s += sum(CODON.get(c[:i] + b + c[i + 1:]) == aa for b in "ACGT" if b != c[i]) / 3.0
    return s


def _diffs(c1, c2):
    pos = [i for i in range(3) if c1[i] != c2[i]]
    if not pos:
        return 0.0, 0.0
    tot_s = tot_n = 0.0
    n_ok = 0
    for path in itertools.permutations(pos):
        cur, s, n, ok = c1, 0, 0, True
        for i in path:
            nxt = cur[:i] + c2[i] + cur[i + 1:]
            if CODON.get(nxt, "*") == "*":
                ok = False
                break
            s += CODON[cur] == CODON[nxt]
            n += CODON[cur] != CODON[nxt]
            cur = nxt
        if ok:
            tot_s += s; tot_n += n; n_ok += 1
    if not n_ok:
        return len(pos) / 2.0, len(pos) / 2.0
    return tot_s / n_ok, tot_n / n_ok


def ng86(a1, a2, c1, c2):
    """Nei-Gojobori 1986 on a protein alignment back-translated with the two CDS.
    Returns pS, pN, dS, dN (JC); dS/dN = nan when p >= 0.75 (saturated)."""
    i1 = i2 = 0
    S = N = Sd = Nd = 0.0
    for x, y in zip(a1, a2):
        k1 = c1[3 * i1:3 * i1 + 3] if x != "-" else None
        k2 = c2[3 * i2:3 * i2 + 3] if y != "-" else None
        i1 += x != "-"
        i2 += y != "-"
        if not k1 or not k2:
            continue
        s1, s2 = _syn_sites(k1), _syn_sites(k2)
        if s1 is None or s2 is None:
            continue
        S += (s1 + s2) / 2.0
        N += 3.0 - (s1 + s2) / 2.0
        sd, nd = _diffs(k1, k2)
        Sd += sd; Nd += nd
    pS, pN = (Sd / S if S else float("nan")), (Nd / N if N else float("nan"))
    jc = lambda p: -0.75 * math.log(1 - 4 * p / 3) if p < 0.75 else float("nan")
    return pS, pN, jc(pS), jc(pN)


def mean_sem(x):
    x = [v for v in x if v == v]
    if not x:
        return float("nan"), float("nan"), 0
    return float(np.mean(x)), (float(ss.sem(x)) if len(x) > 1 else float("nan")), len(x)


def fmt_ms(x, d=3):
    m, s, n = mean_sem(x)
    return f"{m:.{d}f} +/- {s:.{d}f} (n = {n})"


# ----------------------------------------------------------------------------- shared loaders
class Ctx:
    def __init__(self, cfg):
        self.cfg = cfg
        self.W, self.R = Path(cfg["workdir"]), Path(cfg["resultsdir"])
        self.T = self.W / "tkp_3.3"
        self.T.mkdir(parents=True, exist_ok=True)
        (HERE / "tables").mkdir(exist_ok=True)
        (HERE / "figures").mkdir(exist_ok=True)
        self.cat = catalogue(cfg)
        self.panel = {p["asm_id"]: p for p in panel_rows(cfg)}

    def ref_arch(self):
        q = fasta_dict(self.cfg["queries"])
        dt = read_domtbl(self.W / "loci" / "queries.domtbl", set(CLASS3))
        out = {}
        for g in CLASS3:
            kin, acc, span = architecture(dt[g])
            out[g] = {"seq": q[g], "kin": kin, "acc": acc, "acc_span": span}
        return out

    def homologues(self):
        prune = {r["cand_id"] for r in read_tsv(self.R / "S9_paralogue_pruning.tsv") if r["verdict"] == "retained"} \
            if (self.R / "S9_paralogue_pruning.tsv").exists() else None
        return [r for r in read_tsv(self.R / "candidate_loci.tsv")
                if r["best_query"] in CLASS3 and r["anchor_ok"] == "1" and (prune is None or r["cand_id"] in prune)]


def miniprot_models(W, asm_filter=None, query_filter=CLASS3):
    """(asm, seqid, query) -> list of models {start,end,strand,qst,qen,qlen,ident,cds:[(s,e)]}"""
    models = collections.defaultdict(list)
    for gff in sorted(glob.glob(str(W / "miniprot" / "*" / "*.gff"))):
        asm = Path(gff).parent.name
        if asm_filter and asm not in asm_filter:
            continue
        paf, cur = None, None
        for line in open(gff):
            if line.startswith("##PAF\t"):
                p = line.rstrip("\n").split("\t")
                paf = (p[1], int(p[2]), int(p[3]), int(p[4]), int(p[10]) / max(1, int(p[11])))
                cur = None
            elif "\tmRNA\t" in line and paf:
                p = line.rstrip("\n").split("\t")
                if query_filter is None or paf[0] in query_filter:
                    cur = {"asm": asm, "seqid": p[0], "start": int(p[3]), "end": int(p[4]), "strand": p[6],
                           "query": paf[0], "qlen": paf[1], "qst": paf[2], "qen": paf[3], "ident": paf[4], "cds": []}
                    models[(asm, p[0], paf[0])].append(cur)
                else:
                    cur = None
            elif "\tCDS\t" in line and cur is not None:
                p = line.split("\t")
                cur["cds"].append((int(p[3]), int(p[4])))
    return models


def match_model(models, r):
    """Miniprot model that corresponds to a candidate locus (largest overlap)."""
    best, bo = None, 0
    s, e = int(r["start"]), int(r["end"])
    for m in models.get((r["asm_id"], r["seqid"], r["best_query"]), []):
        o = min(e, m["end"]) - max(s, m["start"])
        if o > bo:
            best, bo = m, o
    return best


# ----------------------------------------------------------------------------- 1 domains
def step_domains(ctx):
    log("=== rust33 step domains")
    ref = ctx.ref_arch()
    homs = ctx.homologues()
    cseq = fasta_dict(ctx.W / "loci" / "candidates.faa")
    dt = read_domtbl(ctx.W / "loci" / "candidates.domtbl", {r["cand_id"] for r in homs})
    rows = []

    def row(kind, rid, g, asm, sub, seq, kin, acc, span, extra):
        k1 = kin[0] if kin else None
        k2 = kin[1] if len(kin) > 1 else None
        st1, m1 = catalytic(seq[k1[0] - 1:k1[1]]) if k1 else ("-", "-")
        st2, m2 = catalytic(seq[k2[0] - 1:k2[1]]) if k2 else ("-", "-")
        if len(kin) >= 2:
            arch = "tandem kinase + accessory" if acc else "tandem kinase"
        elif len(kin) == 1:
            arch = "single kinase + accessory" if acc else "single kinase"
        else:
            arch = "accessory only" if acc else "partial kinase"
        d = {"Record": kind, "ID": rid, "Gene": g, "Class": "III", "Set": ctx.cat[g]["Set"],
             "Catalogue_architecture": ctx.cat[g]["Architecture"], "asm_id": asm, "subgenome": sub,
             "aa_len": len(seq), "n_kinase_modules": len(kin),
             "K1_range": f"{k1[0]}-{k1[1]}" if k1 else "-", "K1_status": st1, "K1_motifs": m1,
             "K2_range": f"{k2[0]}-{k2[1]}" if k2 else "-", "K2_status": st2, "K2_motifs": m2,
             "accessory_domains": ",".join(acc) if acc else "-",
             "accessory_range": f"{span[0]}-{span[1]}" if span else "-", "observed_architecture": arch}
        d.update(extra)
        rows.append(d)

    for g in CLASS3:
        a = ref[g]
        row("reference", g, g, "cloned_reference", "-", a["seq"], a["kin"], a["acc"], a["acc_span"],
            {"identity": "1.0", "qcov": "1.0", "intact": "1", "ploidy": "-"})
    for r in homs:
        seq = cseq.get(r["cand_id"], "")
        kin, acc, span = architecture(dt.get(r["cand_id"], []))
        row("homologue", r["cand_id"], r["best_query"], r["asm_id"], r["subgenome"], seq, kin, acc, span,
            {"identity": r["identity"], "qcov": r["qcov"], "intact": r["intact"],
             "ploidy": ctx.panel.get(r["asm_id"], {}).get("ploidy", "-")})
    cols = list(rows[0].keys())
    write_tsv(ctx.R / "S8a_classIII_domains.tsv", rows, cols)
    n_h = sum(r["Record"] == "homologue" for r in rows)
    log(f"domains done: 8 reference proteins + {n_h} retained Class III homologues")
    for g in CLASS3:
        a = ref[g]
        log(f"  {g:5s} kinase modules={a['kin']} accessory={a['acc'] or '-'}")


# ----------------------------------------------------------------------------- 2 pairs (K1 vs K2)
def extract_genomic_cds(genome_gz, wanted):
    """wanted: {key: (seqid, strand, [(s,e),...])} -> {key: CDS nt (5'->3')}. Streams the genome once."""
    need = collections.defaultdict(list)                 # seqid -> [(s, e, key, idx)]
    for key, (sid, strand, cds) in wanted.items():
        for idx, (s, e) in enumerate(cds):
            need[sid].append((s, e, key, idx))
    for sid in need:
        need[sid].sort()
    pieces = collections.defaultdict(dict)
    proc = subprocess.Popen(["gzip", "-dc", str(genome_gz)], stdout=subprocess.PIPE, bufsize=1 << 24)
    sid, pos, todo, ptr = None, 0, [], 0
    for line in proc.stdout:
        if line[:1] == b">":
            sid = line[1:].split()[0].decode()
            pos, todo, ptr = 0, need.get(sid, []), 0
            continue
        if not todo:
            continue
        L = len(line) - 1 if line.endswith(b"\n") else len(line)
        lo, hi = pos + 1, pos + L                          # 1-based inclusive coordinates of this line
        while ptr < len(todo) and todo[ptr][1] < lo:
            ptr += 1
        j = ptr
        while j < len(todo) and todo[j][0] <= hi:
            s, e, key, idx = todo[j]
            a, b = max(s, lo), min(e, hi)
            if a <= b:
                pieces[key].setdefault(idx, []).append(line[a - lo:b - lo + 1].decode())
            j += 1
        pos += L
    proc.wait()
    out = {}
    for key, (sid, strand, cds) in wanted.items():
        segs = ["".join(pieces[key].get(i, [])).upper() for i in range(len(cds))]
        if any(len(sg) != e - s + 1 for sg, (s, e) in zip(segs, cds)):
            continue
        if strand == "-":
            segs = [sg.translate(COMP)[::-1] for sg in segs]
        out[key] = "".join(segs)
    return out


def order_cds(m):
    cds = sorted(m["cds"])
    return cds if m["strand"] == "+" else cds[::-1]


def step_pairs(ctx):
    log("=== rust33 step pairs")
    al = aligner()
    ref = ctx.ref_arch()
    qcds = cds_with_headers(HERE / "queries_3.1.cds.fna")
    dom = read_tsv(ctx.R / "S8a_classIII_domains.tsv")
    cseq = fasta_dict(ctx.W / "loci" / "candidates.faa")
    rows = []

    # (a) reference TKPs: CDS from NCBI
    for g in CLASS3:
        a = ref[g]
        if len(a["kin"]) < 2:
            continue
        (s1, e1), (s2, e2) = a["kin"][0], a["kin"][1]
        p1, p2 = a["seq"][s1 - 1:e1], a["seq"][s2 - 1:e2]
        nt = qcds[g][0]
        if translate(nt).rstrip("*") != a["seq"]:
            log(f"  WARNING {g}: NCBI CDS does not translate to the query protein; dS skipped")
            nt = None
        a1, a2 = align_pair(al, p1, p2)
        ident, ncol = aa_identity(a1, a2)
        pS = pN = dS = dN = float("nan")
        if nt:
            pS, pN, dS, dN = ng86(a1, a2, nt[3 * (s1 - 1):3 * e1], nt[3 * (s2 - 1):3 * e2])
        rows.append({"Record": "reference", "ID": g, "Gene": g, "asm_id": "cloned_reference", "subgenome": "-",
                     "K1_range": f"{s1}-{e1}", "K2_range": f"{s2}-{e2}", "aligned_aa": ncol,
                     "aa_identity_pct": f"{100 * ident:.2f}", "cds_source": "NCBI CDS" if nt else "-",
                     "pS": f"{pS:.4f}", "pN": f"{pN:.4f}", "dS_JC": f"{dS:.4f}", "dN_JC": f"{dN:.4f}",
                     "dS_saturated": int(pS == pS and pS >= 0.75)})

    # (b) pangenome TKP homologues (>= 2 full kinase modules)
    tkp = [r for r in dom if r["Record"] == "homologue" and int(r["n_kinase_modules"]) >= 2]
    genomes = {a: p for a in {r["asm_id"] for r in tkp}
               for p in [HERE / "genomes" / f"{a}.fa.gz", ctx.W / "genomes" / f"{a}.fa.gz"] if p.exists()}
    nt_of = {}
    for cache in sorted(ctx.T.glob("genomic_cds_*.fna")):
        nt_of.update(fasta_dict(cache))
    if genomes:
        cand = {r["cand_id"]: r for r in read_tsv(ctx.R / "candidate_loci.tsv")}
        models = miniprot_models(ctx.W, set(genomes))
        for asm, gz in genomes.items():
            wanted = {}
            for r in tkp:
                if r["asm_id"] != asm:
                    continue
                m = match_model(models, cand[r["ID"]])
                if m:
                    wanted[r["ID"]] = (m["seqid"], m["strand"], order_cds(m))
            cache = ctx.T / f"genomic_cds_{asm}.fna"
            if cache.exists():
                got = fasta_dict(cache)
            else:
                log(f"  streaming {gz.name} once to extract {len(wanted)} TKP CDS (takes a few minutes)")
                got = extract_genomic_cds(gz, wanted)
                with open(cache, "w") as f:
                    for k, v in got.items():
                        f.write(f">{k}\n{v}\n")
            nt_of.update(got)
    if nt_of:
        log(f"  genomic CDS available for {len(nt_of)} TKP homologues from cached/local genome CDS")
    else:
        log("  no genome FASTA or cached genomic CDS on disk - pangenome dS skipped (protein identity still computed)")

    for r in tkp:
        seq = cseq[r["ID"]].rstrip("*")
        s1, e1 = map(int, r["K1_range"].split("-"))
        s2, e2 = map(int, r["K2_range"].split("-"))
        k1_aa, k2_aa = seq[s1 - 1:e1], seq[s2 - 1:e2]
        a1, a2 = align_pair(al, k1_aa, k2_aa)
        ident, ncol = aa_identity(a1, a2)
        pS = pN = dS = dN = float("nan")
        src = "-"
        nt = nt_of.get(r["ID"])
        if nt:
            tr = translate(nt).rstrip("*")
            i1, i2 = tr.find(k1_aa), tr.find(k2_aa)
            # usable when both K1 and K2 translate in-frame without frameshifts
            if i1 >= 0 and i2 >= 0:
                pS, pN, dS, dN = ng86(a1, a2, nt[3 * i1:3 * (i1 + len(k1_aa))], nt[3 * i2:3 * (i2 + len(k2_aa))])
                src = "genomic CDS"
        rows.append({"Record": "homologue", "ID": r["ID"], "Gene": r["Gene"], "asm_id": r["asm_id"],
                     "subgenome": r["subgenome"], "K1_range": r["K1_range"], "K2_range": r["K2_range"],
                     "aligned_aa": ncol, "aa_identity_pct": f"{100 * ident:.2f}", "cds_source": src,
                     "pS": f"{pS:.4f}", "pN": f"{pN:.4f}", "dS_JC": f"{dS:.4f}", "dN_JC": f"{dN:.4f}",
                     "dS_saturated": int(pS == pS and pS >= 0.75)})
    write_tsv(ctx.R / "S8c_kinase_pair_divergence.tsv", rows, list(rows[0].keys()))
    n_ds = sum(r["cds_source"] == "genomic CDS" for r in rows)
    log(f"pairs done: {len(rows)} Kinase1-Kinase2 pairs ({n_ds} pangenome pairs with genomic dS)")


# ----------------------------------------------------------------------------- 3 phylogeny
def _iqtree(T, prefix, fa, threads):
    tree = T / f"{prefix}.treefile"
    if tree.exists():
        return tree
    aln, trim = T / f"{prefix}.aln.faa", T / f"{prefix}.trim.faa"
    with open(aln, "w") as out:
        subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", str(threads), "--quiet", str(fa)],
                       stdout=out, check=True)
    subprocess.run(["trimal", "-in", str(aln), "-out", str(trim), "-automated1"], check=True)
    subprocess.run(["iqtree", "-s", str(trim), "-m", "MFP", "-B", "1000", "-alrt", "1000", "-T", str(threads),
                    "--prefix", str(T / prefix), "-redo", "--quiet"], check=True)
    return tree


def _support(clade):
    lab = clade.name or (str(clade.confidence) if clade.confidence is not None else "")
    parts = [p for p in str(lab).split("/") if p.replace(".", "", 1).isdigit()]
    return (parts[0], parts[-1]) if parts else ("NA", "NA")


def step_phylogeny(ctx, threads):
    log("=== rust33 step phylogeny")
    which_or_die("mafft", "trimal", "iqtree")
    from Bio import Phylo
    ref = ctx.ref_arch()
    dom = read_tsv(ctx.R / "S8a_classIII_domains.tsv")
    cseq = fasta_dict(ctx.W / "loci" / "candidates.faa")
    seqs, k1, k2 = {}, {}, {}

    for g in CLASS3:
        a = ref[g]
        for i, (s, e) in enumerate(a["kin"][:2], start=1):
            seqs[f"{g}_K{i}_ref"] = a["seq"][s - 1:e]
        if len(a["kin"]) >= 2:
            k1[f"{g}_ref"], k2[f"{g}_ref"] = seqs[f"{g}_K1_ref"], seqs[f"{g}_K2_ref"]

    # one best intact homologue per gene x ploidy tier (2x, 4x, 6x); best single-kinase progenitor per gene
    homs = sorted([r for r in dom if r["Record"] == "homologue"], key=lambda r: -float(r["identity"]))
    seen = set()
    for r in homs:
        g, pl = r["Gene"], r["ploidy"]
        nk = int(r["n_kinase_modules"])
        tag = f"{g}_{r['asm_id']}_{r['subgenome']}"
        if nk >= 2 and (g, pl, "tkp") not in seen and ref[g]["kin"] and len(ref[g]["kin"]) >= 2:
            seen.add((g, pl, "tkp"))
            s1, e1 = map(int, r["K1_range"].split("-")); s2, e2 = map(int, r["K2_range"].split("-"))
            seqs[f"{tag}_K1"], seqs[f"{tag}_K2"] = cseq[r["ID"]][s1 - 1:e1], cseq[r["ID"]][s2 - 1:e2]
            k1[tag], k2[tag] = seqs[f"{tag}_K1"], seqs[f"{tag}_K2"]
        elif nk == 1 and r["accessory_domains"] == "-" and (g, "single") not in seen:
            seen.add((g, "single"))
            s1, e1 = map(int, r["K1_range"].split("-"))
            seqs[f"{tag}_single"] = cseq[r["ID"]][s1 - 1:e1]

    og_fa, og_dt = ctx.T / "outgroup_kinases.faa", ctx.T / "outgroup_kinases.domtbl"
    outgroups = []
    if og_fa.exists() and og_dt.exists():
        ogs = fasta_dict(og_fa)
        dt = read_domtbl(og_dt)
        used = set()
        for acc, hits in dt.items():
            kin = [iv for iv in merge_kinase([(s, e) for d, s, e in hits if d in KIN]) if iv[1] - iv[0] + 1 >= MIN_KIN_LEN]
            if not kin:
                continue
            sq = ogs[acc][kin[0][0] - 1:kin[0][1]]
            if sq in used:
                continue
            used.add(sq)
            if acc.startswith(("Bdistachyon_", "Hvulgare_")):
                name = f"OUT_{acc.split('.')[0]}"
            else:
                sp = "Bdistachyon" if acc.startswith("XP_") else "Hvulgare"
                name = f"OUT_{sp}_{acc.split('.')[0]}"
            seqs[name] = sq
            outgroups.append(name)
    else:
        raise RuntimeError("Missing work/tkp_3.3/outgroup_kinases.faa or outgroup_kinases.domtbl; refusing midpoint rooting")

    for name, d in (("kinase_domains", seqs), ("kinase1_tangle", k1), ("kinase2_tangle", k2)):
        with open(ctx.T / f"{name}.faa", "w") as f:
            for k, v in d.items():
                f.write(f">{k}\n{v}\n")
    trees = {n: _iqtree(ctx.T, n, ctx.T / f"{n}.faa", threads) for n in ("kinase_domains", "kinase1_tangle", "kinase2_tangle")}

    tree = Phylo.read(trees["kinase_domains"], "newick")
    tree_outs = [t.name for t in tree.get_terminals() if t.name.startswith("OUT_")]
    if not tree_outs:
        raise RuntimeError("No usable outgroup tips were loaded; rooted origin inference is not allowed")
    basal_outs = [t for t in tree_outs if any(x in t for x in ("XP_044954759", "XP_044957438"))] or tree_outs
    try:
        tree.root_with_outgroup(*basal_outs)
    except Exception as e:
        raise RuntimeError(f"Outgroup rooting failed; refusing midpoint fallback: {e}") from e
    tips = [t.name for t in tree.get_terminals()]
    og_peers = {"Rpg1": {"Rpg1", "WTK4"}, "WTK4": {"Rpg1", "WTK4"}}
    ref_dom = {r["Gene"]: r for r in dom if r["Record"] == "reference"}
    hom_by_gene = collections.defaultdict(list)
    for r in dom:
        if r["Record"] == "homologue":
            hom_by_gene[r["Gene"]].append(r)
    pairs_path = ctx.R / "S8c_kinase_pair_divergence.tsv"
    ref_pairs = {r["Gene"]: r for r in read_tsv(pairs_path) if r["Record"] == "reference"} if pairs_path.exists() else {}
    rows = []
    for g in CLASS3:
        t1 = [t for t in tips if t.startswith(f"{g}_") and (t.endswith("_K1") or t.endswith("_K1_ref"))]
        t2 = [t for t in tips if t.startswith(f"{g}_") and (t.endswith("_K2") or t.endswith("_K2_ref"))]
        t_single = [t for t in tips if t.startswith(f"{g}_") and t.endswith("_single")]
        allowed = og_peers.get(g, {g})
        rr = ref_dom[g]
        arch_parts = [f"K1({rr['K1_status']})"]
        if rr["K2_range"] != "-":
            arch_parts.append(f"K2({rr['K2_status']})")
        if rr["accessory_domains"] != "-":
            arch_parts.append(rr["accessory_domains"])
        arch_str = "+".join(arch_parts)
        homs_g = hom_by_gene[g]
        n_tkp_g = sum(int(h["n_kinase_modules"]) >= 2 for h in homs_g)

        if t1 and t2:
            m1 = tree.common_ancestor(t1)
            m2 = tree.common_ancestor(t2)
            mrca = tree.common_ancestor(t1 + t2)
            m1_terms = [t.name for t in m1.get_terminals()]
            m2_terms = [t.name for t in m2.get_terminals()]
            mrca_terms = [t.name for t in mrca.get_terminals()]
            c1 = all(n.split("_")[0] in allowed and ("_K2" not in n) for n in m1_terms)
            c2 = all(n.split("_")[0] in allowed and ("_K1" not in n) for n in m2_terms)
            foreign = sorted({n.split("_")[0] for n in mrca_terms} - allowed)
            both = len(foreign) == 0
            shalrt, ufb = _support(mrca)
            sat = ref_pairs.get(g, {}).get("dS_saturated") == "1"
            origin = "intragenic duplication" if (not sat and not foreign) else ("intergenic fusion" if (sat and foreign) else "ambiguous")
            rows.append({"Gene": g, "n_K1_tips": len(t1), "n_K2_tips": len(t2),
                         "K1_monophyletic": int(bool(c1)), "K2_monophyletic": int(bool(c2)),
                         "K1_K2_sister_clades": int(bool(both) and bool(c1) and bool(c2)),
                         "MRCA_SH_aLRT": shalrt, "MRCA_UFBoot": ufb,
                         "other_genes_inside_K1K2_MRCA": ",".join(sorted({n.split("_")[0] for n in mrca_terms} - {g})) or "-",
                         "ref_K1_K2_patristic": f"{tree.distance(f'{g}_K1_ref', f'{g}_K2_ref'):.3f}",
                         "Set": ctx.cat[g]["Set"], "Domain_Architecture": arch_str,
                         "n_single_progenitor_tips": len(t_single),
                         "n_pangenome_homologues": len(homs_g), "n_two_kinase_homologues": n_tkp_g,
                         "Origin_Classification": origin})
        elif t1:
            t_k1_all = t1 + t_single
            if len(t_k1_all) > 1:
                m1 = tree.common_ancestor(t_k1_all)
                m1_terms = [t.name for t in m1.get_terminals()]
                c1 = int(all(n.split("_")[0] in allowed and ("_K2" not in n) for n in m1_terms))
                shalrt, ufb = _support(m1)
                other_in = ",".join(sorted({n.split("_")[0] for n in m1_terms} - {g})) or "-"
            else:
                c1, shalrt, ufb, other_in = 1, "NA", "NA", "-"
            rows.append({"Gene": g, "n_K1_tips": len(t1), "n_K2_tips": 0,
                         "K1_monophyletic": c1, "K2_monophyletic": "NA",
                         "K1_K2_sister_clades": "NA",
                         "MRCA_SH_aLRT": shalrt, "MRCA_UFBoot": ufb,
                         "other_genes_inside_K1K2_MRCA": other_in,
                         "ref_K1_K2_patristic": "NA",
                         "Set": ctx.cat[g]["Set"], "Domain_Architecture": arch_str,
                         "n_single_progenitor_tips": len(t_single),
                         "n_pangenome_homologues": len(homs_g), "n_two_kinase_homologues": n_tkp_g,
                         "Origin_Classification": f"single-kinase + accessory fusion ({rr['accessory_domains']})"})
    write_tsv(ctx.R / "S8b_kinase_phylogeny.tsv", rows, list(rows[0].keys()))
    log(f"phylogeny done: {len(seqs)} domains, {len(outgroups)} outgroups; trees in {ctx.T}")
    for r in rows:
        log(f"  {r['Gene']:5s} sister={r['K1_K2_sister_clades']} UFBoot={r['MRCA_UFBoot']} foreign={r['other_genes_inside_K1K2_MRCA']} origin={r['Origin_Classification']}")


# ----------------------------------------------------------------------------- 4 origins
def step_origins(ctx):
    log("=== rust33 step origins")
    ref = ctx.ref_arch()
    dom = [r for r in read_tsv(ctx.R / "S8a_classIII_domains.tsv") if r["Record"] == "homologue"]
    models = miniprot_models(ctx.W)
    by_asm_gene = collections.defaultdict(list)
    for (asm, _, q), ms in models.items():
        by_asm_gene[(asm, q)].extend(ms)
    hid, hcov = ctx.cfg["homolog"]["min_identity"], ctx.cfg["homolog"]["min_qcov"]
    rows = []
    for asm, p in ctx.panel.items():
        for g in CLASS3:
            a = ref[g]
            sub = [r for r in dom if r["asm_id"] == asm and r["Gene"] == g]
            ref_nk, ref_acc = len(a["kin"]), set(a["acc"])
            fused = [r for r in sub if float(r["identity"]) >= hid and float(r["qcov"]) >= hcov
                     and int(r["n_kinase_modules"]) >= ref_nk
                     and (not ref_acc or set(r["accessory_domains"].split(",")) & ref_acc)]
            single = [r for r in sub if int(r["n_kinase_modules"]) == 1 and r["accessory_domains"] == "-"]
            acc_only = "NA"
            if a["acc_span"]:
                k_end = max(e for _, e in a["kin"])
                s_acc, e_acc = a["acc_span"]
                acc_hits = [m for m in by_asm_gene[(asm, g)]
                            if m["qst"] + 1 > k_end - 10                                   # does not include the kinase
                            and min(m["qen"], e_acc) - max(m["qst"] + 1, s_acc) >= 0.6 * (e_acc - s_acc)]
                acc_only = int(bool(acc_hits))
            state = "fused" if fused else ("progenitor_only" if (single or acc_only == 1) else "absent")
            rows.append({"asm_id": asm, "species": p["species"], "ploidy": p["ploidy"], "genomes": p["genomes"],
                         "Gene": g, "fused_architecture_present": int(bool(fused)),
                         "n_fused_copies": len(fused), "best_fused_identity": f"{max([float(r['identity']) for r in fused], default=0):.4f}",
                         "single_kinase_progenitor_present": int(bool(single)), "n_single_kinase": len(single),
                         "standalone_accessory_progenitor_present": acc_only, "state": state})
    write_tsv(ctx.R / "S8d_architecture_distribution.tsv", rows, list(rows[0].keys()))
    log(f"origins done: {len(rows)} assembly x gene cells")


# ----------------------------------------------------------------------------- 5 exons
def phases_from_cds(lengths):
    """intron phase after each exon except the last = cumulative CDS length mod 3; returns [(aa_pos, phase)]"""
    out, cum = [], 0
    for L in lengths[:-1]:
        cum += L
        out.append((cum // 3, cum % 3))
    return out


def step_exons(ctx):
    log("=== rust33 step exons")
    ref = ctx.ref_arch()
    qcds = cds_with_headers(HERE / "queries_3.1.cds.fna")
    models = miniprot_models(ctx.W)
    te = {r["LocusGroup"]: r for r in read_tsv(ctx.R / "S7a_flanking_te_windows.tsv")} \
        if (ctx.R / "S7a_flanking_te_windows.tsv").exists() else {}
    rows = []
    for g in CLASS3:
        a = ref[g]
        bounds = [("K1_end", a["kin"][0][1])]
        if len(a["kin"]) > 1:
            bounds += [("K2_start", a["kin"][1][0]), ("K2_end", a["kin"][1][1])]
        if a["acc_span"]:
            bounds += [("accessory_start", a["acc_span"][0])]
        sources = []
        loc = qcds[g][1]
        ivs = [(int(x), int(y)) for x, y in re.findall(r"(\d+)\.\.(\d+)", loc)]
        if len(ivs) > 1:
            sources.append(("NCBI " + ctx.cat[g]["Protein_acc"], [e - s + 1 for s, e in ivs]))
        full = [m for ms in models.values() for m in ms
                if m["query"] == g and m["qst"] <= 5 and m["qen"] >= m["qlen"] - 5]
        if full:
            m = max(full, key=lambda m: m["ident"])
            cds = order_cds(m)
            sources.append((f"Miniprot {m['asm']} {m['seqid']}:{m['start']}-{m['end']} (id {m['ident']:.3f})",
                            [e - s + 1 for s, e in cds]))
        t = te.get(g, {})
        for src, lens in sources:
            ph = phases_from_cds(lens)
            for name, aa in bounds:
                near = min(ph, key=lambda x: abs(x[0] - aa)) if ph else None
                rows.append({"Gene": g, "source": src, "n_exons": len(lens),
                             "all_intron_phases": ",".join(str(p) for _, p in ph) or "-",
                             "domain_boundary": name, "boundary_aa": aa,
                             "nearest_intron_aa": near[0] if near else "-",
                             "distance_aa": abs(near[0] - aa) if near else "-",
                             "nearest_intron_phase": near[1] if near else "-",
                             "flank_total_TE_pct": t.get("total_te_pct", "NA"),
                             "flank_window": f"{t.get('ref_asm', 'NA')} {t.get('seqid', '')}:{t.get('window_start', '')}-{t.get('window_end', '')}"})
    write_tsv(ctx.R / "S8e_exon_intron_phases.tsv", rows, list(rows[0].keys()))
    log(f"exons done: {len(rows)} boundary records")


# ----------------------------------------------------------------------------- 6 stats
def classify(ctx):
    """Mechanism per reference TKP from data: intragenic duplication when K1/K2 resolve as sister/monophyletic
    clades within the gene family and reference dS is unsaturated; intergenic fusion when K1 and K2 resolve in
    phylogenetically distant kinase subfamilies and dS is saturated."""
    pairs = {r["Gene"]: r for r in read_tsv(ctx.R / "S8c_kinase_pair_divergence.tsv") if r["Record"] == "reference"}
    phy = {r["Gene"]: r for r in read_tsv(ctx.R / "S8b_kinase_phylogeny.tsv")} if (ctx.R / "S8b_kinase_phylogeny.tsv").exists() else {}
    out = {}
    for g, r in pairs.items():
        sat = r["dS_saturated"] == "1"
        raw_foreign = phy.get(g, {}).get("other_genes_inside_K1K2_MRCA", "-").split(",")
        foreign = set(raw_foreign) - {"-", "WTK4", "Rpg1"}
        if not sat and not foreign:
            out[g] = "intragenic duplication"
        elif sat and foreign:
            out[g] = "intergenic fusion"
        else:
            out[g] = "ambiguous"
    return out


def step_stats(ctx):
    log("=== rust33 step stats")
    dom = read_tsv(ctx.R / "S8a_classIII_domains.tsv")
    pairs = read_tsv(ctx.R / "S8c_kinase_pair_divergence.tsv")
    orig = read_tsv(ctx.R / "S8d_architecture_distribution.tsv")
    exo = read_tsv(ctx.R / "S8e_exon_intron_phases.tsv")
    phy = read_tsv(ctx.R / "S8b_kinase_phylogeny.tsv")
    mech = classify(ctx)
    dup = sorted(g for g, m in mech.items() if m == "intragenic duplication")
    fus = sorted(g for g, m in mech.items() if m == "intergenic fusion")
    amb = sorted(g for g, m in mech.items() if m == "ambiguous")
    f = lambda x: float(x)

    homs = [r for r in dom if r["Record"] == "homologue"]
    homs_main = [r for r in homs if r["Set"] == "Main"]
    ref_p = {r["Gene"]: r for r in pairs if r["Record"] == "reference"}
    hp = [r for r in pairs if r["Record"] == "homologue"]
    id_dup = [f(r["aa_identity_pct"]) for r in hp if r["Gene"] in dup]
    id_fus = [f(r["aa_identity_pct"]) for r in hp if r["Gene"] in fus]
    g_dup = [r for r in hp if r["Gene"] in dup and r["cds_source"] == "genomic CDS"]
    g_fus = [r for r in hp if r["Gene"] in fus and r["cds_source"] == "genomic CDS"]
    ds_dup = [f(r["dS_JC"]) for r in g_dup if r["dS_saturated"] == "0" and r["dS_JC"] != "nan"]
    ds_fus = [f(r["dS_JC"]) for r in g_fus if r["dS_saturated"] == "0" and r["dS_JC"] != "nan"]
    kv = []
    add = lambda k, v, n="": kv.append({"Metric": k, "Value": v, "Notes": n})

    add("N_classIII_queries", f"{len(CLASS3)} ({sum(ctx.cat[g]['Set'] == 'Main' for g in CLASS3)} main + "
        f"{sum(ctx.cat[g]['Set'] != 'Main' for g in CLASS3)} comparator)", ", ".join(CLASS3))
    add("N_classIII_homologues", f"{len(homs_main)} main-query homologues; {len(homs)} including comparators",
        "retained (paralogue-pruned) anchor-verified loci, S8a")
    add("N_tkp_homologues_two_kinase_modules", str(len(hp)), f"homologues with >=2 kinase modules >= {MIN_KIN_LEN} aa")
    for g in CLASS3:
        r = [x for x in dom if x["Record"] == "reference" and x["Gene"] == g][0]
        add(f"Ref_architecture_{g}", f"{r['n_kinase_modules']} kinase module(s); K1 {r['K1_status']}; K2 {r['K2_status']}; accessory {r['accessory_domains']}",
            f"K1 {r['K1_range']} ({r['K1_motifs']}), K2 {r['K2_range']} ({r['K2_motifs']})")
    for g, r in ref_p.items():
        add(f"Ref_K1K2_{g}", f"aa identity {f(r['aa_identity_pct']):.1f}%; pS {f(r['pS']):.3f}; dS {r['dS_JC'] if r['dS_saturated'] == '0' else 'saturated (pS >= 0.75)'}",
            f"mechanism call: {mech[g]}")
    for r in phy:
        sis = "NA (single-kinase fusion)" if r["K1_K2_sister_clades"] == "NA" else ("yes" if r["K1_K2_sister_clades"] == "1" else "no")
        add(f"Tree_{r['Gene']}", f"K1/K2 sister clades = {sis}; SH-aLRT/UFBoot {r['MRCA_SH_aLRT']}/{r['MRCA_UFBoot']}; patristic {r['ref_K1_K2_patristic']}",
            f"other genes inside K1+K2 MRCA: {r['other_genes_inside_K1K2_MRCA']}")
    add("Tree_outgroups_used", "cached Brachypodium/Hordeum kinase domains required for rooted tree", "rooted on basal outgroup clade (XP_044954759, XP_044957438)")
    add("Intragenic_duplication_TKPs", ", ".join(dup) or "none", "sister clades + unsaturated dS")
    add("Intergenic_fusion_TKPs", ", ".join(fus) or "none", "non-sister + saturated dS")
    add("Ambiguous_TKPs", ", ".join(amb) or "none", "tree and dS disagree")
    if id_dup and id_fus:
        # Collapse repeated assembly-level observations to one mean per reference locus.
        # Pair-level values remain descriptive because homologues are not independent.
        id_dup_gene = [float(np.mean([f(r["aa_identity_pct"]) for r in hp if r["Gene"] == g])) for g in dup]
        id_fus_gene = [float(np.mean([f(r["aa_identity_pct"]) for r in hp if r["Gene"] == g])) for g in fus]
        u, p = ss.mannwhitneyu(id_dup_gene, id_fus_gene, alternative="two-sided")
        t_stat, p_welch = ss.ttest_ind(id_dup_gene, id_fus_gene, equal_var=False)
        u_pair, p_pair = ss.mannwhitneyu(id_dup, id_fus, alternative="two-sided")
        add("K1K2_aa_identity_duplication", fmt_ms(id_dup_gene, 1) + " %", "mean of pangenome homologues within each reference locus")
        add("K1K2_aa_identity_fusion", fmt_ms(id_fus_gene, 1) + " %", "mean of pangenome homologues within each reference locus")
        add("MWU_aa_identity_duplication_vs_fusion", f"U = {u:.1f}, P = {p:.4f} (Welch t = {t_stat:.2f}, P = {p_welch:.4f})",
            "one value per reference locus; n = %d vs %d" % (len(id_dup_gene), len(id_fus_gene)))
        add("K1K2_aa_identity_pairs_summary",
            f"duplication {fmt_ms(id_dup, 2)}% vs fusion {fmt_ms(id_fus, 2)}% (MWU U = {u_pair:.1f}, P = {p_pair:.2e})",
            "all pangenome TKP homologue pairs")
    ref_ds_dup = [f(ref_p[g]["dS_JC"]) for g in dup if ref_p[g]["dS_saturated"] == "0"]
    ref_ps_dup = [f(ref_p[g]["pS"]) for g in dup]
    ref_ps_fus = [f(ref_p[g]["pS"]) for g in fus]
    u_ref, p_ref = ss.mannwhitneyu(ref_ps_dup, ref_ps_fus, alternative="two-sided")
    t_ref, pt_ref = ss.ttest_ind(ref_ps_dup, ref_ps_fus, equal_var=False)
    add("Ref_dS_duplication", fmt_ms(ref_ds_dup), "; ".join(f"{g}={ref_p[g]['dS_JC']}" for g in dup))
    add("Ref_pS_duplication", fmt_ms(ref_ps_dup), "; ".join(f"{g}={ref_p[g]['pS']}" for g in dup))
    add("Ref_pS_fusion", fmt_ms(ref_ps_fus), "; ".join(f"{g}={ref_p[g]['pS']}" for g in fus))
    add("MWU_Ref_pS_duplication_vs_fusion", f"U = {u_ref:.1f}, P = {p_ref:.4f} (Welch t = {t_ref:.2f}, P = {pt_ref:.4f})", "cloned reference queries (Nei-Gojobori pS)")
    pat_map = {r["Gene"]: f(r["ref_K1_K2_patristic"]) for r in phy if r["ref_K1_K2_patristic"] != "NA"}
    pat_dup = [pat_map[g] for g in dup if g in pat_map]
    pat_fus = [pat_map[g] for g in fus if g in pat_map]
    if pat_dup and pat_fus:
        u_pat, p_pat = ss.mannwhitneyu(pat_dup, pat_fus, alternative="two-sided")
        t_pat, pt_pat = ss.ttest_ind(pat_dup, pat_fus, equal_var=False)
        add("Ref_patristic_duplication_vs_fusion",
            f"duplication {fmt_ms(pat_dup, 3)} vs fusion {fmt_ms(pat_fus, 3)} (MWU U = {u_pat:.1f}, P = {p_pat:.4f}; Welch t = {t_pat:.2f}, P = {pt_pat:.4f})",
            "rooted ML patristic distance between reference K1 and K2")
    if g_dup and g_fus:
        ps_d, ps_f = [f(r["pS"]) for r in g_dup], [f(r["pS"]) for r in g_fus]
        u, p = ss.mannwhitneyu(ps_d, ps_f, alternative="two-sided")
        add("dS_duplication_genomic", fmt_ms(ds_dup), f"{sum(r['dS_saturated'] == '1' for r in g_dup)} of {len(g_dup)} pairs saturated")
        add("dS_fusion_genomic", fmt_ms(ds_fus), f"{sum(r['dS_saturated'] == '1' for r in g_fus)} of {len(g_fus)} pairs saturated (unsaturated subset shown)")
        add("pS_duplication_genomic", fmt_ms(ps_d), f"n = {len(g_dup)} pangenome homologue pairs with genomic CDS")
        add("pS_fusion_genomic", fmt_ms(ps_f), f"n = {len(g_fus)} pangenome homologue pairs with genomic CDS")
        add("MWU_pS_duplication_vs_fusion", f"U = {u:.1f}, P = {p:.2e}", f"n = {len(g_dup)} vs {len(g_fus)} pairs with genomic CDS")

    dip = [r for r in orig if r["ploidy"] == "2"]
    n_dip = len({r["asm_id"] for r in dip})
    for g in CLASS3:
        sub = [r for r in orig if r["Gene"] == g]
        dsub = [r for r in dip if r["Gene"] == g]
        prog = sum(r["single_kinase_progenitor_present"] == "1" for r in dsub)
        acc = [r for r in dsub if r["standalone_accessory_progenitor_present"] != "NA"]
        fused_in = sorted({r["asm_id"] for r in sub if r["fused_architecture_present"] == "1"})
        tiers = sorted({r["ploidy"] + "x" for r in sub if r["fused_architecture_present"] == "1"})
        add(f"Distribution_{g}", f"fused in {len(fused_in)}/{len({r['asm_id'] for r in sub})} assemblies (ploidy {','.join(tiers) or 'none'}); "
            f"single-kinase progenitor in {prog}/{n_dip} diploids"
            + (f"; standalone accessory in {sum(r['standalone_accessory_progenitor_present'] == '1' for r in acc)}/{n_dip} diploids" if acc else ""),
            "fused in: " + (", ".join(fused_in) if fused_in else "none"))
    for g in CLASS3:
        sub = [r for r in exo if r["Gene"] == g and r["domain_boundary"] in ("accessory_start", "K2_start")]
        for r in sub:
            if r["distance_aa"] != "-" and int(r["distance_aa"]) <= 30:
                add(f"Junction_{g}_{r['domain_boundary']}_{r['source'].split()[0]}",
                    f"intron phase {r['nearest_intron_phase']} at aa {r['nearest_intron_aa']} ({r['distance_aa']} aa from boundary aa {r['boundary_aa']})",
                    r["source"])
    te = {r["Gene"]: r["flank_total_TE_pct"] for r in exo}
    add("Flanking_TE_classIII", "; ".join(f"{g} {v}%" for g, v in te.items() if v != "NA"),
        "S7a total soft-masked TE over +/-200 kb (family-level LTR split not available)")

    write_tsv(ctx.R / "section3.3_values.tsv", kv, ["Metric", "Value", "Notes"])
    md = ["# Section 3.3 - computed values", "", "| Metric | Value | Notes |", "| :--- | :--- | :--- |"]
    md += [f"| `{r['Metric']}` | {r['Value']} | {r['Notes']} |" for r in kv]
    (HERE / "tables" / "Section3.3_Summary_Statistics.md").write_text("\n".join(md) + "\n")
    log("stats done -> results/section3.3_values.tsv")
    for r in kv:
        print(f"  {r['Metric']}: {r['Value']}   {r['Notes']}")


# ----------------------------------------------------------------------------- 8 figure
def step_figure(ctx):
    log("=== rust33 step figure")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from Bio import Phylo
    dom = read_tsv(ctx.R / "S8a_classIII_domains.tsv")
    pairs = read_tsv(ctx.R / "S8c_kinase_pair_divergence.tsv")
    orig = read_tsv(ctx.R / "S8d_architecture_distribution.tsv")
    exo = read_tsv(ctx.R / "S8e_exon_intron_phases.tsv")
    mech = classify(ctx)
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35)

    # a - domain diagrams
    ax = fig.add_subplot(gs[0, 0])
    refs = [r for r in dom if r["Record"] == "reference"]
    for i, r in enumerate(refs):
        y = len(refs) - i
        ax.plot([0, int(r["aa_len"])], [y, y], color="grey", lw=1)
        for key, st in (("K1_range", "K1_status"), ("K2_range", "K2_status")):
            if r[key] != "-":
                s, e = map(int, r[key].split("-"))
                ax.add_patch(Rectangle((s, y - 0.3), e - s, 0.6, color="#1f77b4" if r[st] == "active" else "#9ecae1"))
        if r["accessory_range"] != "-":
            s, e = map(int, r["accessory_range"].split("-"))
            ax.add_patch(Rectangle((s, y - 0.3), e - s, 0.6, color="#d95f02"))
            ax.text((s + e) / 2, y + 0.38, r["accessory_domains"], fontsize=6, ha="center")
        ax.text(-30, y, r["Gene"], ha="right", va="center", fontsize=8, style="italic")
    ax.set_xlim(-200, 1450); ax.set_ylim(0.3, len(refs) + 0.9); ax.set_yticks([])
    ax.set_xlabel("amino-acid position")
    ax.set_title("a  Domain organisation (dark = active, light = pseudokinase)", loc="left", fontsize=9, fontweight="bold")

    # b - tree
    ax = fig.add_subplot(gs[0, 1:])
    tf = ctx.T / "kinase_domains.treefile"
    if tf.exists():
        tree = Phylo.read(tf, "newick")
        outs = [t.name for t in tree.get_terminals() if t.name.startswith("OUT_")]
        basal_outs = [t for t in outs if any(x in t for x in ("XP_044954759", "XP_044957438"))] or outs
        try:
            tree.root_with_outgroup(*basal_outs) if basal_outs else tree.root_at_midpoint()
        except Exception:
            tree.root_at_midpoint()
        tree.ladderize()
        col = lambda n: "#1f77b4" if "_K1" in n else ("#d62728" if "_K2" in n else ("grey" if n.startswith("OUT_") else "#2ca02c"))
        Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: c.name if c.is_terminal() else None,
                   label_colors=col, branch_labels=None)
        for t in ax.texts:
            t.set_fontsize(5)
        ax.set_ylabel(""); ax.set_yticks([])
    ax.set_title("b  ML tree of dissected kinase domains (blue K1, red K2, green single kinase, grey outgroup)",
                 loc="left", fontsize=9, fontweight="bold")

    # c - identity / dS
    ax = fig.add_subplot(gs[1, 0])
    hp = [r for r in pairs if r["Record"] == "homologue"]
    groups = [("intragenic\nduplication", [g for g, m in mech.items() if m == "intragenic duplication"]),
              ("intergenic\nfusion", [g for g, m in mech.items() if m == "intergenic fusion"])]
    data = [[float(r["aa_identity_pct"]) for r in hp if r["Gene"] in gl] for _, gl in groups]
    if all(data):
        ax.boxplot(data, widths=0.5, showfliers=False)
        for i, d in enumerate(data, start=1):
            ax.scatter(np.random.normal(i, 0.06, len(d)), d, s=6, alpha=0.4)
        ax.set_xticks([1, 2]); ax.set_xticklabels([f"{n}\n({', '.join(gl)})" for n, gl in groups], fontsize=7)
        u, p = ss.mannwhitneyu(data[0], data[1], alternative="two-sided")
        ax.set_title(f"c  K1 vs K2 aa identity (MWU P = {p:.1e})", loc="left", fontsize=9, fontweight="bold")
    ax.set_ylabel("Kinase 1 vs Kinase 2 identity (%)")

    # d - presence heatmap
    ax = fig.add_subplot(gs[1, 1])
    asms = list(dict.fromkeys(r["asm_id"] for r in orig))
    code = {"absent": 0, "progenitor_only": 1, "fused": 2}
    M = np.array([[code[[r for r in orig if r["asm_id"] == a and r["Gene"] == g][0]["state"]] for g in CLASS3] for a in asms])
    ax.imshow(M, aspect="auto", cmap=matplotlib.colors.ListedColormap(["#f0f0f0", "#fdae6b", "#08519c"]), vmin=0, vmax=2)
    ax.set_xticks(range(len(CLASS3))); ax.set_xticklabels(CLASS3, rotation=45, fontsize=7, style="italic")
    ax.set_yticks(range(len(asms))); ax.set_yticklabels(asms, fontsize=6)
    ax.set_title("d  fused (blue) / progenitor only (orange) / absent", loc="left", fontsize=9, fontweight="bold")

    # e - exon structures with phases
    ax = fig.add_subplot(gs[1, 2])
    y = 0
    for g in CLASS3:
        sub = [r for r in exo if r["Gene"] == g]
        if not sub:
            continue
        r = sub[0]
        ph = r["all_intron_phases"].split(",") if r["all_intron_phases"] != "-" else []
        n = int(r["n_exons"])
        for i in range(n):
            ax.add_patch(Rectangle((i, y - 0.25), 0.8, 0.5, color="#636363"))
            if i < len(ph):
                ax.text(i + 0.9, y + 0.3, ph[i], fontsize=5, ha="center")
        ax.text(-0.3, y, g, ha="right", va="center", fontsize=7, style="italic")
        y -= 1
    ax.set_xlim(-2.5, 18); ax.set_ylim(y + 0.3, 0.8); ax.axis("off")
    ax.set_title("e  exons (boxes) and intron phases", loc="left", fontsize=9, fontweight="bold")

    out = HERE / "figures" / "Fig3_tandem_kinase_origin"
    fig.savefig(f"{out}.pdf", bbox_inches="tight")
    fig.savefig(f"{out}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log(f"figure done -> {out}.pdf/.png")


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["domains", "pairs", "phylogeny", "origins", "exons", "stats", "figure", "all"])
    ap.add_argument("--config", default=str(HERE / "config_3.1.json"))
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    os.chdir(HERE)
    ctx = Ctx(load_cfg(args.config))
    order = ["domains", "pairs", "phylogeny", "origins", "exons", "stats", "figure"]
    for s in (order if args.step == "all" else [args.step]):
        {"domains": lambda: step_domains(ctx), "pairs": lambda: step_pairs(ctx),
         "phylogeny": lambda: step_phylogeny(ctx, args.threads), "origins": lambda: step_origins(ctx),
         "exons": lambda: step_exons(ctx), "stats": lambda: step_stats(ctx),
         "figure": lambda: step_figure(ctx)}[s]()


if __name__ == "__main__":
    main()
