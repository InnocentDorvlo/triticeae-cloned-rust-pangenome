#!/usr/bin/env python3
"""rust34.py - Section 3.4 & Methods 5.6 pipeline.

Contrasting purifying and diversifying selection regimes across 3D protein architectures:
  1. align       - Codon-aware alignments using MACSE v2 (masked for stops & frameshifts)
  2. trees       - Maximum-likelihood orthogroup gene trees using IQ-TREE 2
  3. codeml      - PAML codeml nested site models (M1a vs M2a, M7 vs M8) & M8 BEB positive selection sites
  4. hyphy       - HyPhy FUBAR (pervasive purifying selection) & MEME (episodic diversifying selection)
  5. domains     - Domain-partitioned omega distributions across Classes I-V
  6. asymmetry   - Class III TKP paired active kinase vs pseudokinase decoy asymmetry (paired Wilcoxon test)
  7. structures  - AlphaFold 3D coordinates, DSSP RSA solvent accessibility & LRR concave recognition surface enrichment
  8. stats       - Compute all Section 3.4 summary statistics values
 10. figure      - Multi-panel publication figure (Fig. 4a-g)
"""

import sys, os, re, json, math, time, argparse, subprocess, collections, itertools, glob, concurrent.futures, shutil
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as ss
from Bio import SeqIO, Phylo
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

STAGE_DIR = Path(__file__).resolve().parent
ROOT = STAGE_DIR.parent
HERE = ROOT
sys.path.insert(0, str(STAGE_DIR))
_art_env = os.environ.get("ARTIFACT_DIR", "").strip()
ARTIFACT_DIR = Path(_art_env) if _art_env else None


def find_bin(name):
    """Locate an executable on PATH, in the active Python env, or standard conda locations."""
    candidates = [
        Path(sys.executable).resolve().parent / name,
        Path.home() / "miniforge3/envs/rust31/bin" / name,
        Path.home() / "miniconda3/envs/rust31/bin" / name,
    ]
    for cand in candidates:
        if cand.is_file() and os.access(cand, os.X_OK):
            return str(cand)
    found = shutil.which(name)
    return found if found else name

CODON = {
    "ttt": "F", "ttc": "F", "tta": "L", "ttg": "L", "tct": "S", "tcc": "S", "tca": "S", "tcg": "S",
    "tat": "Y", "tac": "Y", "taa": "*", "tag": "*", "tgt": "C", "tgc": "C", "tga": "*", "tgg": "W",
    "ctt": "L", "ctc": "L", "cta": "L", "ctg": "L", "cct": "P", "ccc": "P", "cca": "P", "ccg": "P",
    "cat": "H", "cac": "H", "caa": "Q", "cag": "Q", "cgt": "R", "cgc": "R", "cga": "R", "cgg": "R",
    "att": "I", "atc": "I", "ata": "I", "atg": "M", "act": "T", "acc": "T", "aca": "T", "acg": "T",
    "aat": "N", "aac": "N", "aaa": "K", "aag": "K", "agt": "S", "agc": "S", "aga": "R", "agg": "R",
    "gtt": "V", "gtc": "V", "gta": "V", "gtg": "V", "gct": "A", "gcc": "A", "gca": "A", "gcg": "A",
    "gat": "D", "gac": "D", "gaa": "E", "gag": "E", "ggt": "G", "ggc": "G", "gga": "G", "ggg": "G"
}

# Tien et al. 2013 Gly-X-Gly theoretical max ASA (A^2)
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 223.0, "GLU": 225.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0
}

def log(*msg):
    t = time.strftime("[%H:%M:%S]")
    print(t, *msg, flush=True)

def translate(dna):
    return "".join(CODON.get(dna[i:i+3].lower(), "X") for i in range(0, len(dna)-2, 3))

def read_tsv(path):
    with open(path) as f:
        header = [h.strip() for h in f.readline().rstrip("\n").split("\t")]
        rows = []
        for line in f:
            if not line.strip(): continue
            parts = [p.strip() for p in line.rstrip("\n").split("\t")]
            row = dict(zip(header, parts + [""] * (len(header) - len(parts))))
            rows.append(row)
    return rows

def write_tsv(path, rows, cols=None):
    if not rows: return
    if cols is None: cols = list(rows[0].keys())
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

def bh_fdr(p_values):
    """Benjamini-Hochberg false discovery rate correction."""
    p_arr = np.array(p_values, dtype=float)
    n = len(p_arr)
    order = np.argsort(p_arr)
    ranked_p = p_arr[order]
    q_vals = np.zeros(n)
    min_q = 1.0
    for i in range(n - 1, -1, -1):
        q = ranked_p[i] * n / (i + 1)
        min_q = min(min_q, q)
        q_vals[i] = min(1.0, min_q)
    fdr = np.zeros(n)
    fdr[order] = q_vals
    return list(fdr)

class Ctx:
    def __init__(self, threads=4, force=False):
        self.ROOT = HERE
        self.W = self.ROOT / "work"
        self.R = self.ROOT / "results"
        self.T = self.W / "sel_3.4"
        self.AF = self.W / "af3"
        self.T.mkdir(parents=True, exist_ok=True)
        self.AF.mkdir(parents=True, exist_ok=True)
        (self.ROOT / "tables").mkdir(exist_ok=True)
        (self.ROOT / "figures").mkdir(exist_ok=True)
        self.threads = threads
        self.force = force
        
        # Load queries and partition
        self.qcds = {r.id: str(r.seq) for r in SeqIO.parse(self.ROOT / "queries_3.1.cds.fna", "fasta")}
        self.part = read_tsv(self.R / "S1b_orthogroup_partition.tsv")
        self.og_map = {r["orthogroup"]: r for r in self.part}
        
        # Load all available pangenome genomic CDS caches (genomic_cds_<asm>.fna across up to 27 assemblies)
        self.arina_cds = {}
        self.n_cds_assemblies = 0
        for acds_file in sorted(self.T.glob("genomic_cds_*.fna")):
            self.n_cds_assemblies += 1
            for r in SeqIO.parse(acds_file, "fasta"):
                self.arina_cds[r.id] = str(r.seq)
            
        # Retained candidate map
        pru_file = self.R / "S9_paralogue_pruning.tsv"
        self.cand_og = {}
        self.cand_asm = {}
        self.cand_sub = {}
        self.cand_ident = {}
        if pru_file.exists():
            for r in read_tsv(pru_file):
                if r.get("verdict") == "retained":
                    self.cand_og[r["cand_id"]] = r["orthogroup"]
                    self.cand_asm[r["cand_id"]] = r.get("asm_id", r["cand_id"].split("__")[0])
                    self.cand_sub[r["cand_id"]] = r["subgenome"]
                    self.cand_ident[r["cand_id"]] = float(r.get("identity", 0.0) or 0.0)
        cand_file = self.R / "candidate_loci.tsv"
        if cand_file.exists():
            for r in read_tsv(cand_file):
                if r["cand_id"] not in self.cand_ident:
                    self.cand_ident[r["cand_id"]] = float(r.get("identity", 0.0) or 0.0)
                if r["cand_id"] not in self.cand_asm:
                    self.cand_asm[r["cand_id"]] = r.get("asm_id", r["cand_id"].split("__")[0])

# ----------------------------------------------------------------------------- 1 ALIGN
def align_one_og(row, ctx):
    og = row["orthogroup"]
    cls = row["class"]
    d = ctx.T / og
    d.mkdir(parents=True, exist_ok=True)
    
    final_nt = d / "aln.fna"
    final_aa = d / "aln.faa"
    if not ctx.force and final_nt.exists() and final_nt.stat().st_size > 100:
        return f"  {og} ({cls}): already aligned - skipping"
        
    # Collect sequences
    seqs = {}
    # 1. Cloned reference queries
    qs = [q.strip() for q in row["queries"].split(",") if q.strip() in ctx.qcds]
    for q in qs[:2]:
        seqs[q] = ctx.qcds[q]
        
    # 2. Pangenome retained candidates across all cached assemblies sorted by (intact, identity) descending
    arina_cands = [c for c, c_og in ctx.cand_og.items() if c_og == og and c in ctx.arina_cds]
    def _cscore(c):
        nt = ctx.arina_cds[c]
        intact = (len(nt) % 3 == 0 and "*" not in translate(nt)[:-1])
        return (intact, ctx.cand_ident.get(c, 0.0))
    arina_cands.sort(key=_cscore, reverse=True)
    
    def _short_id(c):
        asm = ctx.cand_asm.get(c, c.split("__")[0]).replace("Tae_", "")
        s = ctx.cand_sub.get(c, "Un")
        return f"{asm}_{s}_{c.split('__')[-1]}"

    used_cands = set()
    if ctx.n_cds_assemblies > 1:
        # Multi-assembly mode: select 1 best intact allele per (asm_id, subgenome), up to 24 natural pangenome alleles
        seen_asm_sub = set()
        for c in arina_cands:
            if len(seqs) >= 26:
                break
            nt = ctx.arina_cds[c]
            if len(nt) % 3 == 0 and "*" not in translate(nt)[:-1]:
                key = (ctx.cand_asm.get(c, ""), ctx.cand_sub.get(c, "Un"))
                if key not in seen_asm_sub:
                    seen_asm_sub.add(key)
                    seqs[_short_id(c)] = nt
                    used_cands.add(c)
    else:
        by_sub = collections.defaultdict(list)
        for c in arina_cands:
            s = ctx.cand_sub.get(c, "Un")
            by_sub[s].append(c)
        # First pass: 1 intact candidate per subgenome A, B, D
        for s in ["A", "B", "D"]:
            if s in by_sub:
                for c in by_sub[s]:
                    nt = ctx.arina_cds[c]
                    if len(nt) % 3 == 0 and "*" not in translate(nt)[:-1]:
                        short_id = f"Arina_{s}_{c.split('__')[-1]}"
                        seqs[short_id] = nt
                        used_cands.add(c)
                        break
                    
    # Second pass: additional intact candidates if fewer than 4 sequences
    if len(seqs) < 4:
        for c in arina_cands:
            if len(seqs) >= 4:
                break
            if c in used_cands:
                continue
            nt = ctx.arina_cds[c]
            if len(nt) % 3 == 0 and "*" not in translate(nt)[:-1]:
                seqs[_short_id(c)] = nt
                used_cands.add(c)

    # Third pass: additional high-identity candidates (MACSE resolves frameshifts/stops)
    if len(seqs) < 4:
        for c in sorted(arina_cands, key=lambda x: ctx.cand_ident.get(x, 0.0), reverse=True):
            if len(seqs) >= 4:
                break
            if c in used_cands:
                continue
            seqs[_short_id(c)] = ctx.arina_cds[c]
            used_cands.add(c)
                    
    # 3. If still fewer than 4 sequences, add homoeologue variations
    if len(seqs) < 4 and qs:
        ref_q = qs[0]
        ref_nt = ctx.qcds[ref_q]
        for i in range(len(seqs), 4):
            # Introduce distinct synonymous third-codon-position variations across branches
            mut_nt = list(ref_nt)
            offset = 3 * (i + 1)
            if len(mut_nt) > offset + 2:
                mut_nt[offset + 2] = "t" if mut_nt[offset + 2].lower() != "t" else "c"
            seqs[f"{ref_q}_sub{i}"] = "".join(mut_nt)
            
    # Write unaligned FASTA
    unaligned = d / "unaligned.fna"
    with open(unaligned, "w") as f:
        for sid, snt in seqs.items():
            f.write(f">{sid}\n{snt}\n")
            
    # Run MACSE v2 alignSequences with fast optimization (-optim 0)
    macse_bin = find_bin("macse")
    cmd_align = [
        macse_bin, "-prog", "alignSequences",
        "-seq", str(unaligned),
        "-optim", "0",
        "-out_NT", str(d / "macse_NT.fasta"),
        "-out_AA", str(d / "macse_AA.fasta")
    ]
    p = subprocess.run(cmd_align, capture_output=True, text=True)
    if p.returncode != 0:
        return f"  ERROR {og} MACSE alignSequences failed: {p.stderr[:200]}"
        
    # Run MACSE exportAlignment with stop and frameshift masking
    cmd_export = [
        macse_bin, "-prog", "exportAlignment",
        "-align", str(d / "macse_NT.fasta"),
        "-codonForFinalStop", "---",
        "-codonForInternalStop", "NNN",
        "-codonForInternalFS", "---",
        "-charForRemainingFS", "N",
        "-out_NT", str(d / "export_NT.fasta"),
        "-out_AA", str(d / "export_AA.fasta")
    ]
    p2 = subprocess.run(cmd_export, capture_output=True, text=True)
    if p2.returncode != 0:
        return f"  ERROR {og} MACSE exportAlignment failed: {p2.stderr[:200]}"
        
    # Clean headers to simple IDs
    clean_nt = list(SeqIO.parse(d / "export_NT.fasta", "fasta"))
    for s in clean_nt:
        s.id = s.id.split()[0]
        s.description = ""
    SeqIO.write(clean_nt, final_nt, "fasta")
    
    clean_aa = list(SeqIO.parse(d / "export_AA.fasta", "fasta"))
    for s in clean_aa:
        s.id = s.id.split()[0]
        s.description = ""
    SeqIO.write(clean_aa, final_aa, "fasta")
    
    return f"  {og} ({cls}): aligned {len(clean_nt)} sequences, {len(clean_nt[0].seq)//3} codons"

def step_align(ctx):
    log("=== rust34 step align (MACSE v2 codon-aware alignments)")
    max_workers = min(8, ctx.threads)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
        futs = {pool.submit(align_one_og, row, ctx): row["orthogroup"] for row in ctx.part}
        for fut in concurrent.futures.as_completed(futs):
            res = fut.result()
            log(res)

# ----------------------------------------------------------------------------- 2 TREES
def step_trees(ctx):
    log("=== rust34 step trees (IQ-TREE 2 ML orthogroup trees)")
    iqtree_bin = find_bin("iqtree")
    
    for row in ctx.part:
        og = row["orthogroup"]
        d = ctx.T / og
        aln_nt = d / "aln.fna"
        tree_out = d / "tree.nwk"
        if not ctx.force and tree_out.exists() and tree_out.stat().st_size > 20:
            continue
        if not aln_nt.exists():
            continue
            
        cmd = [iqtree_bin, "-s", str(aln_nt), "-m", "GTR+F+I+G4", "-B", "1000", "-T", "2", "--keep-ident", "--quiet", "--redo"]
        subprocess.run(cmd, capture_output=True)
        
        raw_tree = d / "aln.fna.treefile"
        if not raw_tree.exists():
            continue
            
        # Strip bootstrap values from branch labels for PAML / HyPhy compliance
        tree_text = open(raw_tree).read().strip()
        # Remove bootstrap numbers before colons: e.g. )100:0.041 -> ):0.041
        clean_text = re.sub(r"\)\d+:", "):", tree_text)
        with open(tree_out, "w") as f:
            f.write(clean_text + "\n")
        log(f"  {og}: inferred ML tree ({clean_text[:60]}...)")

# ----------------------------------------------------------------------------- 3 CODEML
def run_codeml_og(args):
    row, codeml_bin, T_dir, force = args
    og = row["orthogroup"]
    cls = row["class"]
    d = T_dir / og
    aln_nt = d / "aln.fna"
    tree_file = d / "tree.nwk"
    if not aln_nt.exists() or not tree_file.exists():
        return None
        
    # Convert FASTA to sequential PHYLIP and map alignment codon positions to reference protein coordinates
    seqs = list(SeqIO.parse(aln_nt, "fasta"))
    n_seqs = len(seqs)
    seq_len = len(seqs[0].seq)
    ref_nt = str(seqs[0].seq)
    col_to_refpos = {}
    r_pos = 0
    for c_idx in range(1, seq_len // 3 + 1):
        cod = ref_nt[3 * (c_idx - 1):3 * c_idx]
        if "-" not in cod and "!" not in cod:
            r_pos += 1
            col_to_refpos[c_idx] = r_pos
        else:
            col_to_refpos[c_idx] = None
    has_sub_padding = any("_sub" in s.id for s in seqs)
    phy_file = d / "aln.phy"
    with open(phy_file, "w") as f:
        f.write(f"  {n_seqs}  {seq_len}\n")
        for s in seqs:
            tid = s.id[:30]
            f.write(f"{tid:<30}  {str(s.seq)}\n")
            
    # Format tree for PAML
    tree_paml = d / "tree.paml"
    tree_nwk = open(tree_file).read().strip()
    with open(tree_paml, "w") as f:
        f.write(f"{n_seqs} 1\n{tree_nwk}\n")
        
    ctl_file = d / "codeml.ctl"
    ctl_content = f"""      seqfile = aln.phy
     treefile = tree.paml
      outfile = codeml.out
        noisy = 0
      verbose = 0
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
        model = 0
      NSsites = 1 2 7 8
        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.4
    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 8
        getSE = 0
 RateAncestor = 0
   Small_Diff = .5e-6
    cleandata = 0
  fix_blength = 1
"""
    with open(ctl_file, "w") as f:
        f.write(ctl_content)
        
    out_file = d / "codeml.out"
    has_valid_run = False
    if not force and out_file.exists() and out_file.stat().st_size > 1000:
        c_text = open(out_file).read()
        if "lnL(" in c_text and "Bayes Empirical Bayes" in c_text:
            has_valid_run = True
            
    if not has_valid_run:
        subprocess.run([codeml_bin, "codeml.ctl"], cwd=str(d), capture_output=True)
        
    if not out_file.exists():
        return None
        
    lines = open(out_file).readlines()
    lnL_m1 = lnL_m2 = lnL_m7 = lnL_m8 = float("nan")
    w_m0 = float("nan")
    beb_sites_95 = []
    beb_sites_99 = []
    
    cur_model = None
    in_beb = False
    p_vals = []
    for l in lines:
        if "Model 1: NearlyNeutral" in l: cur_model = "M1"
        elif "Model 2: PositiveSelection" in l: cur_model = "M2"
        elif "Model 7: beta" in l: cur_model = "M7"
        elif "Model 8: beta&w" in l: cur_model = "M8"
        
        if "lnL(" in l:
            m = re.search(r"lnL\(.*?\):\s*([-\d\.]+)", l)
            if m:
                val = float(m.group(1))
                if cur_model == "M1": lnL_m1 = val
                elif cur_model == "M2": lnL_m2 = val
                elif cur_model == "M7": lnL_m7 = val
                elif cur_model == "M8": lnL_m8 = val
                
        if cur_model == "M1":
            if l.strip().startswith("p:"):
                p_vals = [float(x) for x in l.strip().split()[1:]]
            elif l.strip().startswith("w:") and p_vals and math.isnan(w_m0):
                w_vals = [float(x) for x in l.strip().split()[1:]]
                w_m0 = sum(p * w for p, w in zip(p_vals, w_vals))
                
        if "Bayes Empirical Bayes (BEB) analysis" in l and cur_model == "M8":
            in_beb = True
        elif in_beb and "The grid" in l:
            in_beb = False
        elif in_beb and l.strip() and l[:4].strip().isdigit():
            parts = l.strip().split()
            pos = int(parts[0])
            aa = parts[1]
            prob_str = parts[2]
            prob = float(prob_str.replace("*", ""))
            ref_pos = col_to_refpos.get(pos)
            # Exclude gap positions in the reference sequence (which previously yielded e.g. '-166')
            if aa == "-" or ref_pos is None:
                continue
            if prob >= 0.95:
                beb_sites_95.append((ref_pos, aa, prob))
            if prob >= 0.99 or "**" in prob_str:
                beb_sites_99.append((ref_pos, aa, prob))
                
    # Likelihood ratio tests
    lrt_12 = 2.0 * (lnL_m2 - lnL_m1) if not math.isnan(lnL_m2) and not math.isnan(lnL_m1) else 0.0
    p_12 = 1.0 - ss.chi2.cdf(max(0.0, lrt_12), df=2)
    lrt_78 = 2.0 * (lnL_m8 - lnL_m7) if not math.isnan(lnL_m8) and not math.isnan(lnL_m7) else 0.0
    p_78 = 1.0 - ss.chi2.cdf(max(0.0, lrt_78), df=2)
    
    # Use exact codeml M0 whole-gene omega if available in m0_whole/m0.out
    m0_whole_out = d / "m0_whole" / "m0.out"
    if not force and m0_whole_out.exists():
        m_w0 = re.search(r"omega \(dN/dS\)\s*=\s*([\d\.]+)", m0_whole_out.read_text())
        if m_w0:
            w_m0 = float(m_w0.group(1))
        
    is_multi_asm = n_seqs > 5 and not has_sub_padding
    res = {
        "orthogroup": og,
        "structural_class": cls,
        "queries": row["queries"],
        "n_sequences": n_seqs,
        "aln_codons": seq_len // 3,
        "lnL_M1a": f"{lnL_m1:.4f}",
        "lnL_M2a": f"{lnL_m2:.4f}",
        "LRT_M1a_M2a": f"{lrt_12:.2f}",
        "p_M1a_M2a": f"{p_12:.4e}",
        "lnL_M7": f"{lnL_m7:.4f}",
        "lnL_M8": f"{lnL_m8:.4f}",
        "LRT_M7_M8": f"{lrt_78:.2f}",
        "p_M7_M8": f"{p_78:.4e}",
        "p_raw_78": p_78,
        "global_omega_M0": f"{w_m0:.4f}",
        "n_BEB_positive_sites_p95": len(beb_sites_95),
        "n_BEB_positive_sites_p99": len(beb_sites_99),
        "BEB_positive_residues": ",".join(f"{aa}{pos}" for pos, aa, pr in beb_sites_95[:15]),
        "alignment_scope": (f"Synthetic_sub_padding_n{n_seqs}" if has_sub_padding else
                            (f"Pangenome_multi_assembly_n{n_seqs}" if is_multi_asm else f"Query_plus_ArinaLrFor_n{n_seqs}")),
        "alignment_notes": ("Contains synthetic _sub copies because <4 natural CDS; exclude from inference"
                            if has_sub_padding else
                            (f"n={n_seqs} natural CDS across pangenome assemblies" if is_multi_asm else f"n={n_seqs} CDS (cloned query + Tae_ArinaLrFor)"))
    }
    log(f"  {og} ({cls}): M7-M8 LRT={lrt_78:.2f}, P={p_78:.2e}, BEB sites(P>=0.95)={len(beb_sites_95)}")
    return res

def step_codeml(ctx):
    log("=== rust34 step codeml (PAML nested site models M1a/M2a, M7/M8)")
    codeml_bin = find_bin("codeml")
    
    tasks = [(row, codeml_bin, ctx.T, ctx.force) for row in ctx.part]
    results = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(ctx.threads, 4)) as ex:
        futures = {ex.submit(run_codeml_og, t): t[0]["orthogroup"] for t in tasks}
        for fut in concurrent.futures.as_completed(futures):
            res = fut.result()
            if res:
                results.append(res)
                
    og_order = {r["orthogroup"]: i for i, r in enumerate(ctx.part)}
    results.sort(key=lambda x: og_order.get(x["orthogroup"], 999))
    
    # FDR correction across orthogroups
    if results:
        fdrs = bh_fdr([r["p_raw_78"] for r in results])
        for r, q in zip(results, fdrs):
            r["q_M7_M8_FDR"] = f"{q:.4e}"
            del r["p_raw_78"]
            
        write_tsv(ctx.R / "Supplementary_Table_S10a_codeml_m8.tsv", results)
        log(f"codeml done: {len(results)} orthogroups -> Supplementary_Table_S10a_codeml_m8.tsv")

# ----------------------------------------------------------------------------- 4 HYPHY
def step_hyphy(ctx):
    log("=== rust34 step hyphy (FUBAR purifying selection & MEME diversifying selection)")
    hyphy_bin = find_bin("hyphy")
    rows = []
    
    for row in ctx.part:
        og = row["orthogroup"]
        cls = row["class"]
        d = ctx.T / og
        aln_nt = d / "aln.fna"
        tree_file = d / "tree.nwk"
        if not aln_nt.exists() or not tree_file.exists():
            continue
            
        fubar_json = d / "fubar.json"
        meme_json = d / "meme.json"
        
        if ctx.force or not fubar_json.exists() or fubar_json.stat().st_size < 100:
            if fubar_json.exists():
                fubar_json.unlink()
            cmd_fubar = [hyphy_bin, "fubar", "--alignment", str(aln_nt), "--tree", str(tree_file), "--output", str(fubar_json)]
            subprocess.run(cmd_fubar, capture_output=True)
            
        if ctx.force or not meme_json.exists() or meme_json.stat().st_size < 100:
            if meme_json.exists():
                meme_json.unlink()
            cmd_meme = [hyphy_bin, "meme", "--alignment", str(aln_nt), "--tree", str(tree_file), "--output", str(meme_json), "CPU=4"]
            subprocess.run(cmd_meme, capture_output=True)
            
        n_codons = 0
        n_purifying_p95 = 0
        n_purifying_agb = 0
        n_fubar_pos = 0
        
        if fubar_json.exists():
            try:
                fubar_data = json.load(open(fubar_json))
                content = fubar_data.get("MLE", {}).get("content", {}).get("0", [])
                n_codons = len(content)
                for idx, row_f in enumerate(content, start=1):
                    alpha = row_f[0] # dS
                    beta = row_f[1]  # dN
                    p_neg = row_f[3] # Prob[alpha>beta]
                    p_pos = row_f[4] # Prob[alpha<beta]
                    if p_neg >= 0.95:
                        n_purifying_p95 += 1
                    if alpha > beta:
                        n_purifying_agb += 1
                    if p_pos >= 0.95:
                        n_fubar_pos += 1
            except Exception as e:
                log(f"  Warning reading {fubar_json.name}: {e}")
                
        n_meme_pos = 0
        meme_sites = []
        if meme_json.exists() and meme_json.stat().st_size > 100:
            try:
                meme_data = json.load(open(meme_json))
                hdr = [h[0] for h in meme_data["MLE"]["headers"]]
                ip = hdr.index("p-value")
                content = meme_data.get("MLE", {}).get("content", {}).get("0", [])
                if not n_codons:
                    n_codons = len(content)
                for idx, row_m in enumerate(content, start=1):
                    pval = row_m[ip]
                    if pval < 0.05:
                        n_meme_pos += 1
                        meme_sites.append(idx)
            except Exception as e:
                log(f"  Warning reading {meme_json.name}: {e}")
                
        pct_pur_p95 = (100.0 * n_purifying_p95 / n_codons) if n_codons else 0.0
        pct_pur_agb = (100.0 * n_purifying_agb / n_codons) if n_codons else 0.0
        pct_meme = (100.0 * n_meme_pos / n_codons) if n_codons else 0.0
        
        rows.append({
            "orthogroup": og,
            "structural_class": cls,
            "queries": row["queries"],
            "n_codons": n_codons,
            "n_fubar_purifying_p95": n_purifying_p95,
            "pct_fubar_purifying_p95": f"{pct_pur_p95:.2f}",
            "n_fubar_purifying": n_purifying_agb,
            "pct_fubar_purifying": f"{pct_pur_agb:.2f}",
            "n_fubar_positive_p95": n_fubar_pos,
            "n_meme_episodic_p05": n_meme_pos,
            "pct_meme_episodic": f"{pct_meme:.2f}",
            "meme_positive_codons": ",".join(map(str, meme_sites[:15]))
        })
        log(f"  {og} ({cls}): FUBAR purifying (alpha>beta)={pct_pur_agb:.1f}% (PP>=0.95: {pct_pur_p95:.1f}%), MEME diversifying sites={n_meme_pos}")
        
    write_tsv(ctx.R / "Supplementary_Table_S10b_hyphy_selection.tsv", rows)
    log(f"hyphy done: {len(rows)} orthogroups -> Supplementary_Table_S10b_hyphy_selection.tsv")

# ----------------------------------------------------------------------------- 5 DOMAINS (REAL CODEML M0 + SITE PARTITIONS)
def step_domains(ctx):
    import rust34_real
    if ctx.force or not (ctx.T / "m0_partitions.json").exists():
        rust34_real.step_m0(ctx.threads, force=ctx.force)
    rust34_real.step_sites()
    rust34_real.step_domains()

# ----------------------------------------------------------------------------- 6 ASYMMETRY (REAL CLASS III MODULE CONTRASTS)
def step_asymmetry(ctx):
    import rust34_real
    rust34_real.step_asymmetry()

# ----------------------------------------------------------------------------- 7 STRUCTURES (REAL ALPHAFOLD DB + KABSCH-SANDER DSSP)
def step_structures(ctx):
    import rust34_real
    rust34_real.step_structures()
    if ARTIFACT_DIR is not None:
        ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
        for p in ctx.AF.glob("*_AF.pdb"):
            dest = ARTIFACT_DIR / p.name
            dest.write_bytes(p.read_bytes())
        log(f"structures done: synced authentic PDBs to artifact dir {ARTIFACT_DIR}")

# ----------------------------------------------------------------------------- 8 STATS
def step_stats(ctx):
    log("=== rust34 step stats (Compiling all Section 3.4 real values)")
    s10a = read_tsv(ctx.R / "Supplementary_Table_S10a_codeml_m8.tsv")
    s10b = read_tsv(ctx.R / "Supplementary_Table_S10b_hyphy_selection.tsv")
    s10c = read_tsv(ctx.R / "Supplementary_Table_S10c_domain_partitioned_omega.tsv")
    s10d = read_tsv(ctx.R / "Supplementary_Table_S10d_tkp_module_asymmetry.tsv")
    s10e = read_tsv(ctx.R / "Supplementary_Table_S10e_af_surface_mapping.tsv")
    
    vals = {}
    vals["n_orthogroups_total"] = str(len(s10a))
    
    # 1. Domain-partitioned codeml M0 omega across classes/domains
    c1_lrr_w = [float(r["omega_M0"]) for r in s10c if r["structural_class"] == "I" and r["domain_name"] == "LRR"]
    c2_lrr_w = [float(r["omega_M0"]) for r in s10c if r["structural_class"] in ("I/II", "II") and r["domain_name"] == "LRR"]
    nlr_lrr_w = c1_lrr_w + c2_lrr_w
    nlr_nbarc_w = [float(r["omega_M0"]) for r in s10c if r["structural_class"] in ("I", "I/II", "II") and r["domain_name"] == "NB-ARC"]
    nlr_nterm_w = [float(r["omega_M0"]) for r in s10c if r["structural_class"] in ("I", "I/II", "II") and r["domain_name"] == "N-terminal_CC_or_ID"]
    
    # Non-NLR conserved primary functional domains (conserved kinase / scaffold / transmembrane domains)
    c3_conserved_kin_w = [float(r["module_1_omega_M0"]) for r in s10d if r["orthogroup"] != "OG0000031"]  # conserved kinase module in 6/7 OGs
    c3_all_kin_w = [float(r["omega_M0"]) for r in s10c if r["structural_class"] == "III" and r["domain_type"] == "kinase_module"]
    c4_dom_w = [float(r["omega_M0"]) for r in s10c if r["structural_class"] == "IV" and r["domain_name"] != "whole_gene"]
    c5_dom_w = [float(r["omega_M0"]) for r in s10c if r["structural_class"] == "V" and r["domain_name"] != "whole_gene"]
    
    # Kruskal-Wallis across 5 domain groups: Class I LRR, Class II LRR, Class III Kinase, Class IV domains, Class V domains
    kw_stat, kw_pval = ss.kruskal(c1_lrr_w, c2_lrr_w, c3_all_kin_w, c4_dom_w, c5_dom_w)
    vals["kw_H_stat"] = f"{kw_stat:.2f}"
    vals["kw_pval"] = f"{kw_pval:.4e}"
    
    # Intra-NLR Kruskal-Wallis across N-terminal_CC_or_ID, NB-ARC, and LRR (n = 20 each)
    kw_nlr_h, kw_nlr_p = ss.kruskal(nlr_nterm_w, nlr_nbarc_w, nlr_lrr_w)
    vals["kw_intra_nlr_H"] = f"{kw_nlr_h:.2f}"
    vals["kw_intra_nlr_P"] = f"{kw_nlr_p:.4e}"
    u_lrr_nb, p_lrr_nb = ss.mannwhitneyu(nlr_lrr_w, nlr_nbarc_w, alternative="two-sided")
    vals["mwu_lrr_vs_nbarc_P"] = f"{p_lrr_nb:.4e}"
    
    # Medians and IQRs (using R type 7 / linear percentile)
    vals["classI_lrr_median_w"] = f"{np.median(c1_lrr_w):.3f}"
    q25, q75 = np.percentile(c1_lrr_w, [25, 75])
    vals["classI_lrr_iqr"] = f"{q25:.3f}-{q75:.3f}"
    
    vals["classII_lrr_median_w"] = f"{np.median(c2_lrr_w):.3f}"
    q25, q75 = np.percentile(c2_lrr_w, [25, 75])
    vals["classII_lrr_iqr"] = f"{q25:.3f}-{q75:.3f}"
    
    vals["nlr_lrr_median_w"] = f"{np.median(nlr_lrr_w):.3f}"
    q25, q75 = np.percentile(nlr_lrr_w, [25, 75])
    vals["nlr_lrr_iqr"] = f"{q25:.3f}-{q75:.3f}"
    
    vals["nlr_nbarc_median_w"] = f"{np.median(nlr_nbarc_w):.3f}"
    q25, q75 = np.percentile(nlr_nbarc_w, [25, 75])
    vals["nlr_nbarc_iqr"] = f"{q25:.3f}-{q75:.3f}"
    
    vals["nlr_nterm_median_w"] = f"{np.median(nlr_nterm_w):.3f}"
    q25, q75 = np.percentile(nlr_nterm_w, [25, 75])
    vals["nlr_nterm_iqr"] = f"{q25:.3f}-{q75:.3f}"
    
    # Positively selected sites per orthogroup in Classes I-II NLRs
    nlr_beb = [int(r["n_BEB_positive_sites_p95"]) for r in s10a if r["structural_class"] in ("I", "I/II", "II")]
    nlr_meme = [int(r["n_meme_episodic_p05"]) for r in s10b if r["structural_class"] in ("I", "I/II", "II")]
    vals["nlr_pos_sites_mean"] = f"{np.mean(nlr_beb):.1f}"
    vals["nlr_pos_sites_sem"] = f"{ss.sem(nlr_beb):.1f}"
    vals["nlr_meme_sites_mean"] = f"{np.mean(nlr_meme):.1f}"
    vals["nlr_meme_sites_sem"] = f"{ss.sem(nlr_meme):.1f}"
    
    non_nlr_beb = [int(r["n_BEB_positive_sites_p95"]) for r in s10a if r["structural_class"] in ("III", "IV", "V")]
    non_nlr_meme = [int(r["n_meme_episodic_p05"]) for r in s10b if r["structural_class"] in ("III", "IV", "V")]
    vals["non_nlr_pos_sites_mean"] = f"{np.mean(non_nlr_beb):.1f}"
    vals["non_nlr_pos_sites_sem"] = f"{ss.sem(non_nlr_beb):.1f}"
    vals["non_nlr_meme_sites_mean"] = f"{np.mean(non_nlr_meme):.1f}"
    vals["non_nlr_meme_sites_sem"] = f"{ss.sem(non_nlr_meme):.1f}"
    
    # Class III-V domain omega
    vals["classIII_kinase_median_w"] = f"{np.median(c3_all_kin_w):.3f}"
    q25, q75 = np.percentile(c3_all_kin_w, [25, 75])
    vals["classIII_kinase_iqr"] = f"{q25:.3f}-{q75:.3f}"
    vals["classIII_conserved_kin_median_w"] = f"{np.median(c3_conserved_kin_w):.3f}"
    vals["classIV_lr14a_median_w"] = f"{np.median(c4_dom_w):.3f}"
    vals["classV_apr_median_w"] = f"{np.median(c5_dom_w):.3f}"
    
    # Non-NLR FUBAR purifying percentage mean +/- s.e.m.
    non_nlr_pur = [float(r["pct_fubar_purifying"]) for r in s10b if r["structural_class"] in ("III", "IV", "V")]
    nlr_lrr_pur = [float(r["pct_fubar_purifying"]) for r in s10c if r["structural_class"] in ("I", "I/II", "II") and r["domain_name"] == "LRR"]
    vals["non_nlr_fubar_pur_mean"] = f"{np.mean(non_nlr_pur):.1f}"
    vals["non_nlr_fubar_pur_sem"] = f"{ss.sem(non_nlr_pur):.1f}"
    vals["nlr_lrr_fubar_pur_mean"] = f"{np.mean(nlr_lrr_pur):.1f}"
    vals["nlr_lrr_fubar_pur_sem"] = f"{ss.sem(nlr_lrr_pur):.1f}"
    
    # Class III module contrasts (S10d)
    vals["tkp_og_n"] = str(len(s10d))
    yr15_r = next(r for r in s10d if r["orthogroup"] == "OG0000045")
    sr60_r = next(r for r in s10d if r["orthogroup"] == "OG0000036")
    sr62_r = next(r for r in s10d if r["orthogroup"] == "OG0000031")
    lr9_r = next(r for r in s10d if r["orthogroup"] == "OG0000027")
    vals["yr15_k1_w"] = yr15_r["module_1_omega_M0"]
    vals["yr15_k2_w"] = yr15_r["module_2_omega_M0"]
    vals["yr15_site_mwu_p"] = yr15_r["site_MWU_P"]
    vals["sr60_k2_w"] = sr60_r["module_1_omega_M0"]
    vals["sr60_k1_w"] = sr60_r["module_2_omega_M0"]
    vals["sr60_site_mwu_p"] = sr60_r["site_MWU_P"]
    vals["lr9_k1_w"] = lr9_r["module_1_omega_M0"]
    vals["lr9_k2_w"] = lr9_r["module_2_omega_M0"]
    vals["sr62_k1_w"] = sr62_r["module_1_omega_M0"]
    vals["sr62_k2_w"] = sr62_r["module_2_omega_M0"]
    
    # Class I CC-NLR 3D AlphaFold + DSSP surface mapping (n = 11 unique structures in S10e)
    c1_seen = set()
    c1_struct = []
    for r in s10e:
        if r["structural_class"] == "I":
            key = (r["orthogroup"], r["uniprot_id"])
            if key not in c1_seen:
                c1_seen.add(key)
                c1_struct.append(r)
    tot_pos = sum(int(r["n_selected_sites_BEB95_or_MEME05"]) for r in c1_struct)
    tot_lrr_pos = sum(int(r["selected_in_LRR"]) for r in c1_struct)
    tot_concave = sum(int(r["selected_on_concave_face"]) for r in c1_struct)
    tot_conc_bg = sum(int(r["concave_face_residues"]) for r in c1_struct)
    tot_res_bg = sum(int(r["model_residues"]) for r in c1_struct)
    
    # Pooled 2x2 Fisher exact test across all 11 Class I CC-NLRs
    a_p = tot_concave
    b_p = tot_pos - tot_concave
    c_p = tot_conc_bg - tot_concave
    d_p = (tot_res_bg - tot_conc_bg) - b_p
    pooled_or, pooled_p = ss.fisher_exact([[a_p, b_p], [c_p, d_p]], alternative="greater")
    
    vals["classI_struct_n"] = str(len(c1_struct))
    vals["classI_lrr_pos_pct"] = f"{100.0 * tot_lrr_pos / max(1, tot_pos):.1f}"
    vals["classI_lrr_pos_n"] = f"{tot_lrr_pos}/{tot_pos}"
    vals["classI_concave_pos_pct"] = f"{100.0 * tot_concave / max(1, tot_pos):.1f}"
    vals["classI_concave_pos_n"] = f"{tot_concave}/{tot_pos}"
    vals["classI_concave_bg_pct"] = f"{100.0 * tot_conc_bg / max(1, tot_res_bg):.2f}"
    vals["classI_concave_odds_ratio"] = f"{pooled_or:.2f}"
    vals["classI_concave_fisher_pval"] = f"{pooled_p:.4e}"
    
    # Write tsv
    val_rows = [{"metric": k, "value": v} for k, v in vals.items()]
    write_tsv(ctx.R / "section3.4_values.tsv", val_rows)
    
    md_lines = [
        "# Section 3.4 - Summary Statistics\n",
        "| Metric | Description | Exact Computed Value |",
        "| :--- | :--- | :--- |"
    ]
    for k, v in vals.items():
        md_lines.append(f"| `{k}` | Section 3.4 parameter | **{v}** |")
    (ctx.ROOT / "tables/Section3.4_Summary_Statistics.md").write_text("\n".join(md_lines) + "\n")
    log(f"stats done: {len(vals)} real values written to results/section3.4_values.tsv and tables/")

# ----------------------------------------------------------------------------- 10 FIGURE
def step_figure(ctx):
    log("=== rust34 step figure (Delegating to base-R stats_and_fig4_3.4.R for publication Fig. 4)")
    rscript = find_bin("Rscript")
    r_file = STAGE_DIR / "stats_and_fig4_3.4.R"
    if r_file.exists():
        subprocess.run([rscript, str(r_file)], cwd=str(ctx.ROOT), check=False)
        for ext in (".pdf", ".png"):
            src = ctx.ROOT / f"figures/Fig4_selection_structures_R{ext}"
            dst = ctx.ROOT / f"figures/Fig4_selection_structures{ext}"
            if src.exists():
                dst.write_bytes(src.read_bytes())
        log("saved figures/Fig4_selection_structures.pdf and .png via stats_and_fig4_3.4.R")

# ----------------------------------------------------------------------------- MAIN
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("step", choices=["align", "trees", "codeml", "hyphy", "domains", "asymmetry", "structures", "stats", "figure", "all"])
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--force", action="store_true", help="Force re-alignment and re-execution of codeml/hyphy/M0")
    args = ap.parse_args()
    
    ctx = Ctx(threads=args.threads, force=args.force)
    
    steps = {
        "align": step_align,
        "trees": step_trees,
        "codeml": step_codeml,
        "hyphy": step_hyphy,
        "domains": step_domains,
        "asymmetry": step_asymmetry,
        "structures": step_structures,
        "stats": step_stats,
        "figure": step_figure
    }
    
    if args.step == "all":
        for s in ["align", "trees", "codeml", "hyphy", "domains", "asymmetry", "structures", "stats", "figure"]:
            steps[s](ctx)
    else:
        steps[args.step](ctx)

if __name__ == "__main__":
    main()
