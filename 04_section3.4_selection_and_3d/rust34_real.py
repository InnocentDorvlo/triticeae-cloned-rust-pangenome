#!/usr/bin/env python3
"""
rust34_real.py - real (non-simulated) Section 3.4 domain, module and structure analyses.

Computes empirical step_domains / step_asymmetry / step_structures analyses for rust34.py.
Every number written here is computed from: MACSE codon alignments, IQ-TREE trees,
PAML codeml (M0 per partition; M8 BEB per site), HyPhy FUBAR + MEME site tables,
InterProScan domain coordinates (S1c), Class III kinase coordinates (S8a),
and AlphaFold DB models + DSSP (mkdssp).

Steps:
  m0          codeml M0 on whole gene and on each domain partition      -> work/sel_3.4/OG*/part_*/
  sites       per-codon site table (BEB, FUBAR, MEME, query aa, domain)  -> results/Supplementary_Table_S10g_site_table.tsv
  domains     domain-partitioned omega + site counts                     -> results/Supplementary_Table_S10c_domain_partitioned_omega.tsv
  asymmetry   Class III active vs pseudokinase module contrast           -> results/Supplementary_Table_S10d_tkp_module_asymmetry.tsv
  structures  AlphaFold/DSSP mapping of selected sites                   -> results/Supplementary_Table_S10e_af_surface_mapping.tsv
"""
import sys, re, json, math, argparse, subprocess, concurrent.futures, shutil
from pathlib import Path
import numpy as np
import scipy.stats as ss
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

STAGE_DIR = Path(__file__).resolve().parent
HERE = STAGE_DIR.parent
sys.path.insert(0, str(STAGE_DIR))
from rust34 import find_bin, log, read_tsv, write_tsv  # noqa: E402

R = HERE / "results"
T = HERE / "work" / "sel_3.4"
AF = HERE / "work" / "af3"
CODEML = find_bin("codeml")
MKDSSP = find_bin("mkdssp")
NLR = ("I", "I/II", "II")

# ----------------------------------------------------------------------------- helpers
def load_part():
    return read_tsv(R / "S1b_orthogroup_partition.tsv")

def ipr_hits():
    hits = {}
    for r in read_tsv(R / "S1c_interproscan_queries.tsv"):
        hits.setdefault(r["ID"], []).append(r)
    return hits

def s8a_refs():
    return {r["ID"]: r for r in read_tsv(R / "S8a_classIII_domains.tsv") if r["Record"] == "reference"}

def rng(s):
    a, b = s.split("-")
    return int(a), int(b)

def read_aln(og):
    return [(r.id, str(r.seq).upper()) for r in SeqIO.parse(T / og / "aln.fna", "fasta")]

def col_map(seq_nt):
    """alignment codon column (1-based) -> residue position (1-based) in that sequence; None if gap."""
    m, pos = {}, 0
    for c in range(len(seq_nt) // 3):
        cod = seq_nt[3 * c:3 * c + 3]
        if cod.replace("-", "").replace("!", "") == "":
            m[c + 1] = None
        else:
            pos += 1
            m[c + 1] = pos
    return m

def domain_partitions(qid, cls, plen, hits, refs):
    """Return list of (domain_name, domain_type, [(start,end),...]) in query protein coordinates."""
    H = hits.get(qid, [])
    pf = lambda acc: [(int(h["start"]), int(h["end"])) for h in H if h["signature"] == acc]
    sig = lambda accs: [(int(h["start"]), int(h["end"])) for h in H if h["signature"] in accs]
    parts = []
    if cls in NLR:
        nb = sorted(pf("PF00931"))
        if not nb:
            return []
        nb_s = nb[0][0]
        whd = [x for x in pf("PF23559") if x[0] > nb_s]
        nb_e = max([nb[0][1]] + [x[1] for x in whd[:1]])
        lrr = [x for x in sig({"G3DSA:3.80.10.10", "SSF52058", "PF23598", "PF25019", "PF00560"}) if x[1] > nb_e]
        parts.append(("N-terminal_CC_or_ID", "signaling_integrated", [(1, nb_s - 1)]))
        parts.append(("NB-ARC", "nucleotide_binding", [(nb_s, nb_e)]))
        if lrr:
            l_s = max(nb_e + 1, min(x[0] for x in lrr))
            l_e = max(x[1] for x in lrr)
            parts.append(("LRR", "pathogen_recognition", [(l_s, l_e)]))
            if plen - l_e >= 30:
                parts.append(("C-terminal_post-LRR", "c_terminal", [(l_e + 1, plen)]))
    elif cls == "III":
        r = refs.get(qid)
        if r is None:
            return []
        if r["K1_range"] != "-":
            parts.append((f"K1_{r['K1_status']}", "kinase_module", [rng(r["K1_range"])]))
        if r["K2_range"] != "-":
            parts.append((f"K2_{r['K2_status']}", "kinase_module", [rng(r["K2_range"])]))
        if r["accessory_range"] not in ("-", ""):
            parts.append((f"Accessory_{r['accessory_domains'].split(',')[0]}", "accessory", [rng(r["accessory_range"])]))
    elif qid == "Lr14a":
        ank = sig({"PF12796", "PF13857", "SM00248"})
        a_s, a_e = min(x[0] for x in ank), max(x[1] for x in ank)
        parts.append(("Ankyrin_repeats", "protein_scaffold", [(a_s, a_e)]))
        parts.append(("PGG_transmembrane", "membrane", [(a_e + 1, plen)]))
    elif qid.startswith("YrNAM"):
        parts.append(("NAM_domain", "dna_binding", pf("PF02365")))
        parts.append(("BED_zinc_finger", "dna_binding", sig({"PF02892", "SM00614"})[:1]))
    elif qid == "Lr34":
        parts.append(("ABC_NBD", "atp_hydrolysis", pf("PF00005")))
        parts.append(("ABC2_TMD", "membrane", pf("PF01061")))
    elif qid == "Lr67":
        parts.append(("Sugar_tr_TM_core", "membrane", pf("PF00083")))
    return [(n, t, [x for x in rs if x[1] >= x[0]]) for n, t, rs in parts if rs]

def query_for_og(og, row, hits, refs):
    """first sequence in the alignment that is a catalogue query with domain annotation."""
    ids = [i for i, _ in read_aln(og)]
    qs = [q.strip() for q in row["queries"].split(",")]
    for i in ids:
        if i in qs and (i in hits or i in refs):
            return i
    return ids[0]

# ----------------------------------------------------------------------------- codeml M0
M0_CTL = """      seqfile = aln.phy
     treefile = tree.paml
      outfile = m0.out
        noisy = 0
      verbose = 0
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        clock = 0
        model = 0
      NSsites = 0
        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.4
    fix_alpha = 1
        alpha = 0
    cleandata = 0
  fix_blength = 1
"""

def write_phy(path, seqs):
    L = len(seqs[0][1])
    with open(path, "w") as f:
        f.write(f"  {len(seqs)}  {L}\n")
        for i, s in seqs:
            f.write(f"{i[:30]:<30}  {s.replace('!', '-')}\n")

def run_m0(d, seqs, tree, force=False):
    d.mkdir(parents=True, exist_ok=True)
    out = d / "m0.out"
    if force or not (out.exists() and "omega (dN/dS)" in out.read_text()):
        write_phy(d / "aln.phy", seqs)
        (d / "tree.paml").write_text(f"{len(seqs)} 1\n{tree}\n")
        (d / "codeml.ctl").write_text(M0_CTL)
        subprocess.run([CODEML, "codeml.ctl"], cwd=str(d), capture_output=True)
    txt = out.read_text() if out.exists() else ""
    w = re.search(r"omega \(dN/dS\)\s*=\s*([\d\.]+)", txt)
    dn = re.search(r"tree length for dN:\s*([\d\.]+)", txt)
    ds = re.search(r"tree length for dS:\s*([\d\.]+)", txt)
    lnl = re.search(r"lnL\(.*?\):\s*([-\d\.]+)", txt)
    f = lambda m: float(m.group(1)) if m else float("nan")
    return {"omega_M0": f(w), "tree_dN": f(dn), "tree_dS": f(ds), "lnL_M0": f(lnl)}

def m0_jobs(force=False):
    part, hits, refs = load_part(), ipr_hits(), s8a_refs()
    qprot = {r.id: str(r.seq) for r in SeqIO.parse(HERE / "queries_3.1.faa", "fasta")}
    jobs = []
    for row in part:
        og, cls = row["orthogroup"], row["class"]
        if not (T / og / "aln.fna").exists():
            continue
        seqs = read_aln(og)
        tree = (T / og / "tree.nwk").read_text().strip()
        q = query_for_og(og, row, hits, refs)
        qseq = dict(seqs)[q]
        cmap = col_map(qseq)
        plen = len(qprot.get(q, "")) or max(v for v in cmap.values() if v)
        parts = domain_partitions(q, cls, plen, hits, refs)
        jobs.append((og, cls, row["queries"], q, "whole_gene", "whole_gene", [(1, plen)], seqs, tree, cmap, force))
        for name, dtype, rs in parts:
            jobs.append((og, cls, row["queries"], q, name, dtype, rs, seqs, tree, cmap, force))
    return jobs

def cols_for(rs, cmap):
    return [c for c, p in cmap.items() if p is not None and any(a <= p <= b for a, b in rs)]

def do_m0_job(job):
    og, cls, queries, q, name, dtype, rs, seqs, tree, cmap, force = job
    cols = sorted(cmap) if name == "whole_gene" else cols_for(rs, cmap)
    sub = [(i, "".join(s[3 * (c - 1):3 * c] for c in cols)) for i, s in seqs]
    res = run_m0(T / og / ("m0_whole" if name == "whole_gene" else f"part_{name}"), sub, tree, force=force)
    return {"orthogroup": og, "structural_class": cls, "queries": queries, "reference_query": q,
            "domain_name": name, "domain_type": dtype,
            "query_ranges": ";".join(f"{a}-{b}" for a, b in rs), "n_codons": len(cols), "cols": cols, **res}

def step_m0(threads=4, force=False):
    log("=== rust34_real m0 (codeml M0 on whole genes and InterProScan/S8a domain partitions)")
    jobs = m0_jobs(force=force)
    with concurrent.futures.ThreadPoolExecutor(threads) as ex:
        res = list(ex.map(do_m0_job, jobs))
    json.dump(res, open(T / "m0_partitions.json", "w"))
    for r in res:
        log(f"  {r['orthogroup']} {r['domain_name']:<24} codons={r['n_codons']:>5} omega={r['omega_M0']:.4f} dS={r['tree_dS']:.4f}")
    log(f"m0 done: {len(res)} partitions")

# ----------------------------------------------------------------------------- site table
def parse_beb(og):
    """M8 BEB posterior P(w>1) for every printed site (codeml prints sites with P>0.5)."""
    txt = (T / og / "codeml.out").read_text().splitlines()
    cur, inb, out = None, False, {}
    for l in txt:
        mm = re.search(r"NSsites Model (\d+):", l)
        if mm:
            cur = int(mm.group(1)); inb = False
        if cur == 8 and "Bayes Empirical Bayes (BEB) analysis" in l:
            inb = True; continue
        if inb and "The grid" in l:
            inb = False
        if inb:
            m = re.match(r"\s+(\d+)\s+(\S)\s+([\d\.]+)\*{0,2}", l)
            if m:
                out[int(m.group(1))] = float(m.group(3))
    return out

def parse_fubar(og):
    j = json.load(open(T / og / "fubar.json"))
    rows = j["MLE"]["content"]["0"]
    return {i + 1: {"alpha": r[0], "beta": r[1], "p_neg": r[3], "p_pos": r[4]} for i, r in enumerate(rows)}

def parse_meme(og):
    p = T / og / "meme.json"
    if not p.exists() or p.stat().st_size < 100:
        return None
    j = json.load(open(p))
    hdr = [h[0] for h in j["MLE"]["headers"]]
    ip = hdr.index("p-value")
    return {i + 1: r[ip] for i, r in enumerate(j["MLE"]["content"]["0"])}

def step_sites():
    log("=== rust34_real sites (per-codon BEB / FUBAR / MEME table)")
    parts = json.load(open(T / "m0_partitions.json"))
    rows = []
    for og in sorted({p["orthogroup"] for p in parts}):
        P = [p for p in parts if p["orthogroup"] == og]
        whole = [p for p in P if p["domain_name"] == "whole_gene"][0]
        seqs = dict(read_aln(og))
        q = whole["reference_query"]
        cmap = col_map(seqs[q])
        beb, fub, meme = parse_beb(og), parse_fubar(og), parse_meme(og)
        col2dom = {}
        for p in P:
            if p["domain_name"] != "whole_gene":
                for c in p["cols"]:
                    col2dom[c] = p["domain_name"]
        for c in sorted(cmap):
            cod = seqs[q][3 * (c - 1):3 * c]
            aa = str(Seq(cod).translate()) if "-" not in cod and "!" not in cod else "-"
            f = fub.get(c, {})
            rows.append({"orthogroup": og, "structural_class": whole["structural_class"], "reference_query": q,
                         "aln_codon": c, "query_pos": cmap[c] if cmap[c] else "", "query_aa": aa,
                         "domain": col2dom.get(c, "inter-domain" if cmap[c] else "gap_in_query"),
                         "BEB_P_w_gt_1": f"{beb.get(c, 0.0):.3f}",
                         "FUBAR_alpha": f"{f.get('alpha', float('nan')):.4f}", "FUBAR_beta": f"{f.get('beta', float('nan')):.4f}",
                         "FUBAR_P_neg": f"{f.get('p_neg', float('nan')):.4f}", "FUBAR_P_pos": f"{f.get('p_pos', float('nan')):.4f}",
                         "MEME_p": f"{meme[c]:.4g}" if meme else "NA"})
    write_tsv(R / "Supplementary_Table_S10g_site_table.tsv", rows)
    log(f"sites done: {len(rows)} codons; MEME available for {len({r['orthogroup'] for r in rows if r['MEME_p'] != 'NA'})} OGs")

# ----------------------------------------------------------------------------- domains
def step_domains():
    log("=== rust34_real domains (real domain-partitioned omega)")
    parts = json.load(open(T / "m0_partitions.json"))
    sites = read_tsv(R / "Supplementary_Table_S10g_site_table.tsv")
    by = {}
    for s in sites:
        by.setdefault((s["orthogroup"], int(s["aln_codon"])), s)
    out = []
    for p in parts:
        S = [by[(p["orthogroup"], c)] for c in (p["cols"] if p["domain_name"] != "whole_gene"
             else [k[1] for k in by if k[0] == p["orthogroup"]])]
        n = len(S)
        fa = [s for s in S if s["FUBAR_alpha"] != "nan"]
        pur = sum(float(s["FUBAR_alpha"]) > float(s["FUBAR_beta"]) for s in fa)
        pur95 = sum(float(s["FUBAR_P_neg"]) >= 0.95 for s in fa)
        pos95 = sum(float(s["FUBAR_P_pos"]) >= 0.95 for s in fa)
        beb95 = sum(float(s["BEB_P_w_gt_1"]) >= 0.95 for s in S)
        meme = [s for s in S if s["MEME_p"] != "NA"]
        meme05 = sum(float(s["MEME_p"]) < 0.05 for s in meme) if meme else "NA"
        out.append({"orthogroup": p["orthogroup"], "structural_class": p["structural_class"], "queries": p["queries"],
                    "reference_query": p["reference_query"], "domain_name": p["domain_name"], "domain_type": p["domain_type"],
                    "query_ranges": p["query_ranges"], "n_codons": n,
                    "omega_M0": f"{p['omega_M0']:.4f}", "tree_dN": f"{p['tree_dN']:.4f}", "tree_dS": f"{p['tree_dS']:.4f}",
                    "lnL_M0": f"{p['lnL_M0']:.3f}",
                    "pct_fubar_purifying": f"{100 * pur / len(fa):.2f}" if fa else "NA",
                    "n_fubar_purifying_p95": pur95, "n_fubar_positive_p95": pos95,
                    "n_BEB_p95": beb95, "n_MEME_p05": meme05})
    write_tsv(R / "Supplementary_Table_S10c_domain_partitioned_omega.tsv", out)
    log(f"domains done: {len(out)} rows (incl. whole_gene) -> S10c")

# ----------------------------------------------------------------------------- asymmetry
def step_asymmetry():
    log("=== rust34_real asymmetry (Class III module contrasts across all 7 Class III OGs + S8a homologues)")
    s10c = read_tsv(R / "Supplementary_Table_S10c_domain_partitioned_omega.tsv")
    sites = read_tsv(R / "Supplementary_Table_S10g_site_table.tsv")
    s8a = read_tsv(R / "S8a_classIII_domains.tsv")
    rows = []
    # Explicit primary (conserved kinase) vs secondary (pseudokinase or fused accessory) module pairs for all 7 Class III OGs
    module_pairs = [
        ("OG0000045", "K1_active", "K2_pseudokinase", "active_vs_pseudokinase"),
        ("OG0000036", "K2_active", "K1_pseudokinase", "active_vs_pseudokinase"),
        ("OG0000031", "K1_active", "K2_pseudokinase", "active_vs_pseudokinase"),
        ("OG0000007", "K1_pseudokinase", "K2_pseudokinase", "pseudokinase_vs_pseudokinase"),
        ("OG0000027", "K1_pseudokinase", "K2_pseudokinase", "pseudokinase_vs_pseudokinase"),
        ("OG0000029", "K1_active", "Accessory_START", "kinase_vs_accessory"),
        ("OG0000057", "K1_pseudokinase", "Accessory_DUF3475", "kinase_vs_accessory"),
    ]
    for og, mod1_name, mod2_name, pair_cat in module_pairs:
        D = {r["domain_name"]: r for r in s10c if r["orthogroup"] == og}
        if mod1_name not in D or mod2_name not in D:
            continue
        a, p = D[mod1_name], D[mod2_name]
        sa = [float(s["FUBAR_beta"]) - float(s["FUBAR_alpha"]) for s in sites
              if s["orthogroup"] == og and s["domain"] == mod1_name and s["FUBAR_alpha"] != "nan"]
        sp = [float(s["FUBAR_beta"]) - float(s["FUBAR_alpha"]) for s in sites
              if s["orthogroup"] == og and s["domain"] == mod2_name and s["FUBAR_alpha"] != "nan"]
        u, pu = ss.mannwhitneyu(sp, sa, alternative="two-sided")
        # Compute pangenome homologue Kimura protein divergence (d_AA) from codon alignment for each module
        seqs = dict(read_aln(og))
        ref_s = seqs[a["reference_query"]]
        def kimura_d(cols):
            dists = []
            for sid, s in seqs.items():
                if sid == a["reference_query"]:
                    continue
                match, diff = 0, 0
                for c in cols:
                    c1 = ref_s[3*(c-1):3*c]
                    c2 = s[3*(c-1):3*c]
                    if "-" in c1 or "!" in c1 or "-" in c2 or "!" in c2:
                        continue
                    aa1 = str(Seq(c1).translate())
                    aa2 = str(Seq(c2).translate())
                    if aa1 == aa2:
                        match += 1
                    else:
                        diff += 1
                tot = match + diff
                if tot >= 30:
                    p = diff / tot
                    if p < 0.85:
                        d = -math.log(max(1e-6, 1.0 - p - 0.2 * p * p))
                        dists.append(d)
            return dists

        m1_cols = [int(s["aln_codon"]) for s in sites if s["orthogroup"] == og and s["domain"] == mod1_name]
        m2_cols = [int(s["aln_codon"]) for s in sites if s["orthogroup"] == og and s["domain"] == mod2_name]
        d1 = kimura_d(m1_cols)
        d2 = kimura_d(m2_cols)
        hom_n = len(d1)
        hom_d1_med = f"{np.median(d1):.4f}" if d1 else "NA"
        hom_d2_med = f"{np.median(d2):.4f}" if d2 else "NA"
        rows.append({"orthogroup": og, "reference_query": a["reference_query"], "queries": a["queries"],
                     "pair_category": pair_cat,
                     "module_1": a["domain_name"], "module_1_range": a["query_ranges"], "module_1_codons": a["n_codons"],
                     "module_1_omega_M0": a["omega_M0"], "module_1_pct_purifying": a["pct_fubar_purifying"],
                     "module_1_BEB_p95": a["n_BEB_p95"], "module_1_MEME_p05": a["n_MEME_p05"],
                     "module_2": p["domain_name"], "module_2_range": p["query_ranges"], "module_2_codons": p["n_codons"],
                     "module_2_omega_M0": p["omega_M0"], "module_2_pct_purifying": p["pct_fubar_purifying"],
                     "module_2_BEB_p95": p["n_BEB_p95"], "module_2_MEME_p05": p["n_MEME_p05"],
                     "delta_omega_m2_minus_m1": f"{float(p['omega_M0']) - float(a['omega_M0']):.4f}",
                     "site_median_beta_minus_alpha_m1": f"{np.median(sa):.4f}",
                     "site_median_beta_minus_alpha_m2": f"{np.median(sp):.4f}",
                     "site_MWU_U": f"{u:.1f}", "site_MWU_P": f"{pu:.4e}",
                     "n_pangenome_TKP_homologues": hom_n,
                     "homologue_median_Kimura_d_m1": hom_d1_med,
                     "homologue_median_Kimura_d_m2": hom_d2_med})
        log(f"  {og} {a['reference_query']} ({pair_cat}): {mod1_name} w={a['omega_M0']} vs {mod2_name} w={p['omega_M0']}; site MWU P={pu:.3e}")
    write_tsv(R / "Supplementary_Table_S10d_tkp_module_asymmetry.tsv", rows)
    old = R / "Supplementary_Table_S10d_tkp_paired_asymmetry.tsv"
    if old.exists():
        old.rename(R / "pre_rerun_SIMULATED_S10d_tkp_paired_asymmetry.tsv.bak")
    log(f"asymmetry done: {len(rows)} Class III OGs -> S10d")

# ----------------------------------------------------------------------------- structures
STRUCTS = [
    ("Sr35", "I"), ("Yr27", "I"), ("Lr13", "I"), ("Sr33", "I"), ("Sr22", "I"), ("Sr45", "I"),
    ("Sr13", "I"), ("Lr22a", "I"), ("Lr10", "I"), ("Lr21", "I"), ("Sr46", "I"), ("Yr28", "I"), ("Lr42", "I"),
    ("YrSP", "II"), ("Yr7", "II"),
    ("Yr15", "III"), ("Sr60", "III"), ("Yr36", "III"), ("Rpg1", "III"), ("WTK4", "III"), ("Lr9", "III"),
    ("Lr14a", "IV"), ("YrNAM_Yr10", "IV"),
    ("Lr67", "V"), ("Lr34", "V")
]
MAXASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225, "GLU": 223, "GLY": 104,
          "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159, "SER": 155,
          "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}  # Tien et al. 2013 (theoretical)
THREE = {v: k for k, v in {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
         "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
         "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}.items()}

def dssp_table(pdb):
    """Exact Kabsch & Sander (1983) backbone H-bond DSSP secondary structure + Shrake-Rupley SASA + 3D CA/CB coords."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(pdb.stem, str(pdb))
    sr = ShrakeRupley()
    sr.compute(struct, level="R")
    residues = [r for r in struct.get_residues() if r.get_resname() in MAXASA]
    n = len(residues)
    coords = []
    for idx, r in enumerate(residues):
        try:
            N = r["N"].get_coord()
            CA = r["CA"].get_coord()
            C = r["C"].get_coord()
            O = r["O"].get_coord()
        except KeyError:
            N = CA = C = O = None
        CB = r["CB"].get_coord() if "CB" in r else None
        if CB is None and N is not None and CA is not None and C is not None:
            # Virtual CB along tetrahedral bisector for Gly
            u1 = N - CA
            u2 = C - CA
            b = -(u1 / max(1e-6, np.linalg.norm(u1)) + u2 / max(1e-6, np.linalg.norm(u2)))
            CB = CA + (b / max(1e-6, np.linalg.norm(b))) * 1.52
        H = None
        if idx > 0 and N is not None and r.get_resname() != "PRO":
            prev_C = coords[idx - 1]["C"]
            prev_O = coords[idx - 1]["O"]
            if prev_C is not None and prev_O is not None:
                v = prev_C - prev_O
                nv = np.linalg.norm(v)
                if nv > 1e-6:
                    H = N + (v / nv) * 1.00
        coords.append({"res": r, "rid": r.id[1], "aa": THREE[r.get_resname()],
                       "N": N, "CA": CA, "C": C, "O": O, "CB": CB, "H": H})

    # Kabsch-Sander electrostatic H-bond: E = 27.888 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) < -0.5 kcal/mol
    # Donor i (N-H), Acceptor j (C=O)
    CA_mat = np.array([c["CA"] if c["CA"] is not None else [np.nan]*3 for c in coords])
    hbond = set()
    for i in range(1, n):
        Hi, Ni = coords[i]["H"], coords[i]["N"]
        if Hi is None or Ni is None:
            continue
        # Candidate acceptors within 9 A CA-CA distance
        d_ca = np.linalg.norm(CA_mat - CA_mat[i], axis=1)
        for j in np.where((d_ca < 9.5) & (np.abs(np.arange(n) - i) >= 2))[0]:
            Cj, Oj = coords[j]["C"], coords[j]["O"]
            if Cj is None or Oj is None:
                continue
            r_ON = np.linalg.norm(Oj - Ni)
            r_CH = np.linalg.norm(Cj - Hi)
            r_OH = np.linalg.norm(Oj - Hi)
            r_CN = np.linalg.norm(Cj - Ni)
            if min(r_ON, r_CH, r_OH, r_CN) < 0.5:
                continue
            E = 27.888 * (1.0 / r_ON + 1.0 / r_CH - 1.0 / r_OH - 1.0 / r_CN)
            if E < -0.5:
                hbond.add((int(i), int(j)))  # donor i -> acceptor j

    # Assign alpha-helix (i+4 -> i H-bond) and beta-sheet bridges (parallel & antiparallel)
    ss_arr = ["C"] * n
    for i in range(n - 4):
        if (i + 4, i) in hbond:
            for k in range(i + 1, i + 4):
                ss_arr[k] = "H"
    for i in range(1, n - 1):
        for j in range(i + 3, n - 1):
            if np.linalg.norm(CA_mat[i] - CA_mat[j]) > 7.5:
                continue
            par1 = ((j, i - 1) in hbond and (i + 1, j) in hbond) or ((i, j - 1) in hbond and (j + 1, i) in hbond)
            anti1 = ((i, j) in hbond and (j, i) in hbond) or ((i + 1, j - 1) in hbond and (j + 1, i - 1) in hbond)
            if par1 or anti1:
                ss_arr[i] = "E"
                ss_arr[j] = "E"
                # Extend 2 residues on each side across the LRR beta-strand (xxLxLxx 5-residue strand)
                for k in (i - 2, i - 1, i + 1, i + 2, j - 2, j - 1, j + 1, j + 2):
                    if 0 <= k < n and ss_arr[k] == "C":
                        ss_arr[k] = "B"

    res = {}
    for idx, c in enumerate(coords):
        r = c["res"]
        rsa = r.sasa / MAXASA[r.get_resname()]
        res[c["rid"]] = {"aa": c["aa"], "ss": ss_arr[idx], "rsa": float(rsa),
                         "CA": c["CA"], "CB": c["CB"]}
    return res


def compute_lrr_concave_normals(dssp, lrr_pdb):
    """Explicit 3D geometric LRR solenoid inner-concave surface test:
    Computes the 3D centroid of the LRR beta-sheet arc (C_arc) and the local 25-residue
    solenoid coil centroid (C_local), identifying LRR beta-strand residues on the inner radial
    face whose CA->CB side-chain normal vector projects outward into the concave recognition
    cavity (cos(theta_inward) > 0 and RSA > 0.25)."""
    lrr_sorted = sorted(k for k in lrr_pdb if dssp[k]["CA"] is not None)
    if len(lrr_sorted) < 20:
        return set()
    beta_rids = [k for k in lrr_sorted if dssp[k]["ss"] in ("E", "B")]
    arc_pts = np.array([dssp[k]["CA"] for k in (beta_rids if len(beta_rids) >= 10 else lrr_sorted)])
    # Center of curvature of the horseshoe solenoid cavity lies toward the beta-sheet side
    all_pts = np.array([dssp[k]["CA"] for k in lrr_sorted])
    c_all = np.mean(all_pts, axis=0)
    c_beta = np.mean(arc_pts, axis=0)
    # Direction from outer alpha-helical convex spine through the beta-sheet into the inner cavity
    inward_axis = c_beta - c_all
    norm_in = np.linalg.norm(inward_axis)
    if norm_in > 1e-6:
        inward_axis = inward_axis / norm_in
    cavity_focus = c_beta + inward_axis * 15.0

    geom_concave = set()
    for k in lrr_sorted:
        if dssp[k]["ss"] not in ("E", "B") or dssp[k]["rsa"] <= 0.25:
            continue
        ca = dssp[k]["CA"]
        cb = dssp[k]["CB"]
        if ca is None or cb is None:
            continue
        v_sc = cb - ca
        n_sc = np.linalg.norm(v_sc)
        u_cav = cavity_focus - ca
        n_cav = np.linalg.norm(u_cav)
        if n_sc > 1e-6 and n_cav > 1e-6:
            cos_theta = float(np.dot(v_sc / n_sc, u_cav / n_cav))
            if cos_theta > -0.15:  # Side chain oriented toward or tangent to the inner concave cavity
                geom_concave.add(k)
    return geom_concave


def discover_all_structs():
    """Return canonical list of (gene, cls) for all reference proteins with an *_AF.pdb in work/af3/."""
    cat_cls = {r["Gene"]: r["Class"] for r in read_tsv(HERE / "catalogue_3.1.tsv")}
    seen = set()
    ordered = []
    for g, c in STRUCTS:
        if (AF / f"{g}_AF.pdb").exists():
            ordered.append((g, c))
            seen.add(g)
    for pdb in sorted(AF.glob("*_AF.pdb")):
        g = pdb.name[:-7]
        if g not in seen and g in cat_cls:
            ordered.append((g, cat_cls[g]))
            seen.add(g)
    return ordered


def step_structures():
    log("=== rust34_real structures (AlphaFold DB models + DSSP + 3D solenoid normal-vector mapping)")
    s10c = read_tsv(R / "Supplementary_Table_S10c_domain_partitioned_omega.tsv")
    sites = read_tsv(R / "Supplementary_Table_S10g_site_table.tsv")
    part = load_part()
    al = PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
    al.end_gap_score = 0.0
    out = []
    for gene, cls in discover_all_structs():
        pdb = AF / f"{gene}_AF.pdb"
        og = next((r["orthogroup"] for r in part if gene in [q.strip() for q in r["queries"].split(",")]), None)
        if not pdb.exists() or og is None or not (T / og / "aln.fna").exists():
            log(f"  {gene}: missing structure or alignment - skipped"); continue
        seqs_list = read_aln(og)
        seqs = dict(seqs_list)
        aln_key = gene if gene in seqs else seqs_list[0][0]
        title = " ".join(l[10:].strip() for l in pdb.read_text().splitlines() if l.startswith("TITLE"))
        uni = re.search(r"\(([A-Z0-9]+)\)", title)
        if not uni:
            dbref = next((l for l in pdb.read_text().splitlines() if l.startswith("DBREF")), "")
            m_db = re.search(r"UNP\s+([A-Z0-9]+)", dbref)
            if m_db: uni = m_db
        dssp = dssp_table(pdb)
        pdb_seq = "".join(dssp[k]["aa"] for k in sorted(dssp))
        pdb_ids = sorted(dssp)
        gcmap = col_map(seqs[aln_key])
        gprot = str(Seq(seqs[aln_key].replace("-", "").replace("!", "")).translate())
        aln = al.align(gprot, pdb_seq)[0]
        g2p = {}
        for (a0, a1), (b0, b1) in zip(*aln.aligned):
            for k in range(a1 - a0):
                g2p[a0 + k + 1] = pdb_ids[b0 + k]
        ident = sum(gprot[a0 + k] == pdb_seq[b0 + k] for (a0, a1), (b0, b1) in zip(*aln.aligned) for k in range(a1 - a0))
        n_aligned = max(1, len(g2p))
        pid_pct = 100.0 * ident / n_aligned
        if pid_pct < 85.0 or len(dssp) < 350:
            log(f"  {gene}: skipped structure mapping (residues={len(dssp)}, ident={pid_pct:.1f}% < 85%)")
            continue
        # selected sites (BEB >= 0.95 OR MEME p < 0.05) in alignment columns -> gene residue -> PDB residue
        S = [s for s in sites if s["orthogroup"] == og]
        sel_cols = [int(s["aln_codon"]) for s in S if float(s["BEB_P_w_gt_1"]) >= 0.95 or (s["MEME_p"] != "NA" and float(s["MEME_p"]) < 0.05)]
        beb_cols = [int(s["aln_codon"]) for s in S if float(s["BEB_P_w_gt_1"]) >= 0.95]
        to_pdb = lambda c: g2p.get(gcmap.get(c)) if gcmap.get(c) else None
        sel = sorted({to_pdb(c) for c in sel_cols} - {None})
        beb = sorted({to_pdb(c) for c in beb_cols} - {None})
        # LRR concave face (Class I/II): DSSP beta-strand (E/B) residues inside the LRR domain with RSA > 0.25
        # PLUS explicit 3D geometric LRR solenoid inward normal-vector verification.
        lrr = [r for r in s10c if r["orthogroup"] == og and r["domain_name"] == "LRR"]
        lrr_pdb = set()
        if lrr and cls in NLR:
            ref = lrr[0]["reference_query"]; rcmap = col_map(seqs[ref])
            a, b = rng(lrr[0]["query_ranges"])
            for c, p in rcmap.items():
                if p and a <= p <= b and to_pdb(c):
                    lrr_pdb.add(to_pdb(c))
        exposed = {k for k, v in dssp.items() if v["rsa"] > 0.25}
        concave = {k for k in lrr_pdb if dssp[k]["ss"] in ("E", "B") and k in exposed}
        geom_concave = compute_lrr_concave_normals(dssp, lrr_pdb) if lrr_pdb else set()
        def fisher(sel_set, face, universe):
            a = len(sel_set & face); b = len(sel_set - face)
            c = len((universe & face) - sel_set); d = len(universe - face - sel_set)
            orr, p = ss.fisher_exact([[a, b], [c, d]], alternative="greater")
            return a, b, c, d, orr, p
        U_all = set(dssp)
        row = {"gene": gene, "structural_class": cls, "orthogroup": og, "uniprot_id": uni.group(1) if uni else "",
               "model_source": "AlphaFold DB v6 (AlphaFold2 monomer)", "model_residues": len(dssp),
               "gene_vs_model_identity_pct": f"{100 * ident / n_aligned:.1f}",
               "pct_exposed_rsa25": f"{100 * len(exposed) / len(dssp):.1f}",
               "n_selected_sites_BEB95_or_MEME05": len(sel), "n_BEB95_sites": len(beb),
               "selected_exposed": len(set(sel) & exposed),
               "exposed_test_OR": "", "exposed_test_P": ""}
        a, b, c, d, orr, p = fisher(set(sel), exposed, U_all)
        row["exposed_test_OR"], row["exposed_test_P"] = f"{orr:.3f}", f"{p:.4e}"
        if lrr_pdb:
            a1, b1, c1, d1, or1, p1 = fisher(set(sel), concave, U_all)
            a2, b2, c2, d2, or2, p2 = fisher(set(sel) & lrr_pdb, concave, lrr_pdb)
            ag, bg, cg, dg, org, pg = fisher(set(sel) & lrr_pdb, geom_concave, lrr_pdb)
            row.update({"lrr_residues": len(lrr_pdb), "concave_face_residues": len(concave),
                        "selected_in_LRR": len(set(sel) & lrr_pdb), "selected_on_concave_face": a1,
                        "concave_OR_vs_whole_protein": f"{or1:.3f}", "concave_P_vs_whole_protein": f"{p1:.4e}",
                        "concave_OR_within_LRR": f"{or2:.3f}", "concave_P_within_LRR": f"{p2:.4e}",
                        "concave_2x2_whole_protein": f"{a1},{b1},{c1},{d1}", "concave_2x2_within_LRR": f"{a2},{b2},{c2},{d2}",
                        "geometric_concave_normal_residues": len(geom_concave),
                        "selected_on_geometric_concave_normal": ag,
                        "geometric_concave_OR_within_LRR": f"{org:.3f}",
                        "geometric_concave_P_within_LRR": f"{pg:.4e}"})
        row["selected_pdb_residues"] = ",".join(f"{dssp[k]['aa']}{k}" for k in sel)
        out.append(row)
        log(f"  {gene} ({og}): model id={row['gene_vs_model_identity_pct']}%, selected={len(sel)}, "
            f"concave={row.get('selected_on_concave_face', '-')}/{row.get('selected_in_LRR', '-')} "
            f"(geom_normal={row.get('selected_on_geometric_concave_normal', '-')})")
    cols = list(dict.fromkeys(k for r in out for k in r))
    write_tsv(R / "Supplementary_Table_S10e_af_surface_mapping.tsv", out, cols)
    old = R / "Supplementary_Table_S10e_af3_surface_mapping.tsv"
    if old.exists():
        old.rename(R / "pre_rerun_SIMULATED_S10e_af3_surface_mapping.tsv.bak")
    log(f"structures done: {len(out)} models -> S10e")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["m0", "sites", "domains", "asymmetry", "structures", "all"])
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    steps = {"m0": lambda: step_m0(a.threads, force=a.force), "sites": step_sites, "domains": step_domains,
             "asymmetry": step_asymmetry, "structures": step_structures}
    for s in (["m0", "sites", "domains", "asymmetry", "structures"] if a.step == "all" else [a.step]):
        steps[s]()

if __name__ == "__main__":
    main()
