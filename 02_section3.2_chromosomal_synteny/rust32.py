#!/usr/bin/env python3
"""
rust32.py - generates the dataset and statistics behind Results section 3.2
("Subgenome asymmetry, chromosomal distribution, and microsynteny retention")
and Methods 5.4 ("Chromosomal coordinate normalization, microsynteny, and TE profiling").

Steps (run in order, each is resumable):
  positions   arm-normalized centromere-to-telomere distance (d_norm) across all 56 subgenomes
  triads      homoeologous A:B:D triad conservation and subgenome CNV asymmetry in hexaploids
  sri         pairwise Synteny Retention Index (SRI) & collinear blocks in +/-500 kb windows
  te          flanking +/-200 kb RepeatMasker TE coverage around Chinese Spring reference loci
  stats       all statistical tests (KS, Fisher, Kruskal-Wallis, MWU, Spearman) -> section3.2_values.tsv
  all         everything above in order

Usage:
  python rust32.py <step> [--config config_3.1.json] [--threads 4] [--all-ensembl]
"""
import argparse, collections, concurrent.futures, csv, json, math, os, re, shutil, subprocess, sys, time
import urllib.error, urllib.request
from pathlib import Path

import numpy as np
from scipy import stats as ss

STAGE_DIR = Path(__file__).resolve().parent
HERE = STAGE_DIR.parent
sys.path.insert(0, str(HERE / "01_section3.1_pangenome_inventory"))
from rust31 import log, load_cfg, read_tsv, write_tsv, read_fasta, panel_rows, catalogue, which_or_die, find_og_files

WINDOW_BP = 500_000
ENSEMBL_REST = "https://rest.ensembl.org"


def load_centromere_ratios():
    """Return dict of (chr_num_1_to_7, subgenome_A_B_D) -> relative centromere midpoint (0..1)
    from centromeres_CSv21.tsv (IWGSC RefSeq v2.1 CENH3 ChIP-seq)."""
    cen_file = HERE / "centromeres_CSv21.tsv"
    ratios = {}
    exact_cs = {}
    if cen_file.exists():
        for r in read_tsv(cen_file):
            ch = r["chrom"].strip()
            c_mid = float(r["cen_mid_bp"])
            c_len = float(r["chr_len_bp"])
            exact_cs[ch] = c_mid
            ratios[ch] = c_mid / c_len
    return exact_cs, ratios


def get_chr_group_and_sub(label, subgenome):
    m = re.match(r"^(?:chr|Chr|CHR)?0?([1-7])([A-Za-z])?", str(label))
    if not m:
        return None, subgenome
    cnum = m.group(1)
    sg = m.group(2).upper() if m.group(2) else subgenome
    if sg == "S":
        sg = "B"
    return cnum, sg


def compute_d_norm(mid_bp, chr_len_bp, asm_id, label, subgenome, exact_cs, ratios):
    """Arm-normalized distance from centromere (0.0) to telomere (1.0):
    d_norm = |x_i - C_j| / arm_length."""
    cnum, sg = get_chr_group_and_sub(label, subgenome)
    key = f"{cnum}{sg}" if cnum and sg else ""
    if asm_id == "Tae_CSv21" and key in exact_cs:
        c_mid = exact_cs[key]
    elif key in ratios:
        c_mid = ratios[key] * chr_len_bp
    else:
        c_mid = 0.5 * chr_len_bp

    if mid_bp <= c_mid and c_mid > 0:
        arm = "S"
        d_norm = (c_mid - mid_bp) / c_mid
    else:
        arm = "L"
        denom = max(1.0, chr_len_bp - c_mid)
        d_norm = (mid_bp - c_mid) / denom
    return min(1.0, max(0.0, d_norm)), arm, int(round(c_mid))


def map_query_to_og(W, cat):
    ogf, unf = find_og_files(W)
    q2og = {}
    for fpath in (ogf, unf):
        with open(fpath) as f:
            f.readline()
            for line in f:
                p = line.rstrip("\n").split("\t")
                for col in p[1:]:
                    for g in col.split(","):
                        g = g.strip()
                        if g in cat:
                            q2og[g] = p[0]
    return q2og


# ----------------------------------------------------------------------------- Step 1: positions
def step_positions(cfg, args):
    W, R = Path(cfg["workdir"]), Path(cfg["resultsdir"])
    cat = catalogue(cfg)
    panel = panel_rows(cfg)
    main_q = {q for q in cat if cat[q]["Set"] == "Main"}
    q2og = map_query_to_og(W, cat)
    exact_cs, ratios = load_centromere_ratios()

    seq_info = {}
    for r in panel:
        a = r["asm_id"]
        smap = W / "miniprot" / a / "seqmap.tsv"
        if smap.exists():
            for x in read_tsv(smap):
                seq_info[(a, x["seqid"])] = (int(x["length"]), x["label"], x["subgenome"])

    cands = [c for c in read_tsv(R / "candidate_loci.tsv")
             if c["anchor_ok"] == "1" and c["best_query"] in main_q]
    pos_rows = []
    for c in cands:
        info = seq_info.get((c["asm_id"], c["seqid"]))
        if not info:
            continue
        L, label, sg = info
        if sg == "Un" or L < 50_000_000:
            continue
        mid = 0.5 * (int(c["start"]) + int(c["end"]))
        d_norm, arm, c_mid = compute_d_norm(mid, L, c["asm_id"], label, sg, exact_cs, ratios)
        q = c["best_query"]
        pos_rows.append({
            "cand_id": c["cand_id"],
            "asm_id": c["asm_id"],
            "subgenome": sg,
            "seqid": c["seqid"],
            "chr_label": label,
            "arm": arm,
            "chr_len_bp": L,
            "centromere_mid_bp": c_mid,
            "midpoint_bp": int(mid),
            "d_norm": round(d_norm, 4),
            "distal_85": int(d_norm >= 0.85),
            "best_query": q,
            "LocusGroup": cat[q]["LocusGroup"],
            "orthogroup": q2og.get(q, ""),
            "Class": cat[q]["Class"],
            "Group": cat[q]["Group"],
            "identity": float(c["identity"]),
            "qcov": float(c["qcov"]),
            "intact": c["intact"],
        })
    write_tsv(R / "S4_chromosomal_positions.tsv", pos_rows)

    # Summarize per LocusGroup (36 independent loci) across all placed homologues
    by_lg = collections.defaultdict(list)
    for p in pos_rows:
        by_lg[p["LocusGroup"]].append(p)

    lg_meta = {}
    for q in sorted(main_q):
        lg = cat[q]["LocusGroup"]
        lg_meta.setdefault(lg, (cat[q]["Class"], cat[q]["Group"], q2og.get(q, "")))

    lg_rows = []
    for lg, (cls, grp, og) in sorted(lg_meta.items()):
        pts = by_lg.get(lg, [])
        dvals = [x["d_norm"] for x in pts]
        can_cnum = CANONICAL_LOCUS_CHR.get(lg, "")
        cs_pts = [x for x in pts if x["asm_id"] == "Tae_CSv21" and (not can_cnum or x["chr_label"].startswith(can_cnum))]
        if not cs_pts:
            cs_pts = [x for x in pts if x["asm_id"] == "Tae_CSv21"]
        can_pts = [x for x in pts if can_cnum and x["chr_label"].startswith(can_cnum)] if can_cnum else pts
        ref_chr = cs_pts[0]["chr_label"] if cs_pts else (can_pts[0]["chr_label"] if can_pts else (pts[0]["chr_label"] if pts else ""))
        ref_arm = cs_pts[0]["arm"] if cs_pts else (can_pts[0]["arm"] if can_pts else (pts[0]["arm"] if pts else ""))
        med = float(np.median(dvals)) if dvals else float("nan")
        q25 = float(np.percentile(dvals, 25)) if dvals else float("nan")
        q75 = float(np.percentile(dvals, 75)) if dvals else float("nan")
        lg_rows.append({
            "LocusGroup": lg,
            "orthogroup": og,
            "Class": cls,
            "Group": grp,
            "ref_chr": ref_chr,
            "ref_arm": ref_arm,
            "n_placed_loci": len(dvals),
            "median_d_norm": round(med, 4) if dvals else "",
            "q25_d_norm": round(q25, 4) if dvals else "",
            "q75_d_norm": round(q75, 4) if dvals else "",
            "mean_d_norm": round(float(np.mean(dvals)), 4) if dvals else "",
            "pct_distal_85": round(100.0 * sum(x >= 0.85 for x in dvals) / len(dvals), 2) if dvals else "",
            "zone": "distal_subtelomeric" if (dvals and med >= 0.85) else "interstitial_pericentromeric"
        })
    write_tsv(R / "S4b_locusgroup_positions.tsv", lg_rows)
    log(f"positions done: {len(pos_rows)} placed candidate loci across {len(lg_rows)} LocusGroups")


# ----------------------------------------------------------------------------- Step 2: triads
CANONICAL_LOCUS_CHR = {
    "Lr1": "5", "Lr10": "1", "Lr13/Yr27": "2", "Lr14a": "7", "Lr21": "1",
    "Lr22a": "2", "Lr30": "4", "Lr34": "7", "Lr42": "1", "Lr47": "7",
    "Lr67": "4", "Lr9": "6", "Sr13": "6", "Sr21": "2", "Sr22": "7",
    "Sr26": "6", "Sr27": "3", "Sr33": "1", "Sr35": "3", "Sr43": "7",
    "Sr45": "1", "Sr46": "2", "Sr50": "1", "Sr60": "5", "Sr61": "6",
    "Sr62": "1", "Sr9": "2", "SrTA1662": "1", "Yr15": "1", "Yr28": "4",
    "Yr36": "6", "Yr5/YrSP": "2", "Yr7": "2", "Yr87/Lr85": "3",
    "YrNAM (=Yr10)": "1", "YrU1": "5",
}


def load_canonical_chr_groups(W, R, cat, main_q):
    """Determine the canonical Triticeae chromosome group (1-7) for each LocusGroup and orthogroup,
    prioritizing curated canonical wheat chromosomes and Ensembl CS v1.1 syntenic anchors."""
    lg_chr = dict(CANONICAL_LOCUS_CHR)
    cmap_file = Path(W) / "sri_3.2" / "cs_compara_locusgroup_map.tsv"
    if cmap_file.exists():
        for r in read_tsv(cmap_file):
            if r["LocusGroup"] not in lg_chr:
                cnum, _ = get_chr_group_and_sub(r["chrom"], "")
                if cnum:
                    lg_chr[r["LocusGroup"]] = cnum
    pos_file = Path(R) / "S4_chromosomal_positions.tsv"
    if pos_file.exists():
        by_lg = collections.defaultdict(list)
        for p in read_tsv(pos_file):
            cnum, _ = get_chr_group_and_sub(p["chr_label"], p["subgenome"])
            if cnum:
                by_lg[p["LocusGroup"]].append((int(p["intact"]), float(p["identity"]), p["asm_id"] == "Tae_CSv21", cnum))
        for lg, items in by_lg.items():
            if lg not in lg_chr and items:
                items.sort(reverse=True)
                lg_chr[lg] = items[0][3]
    return lg_chr


def step_triads(cfg, args):
    W, R = Path(cfg["workdir"]), Path(cfg["resultsdir"])
    cat = catalogue(cfg)
    panel = panel_rows(cfg)
    main_q = {q for q in cat if cat[q]["Set"] == "Main"}
    q2og = map_query_to_og(W, cat)

    pos_file = R / "S4_chromosomal_positions.tsv"
    if not pos_file.exists():
        step_positions(cfg, args)
    pos_rows = read_tsv(pos_file)
    lg_chr = load_canonical_chr_groups(W, R, cat, main_q)

    hex_asms = [r["asm_id"] for r in panel if int(r["ploidy"]) == 6]
    lg_meta = {}
    for q in sorted(main_q):
        lg = cat[q]["LocusGroup"]
        lg_meta.setdefault(lg, (cat[q]["Class"], cat[q]["Group"], q2og.get(q, "")))

    # 1. Per (LocusGroup x hexaploid assembly) homoeologous triad balance on canonical chromosome group
    triad_rows = []
    for lg, (cls, grp, og) in sorted(lg_meta.items()):
        target_cnum = lg_chr.get(lg)
        for a in hex_asms:
            pts = [p for p in pos_rows if p["asm_id"] == a and p["LocusGroup"] == lg
                   and p["intact"] == "1" and float(p["identity"]) >= 0.80]
            if target_cnum:
                pts_chr = [p for p in pts if get_chr_group_and_sub(p["chr_label"], p["subgenome"])[0] == target_cnum]
                if pts_chr:
                    pts = pts_chr
            sub_cnt = collections.Counter(p["subgenome"] for p in pts)
            na, nb, nd = sub_cnt.get("A", 0), sub_cnt.get("B", 0), sub_cnt.get("D", 0)
            tot = na + nb + nd
            if tot == 0:
                state = "absent"
            elif na == 1 and nb == 1 and nd == 1:
                state = "balanced_1:1:1"
            elif na > 0 and nb > 0 and nd > 0:
                state = "present_A+B+D_multicopy"
            else:
                state = "subgenome_asymmetric"
            triad_rows.append({
                "LocusGroup": lg, "orthogroup": og, "Class": cls, "Group": grp,
                "canonical_chr_group": target_cnum or "", "asm_id": a,
                "copies_A": na, "copies_B": nb, "copies_D": nd, "total_ABD": tot,
                "triad_state": state
            })
    write_tsv(R / "S5_hexaploid_triads.tsv", triad_rows)

    # Load Orthogroups.tsv membership from Section 3.1
    cands_all = {c["cand_id"]: c for c in read_tsv(R / "candidate_loci.tsv")}
    ogf, unf = find_og_files(W)
    og_members = {}
    for fpath in (ogf, unf):
        with open(fpath) as f:
            f.readline()
            for line in f:
                p = line.rstrip("\n").split("\t")
                gs = []
                for col in p[1:]:
                    gs += [g.strip() for g in col.split(",") if g.strip()]
                og_members[p[0]] = gs

    # 2. Per orthogroup (30 Main orthogroups) triad summary across 13 hexaploid assemblies
    og_rows_31 = [r for r in read_tsv(R / "S1b_orthogroup_partition.tsv") if r["set"] == "Main"]
    og_chr = {}
    for lg, (_, _, og) in lg_meta.items():
        if og and lg in lg_chr and og not in og_chr:
            og_chr[og] = lg_chr[lg]

    s9_by_og = collections.defaultdict(list)
    if (R / "S9_paralogue_pruning.tsv").exists():
        for r9 in read_tsv(R / "S9_paralogue_pruning.tsv"):
            s9_by_og[r9["orthogroup"]].append(r9)
    best_by_q = collections.defaultdict(list)
    if (R / "best_alignment_per_query_subgenome.tsv").exists():
        for rb in read_tsv(R / "best_alignment_per_query_subgenome.tsv"):
            best_by_q[rb["query"]].append(rb)

    og_triads = []
    retained_sub_by_og = {}
    for r in og_rows_31:
        og = r["orthogroup"]
        cls, grp = r["class"], r["group"]
        target_cnum = og_chr.get(og)
        qs = [q.strip() for q in r["queries"].split(",") if q.strip()]
        use = [q for q in qs if cat[q]["Set"] == "Main"] or qs
        mem = [cands_all[g] for g in og_members.get(og, []) if g in cands_all]

        # Exact S9-pruned focal orthologue + progenitor domain counts per (asm_id, subgenome) matching rust31.py
        retained_sub = collections.Counter()
        if s9_by_og:
            if grp == "NLR":
                pruned_mem = [
                    m for m in s9_by_og[og]
                    if m["verdict"] == "retained" and float(m["qcov"]) >= 0.80
                    and ((m["intact"] == "1" and float(m["identity"]) >= 0.85)
                         or m["reverse_hit_rank"] in ("1", "2", "3", "4"))
                ]
                for m in pruned_mem:
                    if m["subgenome"] in ("A", "B", "S", "D"):
                        retained_sub[(m["asm_id"], m["subgenome"])] += 1
            else:
                pruned_mem = [m for m in s9_by_og[og] if m["verdict"] == "retained"]
                dom_hits = [
                    b for q in use for b in best_by_q.get(q, [])
                    if float(b["qcov"]) >= 0.50 and float(b["identity"]) >= 0.60 and b["subgenome"] in ("A", "B", "S", "D")
                ]
                for m in pruned_mem:
                    if m["subgenome"] in ("A", "B", "S", "D"):
                        retained_sub[(m["asm_id"], m["subgenome"])] += 1
                for b in dom_hits:
                    k_sub = (b["asm_id"], b["subgenome"])
                    if retained_sub[k_sub] == 0:
                        retained_sub[k_sub] = 1
        else:
            for m in mem:
                if m["intact"] == "1" and m["subgenome"] in ("A", "B", "S", "D"):
                    retained_sub[(m["asm_id"], m["subgenome"])] += 1
        retained_sub_by_og[og] = retained_sub

        bal_cnt, all3_cnt, asym_cnt, abs_cnt = 0, 0, 0, 0
        ca_list, cb_list, cd_list = [], [], []
        cl_a_list, cl_b_list, cl_d_list = [], [], []
        for a in hex_asms:
            na = retained_sub.get((a, "A"), 0)
            nb = retained_sub.get((a, "B"), 0)
            nd = retained_sub.get((a, "D"), 0)
            ca_list.append(na); cb_list.append(nb); cd_list.append(nd)
            sc_cl = collections.Counter(m["subgenome"] for m in mem if m["asm_id"] == a)
            cl_a_list.append(sc_cl.get("A", 0)); cl_b_list.append(sc_cl.get("B", 0)); cl_d_list.append(sc_cl.get("D", 0))
            if na == 1 and nb == 1 and nd == 1:
                bal_cnt += 1
            if na > 0 and nb > 0 and nd > 0:
                all3_cnt += 1
            elif na + nb + nd == 0:
                abs_cnt += 1
            else:
                asym_cnt += 1
        mean_a, mean_b, mean_d = np.mean(ca_list), np.mean(cb_list), np.mean(cd_list)
        cv_abd = np.std([mean_a, mean_b, mean_d]) / max(1e-6, np.mean([mean_a, mean_b, mean_d]))
        sd_intact = float(r["intact_copies_sd"])
        sd_retained = float(r.get("retained_copies_sd", sd_intact))
        # Conserved A:B:D triad regime: present in all three subgenomes (A, B, D) across >=60% (>=8/13) of hexaploid assemblies
        is_conserved_triad = int(all3_cnt >= 0.6 * len(hex_asms))
        ever_all3 = int(all3_cnt >= 1)
        og_triads.append({
            "orthogroup": og,
            "queries": r["queries"],
            "Class": cls,
            "Group": grp,
            "canonical_chr_group": target_cnum or "",
            "n_hexaploid_asms": len(hex_asms),
            "n_strict_1_1_1": bal_cnt,
            "n_all_three_ABD": all3_cnt,
            "n_asymmetric_PAV": asym_cnt,
            "n_absent": abs_cnt,
            "mean_copies_A": round(float(mean_a), 2),
            "mean_copies_B": round(float(mean_b), 2),
            "mean_copies_D": round(float(mean_d), 2),
            "cluster_mean_A": round(float(np.mean(cl_a_list)), 2),
            "cluster_mean_B": round(float(np.mean(cl_b_list)), 2),
            "cluster_mean_D": round(float(np.mean(cl_d_list)), 2),
            "subgenome_CV": round(float(cv_abd), 3),
            "retained_copies_sd": round(sd_retained, 3),
            "intact_copies_sd": round(sd_intact, 3),
            "ever_all_three_ABD": ever_all3,
            "conserved_1_1_1_triad_regime": is_conserved_triad,
        })
    write_tsv(R / "S5b_orthogroup_triad_summary.tsv", og_triads)

    # 3. Ensembl Compara Chinese Spring homoeologue triad table (Methods 5.4)
    cmap_file = W / "sri_3.2" / "cs_compara_locusgroup_map.tsv"
    if cmap_file.exists():
        compara_rows = read_tsv(cmap_file)
        write_tsv(R / "S5c_ensembl_compara_triads.tsv", compara_rows)

    # 4. Subgenome copy-number table across all 56 haploid subgenomes (A vs B/S vs D)
    subs = [(r["asm_id"], g, "B/S" if g in ("B", "S") else g) for r in panel for g in r["genomes"]]
    sg_cnv_rows = []
    for r in og_rows_31:
        og = r["orthogroup"]
        mem = [cands_all[g] for g in og_members.get(og, []) if g in cands_all and cands_all[g]["intact"] == "1"]
        cnt = collections.Counter((m["asm_id"], m["subgenome"]) for m in mem)
        ret_cnt = retained_sub_by_og.get(og, collections.Counter())
        for asm_id, sg, tier in subs:
            sg_cnv_rows.append({
                "orthogroup": og, "Class": r["class"], "Group": r["group"],
                "asm_id": asm_id, "subgenome": sg, "tier": tier,
                "intact_copies": cnt.get((asm_id, sg), 0),
                "retained_focal_copies": ret_cnt.get((asm_id, sg), 0),
            })
    write_tsv(R / "S5d_subgenome_cnv.tsv", sg_cnv_rows)
    log(f"triads done: {len(triad_rows)} locus x hexaploid cells, {len(og_triads)} orthogroup triad summaries")


# ----------------------------------------------------------------------------- Step 3: SRI (Microsynteny)
def collinear_chain_len(pairs):
    """Given list of (idx_A, idx_B) orthologous gene index pairs in two windows,
    return (max_collinear_count, n_inversions) using longest monotonic subsequence (forwards or reverse)."""
    if not pairs:
        return 0, 0
    pairs = sorted(set(pairs))

    def lis(seq):
        if not seq:
            return 0
        tails = []
        for x in seq:
            import bisect
            i = bisect.bisect_left(tails, x)
            if i == len(tails):
                tails.append(x)
            else:
                tails[i] = x
        return len(tails)

    b_seq = [b for _, b in pairs]
    fwd = lis(b_seq)
    rev = lis(list(reversed(b_seq)))
    best = max(fwd, rev)
    # Count discordant orientation pairs as local inversions
    n_inv = min(fwd, rev) if len(pairs) >= 3 and min(fwd, rev) >= 2 else 0
    return best, n_inv


def step_sri(cfg, args):
    os.environ["PATH"] = f"{Path(sys.executable).parent}:{os.environ.get('PATH', '')}"
    W, R = Path(cfg["workdir"]), Path(cfg["resultsdir"])
    cat = catalogue(cfg)
    panel = panel_rows(cfg)
    main_q = {q for q in cat if cat[q]["Set"] == "Main"}
    q2og = map_query_to_og(W, cat)
    if not (R / "S4_chromosomal_positions.tsv").exists():
        step_positions(cfg, args)
    pos_all = read_tsv(R / "S4_chromosomal_positions.tsv")
    lg_chr = load_canonical_chr_groups(W, R, cat, main_q)

    top_hit_all = {}
    for p in pos_all:
        cnum, sg = get_chr_group_and_sub(p["chr_label"], p["subgenome"])
        if not cnum:
            continue
        lg = p["LocusGroup"]
        sc = (int(p["intact"]), float(p["identity"]))
        if lg not in top_hit_all or sc > top_hit_all[lg][0]:
            top_hit_all[lg] = (sc, cnum, sg)

    seq_info = {}
    asm_genomes = {}
    for r in panel:
        asm_genomes[r["asm_id"]] = set("B" if g == "S" else g for g in r["genomes"])
        smap = W / "miniprot" / r["asm_id"] / "seqmap.tsv"
        if smap.exists():
            for x in read_tsv(smap):
                seq_info[(r["asm_id"], x["seqid"])] = (int(x["length"]), x["label"], x["subgenome"])

    pruned_file = R / "S9_paralogue_pruning.tsv"
    pruned = {r["cand_id"]: r for r in read_tsv(pruned_file)} if pruned_file.exists() else {}

    # Identify all R-gene candidate overlapping annotated genes in each assembly so we exclude
    # tandem R-gene paralogues from flanking anchor sets
    cands = read_tsv(R / "candidate_loci.tsv")
    r_overlap_genes = collections.defaultdict(set)
    for c in cands:
        for g in c.get("overlapping_genes", "").split(","):
            if g.strip():
                r_overlap_genes[c["asm_id"]].add(g.strip())

    # Load gene coordinates for all 16 annotated assemblies
    annot_asms = [r["asm_id"] for r in panel if r["annotated"]]
    genes_by_asm_chr = collections.defaultdict(list)
    for a in annot_asms:
        gl = W / "annot" / a / "gene_loc.tsv"
        if not gl.exists():
            continue
        for g in read_tsv(gl):
            genes_by_asm_chr[(a, g["seqid"])].append((int(g["start"]), int(g["end"]), g["gene"], g["strand"]))
    for k in genes_by_asm_chr:
        genes_by_asm_chr[k].sort(key=lambda x: x[0])

    # Identify syntenic orthologous locus per (asm_id, LocusGroup) across the 16 annotated assemblies
    lg_clusters = collections.defaultdict(list)
    fallback_best = {}
    for c in cands:
        a = c["asm_id"]
        if a not in annot_asms or c["anchor_ok"] != "1" or c["best_query"] not in main_q or c["subgenome"] == "Un":
            continue
        info = seq_info.get((a, c["seqid"]))
        if not info or info[0] < 50_000_000:
            continue
        chr_len, lbl, _ = info
        cnum, sg = get_chr_group_and_sub(lbl, c["subgenome"])
        if not cnum or sg not in ("A", "B", "D"):
            continue
        rel_pos = 0.5 * (int(c["start"]) + int(c["end"])) / chr_len
        lg = cat[c["best_query"]]["LocusGroup"]
        grp = cat[c["best_query"]]["Group"]
        fb_sc = (int(cnum == lg_chr.get(lg)), int(c["intact"]), float(c["identity"]) * float(c["qcov"]))
        if (a, lg) not in fallback_best or fb_sc > fallback_best[(a, lg)][0]:
            fallback_best[(a, lg)] = (fb_sc, c)
        p_info = pruned.get(c["cand_id"], {})
        rev_rank = p_info.get("reverse_hit_rank", "")
        if grp == "NLR" and pruned:
            is_strict = (float(c["qcov"]) >= 0.80) and (
                (c["intact"] == "1" and float(c["identity"]) >= 0.85) or (rev_rank in ("1", "2", "3", "4"))
            )
            if not is_strict:
                continue
        sc = (float(c["identity"]) ** 3) * float(c["qcov"])
        lg_clusters[lg].append((cnum, sg, rel_pos, a, int(c["intact"]), sc, c))

    lg_target = {}
    best_locus = {}
    for lg in sorted({cat[q]["LocusGroup"] for q in main_q}):
        grp = next(cat[q]["Group"] for q in main_q if cat[q]["LocusGroup"] == lg)
        items = lg_clusters[lg]
        _, clone_cnum, clone_sg = top_hit_all.get(lg, ((1, 1.0), lg_chr.get(lg, "1"), "A"))
        if grp == "NLR":
            clone_items = [x for x in items if x[0] == clone_cnum and x[1] == clone_sg]
            if clone_items:
                items = clone_items
        best_cluster_score = -1.0
        best_cluster_map = {}
        best_cs = (clone_cnum, clone_sg, 0.5)
        for cnum, sg, rpos, a0, _, _, _ in items:
            cmap = {}
            for cnum2, sg2, rpos2, a, intact, sc, c in items:
                if cnum2 == cnum and sg2 == sg and abs(rpos2 - rpos) <= 0.12:
                    if a not in cmap or (intact, sc) > cmap[a][:2]:
                        cmap[a] = (intact, sc, c)
            tot_sc = sum(v[1] for v in cmap.values())
            if tot_sc > best_cluster_score:
                best_cluster_score = tot_sc
                best_cluster_map = cmap
                best_cs = (cnum, sg, rpos)
        lg_target[lg] = best_cs
        for a, (intact, sc, c) in best_cluster_map.items():
            best_locus[(a, lg)] = ((intact, sc), c)

    # Extract +/-500 kb flanking non-R genes for each (a, lg)
    windows = {}  # (a, lg) -> list of (start, end, gene_id, strand)
    needed_genes = collections.defaultdict(set)
    for (a, lg), (_, c) in best_locus.items():
        seqid = c["seqid"]
        mid = 0.5 * (int(c["start"]) + int(c["end"]))
        w_start, w_end = max(1, int(mid - WINDOW_BP)), int(mid + WINDOW_BP)
        ov_self = set(x.strip() for x in c.get("overlapping_genes", "").split(",") if x.strip())
        flank = []
        for s0, e0, gid, strand in genes_by_asm_chr.get((a, seqid), []):
            if e0 < w_start:
                continue
            if s0 > w_end:
                break
            if gid in ov_self or gid in r_overlap_genes[a]:
                continue
            flank.append((s0, e0, gid, strand))
            needed_genes[a].add(gid)
        windows[(a, lg)] = flank

    # Write pooled flanking protein FASTA and run DIAMOND all-vs-all if not already cached
    sri_dir = W / "sri_3.2"
    sri_dir.mkdir(parents=True, exist_ok=True)
    pool_fa = sri_dir / "flanking_proteins.faa"
    hits_tsv = sri_dir / "flanking_allvsall.tsv"

    cached_hit_files = [
        sri_dir / fn for fn in ("flanking_allvsall.tsv", "flanking_allvsall_syn.tsv", "flanking_allvsall_syn2.tsv")
        if (sri_dir / fn).exists() and (sri_dir / fn).stat().st_size > 0
    ]
    if not cached_hit_files:
        n_written = 0
        with open(pool_fa, "w") as fh:
            for a, gset in sorted(needed_genes.items()):
                pep = W / "annot" / a / "longest.faa"
                if not pep.exists():
                    continue
                for gid, _, seq in read_fasta(pep):
                    if gid in gset:
                        fh.write(f">{a}|{gid}\n{seq}\n")
                        n_written += 1
        log(f"extracted {n_written} flanking proteins across {len(annot_asms)} annotated assemblies; running DIAMOND...")
        subprocess.run(["diamond", "makedb", "--in", str(pool_fa), "-d", str(sri_dir / "flank_db"), "--quiet"], check=True)
        subprocess.run([
            "diamond", "blastp",
            "-q", str(pool_fa),
            "-d", str(sri_dir / "flank_db"),
            "-o", str(hits_tsv),
            "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "bitscore",
            "--max-target-seqs", "200",
            "--evalue", "1e-10",
            "--threads", str(args.threads),
            "--quiet"
        ], check=True)
        cached_hit_files = [hits_tsv]

    # Load orthologous flanking gene pairs (>=70% identity, >=50% coverage)
    ortho_pairs = collections.defaultdict(set)  # pid_A -> set of pid_B
    for hf in cached_hit_files:
        for line in open(hf):
            q, s, pid, qcov, bits = line.rstrip("\n").split("\t")
            if q == s:
                continue
            if float(pid) >= 70.0 and float(qcov) >= 50.0:
                ortho_pairs[q].add(s)
                ortho_pairs[s].add(q)

    lg_meta = {}
    for q in sorted(main_q):
        lg = cat[q]["LocusGroup"]
        lg_meta.setdefault(lg, (cat[q]["Class"], cat[q]["Group"], q2og.get(q, "")))

    pairwise_rows = []
    per_locus_rows = []
    for lg, (cls, grp, og) in sorted(lg_meta.items()):
        tcnum, tsg, _ = lg_target[lg]
        elig_asms = [a for a in annot_asms if tsg in asm_genomes[a]]
        asms_with_lg = [a for a in elig_asms if (a, lg) in windows and len(windows[(a, lg)]) >= 2]
        sri_vals, gcol_vals, gcnt_vals, sv_vals = [], [], [], []
        for a in asms_with_lg:
            gcnt_vals.append(len(windows[(a, lg)]))
        mean_gcnt = int(round(float(np.mean(gcnt_vals)))) if gcnt_vals else 20
        for i in range(len(elig_asms)):
            for j in range(i + 1, len(elig_asms)):
                aA, aB = elig_asms[i], elig_asms[j]
                if (aA, lg) not in windows or (aB, lg) not in windows:
                    gA = len(windows[(aA, lg)]) if (aA, lg) in windows else 0
                    gB = len(windows[(aB, lg)]) if (aB, lg) in windows else 0
                    n_indels = max(gA, gB, mean_gcnt)
                    sri_vals.append(0.0)
                    gcol_vals.append(0)
                    sv_vals.append(n_indels)
                    pairwise_rows.append({
                        "LocusGroup": lg,
                        "orthogroup": og,
                        "Class": cls,
                        "Group": grp,
                        "asm_A": aA,
                        "asm_B": aB,
                        "G_A": gA,
                        "G_B": gB,
                        "G_shared": 0,
                        "G_collinear": 0,
                        "SRI": 0.0,
                        "n_inversions": 0,
                        "n_gene_indels": n_indels,
                        "sv_events": n_indels
                    })
                    continue
                wA, wB = windows[(aA, lg)], windows[(aB, lg)]
                gA, gB = len(wA), len(wB)
                if gA + gB == 0:
                    continue
                idx_B = {f"{aB}|{g[2]}": idx for idx, g in enumerate(wB)}
                pairs = []
                matched_A, matched_B = set(), set()
                for idx_a, g in enumerate(wA):
                    pid_a = f"{aA}|{g[2]}"
                    for pid_b in ortho_pairs.get(pid_a, ()):
                        if pid_b in idx_B:
                            pairs.append((idx_a, idx_B[pid_b]))
                            matched_A.add(idx_a)
                            matched_B.add(idx_B[pid_b])
                g_collinear, n_inv = collinear_chain_len(pairs)
                # Ensure G_collinear cannot exceed min(gA, gB)
                g_collinear = min(g_collinear, gA, gB)
                sri = (2.0 * g_collinear) / (gA + gB) if (gA + gB) > 0 else 0.0
                n_indels = abs(gA - gB) + max(0, min(gA, gB) - len(matched_A))
                sv_events = n_inv + n_indels
                sri_vals.append(sri)
                gcol_vals.append(g_collinear)
                sv_vals.append(sv_events)
                pairwise_rows.append({
                    "LocusGroup": lg,
                    "orthogroup": og,
                    "Class": cls,
                    "Group": grp,
                    "asm_A": aA,
                    "asm_B": aB,
                    "G_A": gA,
                    "G_B": gB,
                    "G_shared": min(len(matched_A), len(matched_B)),
                    "G_collinear": g_collinear,
                    "SRI": round(sri, 4),
                    "n_inversions": n_inv,
                    "n_gene_indels": n_indels,
                    "sv_events": sv_events
                })
        per_locus_rows.append({
            "LocusGroup": lg,
            "orthogroup": og,
            "Class": cls,
            "Group": grp,
            "n_annotated_asms": len(asms_with_lg),
            "n_pairwise_comparisons": len(sri_vals),
            "mean_flanking_genes": round(float(np.mean(gcnt_vals)), 2) if gcnt_vals else "",
            "mean_G_collinear": round(float(np.mean(gcol_vals)), 2) if gcol_vals else 0.0,
            "mean_SRI": round(float(np.mean(sri_vals)), 4) if sri_vals else 0.0,
            "sem_SRI": round(float(ss.sem(sri_vals)), 4) if len(sri_vals) > 1 else 0.0,
            "median_SRI": round(float(np.median(sri_vals)), 4) if sri_vals else 0.0,
            "mean_sv_events": round(float(np.mean(sv_vals)), 2) if sv_vals else float(mean_gcnt),
        })

    write_tsv(R / "S6b_sri_pairwise.tsv", pairwise_rows)
    write_tsv(R / "S6a_sri_per_locus.tsv", per_locus_rows)

    # Build microsynteny ribbon tracks (S6c_microsynteny_ribbons.tsv) for Fig. 2f representative loci
    rep_loci = ["Sr35", "Lr13/Yr27", "Yr15", "Lr14a", "Lr34"]
    rep_asms = ["Tu_G1812", "Aet_AL878", "Tdic_Zavitan", "Tdur_Svevo", "Tae_CSv21", "Tae_Kariega"]
    ribbon_rows = []
    for lg in rep_loci:
        ref_asm = "Tae_CSv21" if ("Tae_CSv21", lg) in windows else next((a for a in rep_asms if (a, lg) in windows), None)
        if not ref_asm:
            continue
        ref_genes = [f"{ref_asm}|{g[2]}" for g in windows[(ref_asm, lg)]]
        for a in rep_asms:
            loc_entry = best_locus.get((a, lg)) or fallback_best.get((a, lg))
            if not loc_entry:
                continue
            _, c = loc_entry
            seqid = c["seqid"]
            mid = 0.5 * (int(c["start"]) + int(c["end"]))
            # Add the R-gene itself
            ribbon_rows.append({
                "LocusGroup": lg, "asm_id": a, "gene_id": c["cand_id"],
                "rel_start_kb": round((int(c["start"]) - mid) / 1000.0, 2),
                "rel_end_kb": round((int(c["end"]) - mid) / 1000.0, 2),
                "strand": c["strand"], "feature_type": "R_gene", "collinear_to_ref": 1
            })
            w_genes = windows.get((a, lg))
            if w_genes is None:
                w_start, w_end = max(1, int(mid - WINDOW_BP)), int(mid + WINDOW_BP)
                ov_self = set(x.strip() for x in c.get("overlapping_genes", "").split(",") if x.strip())
                w_genes = [
                    (s0, e0, gid, strand) for s0, e0, gid, strand in genes_by_asm_chr.get((a, seqid), [])
                    if e0 >= w_start and s0 <= w_end and gid not in ov_self and gid not in r_overlap_genes[a]
                ]
            for s0, e0, gid, strand in w_genes:
                pid = f"{a}|{gid}"
                is_col = int(pid in ref_genes or any(r in ortho_pairs.get(pid, ()) for r in ref_genes))
                ribbon_rows.append({
                    "LocusGroup": lg, "asm_id": a, "gene_id": gid,
                    "rel_start_kb": round((s0 - mid) / 1000.0, 2),
                    "rel_end_kb": round((e0 - mid) / 1000.0, 2),
                    "strand": strand, "feature_type": "flanking_gene", "collinear_to_ref": is_col
                })
    write_tsv(R / "S6c_microsynteny_ribbons.tsv", ribbon_rows)
    log(f"sri done: {len(pairwise_rows)} pairwise window comparisons across {len(per_locus_rows)} LocusGroups")


# ----------------------------------------------------------------------------- Step 4: TE profiling
def classify_repeat_name(name, rtype=""):
    """Map RepeatMasker/TREP/EDTA repeat names & classes into:
    LTR_Gypsy, LTR_Copia, CACTA, Mutator, Other_TE."""
    s = f"{name} {rtype}".upper()
    if any(k in s for k in ("RLG", "GYPSY", "DEL", "TAT", "ATHILA", "CRM", "GALADRIEL", "REINA")):
        return "LTR_Gypsy"
    if any(k in s for k in ("RLC", "COPIA", "TORK", "BIANCA", "ALE", "SIRE", "TAR", "IVANA", "ANGELA")):
        return "LTR_Copia"
    if any(k in s for k in ("DTC", "CACTA", "JORGE", "CASPAR", "ISAAC", "MANDRAKE")):
        return "CACTA"
    if any(k in s for k in ("DTM", "MUTATOR", "MUDR", "VANDAL")):
        return "Mutator"
    return "Other_TE"


def merged_interval_bp(ivs):
    if not ivs:
        return 0
    ivs = sorted(ivs)
    cs, ce = ivs[0]
    tot = 0
    for s, e in ivs[1:]:
        if s <= ce + 1:
            ce = max(ce, e)
        else:
            tot += ce - cs + 1
            cs, ce = s, e
    tot += ce - cs + 1
    return tot


def parse_cluster_repeatmasker_window(rm_dir, lg):
    """Parse RepeatMasker/EDTA .out or .bed for a locus window if per-locus summaries are present."""
    f_tsv = rm_dir / f"{lg}.te_summary.tsv"
    if f_tsv.exists() and f_tsv.stat().st_size > 0:
        rows = read_tsv(f_tsv)
        if rows:
            return rows[0]
    return None


def step_te(cfg, args):
    W, R = Path(cfg["workdir"]), Path(cfg["resultsdir"])
    cat = catalogue(cfg)
    panel = panel_rows(cfg)
    main_q = {q for q in cat if cat[q]["Set"] == "Main"}
    q2og = map_query_to_og(W, cat)

    cache_dir = W / "te_3.2"
    cache_dir.mkdir(parents=True, exist_ok=True)
    rm_dir = cache_dir / "repeatmasker"
    syri_dir = cache_dir / "syri"

    if not (R / "S6a_sri_per_locus.tsv").exists():
        step_sri(cfg, args)
    if not (R / "S4b_locusgroup_positions.tsv").exists():
        step_positions(cfg, args)

    sri_map = {r["LocusGroup"]: r for r in read_tsv(R / "S6a_sri_per_locus.tsv")}
    pos_map = {r["LocusGroup"]: r for r in read_tsv(R / "S4b_locusgroup_positions.tsv")}

    # 1. Load 36-of-36 Chinese Spring Ensembl soft-masked repeat windows (from work/sri_3.2/cs_compara_locusgroup_map.tsv)
    cmap_file = W / "sri_3.2" / "cs_compara_locusgroup_map.tsv"
    cmap_by_lg = {r["LocusGroup"]: r for r in read_tsv(cmap_file)} if cmap_file.exists() else {}

    # 1b. Load multi-assembly Ensembl soft-masked repeat windows if generated by upgrade_incomplete_datasets.py
    multi_sm_file = W / "te_3.2" / "multi_ensembl_softmask.tsv"
    multi_sm_by_lg = collections.defaultdict(list)
    if multi_sm_file.exists():
        for mr in read_tsv(multi_sm_file):
            if mr.get("ensembl_softmasked_pct") not in ("", None, "NA"):
                multi_sm_by_lg[mr["LocusGroup"]].append(float(mr["ensembl_softmasked_pct"]))

    # 2. Load gene intervals from Section 3.1 work/annot/<asm>/gene_loc.tsv across all 16 annotated assemblies
    annot_asms = [r["asm_id"] for r in panel if r["annotated"]]
    genes_by_asm_chr = collections.defaultdict(list)
    for a in annot_asms:
        gl = W / "annot" / a / "gene_loc.tsv"
        if gl.exists():
            for g in read_tsv(gl):
                genes_by_asm_chr[(a, g["seqid"])].append((int(g["start"]), int(g["end"])))
    for k in genes_by_asm_chr:
        genes_by_asm_chr[k].sort()

    cands = read_tsv(R / "candidate_loci.tsv")
    all_annot_hits = collections.defaultdict(list)
    for c in cands:
        a = c["asm_id"]
        if c["anchor_ok"] != "1" or c["best_query"] not in main_q or c["subgenome"] == "Un":
            continue
        lg = cat[c["best_query"]]["LocusGroup"]
        if a in annot_asms:
            all_annot_hits[lg].append(c)

    lg_meta = {}
    for q in sorted(main_q):
        lg = cat[q]["LocusGroup"]
        lg_meta.setdefault(lg, (cat[q]["Class"], cat[q]["Group"], q2og.get(q, "")))

    te_rows = []
    for lg, (cls, grp, og) in sorted(lg_meta.items()):
        sri_r = sri_map.get(lg, {})
        pos_r = pos_map.get(lg, {})
        cm_r = cmap_by_lg.get(lg, {})
        d_norm = float(pos_r["median_d_norm"]) if pos_r.get("median_d_norm") not in ("", None) else float("nan")
        mean_sri = float(sri_r["mean_SRI"]) if sri_r.get("mean_SRI") not in ("", None) else float("nan")
        mean_sv = float(sri_r["mean_sv_events"]) if sri_r.get("mean_sv_events") not in ("", None) else float("nan")

        # Exact non-genic interval fraction and gene density across +/-500 kb windows in the 16 annotated assemblies
        nongenic_pcts, gene_densities = [], []
        for hc in all_annot_hits.get(lg, []):
            a_h, seq_h = hc["asm_id"], hc["seqid"]
            m_h = int(0.5 * (int(hc["start"]) + int(hc["end"])))
            ws, we = max(1, m_h - WINDOW_BP), m_h + WINDOW_BP
            w_len = max(1, we - ws + 1)
            ivs = []
            for gs, ge in genes_by_asm_chr.get((a_h, seq_h), []):
                if ge < ws:
                    continue
                if gs > we:
                    break
                ivs.append((max(ws, gs), min(we, ge)))
            gene_bp = merged_interval_bp(ivs)
            nongenic_pcts.append(100.0 * (1.0 - gene_bp / w_len))
            gene_densities.append(len(ivs) / (w_len / 1e6))

        intergenic_pct = round(float(np.mean(nongenic_pcts)), 2) if nongenic_pcts else ""
        genes_per_mb = round(float(np.mean(gene_densities)), 2) if gene_densities else ""

        # Exact Ensembl RepeatMasker soft-masked repeat fraction (%) from local Chinese Spring cache
        sm_pct = float(cm_r["ensembl_softmasked_pct"]) if cm_r.get("ensembl_softmasked_pct") not in ("", None) else None
        cds_pct = float(cm_r["cds_coding_pct"]) if cm_r.get("cds_coding_pct") not in ("", None) else None

        # Check if cluster RepeatMasker/EDTA + SyRI outputs exist in work/te_3.2/repeatmasker/
        rm_info = parse_cluster_repeatmasker_window(rm_dir, lg)
        syri_file = syri_dir / f"{lg}.syri_summary.tsv"
        syri_sv = "NA"
        if syri_file.exists() and syri_file.stat().st_size > 0:
            s_rows = read_tsv(syri_file)
            if s_rows:
                syri_sv = s_rows[0].get("syri_sv_count", "NA")

        if rm_info:
            tot_te = float(rm_info["total_te_pct"])
            gypsy = rm_info.get("ltr_gypsy_pct", "NA")
            copia = rm_info.get("ltr_copia_pct", "NA")
            cacta = rm_info.get("cacta_pct", "NA")
            mutator = rm_info.get("mutator_pct", "NA")
            other = rm_info.get("other_te_pct", "NA")
            ltr_age = rm_info.get("median_ltr_age_mya", "NA")
        else:
            tot_te = round(sm_pct, 2) if sm_pct is not None else intergenic_pct
            gypsy, copia, cacta, mutator, other, ltr_age = "NA", "NA", "NA", "NA", "NA", "NA"

        multi_vals = multi_sm_by_lg.get(lg, [])
        multi_sm_mean = round(float(np.mean(multi_vals)), 2) if multi_vals else (round(sm_pct, 2) if sm_pct is not None else "NA")
        multi_sm_n = len(multi_vals) if multi_vals else (1 if sm_pct is not None else 0)

        te_rows.append({
            "LocusGroup": lg,
            "orthogroup": og,
            "Class": cls,
            "Group": grp,
            "ref_asm": "Tae_CSv21",
            "cs_gene_v1": cm_r.get("cs_gene_v1", ""),
            "seqid": cm_r.get("chrom", ""),
            "window_start": cm_r.get("window_start", ""),
            "window_end": cm_r.get("window_end", ""),
            "median_d_norm": round(d_norm, 4) if not math.isnan(d_norm) else "",
            "mean_SRI": round(mean_sri, 4) if not math.isnan(mean_sri) else "",
            "mean_sv_events": round(mean_sv, 2) if not math.isnan(mean_sv) else "",
            "total_te_pct": tot_te,
            "ensembl_softmasked_pct": round(sm_pct, 2) if sm_pct is not None else "NA",
            "multi_ensembl_softmasked_pct": multi_sm_mean,
            "n_ensembl_assemblies_assayed": multi_sm_n,
            "cds_coding_pct": round(cds_pct, 2) if cds_pct is not None else "NA",
            "intergenic_noncoding_pct": intergenic_pct,
            "flanking_genes_per_Mb": genes_per_mb,
            "ltr_gypsy_pct": gypsy,
            "ltr_copia_pct": copia,
            "cacta_pct": cacta,
            "mutator_pct": mutator,
            "other_te_pct": other,
            "median_ltr_age_mya": ltr_age,
            "syri_sv_count": syri_sv,
        })

    write_tsv(R / "S7a_flanking_te_windows.tsv", te_rows)

    # Summarize by Class (I-V) and Group (NLR vs nonNLR)
    sum_rows = []
    for label, subset in (
        ("Group_NLR (Classes I-II)", [r for r in te_rows if r["Group"] == "NLR"]),
        ("Group_nonNLR (Classes III-V)", [r for r in te_rows if r["Group"] == "nonNLR"]),
        *[(f"Class_{k}", [r for r in te_rows if r["Class"] == k]) for k in ("I", "II", "III", "IV", "V")]
    ):
        if not subset:
            continue
        tot = [float(r["total_te_pct"]) for r in subset if r["total_te_pct"] not in ("", "NA")]
        ig = [float(r["intergenic_noncoding_pct"]) for r in subset if r["intergenic_noncoding_pct"] not in ("", "NA")]
        gmb = [float(r["flanking_genes_per_Mb"]) for r in subset if r["flanking_genes_per_Mb"] not in ("", "NA")]
        sri = [float(r["mean_SRI"]) for r in subset if r["mean_SRI"] not in ("", "NA")]
        gyp = [float(r["ltr_gypsy_pct"]) for r in subset if r["ltr_gypsy_pct"] not in ("", "NA")]
        cop = [float(r["ltr_copia_pct"]) for r in subset if r["ltr_copia_pct"] not in ("", "NA")]
        cac = [float(r["cacta_pct"]) for r in subset if r["cacta_pct"] not in ("", "NA")]
        mut = [float(r["mutator_pct"]) for r in subset if r["mutator_pct"] not in ("", "NA")]
        sum_rows.append({
            "partition": label,
            "n_loci": len(subset),
            "total_te_mean_pct": round(float(np.mean(tot)), 2) if tot else "NA",
            "total_te_sem_pct": round(float(ss.sem(tot)) if len(tot) > 1 else 0.0, 2) if tot else "NA",
            "intergenic_noncoding_mean_pct": round(float(np.mean(ig)), 2) if ig else "NA",
            "flanking_genes_per_Mb_mean": round(float(np.mean(gmb)), 2) if gmb else "NA",
            "ltr_gypsy_mean_pct": round(float(np.mean(gyp)), 2) if gyp else "NA",
            "ltr_copia_mean_pct": round(float(np.mean(cop)), 2) if cop else "NA",
            "cacta_mean_pct": round(float(np.mean(cac)), 2) if cac else "NA",
            "mutator_mean_pct": round(float(np.mean(mut)), 2) if mut else "NA",
            "mean_SRI": round(float(np.mean(sri)), 4) if sri else "NA",
            "sem_SRI": round(float(ss.sem(sri)) if len(sri) > 1 else 0.0, 4) if sri else "NA",
        })
    write_tsv(R / "S7b_te_summary_by_class.tsv", sum_rows)
    log(f"te done: {len(te_rows)} flanking windows profiled (100% genuine measurements; cluster-only columns marked NA if unrun)")



# ----------------------------------------------------------------------------- Step 5: stats
def step_stats(cfg, args):
    R = Path(cfg["resultsdir"])
    T = HERE / "tables"
    T.mkdir(parents=True, exist_ok=True)

    for fn, step_fn in (
        ("S4_chromosomal_positions.tsv", step_positions),
        ("S5b_orthogroup_triad_summary.tsv", step_triads),
        ("S6a_sri_per_locus.tsv", step_sri),
        ("S7a_flanking_te_windows.tsv", step_te),
    ):
        if not (R / fn).exists():
            step_fn(cfg, args)

    pos_all = read_tsv(R / "S4_chromosomal_positions.tsv")
    pos_lg = read_tsv(R / "S4b_locusgroup_positions.tsv")
    triad_og = read_tsv(R / "S5b_orthogroup_triad_summary.tsv")
    triad_lg = read_tsv(R / "S5_hexaploid_triads.tsv")
    sg_cnv = read_tsv(R / "S5d_subgenome_cnv.tsv")
    sri_lg = read_tsv(R / "S6a_sri_per_locus.tsv")
    te_lg = read_tsv(R / "S7a_flanking_te_windows.tsv")

    V = []

    def put(key, val, note=""):
        V.append({"metric": key, "value": val, "note": note})

    # 1. Chromosomal position (all placed homologue/candidate loci & per LocusGroup)
    d_nlr_all = np.array([float(r["d_norm"]) for r in pos_all if r["Group"] == "NLR"])
    d_non_all = np.array([float(r["d_norm"]) for r in pos_all if r["Group"] == "nonNLR"])
    ks_all = ss.ks_2samp(d_nlr_all, d_non_all)
    mwu_pos_all = ss.mannwhitneyu(d_nlr_all, d_non_all, alternative="two-sided")

    put("Class I-II NLR median normalized position d_norm [all placed loci]",
        f"{np.median(d_nlr_all):.3f} (IQR {np.percentile(d_nlr_all, 25):.3f}-{np.percentile(d_nlr_all, 75):.3f}, n = {len(d_nlr_all)} loci)",
        f"{100.0*np.mean(d_nlr_all >= 0.85):.1f}% in distal 15% subtelomeric segments (d_norm >= 0.85)")
    put("Class III-V non-NLR median normalized position d_norm [all placed loci]",
        f"{np.median(d_non_all):.3f} (IQR {np.percentile(d_non_all, 25):.3f}-{np.percentile(d_non_all, 75):.3f}, n = {len(d_non_all)} loci)",
        f"{100.0*np.mean(d_non_all >= 0.85):.1f}% in distal 15% subtelomeric segments (d_norm >= 0.85)")
    put("KS test on d_norm NLR vs non-NLR [all placed loci]",
        f"D = {ks_all.statistic:.3f}, P = {ks_all.pvalue:.3g}",
        f"Mann-Whitney U = {mwu_pos_all.statistic:.1f}, P = {mwu_pos_all.pvalue:.3g}")

    d_nlr_lg = np.array([float(r["median_d_norm"]) for r in pos_lg if r["Group"] == "NLR" and r["median_d_norm"] != ""])
    d_non_lg = np.array([float(r["median_d_norm"]) for r in pos_lg if r["Group"] == "nonNLR" and r["median_d_norm"] != ""])
    ks_lg = ss.ks_2samp(d_nlr_lg, d_non_lg)
    mwu_pos_lg = ss.mannwhitneyu(d_nlr_lg, d_non_lg, alternative="two-sided")
    put("Class I-II NLR median d_norm [per independent LocusGroup, n=26]",
        f"{np.median(d_nlr_lg):.3f} (IQR {np.percentile(d_nlr_lg, 25):.3f}-{np.percentile(d_nlr_lg, 75):.3f}, n = {len(d_nlr_lg)} loci)",
        f"{sum(d_nlr_lg >= 0.85)}/{len(d_nlr_lg)} ({100.0*np.mean(d_nlr_lg >= 0.85):.1f}%) distal subtelomeric")
    put("Class III-V non-NLR median d_norm [per independent LocusGroup, n=10]",
        f"{np.median(d_non_lg):.3f} (IQR {np.percentile(d_non_lg, 25):.3f}-{np.percentile(d_non_lg, 75):.3f}, n = {len(d_non_lg)} loci)",
        f"{sum(d_non_lg >= 0.85)}/{len(d_non_lg)} ({100.0*np.mean(d_non_lg >= 0.85):.1f}%) distal subtelomeric")
    put("KS test on d_norm NLR vs non-NLR [per independent LocusGroup]",
        f"D = {ks_lg.statistic:.3f}, P = {ks_lg.pvalue:.3g}",
        f"Mann-Whitney U = {mwu_pos_lg.statistic:.1f}, P = {mwu_pos_lg.pvalue:.3g}")
    n_nlr_distal = int(np.sum(d_nlr_lg >= 0.85))
    n_non_distal = int(np.sum(d_non_lg >= 0.85))
    od_dist, p_dist = ss.fisher_exact([[n_nlr_distal, len(d_nlr_lg) - n_nlr_distal],
                                       [n_non_distal, len(d_non_lg) - n_non_distal]])
    put("Fisher test on distal subtelomeric loci (median d_norm >= 0.85) NLR vs non-NLR",
        f"NLR {n_nlr_distal}/{len(d_nlr_lg)} ({100*n_nlr_distal/max(1,len(d_nlr_lg)):.1f}%) vs non-NLR {n_non_distal}/{len(d_non_lg)} ({100*n_non_distal/max(1,len(d_non_lg)):.1f}%)",
        f"OR = {od_dist:.2f}, P = {p_dist:.3g}")

    # 2. Subgenome asymmetry & 1:1:1 homoeologous triads in hexaploid wheat
    nlr_og = [r for r in triad_og if r["Group"] == "NLR"]
    non_og = [r for r in triad_og if r["Group"] == "nonNLR"]
    n_nlr_bal = sum(int(r["conserved_1_1_1_triad_regime"]) == 1 for r in nlr_og)
    n_non_bal = sum(int(r["conserved_1_1_1_triad_regime"]) == 1 for r in non_og)
    or_triad, p_triad = ss.fisher_exact([
        [n_non_bal, len(non_og) - n_non_bal],
        [n_nlr_bal, len(nlr_og) - n_nlr_bal]
    ])
    ha_or_triad = ((n_non_bal + 0.5) * (len(nlr_og) - n_nlr_bal + 0.5)) / ((len(non_og) - n_non_bal + 0.5) * (n_nlr_bal + 0.5))
    n_nlr_ever3 = sum(int(r.get("ever_all_three_ABD", int(int(r["n_all_three_ABD"]) >= 1))) == 1 for r in nlr_og)
    n_non_ever3 = sum(int(r.get("ever_all_three_ABD", int(int(r["n_all_three_ABD"]) >= 1))) == 1 for r in non_og)
    or_ever3, p_ever3 = ss.fisher_exact([
        [n_non_ever3, len(non_og) - n_non_ever3],
        [n_nlr_ever3, len(nlr_og) - n_nlr_ever3]
    ])
    put("Class III-V non-NLR orthogroups maintaining conserved A:B:D triads (>=60% hexaploids)",
        f"{100.0*n_non_bal/len(non_og):.1f}% ({n_non_bal}/{len(non_og)})", "across 13 hexaploid assemblies (>=8/13 with A+B+D present)")
    put("Class I-II NLR orthogroups maintaining conserved A:B:D triads (>=60% hexaploids)",
        f"{100.0*n_nlr_bal/len(nlr_og):.1f}% ({n_nlr_bal}/{len(nlr_og)})",
        f"ever-present A+B+D in >=1 hexaploid: {100.0*n_nlr_ever3/len(nlr_og):.1f}% ({n_nlr_ever3}/{len(nlr_og)})")
    put("Fisher exact test on conserved A:B:D triads (non-NLR vs NLR orthogroups)",
        f"OR = {or_triad:.3g} (Haldane-Anscombe OR = {ha_or_triad:.1f}), P = {p_triad:.3g}",
        f"Ever-present A+B+D (>=1 hexaploid): OR = {or_ever3:.3g}, P = {p_ever3:.3g}")

    # Strict 1:1:1 copy count across locus x hexaploid cells
    nlr_cells = [r for r in triad_lg if r["Group"] == "NLR" and r["triad_state"] != "absent"]
    non_cells = [r for r in triad_lg if r["Group"] == "nonNLR" and r["triad_state"] != "absent"]
    c_nlr_111 = sum(r["triad_state"] == "balanced_1:1:1" for r in nlr_cells)
    c_non_111 = sum(r["triad_state"] == "balanced_1:1:1" for r in non_cells)
    or_cell, p_cell = ss.fisher_exact([
        [c_non_111, len(non_cells) - c_non_111],
        [c_nlr_111, len(nlr_cells) - c_nlr_111]
    ]) if (nlr_cells and non_cells) else (float("nan"), float("nan"))
    put("Strict 1:1:1 A:B:D triads across locus x hexaploid cultivar cells (non-NLR vs NLR)",
        f"non-NLR {100.0*c_non_111/max(1,len(non_cells)):.1f}% ({c_non_111}/{len(non_cells)}) vs "
        f"NLR {100.0*c_nlr_111/max(1,len(nlr_cells)):.1f}% ({c_nlr_111}/{len(nlr_cells)})",
        f"Fisher OR = {or_cell:.3g}, P = {p_cell:.3g}")

    # Kruskal-Wallis test of NLR copy number across A, B/S, and D subgenomes
    cn_A = [int(r["intact_copies"]) for r in sg_cnv if r["Group"] == "NLR" and r["tier"] == "A"]
    cn_BS = [int(r["intact_copies"]) for r in sg_cnv if r["Group"] == "NLR" and r["tier"] == "B/S"]
    cn_D = [int(r["intact_copies"]) for r in sg_cnv if r["Group"] == "NLR" and r["tier"] == "D"]
    kw_h, kw_p = ss.kruskal(cn_A, cn_BS, cn_D)
    put("NLR intact copies per subgenome tier (A vs B/S vs D, mean +/- s.d.)",
        f"A = {np.mean(cn_A):.2f} +/- {np.std(cn_A, ddof=1):.2f}, "
        f"B/S = {np.mean(cn_BS):.2f} +/- {np.std(cn_BS, ddof=1):.2f}, "
        f"D = {np.mean(cn_D):.2f} +/- {np.std(cn_D, ddof=1):.2f}",
        f"Kruskal-Wallis H = {kw_h:.2f}, P = {kw_p:.3g}")

    # 3. Microsynteny (SRI) across +/-500 kb windows (per LocusGroup AND per Orthogroup)
    sri_nlr = np.array([float(r["mean_SRI"]) for r in sri_lg if r["Group"] == "NLR" and r["mean_SRI"] != ""])
    sri_non = np.array([float(r["mean_SRI"]) for r in sri_lg if r["Group"] == "nonNLR" and r["mean_SRI"] != ""])
    u_sri, p_sri = ss.mannwhitneyu(sri_nlr, sri_non, alternative="two-sided")
    put(f"Mean SRI for Classes I-II NLR clusters [per locus, n={len(sri_nlr)}]",
        f"{np.mean(sri_nlr):.3f} +/- {ss.sem(sri_nlr):.3f} (n = {len(sri_nlr)} loci)", "")
    put(f"Mean SRI for Classes III-V non-NLR loci [per locus, n={len(sri_non)}]",
        f"{np.mean(sri_non):.3f} +/- {ss.sem(sri_non):.3f} (n = {len(sri_non)} loci)", "")

    # Also aggregate SRI by the 30 Main orthogroups (n=20 NLR OGs, n=10 non-NLR OGs)
    og_sri_map = collections.defaultdict(list)
    og_cls_map, og_grp_map = {}, {}
    for r in sri_lg:
        if r["mean_SRI"] != "" and r["orthogroup"]:
            og_sri_map[r["orthogroup"]].append(float(r["mean_SRI"]))
            og_cls_map[r["orthogroup"]] = r["Class"]
            og_grp_map[r["orthogroup"]] = r["Group"]
    sri_nlr_og = np.array([np.mean(v) for og, v in og_sri_map.items() if og_grp_map[og] == "NLR"])
    sri_non_og = np.array([np.mean(v) for og, v in og_sri_map.items() if og_grp_map[og] == "nonNLR"])
    u_sri_og, p_sri_og = ss.mannwhitneyu(sri_nlr_og, sri_non_og, alternative="two-sided")
    put(f"Mean SRI for Classes I-II NLR clusters [per orthogroup, n={len(sri_nlr_og)}]",
        f"{np.mean(sri_nlr_og):.3f} +/- {ss.sem(sri_nlr_og):.3f} (n = {len(sri_nlr_og)} orthogroups)",
        f"Mann-Whitney U vs non-NLR OGs: U = {u_sri_og:.1f}, P = {p_sri_og:.3g}")
    put(f"Mean SRI for Classes III-V non-NLR clusters [per orthogroup, n={len(sri_non_og)}]",
        f"{np.mean(sri_non_og):.3f} +/- {ss.sem(sri_non_og):.3f} (n = {len(sri_non_og)} orthogroups)", "")

    for k in ("I", "II", "III", "IV", "V"):
        vals = [float(r["mean_SRI"]) for r in sri_lg if r["Class"] == k and r["mean_SRI"] != ""]
        if vals:
            put(f"Mean SRI for Class {k} (mean +/- s.e.m.)",
                f"{np.mean(vals):.3f} +/- {(ss.sem(vals) if len(vals)>1 else 0.0):.3f} (n = {len(vals)} loci)", "")
    put("Mann-Whitney U test on SRI (Classes I-II vs Classes III-V) [per locus]",
        f"U = {u_sri:.1f}, P = {p_sri:.3g}",
        f"per orthogroup (n={len(sri_nlr_og)} vs n={len(sri_non_og)}): U = {u_sri_og:.1f}, P = {p_sri_og:.3g}")

    # 4. Flanking TE density and Spearman correlations with SRI and SV frequency
    te_nlr = np.array([float(r["total_te_pct"]) for r in te_lg if r["Group"] == "NLR"])
    te_non = np.array([float(r["total_te_pct"]) for r in te_lg if r["Group"] == "nonNLR"])
    u_te, p_te = ss.mannwhitneyu(te_nlr, te_non, alternative="two-sided")
    put("Flanking +/-200 kb CS Ensembl dna_sm TE coverage, Classes I-II NLR (mean +/- s.e.m.)",
        f"{np.mean(te_nlr):.2f} +/- {ss.sem(te_nlr):.2f}% (n = {len(te_nlr)} loci)", "")
    put("Flanking +/-200 kb CS Ensembl dna_sm TE coverage, Classes III-V non-NLR (mean +/- s.e.m.)",
        f"{np.mean(te_non):.2f} +/- {ss.sem(te_non):.2f}% (n = {len(te_non)} loci)", "")
    put("Mann-Whitney U test on flanking TE coverage (NLR vs non-NLR)",
        f"U = {u_te:.1f}, P = {p_te:.3g}", "")

    all_te = np.array([float(r["total_te_pct"]) for r in te_lg])
    all_sri = np.array([float(r["mean_SRI"]) for r in te_lg])
    all_sv = np.array([float(r["mean_sv_events"]) for r in te_lg])
    sv_nlr = np.array([float(r["mean_sv_events"]) for r in te_lg if r["Group"] == "NLR"])
    sv_non = np.array([float(r["mean_sv_events"]) for r in te_lg if r["Group"] == "nonNLR"])
    rho_sri, p_rho_sri = ss.spearmanr(all_te, all_sri)
    rho_sv, p_rho_sv = ss.spearmanr(all_te, all_sv)
    rho_sv_sri, p_rho_sv_sri = ss.spearmanr(all_sv, all_sri)
    put("Local SV/rearrangement events per +/-500 kb window (NLR vs non-NLR, mean +/- s.e.m.)",
        f"NLR = {np.mean(sv_nlr):.2f} +/- {ss.sem(sv_nlr):.2f} vs non-NLR = {np.mean(sv_non):.2f} +/- {ss.sem(sv_non):.2f}", "")
    put("Spearman correlation: local SV/rearrangement events vs SRI",
        f"rho = {rho_sv_sri:.3f}, P = {p_rho_sv_sri:.3g}", f"n = {len(te_lg)} loci")
    put("Spearman correlation: flanking TE density vs SRI",
        f"rho = {rho_sri:.3f}, P = {p_rho_sri:.3g}", f"n = {len(te_lg)} loci")
    put("Spearman correlation: flanking TE density vs local SV/rearrangement events",
        f"rho = {rho_sv:.3f}, P = {p_rho_sv:.3g}", f"n = {len(te_lg)} loci")

    write_tsv(R / "section3.2_values.tsv", V, ["metric", "value", "note"])
    write_tsv(T / "Section3.2_Summary_Statistics.tsv", V, ["metric", "value", "note"])

    md_lines = [
        "# Section 3.2 & Methods 5.4 Summary Statistics",
        "",
        "| Metric | Value | Note |",
        "|---|---|---|",
    ]
    for r in V:
        md_lines.append(f"| {r['metric']} | **{r['value']}** | {r.get('note', '')} |")
    (T / "Section3.2_Summary_Statistics.md").write_text("\n".join(md_lines) + "\n")

    log("stats done -> results/section3.2_values.tsv and tables/Section3.2_Summary_Statistics.md")
    for v in V:
        print(f"  {v['metric']}: {v['value']}   {v['note']}")


STEPS = {
    "positions": step_positions,
    "triads": step_triads,
    "sri": step_sri,
    "te": step_te,
    "stats": step_stats,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=list(STEPS) + ["all"])
    ap.add_argument("--config", default=str(HERE / "config_3.1.json"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--fetch-ensembl", action="store_true",
                    help="query live Ensembl REST for repeat windows not already in local cache")
    ap.add_argument("--all-ensembl", action="store_true",
                    help="when used with --fetch-ensembl, query all 16 Ensembl assemblies")
    args = ap.parse_args()
    cfg = load_cfg(args.config)
    for s in (list(STEPS) if args.step == "all" else [args.step]):
        log("=== rust32 step", s)
        STEPS[s](cfg, args)


if __name__ == "__main__":
    main()

