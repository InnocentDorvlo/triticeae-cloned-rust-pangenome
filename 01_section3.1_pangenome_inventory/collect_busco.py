#!/usr/bin/env python3
"""
collect_busco.py - summarizes assembly gene-space completeness verified
against the poales_odb10 lineage dataset.

Supports three output layouts under work/tree/<asm>/:
  1. compleasm/summary.txt or compleasm/**/short_summary*.txt  (cluster genome mode)
  2. busco_prot/**/short_summary*.txt                          (full BUSCO protein mode)
  3. subset_*.domtbl + poales_odb10/scores_cutoff              (Mac subset_500 protein mode)

Input : work/tree/<asm>/..., results/S2_assembly_panel_qc.tsv
Output: results/S2_assembly_panel_qc_busco.tsv
        results/section3.1_busco_values.tsv

Usage: python collect_busco.py [--config config_3.1.json]
"""
import argparse, collections, math, re, sys
from pathlib import Path

STAGE_DIR = Path(__file__).resolve().parent
HERE = STAGE_DIR.parent
sys.path.insert(0, str(STAGE_DIR))
from rust31 import log, load_cfg, read_tsv, write_tsv

KEYS = ("complete_pct", "single_copy_pct", "duplicated_pct", "fragmented_pct",
        "missing_pct", "subgenome_single_copy_pct", "n_markers", "busco_mode")


def parse_summary(p, mode_label="compleasm"):
    t = Path(p).read_text()
    out = {}
    m = re.search(r"S:\s*([\d.]+)%.*?D:\s*([\d.]+)%.*?F:\s*([\d.]+)%.*?M:\s*([\d.]+)%", t, re.S)
    if m:
        s, d, f, ms = float(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4))
        out.update(complete_pct=round(s + d, 2), single_copy_pct=s, duplicated_pct=d,
                   fragmented_pct=f, missing_pct=ms, subgenome_single_copy_pct="", busco_mode=mode_label)
    else:
        for key, pat in (("single_copy_pct", r"\(S\)[^\d]*([\d.]+)"),
                         ("duplicated_pct", r"\(D\)[^\d]*([\d.]+)"),
                         ("fragmented_pct", r"\(F\)[^\d]*([\d.]+)"),
                         ("missing_pct", r"\(M\)[^\d]*([\d.]+)")):
            mm = re.search(pat, t)
            if mm:
                out[key] = float(mm.group(1))
        if "single_copy_pct" in out and "duplicated_pct" in out:
            out["complete_pct"] = round(out["single_copy_pct"] + out["duplicated_pct"], 2)
        out["subgenome_single_copy_pct"] = ""
        out["busco_mode"] = mode_label
    n = re.search(r"N:\s*(\d+)", t) or re.search(r"Total BUSCO groups searched:\s*(\d+)", t)
    if n:
        out["n_markers"] = int(n.group(1))
    return out


def load_cutoffs(W):
    for c in (HERE / "busco_downloads" / "lineages" / "poales_odb10" / "scores_cutoff",
              W.parent / "busco_downloads" / "lineages" / "poales_odb10" / "scores_cutoff"):
        if c.exists():
            cut = {}
            for line in open(c):
                p = line.split()
                if len(p) >= 2:
                    cut[p[0]] = float(p[1])
            return cut
    return {}


def load_subset_models(W):
    for sf in sorted((W / "tree").glob("subset_*.hmm")):
        models = set()
        for line in open(sf):
            if line.startswith("NAME "):
                models.add(line.split()[1])
        if models:
            return models
    return set()


def parse_domtbl(domtbl, cutoffs, subset_models, W, asm_id, genomes):
    n_total = len(subset_models) if subset_models else len(cutoffs)
    if not n_total:
        return {}
    sub_of = {}
    seqmap = W / "miniprot" / asm_id / "seqmap.tsv"
    if seqmap.exists():
        sub_of = {x["seqid"]: x["subgenome"] for x in read_tsv(seqmap)}
    gene_sub = {}
    gl = W / "annot" / asm_id / "gene_loc.tsv"
    if gl.exists():
        gene_sub = {g["gene"]: sub_of.get(g["seqid"], "Un") for g in read_tsv(gl)}

    hits_genes = collections.defaultdict(set)
    any_hit_models = set()
    hits_by_sg = collections.defaultdict(lambda: collections.defaultdict(set))
    for line in open(domtbl):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 14:
            continue
        gene, model, evalue, score = f[0], f[3], float(f[6]), float(f[7])
        if subset_models and model not in subset_models:
            continue
        if evalue <= 1e-5:
            any_hit_models.add(model)
        if score >= cutoffs.get(model, 0.0):
            hits_genes[model].add(gene)
            sg = gene_sub.get(gene, "Un")
            if sg not in (None, "Un"):
                hits_by_sg[sg][model].add(gene)

    n_single = sum(1 for m, gs in hits_genes.items() if len(gs) == 1)
    n_dup = sum(1 for m, gs in hits_genes.items() if len(gs) >= 2)
    n_comp = n_single + n_dup
    n_frag = len(any_hit_models - set(hits_genes.keys()))
    n_miss = max(0, n_total - n_comp - n_frag)

    sg_single_pcts = []
    for sg in genomes:
        if sg in hits_by_sg:
            sg_sc = sum(1 for m, gs in hits_by_sg[sg].items() if len(gs) == 1)
            sg_single_pcts.append(100.0 * sg_sc / n_total)
    sg_mean = round(sum(sg_single_pcts) / len(sg_single_pcts), 2) if sg_single_pcts else ""

    return {
        "complete_pct": round(100.0 * n_comp / n_total, 2),
        "single_copy_pct": round(100.0 * n_single / n_total, 2),
        "duplicated_pct": round(100.0 * n_dup / n_total, 2),
        "fragmented_pct": round(100.0 * n_frag / n_total, 2),
        "missing_pct": round(100.0 * n_miss / n_total, 2),
        "subgenome_single_copy_pct": sg_mean,
        "n_markers": n_total,
        "busco_mode": f"proteome_subset_{n_total}"
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(HERE / "config_3.1.json"))
    args = ap.parse_args()
    cfg = load_cfg(args.config)
    W, R = Path(cfg["workdir"]), Path(cfg["resultsdir"])
    qc = read_tsv(R / "S2_assembly_panel_qc.tsv")
    cutoffs = load_cutoffs(W)
    subset_models = load_subset_models(W)

    found = 0
    for r in qc:
        a = r["asm_id"]
        cands_comp = list((W / "tree" / a).glob("compleasm/summary.txt")) \
            + list((W / "tree" / a).glob("compleasm/**/short_summary*.txt"))
        cands_busco = list((W / "tree" / a).glob("busco_prot/**/short_summary*.txt"))
        cands_dom = sorted((W / "tree" / a).glob("subset_*.domtbl"))

        if cands_comp:
            r.update(parse_summary(cands_comp[0], "compleasm_genome"))
            r["busco_lineage"] = "poales_odb10"
            found += 1
        elif cands_dom and cutoffs:
            # Prefer consistent subset_500 across all 16 annotated proteomes
            r.update(parse_domtbl(cands_dom[-1], cutoffs, subset_models, W, a, r.get("genomes", "")))
            r["busco_lineage"] = "poales_odb10"
            found += 1
        elif cands_busco:
            r.update(parse_summary(cands_busco[0], "busco_proteins"))
            r["busco_lineage"] = "poales_odb10"
            found += 1
        else:
            for k in KEYS:
                r[k] = ""
            r["busco_lineage"] = ""

    if not found:
        sys.exit("no compleasm/BUSCO/subset_domtbl outputs found under work/tree/*/")

    cols = list(qc[0].keys())
    write_tsv(R / "S2_assembly_panel_qc_busco.tsv", qc, cols)

    comp_vals = [float(r["complete_pct"]) for r in qc if r.get("complete_pct") != ""]
    n50_vals = [int(r["contig_n50_bp"]) / 1e6 for r in qc if r.get("contig_n50_bp")]
    V = [
        {"metric": "N assemblies with BUSCO completeness", "value": f"{found}/{len(qc)}",
         "note": "poales_odb10"},
        {"metric": "BUSCO complete % across annotated proteomes (mean +/- s.d., range)",
         "value": f"{sum(comp_vals)/len(comp_vals):.2f} +/- "
                  f"{(sum((x - sum(comp_vals)/len(comp_vals))**2 for x in comp_vals)/max(1, len(comp_vals)-1))**0.5:.2f}% "
                  f"(range {min(comp_vals):.1f}-{max(comp_vals):.1f}%)" if comp_vals else "NA",
         "note": f"n = {len(comp_vals)} assemblies"},
        {"metric": "Contig N50 across 27 assemblies (median, range Mb)",
         "value": f"median {sorted(n50_vals)[len(n50_vals)//2]:.2f} Mb "
                  f"(range {min(n50_vals):.2f}-{max(n50_vals):.2f} Mb)" if n50_vals else "NA",
         "note": f"n = {len(n50_vals)} assemblies"}
    ]
    write_tsv(R / "section3.1_busco_values.tsv", V, ["metric", "value", "note"])

    log(f"completeness added for {found}/{len(qc)} assemblies -> results/S2_assembly_panel_qc_busco.tsv")
    for r in qc:
        if r.get("complete_pct") != "":
            print(f"  {r['asm_id']:16s}: C={r['complete_pct']:5.1f}% (S={r['single_copy_pct']:5.1f}%, "
                  f"D={r['duplicated_pct']:5.1f}%, M={r['missing_pct']:4.1f}%) "
                  f"[subgenome-S={r.get('subgenome_single_copy_pct', '')}%]")


if __name__ == "__main__":
    main()
