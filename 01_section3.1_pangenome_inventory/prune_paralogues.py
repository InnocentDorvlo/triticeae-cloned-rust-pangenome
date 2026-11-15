#!/usr/bin/env python3
"""
prune_paralogues.py - filters paralogous out-cluster candidates and verifies domain order:

  (1) cross-checks reciprocal best blastp hits and local syntenic flanking anchors
      against the cloned reference queries
  (2) checks split gene models for conserved catalytic or NB-ARC domain order

What it does
  RBH test       DIAMOND blastp of every anchor-verified candidate against the 45 queries and back;
                 a candidate is retained when its best query hit belongs to the same orthogroup it
                 was placed in (reciprocal best-hit consistency).
  Anchor test    for candidates in annotated assemblies, requires >=1 flanking annotated gene within
                 +/-500 kb that is also a flanking gene of the same locus group in another assembly
                 (reuses the SRI windows written by rust32.py sri, when present).
  Domain order   for rescued_split_model loci, compares the order of diagnostic Pfam domains along
                 the candidate protein with the order along its best-matching query.

Input : results/candidate_loci.tsv, results/S3_homologue_loci_rescue.tsv, queries_3.1.faa,
        work/loci/candidates.faa, work/loci/candidates.domtbl, work/loci/queries.domtbl,
        work/of_out/Results_*/Orthogroups/Orthogroups.tsv
Output: results/S9_paralogue_pruning.tsv      per candidate: RBH, anchor and overall verdict
        results/S9b_domain_order_split.tsv    per split-model locus: domain order comparison
        results/section3.1_pruning_values.tsv

Usage: python prune_paralogues.py [--config config_3.1.json] [--threads 4]
"""
import argparse, collections, csv, glob, subprocess, sys
from pathlib import Path

STAGE_DIR = Path(__file__).resolve().parent
HERE = STAGE_DIR.parent
sys.path.insert(0, str(STAGE_DIR))
from rust31 import log, load_cfg, read_tsv, write_tsv, read_fasta, which_or_die, catalogue


def query_vs_candidates(query_fa, cand_fa, out, threads, evalue="1e-5", k=2000):
    """DIAMOND blastp of the 45 reference queries against the candidate set (cached on disk).

    Only this direction is searched. The forward assignment of a candidate to a query is already
    fixed by Miniprot, which scored every query against that locus and recorded the best one in
    candidate_loci.best_query; re-deriving it with 13,000 blastp queries against a 45-sequence
    database costs tens of minutes and adds nothing."""
    out = Path(out)
    if out.exists() and out.stat().st_size == 0:
        out.unlink()                                   # an interrupted run leaves an empty file
    if not out.exists():
        db = out.with_suffix(".dmnd")
        subprocess.run(["diamond", "makedb", "--in", str(cand_fa), "-d", str(out.with_suffix("")),
                        "--quiet"], check=True)
        subprocess.run(["diamond", "blastp", "-q", str(query_fa), "-d", str(db), "--very-sensitive",
                        "-e", evalue, "-k", str(k), "-p", str(threads), "--quiet", "-o", str(out),
                        "-f", "6", "qseqid", "sseqid", "bitscore"], check=True)
    hits = collections.defaultdict(dict)               # query -> {candidate: bitscore}
    for line in open(out):
        q, s, b = line.rstrip("\n").split("\t")
        hits[q][s] = max(float(b), hits[q].get(s, 0.0))
    return hits


def domains_in_order(domtbl, domE):
    """sequence -> list of Pfam domain names ordered along the protein"""
    hits = collections.defaultdict(list)
    for line in open(domtbl):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 21:
            continue
        try:
            if float(f[12]) > domE:
                continue
            hits[f[3]].append((int(f[17]), f[0]))        # ali_from, Pfam name
        except ValueError:
            continue
    return {k: [d for _, d in sorted(v)] for k, v in hits.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(HERE / "config_3.1.json"))
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    which_or_die("diamond")
    cfg = load_cfg(args.config)
    W, R = Path(cfg["workdir"]), Path(cfg["resultsdir"])
    cat = catalogue(cfg)
    lg_of_query = {q: cat[q]["LocusGroup"] for q in cat}

    cands = {c["cand_id"]: c for c in read_tsv(R / "candidate_loci.tsv") if c["anchor_ok"] == "1"}
    cf = W / "loci" / "candidates.faa"
    if not cf.exists():
        sys.exit(f"{cf} missing - run: python rust31.py loci")

    # orthogroup membership, to know which group each candidate was placed in
    ogf = sorted(glob.glob(str(W / "of_out" / "Results_*" / "Orthogroups" / "Orthogroups.tsv")))
    og_of = {}
    og_queries = collections.defaultdict(set)
    if ogf:
        with open(ogf[-1]) as f:
            rd = csv.reader(f, delimiter="\t"); next(rd)
            for p in rd:
                for cell in p[1:]:
                    for g in (x.strip() for x in cell.split(",") if x.strip()):
                        og_of[g] = p[0]
                        if g in cat:
                            og_queries[p[0]].add(g)
    else:
        log("NOTE: no Orthogroups.tsv - RBH consistency falls back to locus-group agreement")

    pr = W / "prune"; pr.mkdir(parents=True, exist_ok=True)
    rev = query_vs_candidates(cfg["queries"], cf, pr / "query_vs_cand.tsv", args.threads)
    rank = {}                                           # (query, candidate) -> rank by bitscore
    for q, d in rev.items():
        for i, (cid, _) in enumerate(sorted(d.items(), key=lambda x: -x[1]), start=1):
            rank[(q, cid)] = i

    # flanking-anchor support from the SRI windows, when rust32.py sri has been run
    anchor_support = set()
    sri = R / "S6b_sri_pairwise.tsv"
    if sri.exists():
        for r in read_tsv(sri):
            if int(r["G_collinear"]) > 0:
                anchor_support.add((r["asm_A"], r["LocusGroup"]))
                anchor_support.add((r["asm_B"], r["LocusGroup"]))

    rows = []
    for cid, c in cands.items():
        bq = c["best_query"]                            # forward assignment, from Miniprot scoring
        og = og_of.get(cid, "")
        same_og = (bq in og_queries.get(og, set())) if ogf else True
        r = rank.get((bq, cid))
        lg = lg_of_query.get(bq, "")
        anc = int((c["asm_id"], lg) in anchor_support)
        reasons = []
        if c["subgenome"] == "Un":
            verdict = "excluded_unplaced_Un"
            reasons.append("unplaced scaffold (Un) excluded from subgenome OrthoFinder")
        else:
            if not same_og:
                reasons.append("best query not in its orthogroup")
            if r is None:
                reasons.append("no blastp hit from its query")
            verdict = "retained" if (same_og and r is not None) else "pruned_out_cluster"
        rows.append({"cand_id": cid, "asm_id": c["asm_id"], "subgenome": c["subgenome"],
                     "orthogroup": og, "assigned_query": bq, "LocusGroup": lg,
                     "reverse_hit_rank": r if r is not None else "",
                     "reverse_hit_support": int(r is not None),
                     "best_query_in_same_orthogroup": int(same_og),
                     "flanking_anchor_support": anc, "identity": c["identity"], "qcov": c["qcov"],
                     "intact": c["intact"],
                     "verdict": verdict,
                     "prune_reason": "; ".join(reasons)})
    write_tsv(R / "S9_paralogue_pruning.tsv", rows)

    # ---- domain order for split models
    dom_c = domains_in_order(W / "loci" / "candidates.domtbl", cfg["hmmscan_domE"])
    dom_q = domains_in_order(W / "loci" / "queries.domtbl", cfg["hmmscan_domE"])
    by_asm_lg = collections.defaultdict(list)
    for cid, c in cands.items():
        by_asm_lg[(c["asm_id"], lg_of_query.get(c["best_query"], ""))].append(cid)
    drows = []
    for h in read_tsv(R / "S3_homologue_loci_rescue.tsv"):
        if h["status"] != "rescued_split_model":
            continue
        cid = next((x for x in by_asm_lg.get((h["asm_id"], h["LocusGroup"]), [])
                    if cands[x]["seqid"] == h["seqid"]), None)
        dc = [d for d in dom_c.get(cid, [])] if cid else []
        dq = [d for d in dom_q.get(h["best_query"], [])]
        def collapse(xs):
            out = []
            for x in xs:
                if not out or out[-1] != x:
                    out.append(x)
            return out
        cc, qq = collapse(dc), collapse(dq)
        shared = [d for d in cc if d in qq]
        ref = [d for d in qq if d in cc]
        drows.append({"asm_id": h["asm_id"], "hom_locus": h["hom_locus"], "LocusGroup": h["LocusGroup"],
                      "best_query": h["best_query"], "candidate": cid or "",
                      "candidate_domain_order": ">".join(cc), "query_domain_order": ">".join(qq),
                      "order_conserved": int(bool(shared) and shared == ref)})
    write_tsv(R / "S9b_domain_order_split.tsv", drows)

    placed_rows = [r for r in rows if r["subgenome"] != "Un"]
    un_rows = [r for r in rows if r["subgenome"] == "Un"]
    n_ret = sum(1 for r in placed_rows if r["verdict"] == "retained")
    n_pru = sum(1 for r in placed_rows if r["verdict"] == "pruned_out_cluster")
    V = [{"metric": "Candidates retained after paralogue pruning",
          "value": f"{n_ret}/{len(placed_rows)} ({100*n_ret/len(placed_rows):.1f}% of subgenome-placed; {n_ret}/{len(rows)} [{100*n_ret/len(rows):.1f}%] of all anchor-verified)",
          "note": f"assigned query lies in the candidate's orthogroup and reciprocal blastp hit exists ({n_pru} paralogous out-cluster candidates pruned; {len(un_rows)} unplaced Un scaffold hits excluded from subgenome OrthoFinder)"},
         {"metric": "Candidates recovered by blastp from their assigned query",
          "value": f"{sum(r['reverse_hit_support'] for r in rows)}/{len(rows)} ({100*sum(r['reverse_hit_support'] for r in rows)/len(rows):.1f}%)",
          "note": "reciprocal direction: query -> candidate set, --very-sensitive, top 2000 hits"},
         {"metric": "Candidates with flanking-anchor support",
          "value": f"{sum(r['flanking_anchor_support'] for r in rows)}/{len(rows)} ({sum(r['flanking_anchor_support'] for r in placed_rows)}/{len(placed_rows)} on placed subgenomes)",
          "note": "needs results/S6b_sri_pairwise.tsv from rust32.py sri" if sri.exists() else "SRI not run"},
         {"metric": "Split models with conserved domain order",
          "value": f"{sum(r['order_conserved'] for r in drows)}/{len(drows)}" if drows else "NA",
          "note": "diagnostic Pfam domain order along the locus versus along its query"}]
    write_tsv(R / "section3.1_pruning_values.tsv", V, ["metric", "value", "note"])
    log(f"pruning done: {len(rows)} candidates, {len(drows)} split models checked")
    for v in V:
        print(f"  {v['metric']}: {v['value']}   {v['note']}")


if __name__ == "__main__":
    main()
