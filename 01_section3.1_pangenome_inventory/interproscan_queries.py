#!/usr/bin/env python3
"""
interproscan_queries.py - verifies domain composition across the 45 cloned reference
query proteins using InterProScan (complementing local 73-model Pfam hmmscan).

Submits the 45 curated query sequences to the EBI InterProScan 5 REST service across
full member databases without requiring local ~50 GB database installation.

Input : queries_3.1.faa, catalogue_3.1.tsv
Output: results/S1c_interproscan_queries.tsv   one row per query x signature
        results/S1d_interproscan_summary.tsv   one row per query: InterPro accessions, Pfam, families
        work/iprscan/<ID>.tsv                  raw service output, so re-runs are free

Usage: python interproscan_queries.py --email you@institute.edu [--appl Pfam,SMART,Gene3D,PANTHER]
       EBI requires a contact e-mail; jobs are polled politely, three at a time.
"""
import argparse, csv, json, sys, time, urllib.parse, urllib.request
from pathlib import Path

STAGE_DIR = Path(__file__).resolve().parent
HERE = STAGE_DIR.parent
BASE = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
W, R = HERE / "work" / "iprscan", HERE / "results"


def post(url, data):
    req = urllib.request.Request(url, data=urllib.parse.urlencode(data).encode(),
                                 headers={"Accept": "text/plain"})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.read().decode().strip()
    except urllib.error.HTTPError as e:            # the body says which parameter was rejected
        sys.exit(f"EBI returned HTTP {e.code}:\n{e.read().decode()[:800]}")


def valid_appl():
    """the service's own list of member databases, so a typo fails before 45 submissions do"""
    try:
        d = json.loads(get(f"{BASE}/parameterdetails/appl", accept="application/json"))
        return [v["value"] for v in d["values"]["values"]]
    except Exception:
        return []


def get(url, accept="text/plain"):
    req = urllib.request.Request(url, headers={"Accept": accept})
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read().decode()


def read_fasta(p):
    name, seq = None, []
    for line in open(p):
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True, help="contact address required by the EBI service")
    ap.add_argument("--appl", default="PfamA,SMART,Gene3d,SuperFamily,Panther,NCBIfam,Coils,MobiDBLite",
                    help="member databases; the service rejects anything outside its own list, which "
                         "this script fetches and validates before submitting")
    ap.add_argument("--concurrent", type=int, default=3)
    args = ap.parse_args()
    W.mkdir(parents=True, exist_ok=True)
    R.mkdir(exist_ok=True)
    ok = valid_appl()
    bad = [a for a in args.appl.split(",") if ok and a not in ok]
    if bad:
        sys.exit(f"--appl value(s) {bad} not accepted by the service.\nValid values: {', '.join(ok)}")

    seqs = dict(read_fasta(HERE / "queries_3.1.faa"))
    todo = [(k, v) for k, v in seqs.items() if not (W / f"{k}.tsv").exists()]
    print(f"{len(seqs)} queries, {len(todo)} to submit, {len(seqs)-len(todo)} cached")

    running = []                                  # (job_id, query_id)
    while todo or running:
        while todo and len(running) < args.concurrent:
            qid, seq = todo.pop(0)
            jid = post(f"{BASE}/run", {"email": args.email, "sequence": f">{qid}\n{seq}",
                                       "stype": "p", "appl": args.appl,
                                       "goterms": "false", "pathways": "false"})
            print(f"  submitted {qid} -> {jid}", flush=True)
            running.append((jid, qid))
        time.sleep(15)
        still = []
        for jid, qid in running:
            st = get(f"{BASE}/status/{jid}")
            if st == "FINISHED":
                (W / f"{qid}.tsv").write_text(get(f"{BASE}/result/{jid}/tsv"))
                print(f"  done {qid}", flush=True)
            elif st in ("RUNNING", "QUEUED", "PENDING"):
                still.append((jid, qid))
            else:
                print(f"  WARNING {qid}: status {st}", flush=True)
        running = still

    # ---- collect: InterProScan TSV columns are documented at
    # https://interproscan-docs.readthedocs.io/en/latest/OutputFormats.html
    cat = {r["ID"]: r for r in csv.DictReader(open(HERE / "catalogue_3.1.tsv"), delimiter="\t")}
    rows, summary = [], []
    for qid in seqs:
        f = W / f"{qid}.tsv"
        if not f.exists():
            continue
        sigs, ipr, pfam = [], [], []
        for line in f.read_text().splitlines():
            p = line.split("\t")
            if len(p) < 9:
                continue
            rec = {"ID": qid, "Class": cat.get(qid, {}).get("Class", ""),
                   "database": p[3], "signature": p[4], "description": p[5],
                   "start": p[6], "end": p[7], "evalue": p[8],
                   "interpro": p[11] if len(p) > 11 else "", "interpro_desc": p[12] if len(p) > 12 else ""}
            rows.append(rec)
            sigs.append(f"{p[3]}:{p[4]}")
            if rec["interpro"].startswith("IPR"):
                ipr.append(rec["interpro"])
            if p[3] == "Pfam":
                pfam.append(p[4])
        summary.append({"ID": qid, "Class": cat.get(qid, {}).get("Class", ""),
                        "Architecture": cat.get(qid, {}).get("Architecture", ""),
                        "n_signatures": len(sigs), "pfam": ",".join(sorted(set(pfam))),
                        "interpro": ",".join(sorted(set(ipr))),
                        "all_signatures": ",".join(sorted(set(sigs)))})
    if not rows:
        sys.exit("no results collected")
    with open(R / "S1c_interproscan_queries.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    with open(R / "S1d_interproscan_summary.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(summary)
    print(f"\n{len(summary)} queries annotated, {len(rows)} signature hits "
          f"-> results/S1c_interproscan_queries.tsv, S1d_interproscan_summary.tsv")


if __name__ == "__main__":
    main()
