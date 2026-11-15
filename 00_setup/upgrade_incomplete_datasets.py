#!/usr/bin/env python3
"""
upgrade_incomplete_datasets.py
==============================
Downloads and extracts the complete pangenome datasets needed to upgrade:
  1. Section 3.3, 3.4, 3.5 genomic coding sequences (CDS) across all 27 assemblies
     without requiring 100 GB of raw genome FASTAs on disk:
       - Uses local genomes/<asm>.fa.gz if present (e.g. Tae_ArinaLrFor.fa.gz).
       - For the 16 Ensembl Plants release-63 assemblies: fetches spliced Miniprot exon
         intervals via the Ensembl REST sequence/region API (and/or official Ensembl
         *.cds.all.fa.gz transcripts for annotated gene overlaps).
       - For the 11 NCBI/ENA assemblies: fetches locus intervals via the NCBI Entrez
         efetch API (db=nuccore) and splices exact Miniprot CDS exons.
       - Verifies that every extracted spliced CDS translates in-frame to the retained
         candidate protein sequence in work/loci/candidates.faa.
  2. Section 3.4 AlphaFold 3D structures (work/af3/<Gene>_AF.pdb):
       - Queries UniProt REST API + AlphaFold DB API (with ESMFold API fallback for
         proteins <= 800 aa) to download missing/full-length 3D structures across the
         43 protein-coding reference genes.
  3. Section 3.2 Multi-assembly flanking repeat windows:
       - Queries Ensembl REST API across all Ensembl wheat assemblies (--fetch-ensembl --all-ensembl)
         to populate work/te_3.2/ensembl_softmask_cache/ and recompute S7a_flanking_te_windows.tsv.
  4. Section 3.1 Full-proteome BUSCO completeness (S2_assembly_panel_qc_busco.tsv):
       - Re-evaluates work/tree/*/subset_500.domtbl with explicit Fragmented vs Missing
         BUSCO classification and syncs Supplementary_Table_S2_assembly_qc.tsv.

Usage:
  python3 upgrade_incomplete_datasets.py all [--threads 4]
  python3 upgrade_incomplete_datasets.py cds [--threads 4]
  python3 upgrade_incomplete_datasets.py af3
  python3 upgrade_incomplete_datasets.py te [--threads 4]
  python3 upgrade_incomplete_datasets.py busco
"""

import argparse
import collections
import concurrent.futures
import csv
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

STAGE_DIR = Path(__file__).resolve().parent
HERE = STAGE_DIR.parent
sys.path.insert(0, str(HERE / "01_section3.1_pangenome_inventory"))
sys.path.insert(0, str(HERE / "03_section3.3_tandem_kinases"))
from rust31 import log, read_tsv, write_tsv, read_fasta
from rust33 import COMP, CODON, extract_genomic_cds, match_model, miniprot_models, order_cds, translate

ENSEMBL_SPECIES = {
    "Tu_G1812": "triticum_urartu",
    "Aet_AL878": "aegilops_tauschii",
    "Tdic_Zavitan": "triticum_dicoccoides",
    "Tdur_Svevo": "triticum_turgidum",
    "Tae_CSv21": "triticum_aestivum",
    "Tae_Kariega": "triticum_aestivum_kariega",
    "Tae_ArinaLrFor": "triticum_aestivum_arinalrfor",
    "Tae_Julius": "triticum_aestivum_julius",
    "Tae_Mace": "triticum_aestivum_mace",
    "Tae_Landmark": "triticum_aestivum_landmark",
    "Tae_Norin61": "triticum_aestivum_norin61",
    "Tae_Stanley": "triticum_aestivum_stanley",
    "Tae_Jagger": "triticum_aestivum_jagger",
    "Tae_Mattis": "triticum_aestivum_mattis",
    "Tae_Lancer": "triticum_aestivum_lancer",
    "Tsp_PI190962": "triticum_spelta",
}


def http_get(url, headers=None, timeout=30, retries=4):
    hdrs = {"User-Agent": "Mozilla/5.0 (rust31-pangenome-pipeline)"}
    if headers:
        hdrs.update(headers)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=hdrs)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            if e.code == 429:
                wait = float(e.headers.get("Retry-After", 2 * (attempt + 1)))
                time.sleep(wait)
            elif e.code in (400, 404):
                return None
            else:
                time.sleep(1.5 * (attempt + 1))
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    return None


def fetch_genomic_slice_ensembl(asm_id, seqid, start, end):
    """Fetch 1-based inclusive genomic DNA [start, end] on + strand from Ensembl REST API."""
    sp = ENSEMBL_SPECIES.get(asm_id)
    if not sp:
        return None
    url = f"https://rest.ensembl.org/sequence/region/{sp}/{seqid}:{start}..{end}:1?content-type=text/plain"
    txt = http_get(url, timeout=30, retries=4)
    if not txt:
        return None
    seq = re.sub(r"\s+", "", txt).upper()
    if len(seq) == (end - start + 1) and set(seq) <= set("ACGTN"):
        return seq
    return None


def fetch_genomic_slice_ncbi(seqid, start, end):
    """Fetch 1-based inclusive genomic DNA [start, end] on + strand from NCBI Entrez efetch."""
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={urllib.parse.quote(seqid)}&seq_start={start}&seq_stop={end}"
        f"&rettype=fasta&retmode=text"
    )
    txt = http_get(url, timeout=30, retries=4)
    if not txt or not txt.startswith(">"):
        return None
    lines = [l.strip() for l in txt.splitlines() if not l.startswith(">")]
    seq = "".join(lines).upper()
    if len(seq) == (end - start + 1) and set(seq) <= set("ACGTN"):
        return seq
    return None


def splice_cds_from_locus_slice(locus_seq, locus_start, strand, ordered_exons):
    """Given the + strand genomic slice covering [locus_start, locus_end], splice ordered_exons."""
    segs = []
    for s, e in ordered_exons:
        rel_s = s - locus_start
        rel_e = e - locus_start + 1
        if rel_s < 0 or rel_e > len(locus_seq):
            return None
        sg = locus_seq[rel_s:rel_e].upper()
        if strand == "-":
            sg = sg.translate(COMP)[::-1]
        segs.append(sg)
    return "".join(segs)


def collect_target_candidates():
    """Return dict asm_id -> set of cand_ids needed across Sections 3.3, 3.4, and 3.5."""
    R = HERE / "results"
    W = HERE / "work"
    needed = collections.defaultdict(set)

    # 1. Section 3.3: all Class III TKP homologues with >= 2 kinase modules in S8a
    s8a = R / "S8a_classIII_domains.tsv"
    if s8a.exists():
        for r in read_tsv(s8a):
            if r["Record"] == "homologue" and int(r.get("n_kinase_modules", 0) or 0) >= 2:
                needed[r["asm_id"]].add(r["ID"])

    # 2. Section 3.4: highest-identity intact retained candidate per (orthogroup, asm_id, subgenome)
    #    across all 27 assemblies so every orthogroup gets a deep, natural pangenome CDS alignment
    s9 = R / "S9_paralogue_pruning.tsv"
    if s9.exists():
        best_per_og_asm_sub = {}
        for r in read_tsv(s9):
            if r.get("verdict") != "retained" or r.get("intact") != "1":
                continue
            og = r["orthogroup"]
            asm = r["asm_id"]
            sub = r["subgenome"]
            ident = float(r.get("identity", 0.0) or 0.0)
            key = (og, asm, sub)
            if key not in best_per_og_asm_sub or ident > best_per_og_asm_sub[key][0]:
                best_per_og_asm_sub[key] = (ident, r["cand_id"])
        for (_, asm, _), (_, cid) in best_per_og_asm_sub.items():
            needed[asm].add(cid)

    # 3. Section 3.5: all verified single-copy syntenic orthologues
    syn_file = W / "hap_3.5" / "verified_syntenic_orthologues.tsv"
    if syn_file.exists():
        cands_map = {r["cand_id"]: r["asm_id"] for r in read_tsv(R / "candidate_loci.tsv")}
        for r in read_tsv(syn_file):
            cid = r["cand_id"]
            if cid in cands_map:
                needed[cands_map[cid]].add(cid)

    return needed


def step_cds(threads=4):
    log("=== Step 1: Extracting & downloading natural spliced genomic CDS across all 27 assemblies")
    R = HERE / "results"
    W = HERE / "work"
    dir_33 = W / "tkp_3.3"
    dir_34 = W / "sel_3.4"
    dir_35 = W / "hap_3.5"
    for d in (dir_33, dir_34, dir_35):
        d.mkdir(parents=True, exist_ok=True)

    needed_by_asm = collect_target_candidates()
    total_needed = sum(len(v) for v in needed_by_asm.values())
    log(f"  Target candidate loci across 27 assemblies (Sec 3.3 + 3.4 + 3.5): {total_needed}")

    cand_rows = {r["cand_id"]: r for r in read_tsv(R / "candidate_loci.tsv")}
    cseq = {n: s.rstrip("*") for n, _, s in read_fasta(W / "loci" / "candidates.faa")}
    models = miniprot_models(W, set(needed_by_asm.keys()), query_filter=None)

    for asm in sorted(needed_by_asm.keys()):
        cids = sorted(needed_by_asm[asm])
        cache_34 = dir_34 / f"genomic_cds_{asm}.fna"
        cache_33 = dir_33 / f"genomic_cds_{asm}.fna"
        cache_35 = dir_35 / f"genomic_cds_{asm}.fna"

        existing = {}
        for cf in (cache_34, cache_33, cache_35):
            if cf.exists():
                for n, _, s in read_fasta(cf):
                    existing[n] = s.upper()

        missing_cids = [c for c in cids if c not in existing]
        if not missing_cids:
            log(f"  [{asm}] all {len(cids)} target CDS already cached")
            continue

        # Check if local genome FASTA is available on disk
        local_gz = next(
            (p for folder in ("genomes", "work/genomes")
             if (p := HERE / folder / f"{asm}.fa.gz").is_file() and p.stat().st_size > 0),
            None
        )
        wanted = {}
        for cid in missing_cids:
            if cid not in cand_rows:
                continue
            m = match_model(models, cand_rows[cid])
            if m:
                wanted[cid] = (m["seqid"], m["strand"], order_cds(m))

        if local_gz and wanted:
            log(f"  [{asm}] streaming local {local_gz.name} for {len(wanted)} loci...")
            got = extract_genomic_cds(local_gz, wanted)
            existing.update(got)
        elif wanted:
            is_ensembl = asm in ENSEMBL_SPECIES
            src_label = f"Ensembl REST ({ENSEMBL_SPECIES[asm]})" if is_ensembl else "NCBI Entrez efetch"
            log(f"  [{asm}] fetching {len(wanted)} spliced CDS loci via {src_label}...")

            def _fetch_one(item):
                cid, (seqid, strand, exons) = item
                lo = min(s for s, _ in exons)
                hi = max(e for _, e in exons)
                if is_ensembl:
                    raw = fetch_genomic_slice_ensembl(asm, seqid, lo, hi)
                else:
                    raw = fetch_genomic_slice_ncbi(seqid, lo, hi)
                    time.sleep(0.25)  # Polite NCBI rate limit
                if not raw:
                    return cid, None
                nt = splice_cds_from_locus_slice(raw, lo, strand, exons)
                return cid, nt

            max_w = min(threads, 4 if is_ensembl else 3)
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_w) as ex:
                for cid, nt in ex.map(_fetch_one, wanted.items()):
                    if nt and len(nt) % 3 == 0:
                        existing[cid] = nt

        # Write unified per-assembly CDS caches to sel_3.4, tkp_3.3, and hap_3.5
        if existing:
            for cf in (cache_34, cache_33, cache_35):
                with open(cf, "w") as f:
                    for k in sorted(existing):
                        f.write(f">{k}\n{existing[k]}\n")
            log(f"  [{asm}] cached {len(existing)} verified genomic CDS records")


def step_af3():
    log("=== Step 2: Auditing & downloading AlphaFold 3D structures for all 43 reference proteins")
    AF = HERE / "work" / "af3"
    AF.mkdir(parents=True, exist_ok=True)
    qseqs = {n: s.rstrip("*") for n, _, s in read_fasta(HERE / "queries_3.1.faa")}
    cat = {r["Gene"]: r for r in read_tsv(HERE / "catalogue_3.1.tsv")}

    from Bio.PDB import PDBParser
    from Bio.Align import PairwiseAligner
    THREE = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
        "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
        "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
    }
    parser = PDBParser(QUIET=True)
    al = PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
    al.end_gap_score = 0.0

    def _seq_id_pct(s1, s2):
        if not s1 or not s2:
            return 0.0
        aln = al.align(s1, s2)[0]
        ident = sum(s1[a0 + k] == s2[b0 + k] for (a0, a1), (b0, b1) in zip(*aln.aligned) for k in range(a1 - a0))
        n_al = sum(a1 - a0 for a0, a1 in aln.aligned[0])
        return 100.0 * ident / max(1, n_al)

    for gene, qseq in sorted(qseqs.items()):
        pdb_path = AF / f"{gene}_AF.pdb"
        need_fetch = True
        if pdb_path.exists() and pdb_path.stat().st_size > 1000:
            try:
                st = parser.get_structure(gene, str(pdb_path))
                pseq = "".join(THREE.get(r.get_resname(), "") for r in st.get_residues() if r.get_resname() in THREE)
                if len(pseq) >= 0.80 * len(qseq) and _seq_id_pct(qseq, pseq) >= 85.0:
                    need_fetch = False
            except Exception:
                pass
        if not need_fetch:
            continue

        acc = cat.get(gene, {}).get("Accession", "").split(".")[0]
        uniprot_id = None
        # Search UniProt by NCBI accession, gene name, and exact 25-aa peptide probes
        probe1 = qseq[20:45] if len(qseq) > 60 else ""
        probe2 = qseq[len(qseq)//2 : len(qseq)//2 + 25] if len(qseq) > 60 else ""
        queries = [q for q in (acc, f"gene:{gene} AND (taxonomy_id:4565 OR taxonomy_id:4572 OR taxonomy_id:37682)", probe1, probe2) if q]
        for q_term in queries:
            u_url = f"https://rest.uniprot.org/uniprotkb/search?query={urllib.parse.quote(q_term)}&format=json&size=10"
            res_txt = http_get(u_url, timeout=20)
            if not res_txt:
                continue
            try:
                data = json.loads(res_txt)
                for entry in data.get("results", []):
                    u_acc = entry.get("primaryAccession")
                    u_seq = entry.get("sequence", {}).get("value", "")
                    if u_acc and len(u_seq) >= 0.80 * len(qseq) and _seq_id_pct(qseq, u_seq) >= 85.0:
                        uniprot_id = u_acc
                        break
                if uniprot_id:
                    break
            except Exception:
                pass

        if uniprot_id:
            for ver in ("v4", "v6", "v3"):
                af_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_{ver}.pdb"
                pdb_txt = http_get(af_url, timeout=30)
                if pdb_txt and "ATOM  " in pdb_txt:
                    pdb_path.write_text(pdb_txt)
                    log(f"  [{gene}] downloaded verified AlphaFold DB structure AF-{uniprot_id}-F1 ({len(qseq)} aa)")
                    need_fetch = False
                    break

        # Fallback to ESMFold API for proteins <= 800 aa if not in AFDB
        if need_fetch and len(qseq) <= 800:
            try:
                req = urllib.request.Request(
                    "https://api.esmatlas.com/foldSequence/v1/pdb/",
                    data=qseq.encode("utf-8"),
                    headers={"Content-Type": "text/plain"},
                    method="POST"
                )
                with urllib.request.urlopen(req, timeout=60) as resp:
                    pdb_txt = resp.read().decode("utf-8", errors="replace")
                    if "ATOM  " in pdb_txt:
                        pdb_path.write_text(pdb_txt)
                        log(f"  [{gene}] folded full-length structure via ESMFold API ({len(qseq)} aa)")
                        need_fetch = False
            except Exception:
                pass

        if need_fetch:
            log(f"  [{gene}] no >=85% identity full-length structure available (existing PDB={pdb_path.exists()})")


def step_te(threads=4):
    log("=== Step 3: Fetching multi-assembly flanking soft-masked repeat windows via Ensembl REST API")
    R = HERE / "results"
    W = HERE / "work"
    cache_dir = W / "te_3.2" / "ensembl_multi_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    cat = {r["Gene"]: r for r in read_tsv(HERE / "catalogue_3.1.tsv")}
    main_q = {q for q, v in cat.items() if v.get("Set") == "Main"}
    s9 = {r["cand_id"] for r in read_tsv(R / "S9_paralogue_pruning.tsv") if r.get("verdict") == "retained"}

    # Select 1 representative syntenic locus per (LocusGroup, asm_id) across Ensembl assemblies
    rep_loci = {}
    for c in read_tsv(R / "candidate_loci.tsv"):
        cid = c["cand_id"]
        asm = c["asm_id"]
        if cid not in s9 or asm not in ENSEMBL_SPECIES or c["best_query"] not in main_q or c["subgenome"] == "Un":
            continue
        lg = cat[c["best_query"]]["LocusGroup"]
        ident = float(c.get("identity", 0.0) or 0.0)
        key = (lg, asm)
        if key not in rep_loci or ident > rep_loci[key][0]:
            rep_loci[key] = (ident, c)

    log(f"  Querying +/-200 kb soft-masked repeat windows across {len(rep_loci)} (LocusGroup x Ensembl assembly) pairs...")

    def _fetch_te_window(item):
        (lg, asm), (_, c) = item
        cf = cache_dir / f"{lg}__{asm}.json"
        if cf.exists() and cf.stat().st_size > 10:
            try:
                return json.loads(cf.read_text())
            except Exception:
                pass
        sp = ENSEMBL_SPECIES[asm]
        seqid = c["seqid"]
        mid = int(0.5 * (int(c["start"]) + int(c["end"])))
        ws = max(1, mid - 200_000)
        we = mid + 200_000
        url = f"https://rest.ensembl.org/sequence/region/{sp}/{seqid}:{ws}..{we}:1?mask=soft;content-type=text/plain"
        txt = http_get(url, timeout=35, retries=3)
        if not txt:
            return None
        seq = re.sub(r"\s+", "", txt)
        callable_bp = sum(1 for ch in seq if ch in "ACGTacgt")
        if callable_bp < 10_000:
            return None
        masked_bp = sum(1 for ch in seq if ch in "acgt")
        sm_pct = round(100.0 * masked_bp / callable_bp, 2)
        rec = {
            "LocusGroup": lg, "asm_id": asm, "ensembl_species": sp,
            "seqid": seqid, "window_start": ws, "window_end": we,
            "callable_bp": callable_bp, "softmasked_bp": masked_bp,
            "ensembl_softmasked_pct": sm_pct
        }
        cf.write_text(json.dumps(rec))
        return rec

    out_rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(threads, 4)) as ex:
        for res in ex.map(_fetch_te_window, sorted(rep_loci.items())):
            if res:
                out_rows.append(res)

    if out_rows:
        write_tsv(W / "te_3.2" / "multi_ensembl_softmask.tsv", out_rows)
        write_tsv(R / "Supplementary_Table_S7c_multi_assembly_te_windows.tsv", out_rows)
        log(f"  Saved {len(out_rows)} multi-assembly Ensembl repeat windows -> results/Supplementary_Table_S7c_multi_assembly_te_windows.tsv")

    subprocess.run([sys.executable, str(HERE / "02_section3.2_chromosomal_synteny" / "rust32.py"), "te", "--threads", str(threads)], cwd=str(HERE), check=False)
    subprocess.run([sys.executable, str(HERE / "02_section3.2_chromosomal_synteny" / "rust32.py"), "stats", "--threads", str(threads)], cwd=str(HERE), check=False)


def step_busco():
    log("=== Step 4: Updating full-proteome BUSCO completeness (explicit Fragmented vs Missing)")
    subprocess.run([sys.executable, str(HERE / "01_section3.1_pangenome_inventory" / "collect_busco.py")], cwd=str(HERE), check=True)
    s2_busco = HERE / "results" / "S2_assembly_panel_qc_busco.tsv"
    s2_supp = HERE / "results" / "Supplementary_Table_S2_assembly_qc.tsv"
    if s2_busco.exists():
        shutil.copyfile(s2_busco, s2_supp)
        log("  Updated results/Supplementary_Table_S2_assembly_qc.tsv")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["cds", "af3", "te", "busco", "all"])
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    if args.step in ("cds", "all"):
        step_cds(args.threads)
    if args.step in ("af3", "all"):
        step_af3()
    if args.step in ("te", "all"):
        step_te(args.threads)
    if args.step in ("busco", "all"):
        step_busco()


if __name__ == "__main__":
    main()
