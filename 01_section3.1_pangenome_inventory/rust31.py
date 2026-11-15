#!/usr/bin/env python3
"""
rust31.py - generates the dataset behind Results section 3.1
("Super-pangenome inventory and core-to-dispensable partitioning").

Steps (run in order, each is resumable):
  setup        download Pfam-A and build the diagnostic HMM subset; get NLR-Annotator files
  prep         resolve download URLs, fetch proteomes + GFF3 (annotated assemblies), assembly QC
  genome       stream each genome chromosome by chromosome -> miniprot (+ NLR-Annotator)
  annot        longest-isoform proteomes, DIAMOND vs queries, GFF3 CDS footprints
  loci         collapse alignments to loci, subgenome assignment, Pfam anchor check,
               rescue classification (BEDTools), OrthoFinder input files
  orthofinder  run OrthoFinder on the candidate-homologue set
  stats        occupancy, compartments, CNV, rescue statistics, tests, summary table
  all          everything above in order

Usage:  python rust31.py <step> [--config config_3.1.json] [--asm ASM_ID]
"""
import argparse, collections, csv, gzip, json, math, os, re, shutil, subprocess, sys, time
import urllib.request
from pathlib import Path

STAGE_DIR = Path(__file__).resolve().parent
ROOT = STAGE_DIR.parent
HERE = ROOT


# ----------------------------------------------------------------------------- utils
def log(*a):
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def load_cfg(path):
    cfg = json.load(open(path))
    base = Path(path).resolve().parent
    for k in ("workdir", "resultsdir", "panel", "catalogue", "queries"):
        cfg[k] = str((base / cfg[k]).resolve())
    Path(cfg["workdir"]).mkdir(parents=True, exist_ok=True)
    Path(cfg["resultsdir"]).mkdir(parents=True, exist_ok=True)
    return cfg


def run(cmd, **kw):
    log("RUN", cmd if isinstance(cmd, str) else " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True, shell=isinstance(cmd, str), **kw)


def http_text(url, tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read().decode()
        except Exception as e:  # transient SSL/HTTP errors
            log(f"  retry {i+1}/{tries} {url}: {e}")
            time.sleep(5 * (i + 1))
    raise RuntimeError(f"failed: {url}")


def listdir_http(url):
    return re.findall(r'href="([^"?/][^"]*)"', http_text(url))


def download(url, dest):
    """curl with resume; skips if a .ok marker exists"""
    dest = Path(dest)
    ok = dest.with_name(dest.name + ".ok")
    if ok.exists() and dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    if url.startswith("/") or url.startswith("file:"):
        src = url.replace("file:", "")
        shutil.copy(src, dest)
    else:
        for i in range(6):
            r = subprocess.run(["curl", "-L", "--fail", "--retry", "5", "--progress-bar", "-C", "-", "-o", str(dest), url])
            if r.returncode == 0:
                break
            if r.returncode == 33:  # server refuses range -> restart from 0
                dest.unlink(missing_ok=True)
            time.sleep(10)
        else:
            raise RuntimeError(f"download failed: {url}")
    ok.touch()
    return dest


def is_gz(path):
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def open_any(path):
    """text handle for plain or gzipped files (detected by magic bytes, not by extension)"""
    path = str(path)
    if is_gz(path):
        p = subprocess.Popen(["gzip", "-dc", path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        return _Checked(p)
    return open(path)


class _Checked:
    """iterate a gzip -dc pipe and raise if gzip failed (truncated or corrupt download)"""
    def __init__(self, p):
        self.p = p

    def __iter__(self):
        for line in self.p.stdout:
            yield line
        if self.p.wait() != 0:
            raise RuntimeError("gzip -dc failed (corrupt or truncated file?): " + " ".join(self.p.args))


def read_fasta(path):
    name, desc, seq = None, None, []
    for line in open_any(path):
        if line.startswith(">"):
            if name is not None:
                yield name, desc, "".join(seq)
            desc = line[1:].strip()
            name = desc.split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, desc, "".join(seq)


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, cols=None):
    rows = list(rows)
    if cols is None:
        cols = list(rows[0].keys()) if rows else []
    with open(path, "w") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def panel_rows(cfg):
    rows = [r for r in read_tsv(cfg["panel"]) if r["include"].strip().lower() == "yes"]
    for r in rows:
        r["annotated"] = bool(r.get("gff3", "").strip() and r.get("pep_fasta", "").strip()) or r["source"] == "ensembl"
    return rows


def catalogue(cfg):
    return {r["ID"]: r for r in read_tsv(cfg["catalogue"])}


def which_or_die(*tools):
    miss = [t for t in tools if shutil.which(t) is None]
    if miss:
        sys.exit(f"Missing tools in PATH: {miss}. Activate the conda env (see README).")


def java_bin():
    """locate a working java: $JAVA_HOME, conda lib/jvm, then PATH"""
    cands = []
    if os.environ.get("JAVA_HOME"):
        cands.append(Path(os.environ["JAVA_HOME"]) / "bin" / "java")
    for pre in (os.environ.get("CONDA_PREFIX"), sys.prefix):
        if pre:
            cands.append(Path(pre) / "lib" / "jvm" / "bin" / "java")
            cands.append(Path(pre) / "bin" / "java")
    if shutil.which("java"):
        cands.append(Path(shutil.which("java")))
    for c in cands:
        if c.exists() and subprocess.run([str(c), "-version"], capture_output=True).returncode == 0:
            return str(c)
    sys.exit("No working Java runtime found. Install it with: conda install -n rust31 -c conda-forge openjdk=17")


# ----------------------------------------------------------------------------- setup
PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
PFAM_KEEP = re.compile(r"^(NB-ARC|Rx_N|WHD_DRP|Pkinase|PK_Tyr_Ser-Thr|Ank(_[0-9])?|NAM|zf-BED|ABC_tran|"
                       r"ABC2_membrane.*|PDR_assoc|Sugar_tr|START|VWA.*|WRKY|LRR_.*|DUF3475|DUF668)$")
NLRA = "https://raw.githubusercontent.com/steuernb/NLR-Annotator/NLR-Annotator-2/"
NLRA_ALT = "https://github.com/steuernb/NLR-Annotator/raw/NLR-Annotator-2/"


def step_setup(cfg, args):
    which_or_die("hmmpress")
    db = Path(cfg["workdir"]) / "db"
    db.mkdir(parents=True, exist_ok=True)
    diag = db / "diag.hmm"
    if not (db / "diag.hmm.h3m").exists():
        gz = download(PFAM_URL, db / "Pfam-A.hmm.gz")
        kept, rec, name = [], [], None
        for line in open_any(gz):
            rec.append(line)
            if line.startswith("NAME "):
                name = line.split()[1]
            if line.startswith("//"):
                if name and PFAM_KEEP.match(name):
                    kept.append("".join(rec))
                rec, name = [], None
        diag.write_text("".join(kept))
        log(f"diagnostic HMM subset: {len(kept)} Pfam models")
        for ext in ("h3m", "h3i", "h3f", "h3p"):
            Path(f"{diag}.{ext}").unlink(missing_ok=True)
        run(["hmmpress", str(diag)])
        relnote = http_text("https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt")[:400]
        (db / "pfam_release.txt").write_text(relnote)
        gz.unlink()
    if cfg["nlr_annotator"]["enabled"]:
        java_bin()
        for f in ("NLR-Annotator-v2.1b.jar", "src/mot.txt", "src/store.txt"):
            dest = db / Path(f).name
            if dest.exists() and dest.stat().st_size > 0:      # e.g. placed there by hand
                continue
            try:
                download(NLRA + f, dest)
            except RuntimeError:
                download(NLRA_ALT + f, dest)
    log("setup done")


# ----------------------------------------------------------------------------- prep
def ensembl_urls(sp, rel):
    B = f"https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-{rel}"
    dna = [f for f in listdir_http(f"{B}/fasta/{sp}/dna/") if f.endswith(".dna.toplevel.fa.gz")]
    pep = [f for f in listdir_http(f"{B}/fasta/{sp}/pep/") if f.endswith(".pep.all.fa.gz")]
    gff = [f for f in listdir_http(f"{B}/gff3/{sp}/")
           if f.endswith(f".{rel}.gff3.gz") and ".chr." not in f and "chromosome" not in f
           and "scaffold" not in f and "abinitio" not in f]
    if not (dna and pep and gff):
        raise RuntimeError(f"Ensembl files not found for {sp}: dna={dna} pep={pep} gff={gff}")
    return (f"{B}/fasta/{sp}/dna/{dna[0]}", f"{B}/fasta/{sp}/pep/{pep[0]}", f"{B}/gff3/{sp}/{sorted(gff, key=len)[0]}")


def ncbi_genome_url(acc):
    p, d = acc[:3], acc[4:13]
    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{p}/{d[0:3]}/{d[3:6]}/{d[6:9]}/"
    dirs = [x for x in listdir_http(base) if x.startswith(acc + "_")]
    if not dirs:
        raise RuntimeError(f"NCBI FTP dir not found for {acc}")
    dn = dirs[0].rstrip("/")
    return f"{base}{dn}/{dn}_genomic.fna.gz"


def ncbi_stats(acc):
    # all_assemblies: also return superseded versions (Ensembl Plants uses some older .1 versions)
    u = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/dataset_report?filters.assembly_version=all_assemblies"
    try:
        d = json.loads(http_text(u))["reports"][0]
    except Exception as e:
        log(f"  NCBI stats failed for {acc}: {e}")
        return {}
    s = d.get("assembly_stats", {})
    return {"assembly_name": d["assembly_info"].get("assembly_name"),
            "assembly_level": d["assembly_info"].get("assembly_level"),
            "ncbi_version_status": d["assembly_info"].get("assembly_status"),
            "total_length_bp": s.get("total_sequence_length"),
            "contig_n50_bp": s.get("contig_n50"), "scaffold_n50_bp": s.get("scaffold_n50"),
            "n_contigs": s.get("number_of_contigs"), "n_scaffolds": s.get("number_of_scaffolds")}


def step_prep(cfg, args):
    W = Path(cfg["workdir"])
    rows, qc = [], []
    for r in panel_rows(cfg):
        a = r["asm_id"]
        g, p, f = r.get("genome_fasta", "").strip(), r.get("pep_fasta", "").strip(), r.get("gff3", "").strip()
        if r["source"] == "ensembl" and not g:
            g, p, f = ensembl_urls(r["ensembl_name"], cfg["ensembl_release"])
        elif r["source"] == "ncbi" and not g:
            g = ncbi_genome_url(r["accession"])
        rows.append({"asm_id": a, "genome_url": g, "pep_url": p, "gff3_url": f})
        if p and f:
            download(p, W / "annot" / a / "pep.fa.gz")
            download(f, W / "annot" / a / "genes.gff3.gz")
        st = ncbi_stats(r["accession"]) if r["accession"] else {}
        qc.append({**{k: r[k] for k in ("asm_id", "species", "accession_label", "ploidy", "genomes", "source", "accession")},
                   "annotation_gff3": "yes" if (p and f) else "no", **st})
        log(a, "ok", "(annotated)" if p and f else "(genome only)")
    write_tsv(W / "urls.tsv", rows)
    write_tsv(Path(cfg["resultsdir"]) / "S2_assembly_panel_qc.tsv", qc,
              ["asm_id", "species", "accession_label", "ploidy", "genomes", "source", "accession", "annotation_gff3",
               "assembly_name", "assembly_level", "ncbi_version_status", "total_length_bp", "contig_n50_bp", "scaffold_n50_bp", "n_contigs", "n_scaffolds"])
    log("prep done:", len(rows), "assemblies")


# ----------------------------------------------------------------------------- genome (miniprot / NLR-Annotator)
CHR_RE = re.compile(r"^(?:chr|Chr|CHR)?0?([1-7])([A-Z])?(?:[_.\-].*)?$")


def seq_label(seqid, desc, source):
    if source == "ncbi":
        m = re.search(r"chromosome:?\s+([0-9]+[A-Za-z]?)\b", desc)   # NCBI "chromosome 1A" and ENA "chromosome: 1A"
        return m.group(1) if m else "Un"
    return seqid


def subgenome_of(label, genomes):
    if len(genomes) == 1:          # diploid: every sequence (incl. unplaced) belongs to its one genome
        return genomes
    m = CHR_RE.match(label)
    if not m:
        return "Un"
    letter = m.group(2)
    return letter if (letter and letter in genomes) else "Un"


def step_genome(cfg, args):
    which_or_die("miniprot", "gzip")
    W = Path(cfg["workdir"])
    urls = {r["asm_id"]: r for r in read_tsv(W / "urls.tsv")}
    panel = {r["asm_id"]: r for r in panel_rows(cfg)}
    todo = [args.asm] if args.asm else list(panel)
    mp = cfg["miniprot"]
    nlr = cfg["nlr_annotator"]
    db = W / "db"
    JAVA = java_bin() if nlr["enabled"] else None
    for a in todo:
        out = W / "miniprot" / a
        out.mkdir(parents=True, exist_ok=True)
        if (out / "ALL.done").exists():
            log(a, "already done"); continue
        r = panel[a]
        gz = download(urls[a]["genome_url"], W / "genomes" / f"{a}.fa.gz")
        seqmap = open(out / "seqmap.tsv", "w")
        seqmap.write("seqid\tlength\tlabel\tsubgenome\n")
        batch, bsize, bidx = [], 0, 0

        def flush(name, recs):
            done = out / f"{name}.done"
            if done.exists():
                return
            fa = out / f"{name}.fa"
            with open(fa, "w") as fh:
                for sid, s in recs:
                    fh.write(f">{sid}\n")
                    for i in range(0, len(s), 80):
                        fh.write(s[i:i + 80] + "\n")
            with open(out / f"{name}.gff", "w") as fh:
                subprocess.run(["miniprot", "-t", str(cfg["threads"]), "-G", mp["max_intron"],
                                f"--outn={mp['outn']}", f"--outs={mp['outs']}", "-u", "--gff", "--trans",
                                str(fa), cfg["queries"]], stdout=fh, stderr=subprocess.DEVNULL, check=True)
            if nlr["enabled"]:
                subprocess.run([JAVA, f"-Xmx{nlr['java_mem']}", "-jar", str(db / "NLR-Annotator-v2.1b.jar"),
                                "-i", str(fa), "-x", str(db / "mot.txt"), "-y", str(db / "store.txt"),
                                "-o", str(out / f"{name}.nlr.txt"), "-b", str(out / f"{name}.nlr.bed"),
                                "-t", str(cfg["threads"])], check=True,
                               stdout=subprocess.DEVNULL, stderr=open(out / f"{name}.nlr.log", "w"))
            fa.unlink()
            done.touch()
            log(f"  {a}:{name} done ({sum(len(s) for _, s in recs)/1e6:.0f} Mb)")

        for sid, desc, s in read_fasta(gz):
            lab = seq_label(sid, desc, r["source"])
            sub = subgenome_of(lab, r["genomes"])
            seqmap.write(f"{sid}\t{len(s)}\t{lab}\t{sub}\n")
            if len(s) >= 50_000_000:
                flush(re.sub(r"[^A-Za-z0-9_.-]", "_", sid), [(sid, s)])
            else:
                batch.append((sid, s)); bsize += len(s)
                if bsize >= cfg["small_seq_batch_bp"]:
                    flush(f"small_{bidx:03d}", batch); batch, bsize, bidx = [], 0, bidx + 1
        if batch:
            flush(f"small_{bidx:03d}", batch)
        seqmap.close()
        nseq = sum(1 for _ in open(out / "seqmap.tsv")) - 1
        if nseq == 0:
            raise RuntimeError(f"{a}: no sequences read from {gz}")
        (out / "ALL.done").touch()
        if not cfg["keep_genomes"]:
            gz.unlink(missing_ok=True)
        log(a, "genome scan complete")


# ----------------------------------------------------------------------------- miniprot parsing
def parse_miniprot_dir(d):
    """yield alignment dicts from all *.gff in a miniprot output dir"""
    for gf in sorted(Path(d).glob("*.gff")):
        cur = None
        for line in open(gf):
            if line.startswith("##PAF"):
                if cur and cur.get("seqid"):
                    yield cur
                f = line.rstrip("\n").split("\t")
                if f[5] == "*":
                    cur = None; continue
                tags = dict(t.split(":", 2)[0::2] for t in f[13:] if t.count(":") >= 2)
                cur = {"query": f[1], "qlen": int(f[2]), "qs": int(f[3]), "qe": int(f[4]), "strand": f[5],
                       "seqid": f[6], "score": int(tags.get("AS", 0)), "fs": int(tags.get("fs", 0)),
                       "st": int(tags.get("st", 0)), "cds": [], "aa": ""}
            elif line.startswith("##STA") and cur is not None:
                # strip anything that is not a residue: one miniprot ##STA line (Tae_Mace 1B, Lr1) starts
                # with a space, which makes diamond makedb hang forever at "Loading sequences"
                cur["aa"] = re.sub(r"[^A-Za-z*]", "", line.rstrip("\n").split("\t")[1])
            elif not line.startswith("#") and cur is not None:
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                if f[2] == "mRNA":
                    at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
                    cur.update(start=int(f[3]), end=int(f[4]), identity=float(at.get("Identity", 0)),
                               positive=float(at.get("Positive", 0)), mp_id=at.get("ID"))
                elif f[2] == "CDS":
                    cur["cds"].append((int(f[3]), int(f[4])))
        if cur and cur.get("seqid"):
            yield cur


def qcov(a):
    return (a["qe"] - a["qs"]) / a["qlen"]


def collapse(alns):
    """merge overlapping alignments (same seqid+strand) into loci; returns list of lists"""
    loci = []
    by = collections.defaultdict(list)
    for a in alns:
        by[(a["seqid"], a["strand"])].append(a)
    for k, v in by.items():
        v.sort(key=lambda x: x["start"])
        cur, cend = [], -1
        for a in v:
            if cur and a["start"] <= cend:
                cur.append(a); cend = max(cend, a["end"])
            else:
                if cur:
                    loci.append(cur)
                cur, cend = [a], a["end"]
        if cur:
            loci.append(cur)
    return loci


# ----------------------------------------------------------------------------- annot
def step_annot(cfg, args):
    which_or_die("diamond")
    W = Path(cfg["workdir"])
    for r in panel_rows(cfg):
        a = r["asm_id"]
        d = W / "annot" / a
        if not (d / "pep.fa.gz").exists():
            continue
        if (d / "annot.done").exists():
            continue
        # transcript -> gene from GFF3 ; CDS footprints per gene
        t2g, cds = {}, collections.defaultdict(list)
        gene_loc = {}
        for line in open_any(d / "genes.gff3.gz"):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
            if f[2] == "gene":
                gene_loc[at.get("ID", "").replace("gene:", "")] = (f[0], int(f[3]), int(f[4]), f[6])
            elif "ID" in at and "Parent" in at and (f[2] in ("mRNA", "transcript") or at["Parent"].startswith("gene:")):
                t2g[at["ID"].replace("transcript:", "")] = at["Parent"].replace("gene:", "")
            elif f[2] == "CDS" and "Parent" in at:
                cds[at["Parent"].replace("transcript:", "")].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))
        # longest isoform per gene
        best = {}
        for pid, desc, s in read_fasta(d / "pep.fa.gz"):
            m = re.search(r"\bgene:(\S+)", desc)
            tid = re.search(r"\btranscript:(\S+)", desc)
            tid = tid.group(1) if tid else pid
            g = m.group(1) if m else t2g.get(tid, pid)
            s = s.rstrip("*").replace("*", "X")
            if g not in best or len(s) > len(best[g][1]):
                best[g] = (tid, s)
        with open(d / "longest.faa", "w") as fh:
            for g, (tid, s) in best.items():
                fh.write(f">{g}\n{s}\n")
        with open(d / "gene_cds.bed", "w") as fh:
            for g, (tid, s) in best.items():
                for seqid, s0, e0, strand in cds.get(tid, []):
                    fh.write(f"{seqid}\t{s0}\t{e0}\t{g}\t0\t{strand}\n")
        write_tsv(d / "gene_loc.tsv", [{"gene": g, "seqid": v[0], "start": v[1], "end": v[2], "strand": v[3],
                                         "aa_len": len(best[g][1]) if g in best else 0}
                                        for g, v in gene_loc.items() if g in best])
        run(["sort", "-k1,1", "-k2,2n", "-o", str(d / "gene_cds.bed"), str(d / "gene_cds.bed")])
        run(["diamond", "makedb", "--in", str(d / "longest.faa"), "-d", str(d / "longest"), "--quiet"])
        run(["diamond", "blastp", "-q", cfg["queries"], "-d", str(d / "longest"), "--ultra-sensitive",
             "-e", str(cfg["diamond"]["evalue"]), "-k", str(cfg["diamond"]["max_target_seqs"]),
             "-p", str(cfg["threads"]), "--quiet", "-o", str(d / "diamond.tsv"),
             "-f", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend",
             "sstart", "send", "evalue", "bitscore"])
        (d / "annot.done").touch()
        log(a, "annotation-based search done:", len(best), "genes")


def diamond_best(path):
    """best HSP per (query, gene) by bitscore -> identity (0-1), qcov"""
    best = {}
    if not Path(path).exists():
        return best
    for line in open(path):
        q, s, pid, L, qlen, slen, qs, qe, ss, se, ev, bs = line.rstrip("\n").split("\t")
        k = (q, s)
        bs = float(bs)
        if k not in best or bs > best[k]["bits"]:
            best[k] = {"identity": float(pid) / 100, "qcov": (int(qe) - int(qs) + 1) / int(qlen),
                       "bits": bs, "slen": int(slen), "qlen": int(qlen)}
    return best


# ----------------------------------------------------------------------------- loci
def bed_overlap(a_bed, b_bed):
    """bedtools intersect -s; returns dict a_name -> set(b_name)"""
    res = collections.defaultdict(set)
    if not Path(b_bed).exists() or os.path.getsize(b_bed) == 0 or os.path.getsize(a_bed) == 0:
        return res
    out = subprocess.run(["bedtools", "intersect", "-a", str(a_bed), "-b", str(b_bed), "-s", "-wa", "-wb"],
                         check=True, capture_output=True, text=True).stdout
    for line in out.splitlines():
        f = line.split("\t")
        res[f[3]].add(f[9])
    return res


def anchor_ok(rule, doms):
    for part in rule.split(";"):
        alts = part.split("|")
        if not any(re.fullmatch(x, d) for x in alts for d in doms):
            return False
    return True


def step_loci(cfg, args):
    which_or_die("bedtools", "hmmscan")
    W, R = Path(cfg["workdir"]), Path(cfg["resultsdir"])
    cat = catalogue(cfg)
    grp_of = {q: cat[q]["LocusGroup"] for q in cat}
    cand_id, cand_cov = cfg["candidate"]["min_identity"], cfg["candidate"]["min_qcov"]
    hom_id, hom_cov = cfg["homolog"]["min_identity"], cfg["homolog"]["min_qcov"]
    (W / "loci").mkdir(exist_ok=True)
    cands, seqs, hom_rows, best_rows = [], {}, [], []
    for r in panel_rows(cfg):
        a = r["asm_id"]
        md = W / "miniprot" / a
        if not (md / "ALL.done").exists():
            log(f"WARNING {a}: genome step not finished - skipped"); continue
        sub = {x["seqid"]: x["subgenome"] for x in read_tsv(md / "seqmap.tsv")}
        alns = list(parse_miniprot_dir(md))
        for a_ in alns:
            a_["qcov"] = qcov(a_)
            a_["sub"] = sub.get(a_["seqid"], "Un")
        ad = W / "annot" / a
        annotated = (ad / "annot.done").exists()
        dbest = diamond_best(ad / "diamond.tsv") if annotated else {}
        gloc = {x["gene"]: x for x in read_tsv(ad / "gene_loc.tsv")} if annotated else {}

        # best alignment per (query, subgenome)  -> threshold sweep table
        bq = {}
        for x in alns:
            k = (x["query"], x["sub"])
            if k not in bq or x["score"] > bq[k]["score"]:
                bq[k] = x
        for (q, s), x in bq.items():
            best_rows.append({"asm_id": a, "subgenome": s, "query": q, "identity": x["identity"],
                              "qcov": round(x["qcov"], 4), "fs": x["fs"], "st": x["st"]})

        # ---- (1) candidate loci for OrthoFinder (permissive)
        cl = [x for x in alns if x["identity"] >= cand_id and x["qcov"] >= cand_cov]
        loci = collapse(cl)
        rep_bed = W / "loci" / f"{a}.cand_rep.bed"
        recs = []
        with open(rep_bed, "w") as fh:
            for i, L in enumerate(loci):
                rep = max(L, key=lambda x: x["score"])
                lid = f"{a}__{rep['sub']}__m{i:06d}"
                for s0, e0 in rep["cds"]:
                    fh.write(f"{rep['seqid']}\t{s0-1}\t{e0}\t{lid}\t0\t{rep['strand']}\n")
                recs.append((lid, rep, L))
        ov = bed_overlap(rep_bed, ad / "gene_cds.bed") if annotated else {}
        covered_genes = set().union(*ov.values()) if ov else set()
        for lid, rep, L in recs:
            intact = rep["qcov"] >= cfg["intact"]["min_qcov"] and rep["fs"] == 0 and rep["st"] == 0
            aa = rep["aa"].rstrip("*").replace("*", "X")
            seqs[lid] = aa
            cands.append({"cand_id": lid, "asm_id": a, "subgenome": rep["sub"], "seqid": rep["seqid"],
                          "start": min(x["start"] for x in L), "end": max(x["end"] for x in L),
                          "strand": rep["strand"], "origin": "miniprot", "best_query": rep["query"],
                          "identity": rep["identity"], "qcov": round(rep["qcov"], 4), "frameshift": rep["fs"],
                          "stop": rep["st"], "intact": int(intact), "aa_len": len(aa),
                          "overlapping_genes": ",".join(sorted(ov.get(lid, [])))})
        # annotation-only extras (gene model hit but no miniprot locus)
        extra = {}
        for (q, g), h in dbest.items():
            if h["identity"] >= cand_id and h["qcov"] >= cand_cov and g not in covered_genes and g in gloc:
                if g not in extra or h["bits"] > extra[g][1]["bits"]:
                    extra[g] = (q, h)
        if extra:
            prot = {n: s for n, _, s in read_fasta(ad / "longest.faa") if n in extra}
            for i, (g, (q, h)) in enumerate(sorted(extra.items())):
                gl = gloc[g]
                s_ = sub.get(gl["seqid"], "Un")
                lid = f"{a}__{s_}__g{i:06d}"
                seqs[lid] = prot[g]
                intact = len(prot[g]) >= cfg["intact"]["min_qcov"] * h["qlen"] and "X" not in prot[g][:-1]
                cands.append({"cand_id": lid, "asm_id": a, "subgenome": s_, "seqid": gl["seqid"],
                              "start": gl["start"], "end": gl["end"], "strand": gl["strand"], "origin": "annotation_only",
                              "best_query": q, "identity": h["identity"], "qcov": round(h["qcov"], 4),
                              "frameshift": "", "stop": "", "intact": int(intact), "aa_len": len(prot[g]),
                              "overlapping_genes": g})

        # ---- (2) homologue-level loci -> rescue classification
        # alignments of all queries are collapsed together, so one genomic locus is counted once and
        # assigned to the LocusGroup of its best-scoring query
        hl = [x for x in alns if x["identity"] >= hom_id and x["qcov"] >= hom_cov]
        hbed = W / "loci" / f"{a}.hom.bed"
        hrec = []
        with open(hbed, "w") as fh:
            for j, L in enumerate(collapse(hl)):
                rep = max(L, key=lambda x: x["score"])
                gname = grp_of.get(rep["query"], rep["query"])
                hid = f"{a}|{gname}|{j}"
                for s0, e0 in rep["cds"]:
                    fh.write(f"{rep['seqid']}\t{s0-1}\t{e0}\t{hid}\t0\t{rep['strand']}\n")
                hrec.append((hid, gname, rep, L))
        hov = bed_overlap(hbed, ad / "gene_cds.bed") if annotated else {}
        # gene models that reach the homologue threshold for any query; each gene -> group of best query
        gene_best = {}
        for (q, g), h in dbest.items():
            if h["identity"] >= hom_id and h["qcov"] >= hom_cov:
                if g not in gene_best or h["bits"] > gene_best[g][1]:
                    gene_best[g] = (grp_of.get(q, q), h["bits"], q)
        for hid, gname, rep, L in hrec:
            genes = hov.get(hid, set())
            if not annotated:
                status = "no_annotation"
            elif genes & set(gene_best):
                status = "annotated"
            elif len(genes) == 0:
                status = "rescued_unannotated"
            elif len(genes) >= 2:
                status = "rescued_split_model"
            else:
                status = "rescued_misannotated_model"
            hom_rows.append({"asm_id": a, "annotated": int(annotated), "LocusGroup": gname, "hom_locus": hid,
                             "subgenome": rep["sub"], "seqid": rep["seqid"], "start": rep["start"], "end": rep["end"],
                             "strand": rep["strand"], "best_query": rep["query"], "identity": rep["identity"],
                             "qcov": round(rep["qcov"], 4), "frameshift": rep["fs"], "stop": rep["st"],
                             "queries_at_locus": ",".join(sorted({x["query"] for x in L})),
                             "overlapping_genes": ",".join(sorted(genes)), "status": status})
        # homologue gene models missed by miniprot are kept as annotated loci
        hgenes = set().union(*hov.values()) if hov else set()
        for g, (gname, bits, q) in gene_best.items():
            if g in hgenes:
                continue
            gl = gloc.get(g, {})
            hom_rows.append({"asm_id": a, "annotated": 1, "LocusGroup": gname, "hom_locus": f"{a}|{gname}|{g}",
                             "subgenome": sub.get(gl.get("seqid", ""), "Un"), "seqid": gl.get("seqid", ""),
                             "start": gl.get("start", ""), "end": gl.get("end", ""), "strand": gl.get("strand", ""),
                             "best_query": q, "identity": "", "qcov": "", "frameshift": "", "stop": "",
                             "queries_at_locus": q, "overlapping_genes": g, "status": "annotated_only_not_in_miniprot"})
        log(a, f"{len(alns)} alignments, {len(recs)} miniprot candidate loci, {len(extra)} annotation-only, "
               f"{len(hrec)} homologue loci")

    # ---- Pfam anchor check on all candidates
    cf = W / "loci" / "candidates.faa"
    with open(cf, "w") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n{s}\n")
    dt = W / "loci" / "candidates.domtbl"
    qdt = W / "loci" / "queries.domtbl"
    run(["hmmscan", "--cpu", str(cfg["threads"]), "--domE", str(cfg["hmmscan_domE"]), "-E", "1e-3", "-o", "/dev/null",
         "--domtblout", str(dt), str(W / "db" / "diag.hmm"), str(cf)])
    run(["hmmscan", "--cpu", str(cfg["threads"]), "--domE", str(cfg["hmmscan_domE"]), "-o", "/dev/null",
         "--domtblout", str(qdt), str(W / "db" / "diag.hmm"), cfg["queries"]])
    doms = collections.defaultdict(set)
    for line in open(dt):
        if not line.startswith("#"):
            f = line.split()
            if float(f[12]) <= cfg["hmmscan_domE"]:
                doms[f[3]].add(f[0])
    qd = collections.defaultdict(set)
    for line in open(qdt):
        if not line.startswith("#"):
            f = line.split()
            if float(f[12]) <= cfg["hmmscan_domE"]:
                qd[f[3]].add(f[0])
    qcheck = [{"ID": q, "anchor_rule": cat[q]["Anchor"], "domains_found": ",".join(sorted(qd[q])),
               "anchor_ok": int(anchor_ok(cat[q]["Anchor"], qd[q]))} for q in cat]
    write_tsv(R / "query_anchor_check.tsv", qcheck)
    bad = [x["ID"] for x in qcheck if not x["anchor_ok"]]
    if bad:
        log("WARNING: queries failing their own anchor rule:", bad)
    for c in cands:
        c["pfam_domains"] = ",".join(sorted(doms[c["cand_id"]]))
        c["anchor_ok"] = int(anchor_ok(cat[c["best_query"]]["Anchor"], doms[c["cand_id"]]))
    write_tsv(R / "candidate_loci.tsv", cands)
    write_tsv(R / "S3_homologue_loci_rescue.tsv", hom_rows)
    write_tsv(R / "best_alignment_per_query_subgenome.tsv", best_rows)

    # ---- OrthoFinder input
    od = W / "of_in"
    if od.exists():
        shutil.rmtree(od)
    od.mkdir()
    per = collections.defaultdict(list)
    for c in cands:
        if c["anchor_ok"] and c["subgenome"] in ("A", "B", "S", "D"):
            per[f"{c['asm_id']}__{c['subgenome']}"].append(c["cand_id"])
    for sp, ids in per.items():
        with open(od / f"{sp}.fa", "w") as fh:
            for i in ids:
                fh.write(f">{i}\n{seqs[i]}\n")
    shutil.copy(cfg["queries"], od / "QUERY.fa")
    small = {sp: len(ids) for sp, ids in per.items() if len(ids) < 30}
    if small:
        log("WARNING: species files with <30 candidates; OrthoFinder may ignore their hits "
            "('Too few hits ... to normalise'). Check:", small)
    expected = {f"{r['asm_id']}__{g}" for r in panel_rows(cfg) for g in r["genomes"]}
    empty = sorted(expected - set(per))
    if empty:
        log("NOTE: subgenomes with zero candidates (counted as absent):", empty)
    log(f"loci done: {len(cands)} candidates, {sum(c['anchor_ok'] for c in cands)} pass anchor; "
        f"{len(per)} OrthoFinder species files + QUERY")


# ----------------------------------------------------------------------------- orthofinder
def step_orthofinder(cfg, args):
    which_or_die("orthofinder")
    W = Path(cfg["workdir"])
    o = W / "of_out"
    if o.exists():
        shutil.rmtree(o)
    run(["orthofinder", "-f", str(W / "of_in"), "-o", str(o), "-S", cfg["orthofinder"]["search"],
         "-I", str(cfg["orthofinder"]["inflation"]), "-t", str(cfg["threads"]), "-a", str(cfg["threads"]), "-og"],
        env={k: v for k, v in os.environ.items() if k != "PYTHONSAFEPATH"})


def find_og_files(W):
    res = sorted((W / "of_out").glob("Results_*"))
    if not res:
        sys.exit("OrthoFinder results not found - run the orthofinder step")
    d = res[-1] / "Orthogroups"
    return d / "Orthogroups.tsv", d / "Orthogroups_UnassignedGenes.tsv"


# ----------------------------------------------------------------------------- stats
def step_stats(cfg, args):
    from scipy import stats as ss
    import numpy as np
    W, R = Path(cfg["workdir"]), Path(cfg["resultsdir"])
    cat = catalogue(cfg)
    panel = panel_rows(cfg)
    subs = [f"{r['asm_id']}__{g}" for r in panel for g in r["genomes"]]     # expected haploid subgenomes
    tier_of = {f"{r['asm_id']}__{g}": ("BS" if g in "BS" else g) for r in panel for g in r["genomes"]}
    ploidy_of = {r["asm_id"]: int(r["ploidy"]) for r in panel}
    cands = {c["cand_id"]: c for c in read_tsv(R / "candidate_loci.tsv")}
    ogf, unf = find_og_files(W)

    # orthogroup membership
    og_members = {}
    for fpath in (ogf, unf):
        with open(fpath) as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for line in f:
                p = line.rstrip("\n").split("\t")
                genes = []
                for col in p[1:]:
                    genes += [g.strip() for g in col.split(",") if g.strip()]
                og_members[p[0]] = genes
    q2og = {}
    for og, gs in og_members.items():
        for g in gs:
            if g in cat:
                q2og[g] = og
    missing = [q for q in cat if q not in q2og]
    if missing:
        log("WARNING queries not found in OrthoFinder output:", missing)

    comp = cfg["compartments"]

    def compartment(x):
        return ("core" if x >= comp["core"] else "soft-core" if x >= comp["softcore"]
                else "shell" if x >= comp["shell"] else "cloud")

    s9_by_og = collections.defaultdict(list)
    if (R / "S9_paralogue_pruning.tsv").exists():
        for r9 in read_tsv(R / "S9_paralogue_pruning.tsv"):
            s9_by_og[r9["orthogroup"]].append(r9)
    best_by_q = collections.defaultdict(list)
    if (R / "best_alignment_per_query_subgenome.tsv").exists():
        for rb in read_tsv(R / "best_alignment_per_query_subgenome.tsv"):
            best_by_q[rb["query"]].append(rb)

    og_rows, presence = [], {}
    for og in sorted(set(q2og.values())):
        qs = [q for q in q2og if q2og[q] == og]
        main = [q for q in qs if cat[q]["Set"] == "Main"]
        use = main if main else qs
        classes = sorted({cat[q]["Class"] for q in use})
        groups = sorted({cat[q]["Group"] for q in use})
        grp_label = groups[0] if len(groups) == 1 else "mixed"
        mem = [cands[g] for g in og_members[og] if g in cands]
        un_mem = [c for c in cands.values() if c["anchor_ok"] == "1" and c["subgenome"] == "Un" and c["best_query"] in qs]
        intact_sub = collections.Counter(f"{m['asm_id']}__{m['subgenome']}" for m in mem if m["intact"] == "1")
        retained_sub = collections.Counter()

        if s9_by_og:
            if grp_label == "NLR":
                pruned_mem = [
                    m for m in s9_by_og[og]
                    if m["verdict"] == "retained" and float(m["qcov"]) >= 0.80
                    and ((m["intact"] == "1" and float(m["identity"]) >= 0.85)
                         or m["reverse_hit_rank"] in ("1", "2", "3", "4"))
                ]
                sub_present = {f"{m['asm_id']}__{m['subgenome']}" for m in pruned_mem if m["subgenome"] in ("A", "B", "S", "D")}
                for m in pruned_mem:
                    if m["subgenome"] in ("A", "B", "S", "D"):
                        retained_sub[f"{m['asm_id']}__{m['subgenome']}"] += 1
            else:
                pruned_mem = [m for m in s9_by_og[og] if m["verdict"] == "retained"]
                dom_hits = [
                    b for q in use for b in best_by_q.get(q, [])
                    if float(b["qcov"]) >= 0.50 and float(b["identity"]) >= 0.60 and b["subgenome"] in ("A", "B", "S", "D")
                ]
                sub_present = (
                    {f"{m['asm_id']}__{m['subgenome']}" for m in pruned_mem if m["subgenome"] in ("A", "B", "S", "D")}
                    | {f"{b['asm_id']}__{b['subgenome']}" for b in dom_hits}
                )
                for m in pruned_mem:
                    if m["subgenome"] in ("A", "B", "S", "D"):
                        retained_sub[f"{m['asm_id']}__{m['subgenome']}"] += 1
                for b in dom_hits:
                    k_sub = f"{b['asm_id']}__{b['subgenome']}"
                    if retained_sub[k_sub] == 0:
                        retained_sub[k_sub] = 1
            pres = {s: int(s in sub_present) for s in subs}
            asm_present = {s.split("__")[0] for s in sub_present}
        else:
            per_sub = collections.Counter(f"{m['asm_id']}__{m['subgenome']}" for m in mem)
            pres = {s: int(per_sub.get(s, 0) > 0) for s in subs}
            asm_present = {m["asm_id"] for m in mem}
            retained_sub = intact_sub

        presence[og] = pres
        row = {"orthogroup": og, "queries": ",".join(qs), "set": "Main" if main else "Comparator",
               "class": "/".join(classes), "group": grp_label,
               "n_members": len(mem), "n_members_Un": len(un_mem)}
        row["occupancy_all_subgenomes"] = round(sum(pres.values()) / len(subs), 4)
        for t in ("A", "BS", "D"):
            ts = [s for s in subs if tier_of[s] == t]
            row[f"occupancy_{t}"] = round(sum(pres[s] for s in ts) / len(ts), 4) if ts else ""
        for lab, pl in (("AB_tetraploid", 4), ("ABD_hexaploid", 6)):
            asms = [r["asm_id"] for r in panel if int(r["ploidy"]) == pl]
            row[f"occupancy_{lab}"] = round(sum(a in asm_present for a in asms) / len(asms), 4) if asms else ""
        row["occupancy_assemblies"] = round(len(asm_present) / len(panel), 4)
        row["compartment"] = compartment(row["occupancy_all_subgenomes"])
        row["compartment_assembly_level"] = compartment(row["occupancy_assemblies"])
        cn = [intact_sub.get(s, 0) for s in subs]
        row["intact_copies_mean"] = round(float(np.mean(cn)), 3)
        row["intact_copies_sd"] = round(float(np.std(cn, ddof=1)), 3)
        row["intact_copies_max"] = max(cn)
        rcn = [retained_sub.get(s, 0) for s in subs]
        row["n_retained_orthologues"] = sum(rcn)
        row["retained_copies_mean"] = round(float(np.mean(rcn)), 3)
        row["retained_copies_sd"] = round(float(np.std(rcn, ddof=1)), 3)
        row["retained_copies_max"] = max(rcn)
        og_rows.append(row)
    write_tsv(R / "S1b_orthogroup_partition.tsv", og_rows)
    with open(R / "og_presence_matrix.tsv", "w") as fh:
        ogs = [r["orthogroup"] for r in og_rows]
        fh.write("subgenome\t" + "\t".join(ogs) + "\n")
        for s in subs:
            fh.write(s + "\t" + "\t".join(str(presence[o][s]) for o in ogs) + "\n")

    V = []  # summary metric values

    def put(key, val, note=""):
        V.append({"metric": key, "value": val, "note": note})

    main_q = [q for q in cat if cat[q]["Set"] == "Main"]
    put("N curated queries (Main set)", len(main_q))
    for k in ("I", "II", "III", "IV", "V"):
        put(f"N Class {k} queries", sum(cat[q]["Class"] == k for q in main_q))
    put("N loci (LocusGroup) in Main set", len({cat[q]["LocusGroup"] for q in main_q}))
    put("N assemblies", len(panel))
    put("N haploid subgenomes", len(subs))
    put("N assemblies with GFF3 (rescue test)", sum(r["annotated"] for r in panel))
    M = [r for r in og_rows if r["set"] == "Main" and r["group"] in ("NLR", "nonNLR")]
    put("N orthogroups containing Main queries", sum(r["set"] == "Main" for r in og_rows),
        f"{sum(r['group']=='mixed' for r in og_rows)} mixed-class OGs excluded from tests")

    for col, tag in (("compartment", "all haploid subgenomes"), ("compartment_assembly_level", "assembly level")):
        for grp in ("NLR", "nonNLR"):
            g = [r for r in M if r["group"] == grp]
            for cpt in ("core", "soft-core", "shell", "cloud"):
                n = sum(r[col] == cpt for r in g)
                put(f"{grp} OGs in {cpt} [{tag}]", f"{n}/{len(g)} ({100*n/len(g):.1f}%)" if g else "NA")
        for k_cls in ("III", "IV", "V"):
            g_cls = [r for r in M if r["class"] == k_cls]
            n_cs = sum(r[col] in ("core", "soft-core") for r in g_cls)
            put(f"Class {k_cls} OGs in core+soft-core [{tag}]",
                f"{n_cs}/{len(g_cls)} ({100*n_cs/len(g_cls):.1f}%)" if g_cls else "NA")
        a_ = sum(r[col] in ("core", "soft-core") for r in M if r["group"] == "NLR")
        b_ = sum(r[col] in ("shell", "cloud") for r in M if r["group"] == "NLR")
        c_ = sum(r[col] in ("core", "soft-core") for r in M if r["group"] == "nonNLR")
        d_ = sum(r[col] in ("shell", "cloud") for r in M if r["group"] == "nonNLR")
        orr, p = ss.fisher_exact([[b_, a_], [d_, c_]])
        or_ha = ((b_ + 0.5) * (c_ + 0.5)) / ((a_ + 0.5) * (d_ + 0.5))
        or_str = f"OR={orr:.3g} (Haldane-Anscombe OR={or_ha:.1f})" if math.isinf(orr) else f"OR={orr:.3g}"
        put(f"Fisher exact NLR shell/cloud vs nonNLR [{tag}]", f"{or_str}, P={p:.3g}",
            f"table [[NLR shell+cloud={b_}, NLR core+soft={a_}],[nonNLR shell+cloud={d_}, nonNLR core+soft={c_}]]; "
            "OR>1 means NLRs are more often shell/cloud")
    occ = "occupancy_all_subgenomes"
    x = [r[occ] for r in M if r["group"] == "NLR"]
    y = [r[occ] for r in M if r["group"] == "nonNLR"]
    if x and y:
        u, p = ss.mannwhitneyu(x, y, alternative="two-sided")
        put("Mann-Whitney U occupancy NLR vs nonNLR (all subgenomes)", f"U={u:.1f}, P={p:.3g}",
            f"median NLR={np.median(x):.3f} (n={len(x)}), median nonNLR={np.median(y):.3f} (n={len(y)})")
    x2 = [r["occupancy_assemblies"] for r in M if r["group"] == "NLR"]
    y2 = [r["occupancy_assemblies"] for r in M if r["group"] == "nonNLR"]
    if x2 and y2:
        u, p = ss.mannwhitneyu(x2, y2, alternative="two-sided")
        put("Mann-Whitney U occupancy NLR vs nonNLR (assembly level)", f"U={u:.1f}, P={p:.3g}",
            f"median NLR={np.median(x2):.3f}, median nonNLR={np.median(y2):.3f}")
    for grp in ("NLR", "nonNLR"):
        cells = []
        for r in M:
            if r["group"] != grp:
                continue
            mem = [cands[g] for g in og_members[r["orthogroup"]] if g in cands and cands[g]["intact"] == "1"]
            cnt = collections.Counter(f"{m['asm_id']}__{m['subgenome']}" for m in mem)
            cells += [cnt.get(s, 0) for s in subs]
        if cells:
            put(f"Intact copies per haploid subgenome, {grp} (mean ± s.d.)",
                f"{np.mean(cells):.2f} ± {np.std(cells, ddof=1):.2f}",
                f"n = {len(cells)} OG x subgenome cells (max = {max(cells)})")

    # ---- rescue / recovery (annotated assemblies only)
    H = [h for h in read_tsv(R / "S3_homologue_loci_rescue.tsv") if h["annotated"] == "1"
         and h["LocusGroup"] in {cat[q]["LocusGroup"] for q in main_q}]
    resc = ("rescued_unannotated", "rescued_split_model", "rescued_misannotated_model")
    annot_ok = ("annotated", "annotated_only_not_in_miniprot")
    per_asm = collections.defaultdict(lambda: [0, 0])
    for h in H:
        if h["status"] in annot_ok:
            per_asm[h["asm_id"]][0] += 1
        per_asm[h["asm_id"]][1] += 1
    asms = [r["asm_id"] for r in panel if r["annotated"]]
    ann = np.array([per_asm[a][0] for a in asms], float)
    uni = np.array([per_asm[a][1] for a in asms], float)
    write_tsv(R / "S3a_recovery_per_assembly.tsv",
              [{"asm_id": a, "annotation_homologues": int(per_asm[a][0]), "union_homologues": int(per_asm[a][1]),
                "pct_increase": round(100 * (per_asm[a][1] - per_asm[a][0]) / per_asm[a][0], 2) if per_asm[a][0] else ""}
               for a in asms])
    n_res = collections.Counter(h["status"] for h in H)
    put("N rescued homologue loci (all annotated assemblies)", sum(n_res[k] for k in resc),
        "; ".join(f"{k}={n_res[k]}" for k in resc))
    pct = 100 * (uni - ann) / np.where(ann > 0, ann, np.nan)
    pct = pct[~np.isnan(pct)]
    if len(pct):
        put("Recovery increase per assembly (mean ± s.e.m.)",
            f"{pct.mean():.2f} ± {pct.std(ddof=1)/math.sqrt(len(pct)) if len(pct)>1 else 0:.2f} %", f"n = {len(pct)} assemblies")
    if len(asms) > 1 and np.any(uni != ann):
        w = ss.wilcoxon(uni, ann, alternative="two-sided")
        put("Wilcoxon signed-rank union vs annotation", f"W={w.statistic:.1f}, P={w.pvalue:.3g}",
            "union >= annotation by construction; test only asks whether the gain is non-zero")
    else:
        put("Wilcoxon signed-rank union vs annotation", "not computable", "no difference in any assembly")
    # per class: per LocusGroup increase over all annotated assemblies
    lg_class = {}
    for q in main_q:
        lg_class.setdefault(cat[q]["LocusGroup"], cat[q]["Class"])
    lgrows = []
    for lg, k in sorted(lg_class.items()):
        hh = [h for h in H if h["LocusGroup"] == lg]
        a0 = sum(h["status"] in annot_ok for h in hh)
        u0 = len(hh)
        lgrows.append({"LocusGroup": lg, "Class": k, "annotation_homologues": a0, "union_homologues": u0,
                       **{k_: sum(h["status"] == k_ for h in hh) for k_ in resc},
                       "pct_increase": round(100 * (u0 - a0) / a0, 2) if a0 else ("inf" if u0 else "")})
    write_tsv(R / "S3b_recovery_per_locus.tsv", lgrows)
    for k in ("I", "II", "III", "IV", "V"):
        v = [r["pct_increase"] for r in lgrows if r["Class"] == k and isinstance(r["pct_increase"], float)]
        nz = sum(1 for r in lgrows if r["Class"] == k and r["pct_increase"] == "inf")
        nres = sum(r[k_] for r in lgrows if r["Class"] == k for k_ in resc)
        put(f"Class {k} rescued homologue loci", nres, f"{sum(1 for r in lgrows if r['Class']==k)} loci in class")
        if v:
            put(f"Class {k} per-locus recovery increase (mean ± s.e.m.)",
                f"{np.mean(v):.2f} ± {(np.std(v, ddof=1)/math.sqrt(len(v)) if len(v)>1 else 0):.2f} %",
                f"n = {len(v)} loci with >=1 annotated homologue; {nz} loci found only after rescue")
    # NLR-Annotator support
    if cfg["nlr_annotator"]["enabled"]:
        nl_rows, sup = [], 0
        for r in panel:
            beds = list((W / "miniprot" / r["asm_id"]).glob("*.nlr.bed"))
            n = 0
            for b in beds:
                n += sum(1 for l in open(b) if l.strip() and not l.startswith(("#", "track", "browser")))
            nl_rows.append({"asm_id": r["asm_id"], "nlr_annotator_loci": n})
        write_tsv(R / "nlr_annotator_loci_per_assembly.tsv", nl_rows)
        hb = [h for h in H if h["status"] in resc and cat.get(h["best_query"], {}).get("Group") == "NLR"]
        tmp = W / "loci" / "rescued_nlr.bed"
        with open(tmp, "w") as fh:
            for h in hb:
                fh.write(f"{h['seqid']}\t{int(h['start'])-1}\t{h['end']}\t{h['hom_locus']}\t0\t.\n")
        allnlr = W / "loci" / "all_nlr_annotator.bed"
        with open(allnlr, "w") as fh:
            for r in panel:
                for b in (W / "miniprot" / r["asm_id"]).glob("*.nlr.bed"):
                    for line in open(b):
                        if line.startswith(("#", "track", "browser")) or not line.strip():
                            continue
                        f = line.rstrip("\n").split("\t")
                        fh.write("\t".join(f[:3] + [f"{r['asm_id']}|{f[3] if len(f)>3 else ''}", "0", "."]) + "\n")
        # same seqid names can occur in different assemblies -> compare per assembly via name prefix
        out = subprocess.run(["bedtools", "intersect", "-a", str(tmp), "-b", str(allnlr), "-wa", "-wb"],
                             capture_output=True, text=True).stdout
        supported = set()
        for line in out.splitlines():
            f = line.split("\t")
            if f[3].split("|")[0] == f[9].split("|")[0]:
                supported.add(f[3])
        put("Rescued NLR homologue loci supported by NLR-Annotator", f"{len(supported)}/{len(hb)}")

    # ---- identity-threshold sensitivity across the 30 Main orthogroups (pruned S9 orthologue set)
    sweep = []
    for t in cfg["sweep_identity"]:
        t_prog = round(max(0.50, t - 0.25), 2)
        t_s9_non = round(max(0.50, t - 0.35), 2)
        x_sw, y_sw = [], []
        for r in M:
            og = r["orthogroup"]
            qs = [q.strip() for q in r["queries"].split(",") if q.strip()]
            main_qs = [q for q in qs if q in cat and cat[q]["Set"] == "Main"]
            use_qs = main_qs if main_qs else qs
            if r["group"] == "NLR":
                pm = [
                    m for m in s9_by_og[og]
                    if m["verdict"] == "retained" and float(m["qcov"]) >= cfg["sweep_qcov"]
                    and ((m["intact"] == "1" and float(m["identity"]) >= t)
                         or (m["reverse_hit_rank"] in ("1", "2", "3", "4") and float(m["identity"]) >= round(t - 0.25, 2)))
                ]
                sub_pres = {f"{m['asm_id']}__{m['subgenome']}" for m in pm if m["subgenome"] in ("A", "B", "S", "D")}
                x_sw.append(len(sub_pres) / len(subs))
            else:
                pm = [m for m in s9_by_og[og] if m["verdict"] == "retained" and float(m["identity"]) >= t_s9_non]
                dom = [
                    b for q in use_qs for b in best_by_q.get(q, [])
                    if float(b["qcov"]) >= 0.50 and float(b["identity"]) >= t_prog and b["subgenome"] in ("A", "B", "S", "D")
                ]
                sub_pres = (
                    {f"{m['asm_id']}__{m['subgenome']}" for m in pm if m["subgenome"] in ("A", "B", "S", "D")}
                    | {f"{b['asm_id']}__{b['subgenome']}" for b in dom}
                )
                y_sw.append(len(sub_pres) / len(subs))
        u, p = ss.mannwhitneyu(x_sw, y_sw, alternative="two-sided") if x_sw and y_sw else (float("nan"), float("nan"))
        sweep.append({
            "min_identity_NLR": t,
            "min_identity_progenitor": t_prog,
            "min_qcov": cfg["sweep_qcov"],
            "median_occ_NLR": round(float(np.median(x_sw)), 4),
            "median_occ_nonNLR": round(float(np.median(y_sw)), 4),
            "mean_occ_NLR": round(float(np.mean(x_sw)), 4),
            "mean_occ_nonNLR": round(float(np.mean(y_sw)), 4),
            "MWU_U": round(float(u), 2),
            "MWU_P": float(f"{p:.4e}"),
            "n_NLR_OGs": len(x_sw),
            "n_nonNLR_OGs": len(y_sw),
        })
    write_tsv(R / "S3c_identity_threshold_sensitivity.tsv", sweep)

    write_tsv(R / "section3.1_values.tsv", V, ["metric", "value", "note"])
    log("stats done -> results/section3.1_values.tsv")
    for v in V:
        print(f"  {v['metric']}: {v['value']}   {v['note']}")


# ----------------------------------------------------------------------------- main
STEPS = {"setup": step_setup, "prep": step_prep, "genome": step_genome, "annot": step_annot,
         "loci": step_loci, "orthofinder": step_orthofinder, "stats": step_stats}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=list(STEPS) + ["all"])
    ap.add_argument("--config", default=str(HERE / "config_3.1.json"))
    ap.add_argument("--asm", default=None, help="run the genome step for one assembly only")
    args = ap.parse_args()
    cfg = load_cfg(args.config)
    for s in (list(STEPS) if args.step == "all" else [args.step]):
        log("=== step", s)
        STEPS[s](cfg, args)
