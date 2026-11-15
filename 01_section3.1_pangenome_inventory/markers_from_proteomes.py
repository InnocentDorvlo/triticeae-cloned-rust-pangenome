#!/usr/bin/env python3
"""
markers_from_proteomes.py - Mac-feasible route to the Fig. 1b tree markers.

extracts the same markers from the proteomes already on disk instead: BUSCO in protein mode is
hmmsearch only, no genome, no download.

  for each of the 16 annotated assemblies:
      busco -m proteins -l poales_odb10 on work/annot/<asm>/longest.faa
      keep complete single-copy markers
      assign each marker's gene to a subgenome via gene_loc.tsv + seqmap.tsv
      write work/tree/markers/<asm>__<subgenome>/<busco_id>.faa

Output follows the canonical per-subgenome marker directory layout so build_species_tree.py runs
afterwards unchanged.

Cost and limits
  - 16 BUSCO runs, roughly 10-40 min each on 4 cores: a long evening, not weeks.
  - Only annotated assemblies have proteomes, so the tree has ~42 tips (the subgenomes of 16
    assemblies) rather than 56. The 11 unannotated assemblies cannot appear; say so in the legend.
  - Markers come from the published annotation, so a subgenome's placement depends on that
    annotation's completeness. The genome-based scan does not have that dependency.

Usage: python markers_from_proteomes.py [--threads 4] [--lineage poales_odb10] [--max-markers 0]
       (needs busco in PATH: conda install -n rust31 -c conda-forge -c bioconda busco)
"""
import argparse, collections, csv, os, shutil, subprocess, sys
from pathlib import Path

STAGE_DIR = Path(__file__).resolve().parent
HERE = STAGE_DIR.parent
sys.path.insert(0, str(STAGE_DIR))
from rust31 import log, load_cfg, read_tsv, read_fasta, panel_rows


def busco_table(d):
    """locate full_table.tsv from a BUSCO protein-mode run directory"""
    for pat in ("run_*/full_table.tsv", "*/run_*/full_table.tsv", "full_table.tsv"):
        hits = list(Path(d).glob(pat))
        if hits:
            return hits[0]
    return None


def subset_mode(args):
    """hmmsearch N lineage HMMs directly against each proteome - BUSCO protein mode in miniature.

    BUSCO runs all 4,896 poales_odb10 models, which costs ~2.3 h per hexaploid proteome. A species
    tree does not need 4,896 markers; a few hundred shared, single-copy ones give the same topology
    at a fraction of the cost. Single-copy here means exactly one protein scoring above the
    lineage's own score cutoff for that model, which is the criterion BUSCO applies."""
    import random
    if shutil.which("hmmsearch") is None:
        sys.exit("hmmsearch not in PATH - conda install -n rust31 -c conda-forge -c bioconda hmmer")
    cfg = load_cfg(args.config)
    W = Path(cfg["workdir"])
    lin = Path(args.lineage_dir) if args.lineage_dir else None
    if lin is None:
        for c in (HERE / "busco_downloads" / "lineages" / args.lineage,
                  Path.cwd() / "busco_downloads" / "lineages" / args.lineage):
            if c.is_dir():
                lin = c; break
    if lin is None or not (lin / "hmms").is_dir():
        sys.exit(f"lineage HMMs not found - pass --lineage-dir .../busco_downloads/lineages/{args.lineage}")
    cut = {}
    for line in open(lin / "scores_cutoff"):
        p = line.split()
        if len(p) >= 2:
            cut[p[0]] = float(p[1])
    ids = sorted(p.stem for p in (lin / "hmms").glob("*.hmm") if p.stem in cut)
    random.Random(20261030).shuffle(ids)
    ids = sorted(ids[:args.subset])
    sub = W / "tree" / f"subset_{len(ids)}.hmm"
    sub.parent.mkdir(parents=True, exist_ok=True)
    if not sub.exists():
        with open(sub, "w") as fh:
            for i in ids:
                fh.write((lin / "hmms" / f"{i}.hmm").read_text())
    log(f"subset mode: {len(ids)} of {len(cut)} lineage HMMs -> {sub}")

    md = W / "tree" / "markers"
    md.mkdir(parents=True, exist_ok=True)
    total = collections.Counter()
    for r in panel_rows(cfg):
        a = r["asm_id"]
        pep, gl = W / "annot" / a / "longest.faa", W / "annot" / a / "gene_loc.tsv"
        if not (pep.exists() and gl.exists()):
            log(f"  {a}: no proteome - skipped (unannotated assembly)")
            continue
        tbl = W / "tree" / a / f"subset_{len(ids)}.domtbl"
        tbl.parent.mkdir(parents=True, exist_ok=True)
        if not tbl.exists() or tbl.stat().st_size == 0:
            log(f"  {a}: hmmsearch {len(ids)} models")
            subprocess.run(["hmmsearch", "--cpu", str(args.threads), "--domtblout", str(tbl),
                            "-o", os.devnull, str(sub), str(pep)], check=True)
        best = collections.defaultdict(list)          # busco id -> [(score, gene)]
        for line in open(tbl):
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 14:
                continue
            gene, model, score = f[0], f[3], float(f[7])
            if score >= cut.get(model, 0):
                best[model].append((score, gene))
        sub_of = {x["seqid"]: x["subgenome"] for x in read_tsv(W / "miniprot" / a / "seqmap.tsv")}
        gene_sub = {g["gene"]: sub_of.get(g["seqid"], "Un") for g in read_tsv(gl)}
        prot = {n: s for n, _, s in read_fasta(pep)}
        kept = 0
        for model, hits in best.items():
            # One copy PER SUBGENOME, not per assembly. A BUSCO marker in a hexaploid normally has
            # three homoeologous copies, so the whole-proteome "single copy" test would discard
            # almost every marker in a polyploid; grouping by subgenome is the criterion that makes
            # sense for subgenome-resolution tips.
            by_sg = collections.defaultdict(set)
            for _, g in hits:
                sg = gene_sub.get(g)
                if sg not in (None, "Un"):
                    by_sg[sg].add(g)
            for sg, genes in by_sg.items():
                if len(genes) != 1:
                    continue
                gene = genes.pop()
                if gene not in prot:
                    continue
                d = md / f"{a}__{sg}"
                d.mkdir(parents=True, exist_ok=True)
                (d / f"{model}.faa").write_text(f">{a}__{sg}\n{prot[gene]}\n")
                kept += 1
                total[f"{a}__{sg}"] += 1
        log(f"  {a}: {kept} single-copy markers assigned")
    if not total:
        sys.exit("no markers written")
    log(f"{len(total)} subgenome tips, {sum(total.values())} marker sequences -> {md}")
    for t in sorted(total):
        print(f"  {t}: {total[t]}")
    print("\nnext: python build_species_tree.py --threads 4 --min-occupancy 0.9")
    return


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--lineage", default="poales_odb10")
    ap.add_argument("--config", default=str(HERE / "config_3.1.json"))
    ap.add_argument("--max-markers", type=int, default=0,
                    help="0 = keep every complete single-copy marker; set e.g. 500 to cap the matrix")
    ap.add_argument("--subset", type=int, default=0,
                    help="FAST PATH: skip BUSCO and hmmsearch only N of the 4,896 lineage HMMs "
                         "(e.g. 500). Runtime falls roughly N/4896-fold: ~3 h for all 16 assemblies "
                         "instead of ~30 h, and a few hundred markers is ample for a species tree.")
    ap.add_argument("--lineage-dir", default="",
                    help="path to busco_downloads/lineages/<lineage>; found automatically if BUSCO "
                         "has already downloaded it")
    args = ap.parse_args()
    if args.subset:
        return subset_mode(args)
    if shutil.which("busco") is None:
        sys.exit("busco not in PATH - conda install -n rust31 -c conda-forge -c bioconda busco")
    cfg = load_cfg(args.config)
    W = Path(cfg["workdir"])
    md = W / "tree" / "markers"
    md.mkdir(parents=True, exist_ok=True)

    total = collections.Counter()
    for r in panel_rows(cfg):
        a = r["asm_id"]
        pep = W / "annot" / a / "longest.faa"
        gl = W / "annot" / a / "gene_loc.tsv"
        if not (pep.exists() and gl.exists()):
            log(f"  {a}: no proteome - skipped (unannotated assembly)")
            continue
        out = W / "tree" / a / "busco_prot"
        if not busco_table(out):
            out.parent.mkdir(parents=True, exist_ok=True)
            log(f"  {a}: BUSCO protein mode")
            subprocess.run(["busco", "-i", str(pep), "-m", "proteins", "-l", args.lineage,
                            "-o", "busco_prot", "--out_path", str(out.parent),
                            "-c", str(args.threads), "-f", "--quiet"], check=True)
        ft = busco_table(out)
        if ft is None:
            log(f"  {a}: BUSCO produced no full_table.tsv - skipped")
            continue

        sub = {x["seqid"]: x["subgenome"] for x in read_tsv(W / "miniprot" / a / "seqmap.tsv")}
        gene_sub = {g["gene"]: sub.get(g["seqid"], "Un") for g in read_tsv(gl)}
        prot = {n: s for n, _, s in read_fasta(pep)}

        # Complete AND Duplicated: in a polyploid the homoeologues of one marker are reported as
        # Duplicated, and they are exactly the copies the subgenome tips need. The single-copy test
        # is then applied within each subgenome.
        hits = collections.defaultdict(set)            # (busco id, subgenome) -> genes
        for line in open(ft):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[1] not in ("Complete", "Duplicated"):
                continue
            bid, gene = f[0], f[2]
            sg = gene_sub.get(gene)
            if sg in (None, "Un") or gene not in prot:
                continue
            hits[(bid, sg)].add(gene)
        kept = 0
        for (bid, sg), genes in hits.items():
            if len(genes) != 1:
                continue
            gene = genes.pop()
            d = md / f"{a}__{sg}"
            d.mkdir(parents=True, exist_ok=True)
            (d / f"{bid}.faa").write_text(f">{a}__{sg}\n{prot[gene]}\n")
            kept += 1
            total[f"{a}__{sg}"] += 1
            if args.max_markers and kept >= args.max_markers:
                break
        log(f"  {a}: {kept} single-copy markers assigned")

    if not total:
        sys.exit("no markers written - check that BUSCO ran and that gene_loc.tsv seqids match seqmap.tsv")
    log(f"{len(total)} subgenome tips, {sum(total.values())} marker sequences -> {md}")
    for t in sorted(total):
        print(f"  {t}: {total[t]}")
    print("\nnext: python build_species_tree.py --threads 4 --min-occupancy 0.9")


if __name__ == "__main__":
    main()
