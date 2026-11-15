#!/usr/bin/env python3
"""
build_species_tree.py - builds the subgenome-resolved species tree for Fig. 1b.

Reads single-copy BUSCO marker protein sequences from:
  work/tree/markers/<asm>__<subgenome>/<busco_id>.faa
(written by markers_from_proteomes.py from the annotated proteomes).

Pipeline:
  1. Keep markers present across >= --min-occupancy of subgenome tips (default 0.90).
  2. Align each marker with MAFFT (--auto --quiet) and trim with trimAl (-automated1),
     caching trimmed alignments in work/tree/aln/<marker>.trim.faa.
  3. Concatenate into work/tree/supermatrix.faa (+ work/tree/partitions.nex).
  4. Infer maximum-likelihood phylogeny with IQ-TREE (1,000 UFBoot + 1,000 SH-aLRT replicates).
  5. Write results/species_tree.treefile, results/S2b_species_tree_tips.tsv, and
     results/section3.1_tree_values.tsv.

Usage:
  python build_species_tree.py [--threads 4] [--min-occupancy 0.90] [--max-markers 250] [--mfp]
"""
import argparse, collections, concurrent.futures, os, re, shutil, subprocess, sys
from pathlib import Path

STAGE_DIR = Path(__file__).resolve().parent
HERE = STAGE_DIR.parent
sys.path.insert(0, str(STAGE_DIR))
from rust31 import log, load_cfg, read_fasta, write_tsv, panel_rows


def find_iqtree():
    for exe in ("iqtree2", "iqtree3", "iqtree"):
        p = shutil.which(exe)
        if p:
            return p
    return None


def clean_seq(s):
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWYBZX*-]", "X", s.upper()).replace("*", "X")


def py_trim_gappy(aln_dict, max_gap_frac=0.5):
    """Fallback column trimmer if trimAl removes everything or is unavailable."""
    if not aln_dict:
        return aln_dict
    seqs = list(aln_dict.values())
    L = len(seqs[0])
    keep = [i for i in range(L) if sum(s[i] in "-X?" for s in seqs) / len(seqs) <= max_gap_frac]
    if not keep:
        keep = list(range(L))
    return {k: "".join(v[i] for i in keep) for k, v in aln_dict.items()}


def align_and_trim_one(task):
    marker, tip_seqs, aln_dir, has_trimal = task
    out_trim = aln_dir / f"{marker}.trim.faa"
    if out_trim.exists() and out_trim.stat().st_size > 0:
        return marker, {n: s for n, _, s in read_fasta(out_trim)}

    raw_fa = aln_dir / f"{marker}.raw.faa"
    aln_fa = aln_dir / f"{marker}.aln.faa"
    with open(raw_fa, "w") as fh:
        for tip, seq in sorted(tip_seqs.items()):
            fh.write(f">{tip}\n{clean_seq(seq)}\n")

    with open(aln_fa, "w") as fh:
        subprocess.run(["mafft", "--auto", "--quiet", "--thread", "1", str(raw_fa)],
                       stdout=fh, stderr=subprocess.DEVNULL, check=True)

    trimmed = {}
    if has_trimal:
        subprocess.run(["trimal", "-in", str(aln_fa), "-out", str(out_trim), "-automated1"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        if out_trim.exists() and out_trim.stat().st_size > 0:
            trimmed = {n: s for n, _, s in read_fasta(out_trim)}

    if not trimmed or len(next(iter(trimmed.values()), "")) < 20:
        aln = {n: s for n, _, s in read_fasta(aln_fa)}
        trimmed = py_trim_gappy(aln, max_gap_frac=0.5)
        with open(out_trim, "w") as fh:
            for k, v in sorted(trimmed.items()):
                fh.write(f">{k}\n{v}\n")

    raw_fa.unlink(missing_ok=True)
    aln_fa.unlink(missing_ok=True)
    return marker, trimmed


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(HERE / "config_3.1.json"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-occupancy", type=float, default=0.90,
                    help="minimum fraction of subgenome tips carrying a marker (default: 0.90)")
    ap.add_argument("--max-markers", type=int, default=250,
                    help="cap number of markers in the supermatrix (0 = use all passing --min-occupancy)")
    ap.add_argument("--model", default="LG+F+G4",
                    help="IQ-TREE amino-acid model (default: LG+F+G4; pass MFP or --mfp for full ModelFinder)")
    ap.add_argument("--mfp", action="store_true",
                    help="run full ModelFinder (-m MFP) instead of fixed LG+F+G4")
    args = ap.parse_args()

    if shutil.which("mafft") is None:
        sys.exit("mafft not found in PATH - conda install -n rust31 -c conda-forge -c bioconda mafft trimal iqtree")
    iqtree_bin = find_iqtree()
    if iqtree_bin is None:
        sys.exit("iqtree/iqtree2/iqtree3 not found in PATH - conda install -n rust31 -c conda-forge -c bioconda iqtree")
    has_trimal = shutil.which("trimal") is not None

    cfg = load_cfg(args.config)
    W, R = Path(cfg["workdir"]), Path(cfg["resultsdir"])
    md = W / "tree" / "markers"
    if not md.is_dir():
        sys.exit(f"{md} not found - run markers_from_proteomes.py --subset 500 first")

    tip_dirs = sorted(d for d in md.iterdir() if d.is_dir() and not d.name.startswith("."))
    if not tip_dirs:
        sys.exit(f"no subgenome directories inside {md}")
    tips = [d.name for d in tip_dirs]

    # If subset_500.hmm exists, prefer markers belonging to that consistent 500-model set
    subset_models = set()
    for sf in sorted((W / "tree").glob("subset_*.hmm")):
        for line in open(sf):
            if line.startswith("NAME "):
                subset_models.add(line.split()[1])

    marker_seqs = collections.defaultdict(dict)  # marker -> {tip: seq}
    for d in tip_dirs:
        tip = d.name
        for fa in d.glob("*.faa"):
            if subset_models and fa.stem not in subset_models:
                continue
            for _, _, seq in read_fasta(fa):
                if seq:
                    marker_seqs[fa.stem][tip] = seq
                    break

    min_tips = max(4, int(round(args.min_occupancy * len(tips))))
    passing = sorted(
        (m for m, d in marker_seqs.items() if len(d) >= min_tips),
        key=lambda m: (-len(marker_seqs[m]), m)
    )
    if len(passing) < 10:
        # fallback to 0.75 occupancy if strict cutoff left too few markers
        min_tips = max(4, int(round(0.75 * len(tips))))
        passing = sorted(
            (m for m, d in marker_seqs.items() if len(d) >= min_tips),
            key=lambda m: (-len(marker_seqs[m]), m)
        )
    if args.max_markers and len(passing) > args.max_markers:
        passing = passing[:args.max_markers]

    log(f"subgenome tips: {len(tips)}; markers passing occupancy (>={min_tips}/{len(tips)} tips): "
        f"{len(passing)} selected")

    aln_dir = W / "tree" / "aln"
    aln_dir.mkdir(parents=True, exist_ok=True)
    tasks = [(m, marker_seqs[m], aln_dir, has_trimal) for m in passing]

    aligned_markers = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
        futures = [ex.submit(align_and_trim_one, t) for t in tasks]
        for i, fut in enumerate(concurrent.futures.as_completed(futures), 1):
            m, trimmed = fut.result()
            if trimmed:
                aligned_markers[m] = trimmed
            if i % 25 == 0 or i == len(tasks):
                log(f"  aligned & trimmed {i}/{len(tasks)} markers")

    ordered_markers = [m for m in passing if m in aligned_markers]
    super_seqs = {t: [] for t in tips}
    partitions = []
    pos = 1
    for m in ordered_markers:
        d = aligned_markers[m]
        L = len(next(iter(d.values())))
        if L == 0:
            continue
        for t in tips:
            super_seqs[t].append(d.get(t, "-" * L))
        partitions.append((m, pos, pos + L - 1))
        pos += L

    total_cols = pos - 1
    super_fa = W / "tree" / "supermatrix.faa"
    part_nex = W / "tree" / "partitions.nex"
    with open(super_fa, "w") as fh:
        for t in tips:
            seq = "".join(super_seqs[t])
            fh.write(f">{t}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")

    with open(part_nex, "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for m, s0, e0 in partitions:
            fh.write(f"  charset {m} = {s0}-{e0};\n")
        fh.write("end;\n")

    log(f"supermatrix written: {len(tips)} tips x {total_cols} AA columns across {len(partitions)} markers")

    # Per-tip summary table
    panel_map = {r["asm_id"]: r for r in panel_rows(cfg)}
    tip_rows = []
    for t in tips:
        asm, sg = t.split("__", 1) if "__" in t else (t, "")
        pr = panel_map.get(asm, {})
        n_pres = sum(1 for m in ordered_markers if t in aligned_markers[m])
        seq_str = "".join(super_seqs[t])
        ungapped = sum(1 for c in seq_str if c not in "-X?")
        tip_rows.append({
            "tip": t,
            "asm_id": asm,
            "subgenome": sg,
            "species": pr.get("species", ""),
            "ploidy": pr.get("ploidy", ""),
            "markers_present": n_pres,
            "markers_total": len(ordered_markers),
            "marker_occupancy_pct": round(100.0 * n_pres / max(1, len(ordered_markers)), 2),
            "ungapped_aa_sites": ungapped,
            "supermatrix_aa_sites": total_cols
        })
    write_tsv(R / "S2b_species_tree_tips.tsv", tip_rows)

    # Run IQ-TREE
    prefix = W / "tree" / "species_tree"
    model = "MFP" if args.mfp else args.model
    cmd = [iqtree_bin, "-s", str(super_fa), "-m", model, "-B", "1000", "-alrt", "1000",
           "-T", str(args.threads), "--prefix", str(prefix), "-redo"]
    log("running IQ-TREE:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    tf = Path(f"{prefix}.treefile")
    if not tf.exists():
        sys.exit(f"IQ-TREE did not produce {tf}")
    shutil.copy(tf, R / "species_tree.treefile")
    if Path(f"{prefix}.contree").exists():
        shutil.copy(f"{prefix}.contree", R / "species_tree.contree")

    V = [
        {"metric": "N subgenome tips in species tree", "value": str(len(tips)),
         "note": f"from {len({t.split('__')[0] for t in tips})} annotated assemblies"},
        {"metric": "N single-copy BUSCO markers in supermatrix", "value": str(len(partitions)),
         "note": f">={args.min_occupancy*100:.0f}% tip occupancy"},
        {"metric": "Supermatrix alignment length (AA sites)", "value": str(total_cols),
         "note": f"MAFFT + {'trimAl -automated1' if has_trimal else 'gap-trimmed'}, IQ-TREE ({model}, 1000 UFBoot)"}
    ]
    write_tsv(R / "section3.1_tree_values.tsv", V, ["metric", "value", "note"])
    log(f"species tree written -> {R / 'species_tree.treefile'}")
    for v in V:
        print(f"  {v['metric']}: {v['value']}   ({v['note']})")


if __name__ == "__main__":
    main()
