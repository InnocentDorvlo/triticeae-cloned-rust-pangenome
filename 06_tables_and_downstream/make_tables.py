#!/usr/bin/env python3
"""
make_tables.py - consolidates all Section 3.1 results into publication-ready tables in tables/.

Called by `bash run_missing.sh figures`.

Outputs in tables/:
  Table1_Cloned_Rust_Resistance_Catalogue.tsv
  Supplementary_Table_S1_Queries_and_Domains.tsv
  Supplementary_Table_S1b_Orthogroup_Partition.tsv
  Supplementary_Table_S2_Assembly_Panel_QC.tsv
  Supplementary_Table_S3_Homologue_Loci_Rescue.tsv
  Supplementary_Table_S3a_Recovery_per_Assembly.tsv
  Supplementary_Table_S3b_Recovery_per_Locus.tsv
  Supplementary_Table_S3c_Identity_Sensitivity.tsv
  Supplementary_Table_S3d_Heaps_Law.tsv
  Supplementary_Table_S9_Paralogue_Pruning.tsv
  Supplementary_Table_S9b_Split_Model_Domain_Order.tsv
  Section3.1_Summary_Statistics.tsv
  Section3.1_Summary_Statistics.md

Usage: python make_tables.py [--config config_3.1.json]
"""
import argparse, shutil, sys
from pathlib import Path

STAGE_DIR = Path(__file__).resolve().parent
HERE = STAGE_DIR.parent
sys.path.insert(0, str(HERE / "01_section3.1_pangenome_inventory"))
from rust31 import log, load_cfg, read_tsv, write_tsv, catalogue


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(HERE / "config_3.1.json"))
    args = ap.parse_args()
    cfg = load_cfg(args.config)
    R = Path(cfg["resultsdir"])
    T = HERE / "tables"
    T.mkdir(parents=True, exist_ok=True)

    cat_rows = read_tsv(cfg["catalogue"])
    qcheck = {r["ID"]: r for r in read_tsv(R / "query_anchor_check.tsv")} if (R / "query_anchor_check.tsv").exists() else {}
    ipr = {r["ID"]: r for r in read_tsv(R / "S1d_interproscan_summary.tsv")} if (R / "S1d_interproscan_summary.tsv").exists() else {}

    q2og = {}
    if (R / "S1b_orthogroup_partition.tsv").exists():
        for r in read_tsv(R / "S1b_orthogroup_partition.tsv"):
            for q in r["queries"].split(","):
                q2og[q.strip()] = (r["orthogroup"], r["compartment"], r["occupancy_all_subgenomes"])

    DONOR_ORIGIN_MAP = {
        "Yr15": "Introg. (T. dicoccoides)",
        "Yr36": "Introg. (T. dicoccoides)",
        "Sr60": "Introg. (T. monococcum)",
        "Sr62": "Introg. (Ae. sharonensis)",
        "Lr9": "Introg. (Ae. umbellulata)",
        "Sr43": "Introg. (Th. ponticum)",
        "Lr14a": "Introg. (T. dicoccum)",
        "Sr50": "Introg. (S. cereale)",
        "Sr27": "Introg. (S. cereale)",
        "Sr26": "Introg. (Th. ponticum)",
        "Sr61": "Introg. (Th. ponticum)",
        "Lr47": "Introg. (Ae. speltoides)",
    }

    t1_rows = []
    for r in cat_rows:
        qid = r["ID"]
        og, comp, occ = q2og.get(qid, ("", "", ""))
        ip = ipr.get(qid, {})
        qc = qcheck.get(qid, {})
        pacc = "OU248342.1" if str(r["Protein_acc"]).startswith("derived:miniprot") else r["Protein_acc"]
        t1_rows.append({
            "ID": qid,
            "Gene": r["Gene"],
            "LocusGroup": r["LocusGroup"],
            "Set": r["Set"],
            "Class": r["Class"],
            "Group": r["Group"],
            "Architecture": r["Architecture"],
            "Donor_or_Locus_Origin": DONOR_ORIGIN_MAP.get(str(r["Gene"]), "Native / Triticeae OG"),
            "Protein_acc": pacc,
            "Protein_len_aa": r["Protein_len"],
            "Anchor_rule": r["Anchor"],
            "Diagnostic_Pfam_domains": qc.get("domains_found", ""),
            "InterProScan_Pfam": ip.get("pfam", ""),
            "InterPro_accessions": ip.get("interpro", ""),
            "Orthogroup": og,
            "Subgenome_occupancy": occ,
            "Pangenome_compartment": comp,
        })
    write_tsv(T / "Table1_Cloned_Rust_Resistance_Catalogue.tsv",
              [r for r in t1_rows if r["Set"] == "Main"])
    write_tsv(T / "Supplementary_Table_S1_Queries_and_Domains.tsv", t1_rows)

    # Copy / standardize supplementary tables
    s2_src = R / "S2_assembly_panel_qc_busco.tsv" if (R / "S2_assembly_panel_qc_busco.tsv").exists() else R / "S2_assembly_panel_qc.tsv"
    copy_map = [
        (R / "S1b_orthogroup_partition.tsv", T / "Supplementary_Table_S1b_Orthogroup_Partition.tsv"),
        (s2_src, T / "Supplementary_Table_S2_Assembly_Panel_QC.tsv"),
        (R / "S3_homologue_loci_rescue.tsv", T / "Supplementary_Table_S3_Homologue_Loci_Rescue.tsv"),
        (R / "S3a_recovery_per_assembly.tsv", T / "Supplementary_Table_S3a_Recovery_per_Assembly.tsv"),
        (R / "S3b_recovery_per_locus.tsv", T / "Supplementary_Table_S3b_Recovery_per_Locus.tsv"),
        (R / "S3c_identity_threshold_sensitivity.tsv", T / "Supplementary_Table_S3c_Identity_Sensitivity.tsv"),
        (R / "S3d_heaps_law.tsv", T / "Supplementary_Table_S3d_Heaps_Law.tsv"),
        (R / "S3e_R_statistical_audit_3.1.tsv", T / "Supplementary_Table_S3e_R_Statistical_Audit_3.1.tsv"),
        (R / "S4_chromosomal_positions.tsv", T / "Supplementary_Table_S4_All_Chromosomal_Positions.tsv"),
        (R / "S4b_locusgroup_positions.tsv", T / "Supplementary_Table_S4_Chromosomal_Positions.tsv"),
        (R / "S5_hexaploid_triads.tsv", T / "Supplementary_Table_S5a_Hexaploid_Locus_Triads.tsv"),
        (R / "S5b_orthogroup_triad_summary.tsv", T / "Supplementary_Table_S5_Hexaploid_Triads.tsv"),
        (R / "S5b_orthogroup_triad_summary.tsv", T / "Supplementary_Table_S5b_Orthogroup_Triad_Summary.tsv"),
        (R / "S5c_ensembl_compara_triads.tsv", T / "Supplementary_Table_S5c_Ensembl_Compara_Triads.tsv"),
        (R / "S5d_subgenome_cnv.tsv", T / "Supplementary_Table_S5d_Subgenome_CNV.tsv"),
        (R / "S6a_sri_per_locus.tsv", T / "Supplementary_Table_S6a_SRI_per_Locus.tsv"),
        (R / "S6b_sri_pairwise.tsv", T / "Supplementary_Table_S6b_SRI_Pairwise.tsv"),
        (R / "S7a_flanking_te_windows.tsv", T / "Supplementary_Table_S7a_Flanking_TE_Windows.tsv"),
        (R / "S7b_te_summary_by_class.tsv", T / "Supplementary_Table_S7b_TE_Summary_by_Class.tsv"),
        (R / "S7c_R_statistical_audit_3.2.tsv", T / "Supplementary_Table_S7c_R_Statistical_Audit_3.2.tsv"),
        (R / "S8a_classIII_domains.tsv", T / "Supplementary_Table_S8a_ClassIII_Domains.tsv"),
        (R / "S8b_kinase_phylogeny.tsv", T / "Supplementary_Table_S8b_Kinase_Phylogeny.tsv"),
        (R / "S8c_kinase_pair_divergence.tsv", T / "Supplementary_Table_S8c_Kinase_Pair_Divergence.tsv"),
        (R / "S8d_architecture_distribution.tsv", T / "Supplementary_Table_S8d_Architecture_Distribution.tsv"),
        (R / "S8e_exon_intron_phases.tsv", T / "Supplementary_Table_S8e_Exon_Intron_Phases.tsv"),
        (R / "S8f_R_statistical_audit_3.3.tsv", T / "Supplementary_Table_S8f_R_Statistical_Audit_3.3.tsv"),
        (R / "S9_paralogue_pruning.tsv", T / "Supplementary_Table_S9_Paralogue_Pruning.tsv"),
        (R / "S9b_domain_order_split.tsv", T / "Supplementary_Table_S9b_Split_Model_Domain_Order.tsv"),
        (R / "Supplementary_Table_S10a_codeml_m8.tsv", T / "Supplementary_Table_S10a_Codeml_M8.tsv"),
        (R / "Supplementary_Table_S10b_hyphy_selection.tsv", T / "Supplementary_Table_S10b_HyPhy_Selection.tsv"),
        (R / "Supplementary_Table_S10c_domain_partitioned_omega.tsv", T / "Supplementary_Table_S10c_Domain_Partitioned_Omega.tsv"),
        (R / "Supplementary_Table_S10d_tkp_module_asymmetry.tsv", T / "Supplementary_Table_S10d_TKP_Paired_Asymmetry.tsv"),
        (R / "Supplementary_Table_S10e_af_surface_mapping.tsv", T / "Supplementary_Table_S10e_AF3_Surface_Mapping.tsv"),
        (R / "Supplementary_Table_S10f_R_statistical_audit_3.4.tsv", T / "Supplementary_Table_S10f_R_Statistical_Audit_3.4.tsv"),
        (R / "Supplementary_Table_S10g_site_table.tsv", T / "Supplementary_Table_S10g_Site_Table.tsv"),
        (R / "Supplementary_Table_S11a_population_genetics_tiers.tsv", T / "Supplementary_Table_S11a_population_genetics_tiers.tsv"),
        (R / "Supplementary_Table_S11b_haplotype_diversity_summary.tsv", T / "Supplementary_Table_S11b_haplotype_diversity_summary.tsv"),
        (R / "Supplementary_Table_S11c_candidate_functional_alleles.tsv", T / "Supplementary_Table_S11c_candidate_functional_alleles.tsv"),
        (R / "Supplementary_Table_S11d_watkins_geographic_distribution.tsv", T / "Supplementary_Table_S11d_watkins_geographic_distribution.tsv"),
        (R / "Supplementary_Table_S11e_bottleneck_tests.tsv", T / "Supplementary_Table_S11e_bottleneck_tests.tsv"),
        (R / "Supplementary_Table_S11f_R_statistical_audit_3.5.tsv", T / "Supplementary_Table_S11f_R_Statistical_Audit_3.5.tsv"),
    ]
    for src, dst in copy_map:
        if src.exists():
            shutil.copy(src, dst)

    # Consolidate summary metric values into one master table
    all_vals = []
    for vf in ("section3.1_values.tsv", "section3.1_pruning_values.tsv",
               "section3.1_busco_values.tsv", "section3.1_tree_values.tsv",
               "section3.2_values.tsv", "section3.3_values.tsv",
               "section3.4_values.tsv"):
        p = R / vf
        if p.exists():
            for row in read_tsv(p):
                all_vals.append({
                    "source": vf,
                    "metric": (
                        row.get("metric")
                        or row.get("Metric")
                        or next((v for k, v in row.items() if k not in ("value", "Value", "note", "Notes")), "")
                    ),
                    "value": row.get("value", row.get("Value", "")),
                    "note": row.get("note", row.get("Notes", ""))
                })

    if (R / "S3d_heaps_law.tsv").exists():
        for row in read_tsv(R / "S3d_heaps_law.tsv"):
            grp = row["group"]
            all_vals.append({
                "source": "S3d_heaps_law.tsv",
                "metric": f"Heaps gamma ({grp})",
                "value": f"{float(row['gamma_meancurve']):.4f} +/- {float(row['gamma_meancurve_se']):.4f}",
                "note": f"n_orthogroups={row['n_orthogroups']}, P_label_perm={float(row['P_label_permutation']):.4f}"
            })

    write_tsv(T / "Section3.1_Summary_Statistics.tsv", all_vals,
              ["source", "metric", "value", "note"])

    # Also write a clean Markdown summary report
    md_lines = [
        "# Section 3.1 & Methods 5.1-5.3 Summary Statistics",
        "",
        "| Source | Metric | Value | Note |",
        "|---|---|---|---|",
    ]
    for r in all_vals:
        md_lines.append(f"| `{r['source']}` | {r['metric']} | **{r['value']}** | {r.get('note', '')} |")
    (T / "Section3.1_Summary_Statistics.md").write_text("\n".join(md_lines) + "\n")

    log(f"tables written to {T}/ ({len(all_vals)} total summary metric entries)")


if __name__ == "__main__":
    main()
