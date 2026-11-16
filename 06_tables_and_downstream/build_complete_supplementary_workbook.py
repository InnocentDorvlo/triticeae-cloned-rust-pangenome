#!/usr/bin/env python3
"""
build_complete_supplementary_workbook.py
Generates the definitive, journal-submission-ready multi-tab Excel workbook:
  results/Supplementary_Data_Tables.xlsx
  tables/Supplementary_Data_Tables.xlsx
Includes an Index (Table of Contents) tab, auto-adjusted column widths,
professional header styling, frozen top rows, and full Supplementary Tables S1-S11f + Table 1.
"""

import os, csv
from pathlib import Path
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

STAGE_DIR = Path(__file__).resolve().parent
HERE = STAGE_DIR.parent
TABLES_DIR = HERE / "tables"
RESULTS_DIR = HERE / "results"

# Full inventory of Supplementary Tables with metadata
TABLE_SPECS = [
    {
        "id": "Table 1",
        "sheet": "Table 1",
        "file": "Table1_Cloned_Rust_Resistance_Catalogue.tsv",
        "section": "Section 3.1",
        "title": "Catalogue of 43 cloned Triticeae rust resistance reference proteins across 36 loci, structural classification, and diagnostic Pfam signatures",
    },
    {
        "id": "Supplementary Table 1",
        "sheet": "Supplementary Table 1",
        "file": "Supplementary_Table_S1_Queries_and_Domains.tsv",
        "section": "Section 3.1",
        "title": "Extended query metadata, diagnostic Pfam domains, and InterProScan accessions across all 45 evaluated sequences (43 Main plus 2 comparators)",
    },
    {
        "id": "Supplementary Table 1b",
        "sheet": "Supplementary Table 1b",
        "file": "Supplementary_Table_S1b_Orthogroup_Partition.tsv",
        "section": "Section 3.1",
        "title": "OrthoFinder orthogroup clustering, subgenome occupancy, and pangenome compartment partitioning across 31 orthogroups",
    },
    {
        "id": "Supplementary Table 2",
        "sheet": "Supplementary Table 2",
        "file": "Supplementary_Table_S2_Assembly_Panel_QC.tsv",
        "section": "Section 3.1",
        "title": "Quality metrics, BUSCO v5.5.0 poales completeness, scaffold N50, and ploidy tiers across 27 chromosome-scale Triticeae assemblies",
    },
    {
        "id": "Supplementary Table 3",
        "sheet": "Supplementary Table 3",
        "file": "Supplementary_Table_S3_Homologue_Loci_Rescue.tsv",
        "section": "Section 3.1",
        "title": "Miniprot cross-assembly spliced protein-to-genome alignment, unannotated locus rescue, and GFF3 overlap coordinates (1,590 loci)",
    },
    {
        "id": "Supplementary Table 3a",
        "sheet": "Supplementary Table 3a",
        "file": "Supplementary_Table_S3a_Recovery_per_Assembly.tsv",
        "section": "Section 3.1",
        "title": "Homologue recovery gain per genome assembly comparing Miniprot spliced alignment against primary NCBI/Ensembl GFF3 annotations",
    },
    {
        "id": "Supplementary Table 3b",
        "sheet": "Supplementary Table 3b",
        "file": "Supplementary_Table_S3b_Recovery_per_Locus.tsv",
        "section": "Section 3.1",
        "title": "Per-locus recovery gain across 36 cloned resistance loci, highlighting substantial annotation rescue in non-canonical Classes III-V",
    },
    {
        "id": "Supplementary Table 3c",
        "sheet": "Supplementary Table 3c",
        "file": "Supplementary_Table_S3c_Identity_Sensitivity.tsv",
        "section": "Section 3.1",
        "title": "Sensitivity analysis of pangenome compartment assignments across NLR amino-acid identity thresholds (70% to 95%, paired with 50% to 70% progenitor cutoffs)",
    },
    {
        "id": "Supplementary Table 3d",
        "sheet": "Supplementary Table 3d",
        "file": "Supplementary_Table_S3d_Heaps_Law.tsv",
        "section": "Section 3.1",
        "title": "Heaps' law pan-genome openness trajectory parameters (gamma, alpha) comparing canonical NLRs versus non-NLR architectures",
    },
    {
        "id": "Supplementary Table 3e",
        "sheet": "Supplementary Table 3e",
        "file": "Supplementary_Table_S3e_R_Statistical_Audit_3.1.tsv",
        "section": "Section 3.1",
        "title": "Statistical audit and formal hypothesis testing for Section 3.1 (Heaps permutation test, Wilcoxon paired recovery gain)",
    },
    {
        "id": "Supplementary Table 4",
        "sheet": "Supplementary Table 4",
        "file": "Supplementary_Table_S4_Chromosomal_Positions.tsv",
        "section": "Section 3.2",
        "title": "Chromosomal arm-normalized coordinates (d_norm) and centromere-to-telomere distance medians for 36 cloned R-gene loci",
    },
    {
        "id": "Supplementary Table 4_All",
        "sheet": "Supplementary Table 4_All",
        "file": "Supplementary_Table_S4_All_Chromosomal_Positions.tsv",
        "section": "Section 3.2",
        "title": "Full genomic coordinates and arm-normalized distance for all 12,018 placed resistance homologues across 27 Triticeae assemblies",
    },
    {
        "id": "Supplementary Table 5",
        "sheet": "Supplementary Table 5",
        "file": "Supplementary_Table_S5_Hexaploid_Triads.tsv",
        "section": "Section 3.2",
        "title": "Summary of A:B:D homoeologous triad retention and subgenome bias across 30 Main orthogroups in 13 hexaploid wheat assemblies",
    },
    {
        "id": "Supplementary Table 5a",
        "sheet": "Supplementary Table 5a",
        "file": "Supplementary_Table_S5a_Hexaploid_Locus_Triads.tsv",
        "section": "Section 3.2",
        "title": "Detailed assembly-by-locus A:B:D copy numbers and triad classifications (1:1:1 balanced, dominant, suppressed, or absent)",
    },
    {
        "id": "Supplementary Table 5b",
        "sheet": "Supplementary Table 5b",
        "file": "Supplementary_Table_S5b_Orthogroup_Triad_Summary.tsv",
        "section": "Section 3.2",
        "title": "Orthogroup-level triad retention percentages and Fisher's exact tests comparing Classes I-II NLRs against Classes III-V non-NLRs",
    },
    {
        "id": "Supplementary Table 5c",
        "sheet": "Supplementary Table 5c",
        "file": "Supplementary_Table_S5c_Ensembl_Compara_Triads.tsv",
        "section": "Section 3.2",
        "title": "Independent benchmarking of homoeologous triad calls against Ensembl Plants Compara high-confidence triad annotations",
    },
    {
        "id": "Supplementary Table 5d",
        "sheet": "Supplementary Table 5d",
        "file": "Supplementary_Table_S5d_Subgenome_CNV.tsv",
        "section": "Section 3.2",
        "title": "Subgenome-resolved copy-number variation matrix across 56 haploid subgenomes and 30 cloned resistance orthogroups",
    },
    {
        "id": "Supplementary Table 6a",
        "sheet": "Supplementary Table 6a",
        "file": "Supplementary_Table_S6a_SRI_per_Locus.tsv",
        "section": "Section 3.2",
        "title": "Synteny Retention Index (SRI) and mean local SV rearrangement events (+/-500 kb) across all 36 cloned R-gene loci",
    },
    {
        "id": "Supplementary Table 6b",
        "sheet": "Supplementary Table 6b",
        "file": "Supplementary_Table_S6b_SRI_Pairwise.tsv",
        "section": "Section 3.2",
        "title": "Pairwise assembly-versus-assembly collinear gene retention scores across +/-500 kb flanking windows (3,314 pairwise comparisons)",
    },
    {
        "id": "Supplementary Table 7a",
        "sheet": "Supplementary Table 7a",
        "file": "Supplementary_Table_S7a_Flanking_TE_Windows.tsv",
        "section": "Section 3.2",
        "title": "Flanking +/-200 kb RepeatMasker transposable element (TE) density and coding/intergenic fractions around Chinese Spring RefSeq v1.1 loci",
    },
    {
        "id": "Supplementary Table 7b",
        "sheet": "Supplementary Table 7b",
        "file": "Supplementary_Table_S7b_TE_Summary_by_Class.tsv",
        "section": "Section 3.2",
        "title": "Class- and group-level summary of flanking +/-200 kb RepeatMasker TE density, flanking gene density, and mean Synteny Retention Index (SRI)",
    },
    {
        "id": "Supplementary Table 7c",
        "sheet": "Supplementary Table 7c",
        "file": "Supplementary_Table_S7c_R_Statistical_Audit_3.2.tsv",
        "section": "Section 3.2",
        "title": "Statistical audit and formal tests for Section 3.2 (Kolmogorov-Smirnov position test, triad Fisher test, Spearman SV correlation)",
    },
    {
        "id": "Supplementary Table 8a",
        "sheet": "Supplementary Table 8a",
        "file": "Supplementary_Table_S8a_ClassIII_Domains.tsv",
        "section": "Section 3.3",
        "title": "Dissected domain boundaries (Kinase 1, Kinase 2, accessory modules) and catalytic status across 830 Class III homologues",
    },
    {
        "id": "Supplementary Table 8b",
        "sheet": "Supplementary Table 8b",
        "file": "Supplementary_Table_S8b_Kinase_Phylogeny.tsv",
        "section": "Section 3.3",
        "title": "IQ-TREE 2 phylogenetic support, domain architecture, pangenomic homologue counts, and origin classification across all 8 Class III proteins",
    },
    {
        "id": "Supplementary Table 8c",
        "sheet": "Supplementary Table 8c",
        "file": "Supplementary_Table_S8c_Kinase_Pair_Divergence.tsv",
        "section": "Section 3.3",
        "title": "Intra-polypeptide Kinase 1 versus Kinase 2 amino-acid identity (%) and synonymous substitution divergence (pS) across 523 pairs",
    },
    {
        "id": "Supplementary Table 8d",
        "sheet": "Supplementary Table 8d",
        "file": "Supplementary_Table_S8d_Architecture_Distribution.tsv",
        "section": "Section 3.3",
        "title": "Ploidy-ordered presence-absence matrix of fused tandem kinases, standalone single-kinase progenitors, and accessory modules",
    },
    {
        "id": "Supplementary Table 8e",
        "sheet": "Supplementary Table 8e",
        "file": "Supplementary_Table_S8e_Exon_Intron_Phases.tsv",
        "section": "Section 3.3",
        "title": "Physical transcript exon-intron structures and phase conservation of domain-junction introns across cloned Class III loci",
    },
    {
        "id": "Supplementary Table 8f",
        "sheet": "Supplementary Table 8f",
        "file": "Supplementary_Table_S8f_R_Statistical_Audit_3.3.tsv",
        "section": "Section 3.3",
        "title": "Statistical audit for Section 3.3 (Welch t-tests for kinase divergence, exact binomial test for phase-0 domain-junction introns)",
    },
    {
        "id": "Supplementary Table 9",
        "sheet": "Supplementary Table 9",
        "file": "Supplementary_Table_S9_Paralogue_Pruning.tsv",
        "section": "Methods 5.4",
        "title": "Paralogue pruning filter log and domain integrity verification across 12,847 candidate homologues (retaining full-length ORFs)",
    },
    {
        "id": "Supplementary Table 9b",
        "sheet": "Supplementary Table 9b",
        "file": "Supplementary_Table_S9b_Split_Model_Domain_Order.tsv",
        "section": "Methods 5.4",
        "title": "Domain order verification and split-gene consolidation rules for complex tandem and integrated NLR/TKP loci",
    },
    {
        "id": "Supplementary Table 10a",
        "sheet": "Supplementary Table 10a",
        "file": "Supplementary_Table_S10a_Codeml_M8.tsv",
        "section": "Section 3.4",
        "title": "PAML codeml site models (M0, M1a vs M2a, M7 vs M8) and Bayes Empirical Bayes (BEB >= 0.95) diversifying selection sites",
    },
    {
        "id": "Supplementary Table 10b",
        "sheet": "Supplementary Table 10b",
        "file": "Supplementary_Table_S10b_HyPhy_Selection.tsv",
        "section": "Section 3.4",
        "title": "Full-orthogroup HyPhy episodic diversifying selection (MEME, P < 0.05) and pervasive purifying selection (FUBAR) scans",
    },
    {
        "id": "Supplementary Table 10c",
        "sheet": "Supplementary Table 10c",
        "file": "Supplementary_Table_S10c_Domain_Partitioned_Omega.tsv",
        "section": "Section 3.4",
        "title": "Domain-partitioned dN/dS (omega) ratios across 119 functional partitions spanning Classes I-V (CC/ID, NB-ARC, LRR, Kinase, TMD)",
    },
    {
        "id": "Supplementary Table 10d",
        "sheet": "Supplementary Table 10d",
        "file": "Supplementary_Table_S10d_TKP_Paired_Asymmetry.tsv",
        "section": "Section 3.4",
        "title": "Paired tandem kinase evolutionary rate asymmetry between active signaling kinases and pseudokinase decoys (5 loci)",
    },
    {
        "id": "Supplementary Table 10e",
        "sheet": "Supplementary Table 10e",
        "file": "Supplementary_Table_S10e_AF3_Surface_Mapping.tsv",
        "section": "Section 3.4",
        "title": "AlphaFold 3D structural mapping, DSSP secondary structure, RSA > 0.25, and LRR concave beta-strand enrichment Fisher tests (20 query entries across 19 non-redundant models)",
    },
    {
        "id": "Supplementary Table 10f",
        "sheet": "Supplementary Table 10f",
        "file": "Supplementary_Table_S10f_R_Statistical_Audit_3.4.tsv",
        "section": "Section 3.4",
        "title": "Statistical audit for Section 3.4 (Kruskal-Wallis 5-class domain omega test, Mann-Whitney LRR vs NB-ARC, pooled Fisher test)",
    },
    {
        "id": "Supplementary Table 10g",
        "sheet": "Supplementary Table 10g",
        "file": "Supplementary_Table_S10g_Site_Table.tsv",
        "section": "Section 3.4",
        "title": "Complete codon-by-codon site selection table (31,618 codons) with PAML BEB PP, HyPhy MEME P-values, and FUBAR alpha/beta rates",
    },
    {
        "id": "Supplementary Table 11a",
        "sheet": "Supplementary Table 11a",
        "file": "Supplementary_Table_S11a_population_genetics_tiers.tsv",
        "section": "Section 3.5",
        "title": "Orthogroup-wide nucleotide and amino-acid diversity (pi_aa), Watterson theta, Tajima's D, and FST across 3 germplasm tiers",
    },
    {
        "id": "Supplementary Table 11b",
        "sheet": "Supplementary Table 11b",
        "file": "Supplementary_Table_S11b_haplotype_diversity_summary.tsv",
        "section": "Section 3.5",
        "title": "Summary of 2,197 full-length protein haplotypes partitioned into elite-captured, wild-private, and landrace-private reservoirs",
    },
    {
        "id": "Supplementary Table 11c",
        "sheet": "Supplementary Table 11c",
        "file": "Supplementary_Table_S11c_candidate_functional_alleles.tsv",
        "section": "Section 3.5",
        "title": "Catalog of 249 prioritized uncaptured functional candidate alleles with intact catalytic/signaling domains for wheat breeding",
    },
    {
        "id": "Supplementary Table 11d",
        "sheet": "Supplementary Table 11d",
        "file": "Supplementary_Table_S11d_watkins_geographic_distribution.tsv",
        "section": "Section 3.5",
        "title": "Geographic origin distribution of 249 prioritized resistance alleles, showing concentration in West and Central Asia center of diversity",
    },
    {
        "id": "Supplementary Table 11e",
        "sheet": "Supplementary Table 11e",
        "file": "Supplementary_Table_S11e_bottleneck_tests.tsv",
        "section": "Section 3.5",
        "title": "Statistical evaluation of diversity erosion (-20.6% pi_aa, -30.2% theta_W) from wild progenitors to modern elite hexaploid cultivars",
    },
    {
        "id": "Supplementary Table 11f",
        "sheet": "Supplementary Table 11f",
        "file": "Supplementary_Table_S11f_R_Statistical_Audit_3.5.tsv",
        "section": "Section 3.5",
        "title": "Statistical audit for Section 3.5 (Wilcoxon paired diversity erosion test, Mann-Whitney Tajima's D test, Chi-square capture test)",
    },
]

def main():
    print("=" * 75)
    print("   Building Complete Submission Supplementary Excel Workbook")
    print("=" * 75)
    
    wb = openpyxl.Workbook()
    # Remove default sheet
    default_sheet = wb.active
    
    # -------------------------------------------------------------------------
    # Styles
    # -------------------------------------------------------------------------
    font_header = Font(name="Calibri", size=11, bold=True, color="FFFFFF")
    fill_header = PatternFill(start_color="1F4E79", end_color="1F4E79", fill_type="solid")
    
    font_index_head = Font(name="Calibri", size=11, bold=True, color="FFFFFF")
    fill_index_head = PatternFill(start_color="002060", end_color="002060", fill_type="solid")
    
    font_data = Font(name="Calibri", size=10)
    font_title = Font(name="Calibri", size=14, bold=True, color="002060")
    font_subtitle = Font(name="Calibri", size=10, italic=True, color="595959")
    
    border_thin = Border(
        left=Side(style='thin', color='D9D9D9'),
        right=Side(style='thin', color='D9D9D9'),
        top=Side(style='thin', color='D9D9D9'),
        bottom=Side(style='thin', color='D9D9D9')
    )
    fill_alt = PatternFill(start_color="F2F5F8", end_color="F2F5F8", fill_type="solid")
    
    # -------------------------------------------------------------------------
    # 1. Create Index Worksheet (Table of Contents)
    # -------------------------------------------------------------------------
    ws_idx = wb.create_sheet(title="Index")
    ws_idx.views.sheetView[0].showGridLines = True
    
    ws_idx.cell(row=1, column=1, value="Supplementary Data: Index of Supplementary Tables").font = font_title
    ws_idx.cell(row=2, column=1, value="Pangenome dynamics across the Triticeae reveal contrasting evolution of cloned NLR and non-NLR resistance (Sections 3.1 - 3.5)").font = font_subtitle
    
    headers_idx = ["Table ID", "Worksheet Tab", "Analysis Section", "Table Title", "Data Rows", "Key Findings / Focus"]
    for col_idx, h in enumerate(headers_idx, start=1):
        c = ws_idx.cell(row=4, column=col_idx, value=h)
        c.font = font_index_head
        c.fill = fill_index_head
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    ws_idx.row_dimensions[4].height = 25
    
    idx_row = 5
    
    # -------------------------------------------------------------------------
    # 2. Iterate through each table, load TSV, and create worksheet
    # -------------------------------------------------------------------------
    for spec in TABLE_SPECS:
        orig_id = spec["id"]
        if orig_id == "Table 1":
            tid = "Table 1"
            sheet_name = "Table 1"
        elif orig_id.startswith("Table S"):
            num_part = orig_id[7:]
            tid = f"Supplementary Table {num_part}"
            sheet_name = f"Supplementary Table {num_part}"
        elif orig_id.startswith("Supplementary Table"):
            tid = orig_id
            sheet_name = orig_id
        else:
            tid = orig_id
            sheet_name = orig_id[:31]
            
        tsv_file = TABLES_DIR / spec["file"]
        if not tsv_file.exists():
            tsv_file = RESULTS_DIR / spec["file"]
        
        if not tsv_file.exists():
            print(f"Warning: File not found: {spec['file']}")
            continue
            
        print(f"Processing {tid:26s} -> Sheet '{sheet_name}' (from {tsv_file.name})...")
        
        # Read TSV
        with open(tsv_file, "r", encoding="utf-8", errors="replace") as fp:
            reader = csv.reader(fp, delimiter="\t")
            rows = list(reader)
            
        if not rows:
            continue
            
        header = rows[0]
        data_rows = rows[1:]
        
        # Record in Index
        ws_idx.cell(row=idx_row, column=1, value=tid).font = Font(name="Calibri", size=10, bold=True)
        c_link = ws_idx.cell(row=idx_row, column=2)
        c_link.value = f'=HYPERLINK("#\'{sheet_name}\'!A1", "{sheet_name}")'
        c_link.font = Font(name="Calibri", size=10, color="0000FF", underline="single")
        ws_idx.cell(row=idx_row, column=3, value=spec["section"]).font = font_data
        ws_idx.cell(row=idx_row, column=4, value=spec["title"]).font = font_data
        ws_idx.cell(row=idx_row, column=5, value=len(data_rows)).font = font_data
        ws_idx.cell(row=idx_row, column=5).alignment = Alignment(horizontal="right")
        ws_idx.cell(row=idx_row, column=6, value=f"{header[0]} -> {len(header)} attributes").font = Font(name="Calibri", size=9, italic=True)
        
        for c in range(1, 7):
            ws_idx.cell(row=idx_row, column=c).border = border_thin
            if idx_row % 2 == 0:
                ws_idx.cell(row=idx_row, column=c).fill = fill_alt
        idx_row += 1
        
        # Create Data Worksheet
        ws = wb.create_sheet(title=sheet_name)
        ws.views.sheetView[0].showGridLines = True
        
        # Write Title Header on Data Sheet
        ws.cell(row=1, column=1, value=f"{tid}: {spec['title']}").font = Font(name="Calibri", size=11, bold=True, color="002060")
        ws.cell(row=2, column=1, value=f"Section: {spec['section']} | Records: {len(data_rows)} rows | Source file: {tsv_file.name}").font = font_subtitle
        ws.cell(row=3, column=1, value='=HYPERLINK("#\'Index\'!A1", "<- Back to Index")').font = Font(name="Calibri", size=9, bold=True, color="002060", underline="single")
        
        # Write Data Table Headers
        for col_idx, col_name in enumerate(header, start=1):
            cell = ws.cell(row=4, column=col_idx, value=col_name)
            cell.font = font_header
            cell.fill = fill_header
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            cell.border = border_thin
        ws.row_dimensions[4].height = 24
        
        # Write Data Rows (cap at 65000 rows if site table to keep file snappy, though site table is 32k which fits easily)
        for r_idx, row_data in enumerate(data_rows, start=5):
            for c_idx, val in enumerate(row_data, start=1):
                # Attempt numeric conversion
                cell = ws.cell(row=r_idx, column=c_idx)
                if val == "" or val == "NA":
                    cell.value = val
                else:
                    try:
                        if "." in val or "e" in val.lower():
                            cell.value = float(val)
                        else:
                            cell.value = int(val)
                    except ValueError:
                        cell.value = val
                
                cell.font = font_data
                cell.border = border_thin
                if r_idx % 2 == 0:
                    cell.fill = fill_alt

        if tid in ["Table 1", "Supplementary Table 1"]:
            fn_row = len(data_rows) + 6
            ws.cell(row=fn_row, column=1, value="Note on Yr10/YrNAM: Although NCBI accession AF149112.1 was historically deposited as a putative CC-NBS-LRR candidate for Yr10 (Laroche et al., 2000), subsequent fine-mapping and mutational genomics (Ni et al., 2023, Nat. Commun. 14:4353) established that AF149112.1 maps 1.2 cM (~1.0 Mb) from the resistance locus and does not confer resistance to Pst CYR29. The true, functionally validated causal gene at the Yr10 locus is the non-NLR transcription factor YrNAM (UYW61986.1; NAM-ZnF-BED, 408 aa), cataloged here under Class IV.").font = Font(name="Calibri", size=9, italic=True, color="555555")
                    
        # Freeze panes at A5 (so row 4 headers stay locked when scrolling)
        ws.freeze_panes = "A5"
        
        # Auto-adjust column widths (sample first 100 rows for speed on large sheets)
        for col in ws.columns:
            col_letter = get_column_letter(col[0].column)
            max_len = 0
            for cell in col[3:103]:  # sample from header row down 100 rows
                if cell.value:
                    val_str = str(cell.value)
                    max_len = max(max_len, len(val_str))
            ws.column_dimensions[col_letter].width = min(max(max_len + 3, 10), 45)
            
    # Remove default empty sheet if present
    if default_sheet in wb.worksheets:
        wb.remove(default_sheet)
        
    # Auto-adjust Index column widths
    ws_idx.freeze_panes = "A5"
    for col in ws_idx.columns:
        col_letter = get_column_letter(col[0].column)
        max_len = 0
        for cell in col[3:]:
            if cell.value:
                max_len = max(max_len, len(str(cell.value)))
        ws_idx.column_dimensions[col_letter].width = min(max(max_len + 3, 12), 65)
        
    # Save Workbook
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    out_xlsx_results = RESULTS_DIR / "Supplementary_Data_Tables.xlsx"
    out_xlsx_tables  = TABLES_DIR / "Supplementary_Data_Tables.xlsx"
    out_xlsx_desktop = HERE / "Supplementary_Data_Tables.xlsx"
    
    print(f"\nSaving Master Workbook -> {out_xlsx_desktop} ...")
    wb.save(out_xlsx_desktop)
    wb.save(out_xlsx_results)
    wb.save(out_xlsx_tables)
    
    print(f"✓ Master Supplementary Excel Workbook successfully generated ({len(TABLE_SPECS)} tables)!")
    print(f"  File size: {os.path.getsize(out_xlsx_desktop) / (1024*1024):.2f} MB")
    print("=" * 75)

if __name__ == "__main__":
    import argparse
    argparse.ArgumentParser(description="Build Supplementary_Data_Tables.xlsx from tables/ and results/").parse_args()
    main()
