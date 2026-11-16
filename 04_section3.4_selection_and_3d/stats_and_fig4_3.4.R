#!/usr/bin/env Rscript
# ==============================================================================
# stats_and_fig4_3.4.R
# Standalone R statistical audit and publication figure generation for Result 3.4:
#   - Supplementary_Table_S10a_codeml_m8.tsv (PAML codeml M1a/M2a, M7/M8 & BEB sites)
#   - Supplementary_Table_S10b_hyphy_selection.tsv (HyPhy FUBAR & MEME selection)
#   - Supplementary_Table_S10c_domain_partitioned_omega.tsv (Domain-partitioned codeml M0 dN/dS)
#   - Supplementary_Table_S10d_tkp_module_asymmetry.tsv (Class III TKP module asymmetry)
#   - Supplementary_Table_S10e_af_surface_mapping.tsv (AlphaFold 3D surface mapping)
#   - Outputs:
#       results/Supplementary_Table_S10f_R_statistical_audit_3.4.tsv
#       tables/Supplementary_Table_S10f_R_Statistical_Audit_3.4.tsv
#       figures/Fig4_selection_structures_R.pdf and .png
# ==============================================================================


.cmd_args <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", .cmd_args[grep("^--file=", .cmd_args)])
.root_dir <- if (length(.file_arg) > 0) {
  dirname(dirname(normalizePath(.file_arg[1])))
} else if (file.exists("config_3.1.json")) {
  getwd()
} else if (file.exists("../config_3.1.json")) {
  normalizePath("..")
} else {
  getwd()
}
setwd(.root_dir)

options(stringsAsFactors = FALSE, width = 120, scipen = 999)

cat("======================================================================\n")
cat("  Section 3.4 R Statistical Audit & 3D Selection Regime Analysis\n")
cat("======================================================================\n\n")

# ------------------------------------------------------------------------------
# Load Section 3.4 Datasets
# ------------------------------------------------------------------------------
s10a <- read.delim("results/Supplementary_Table_S10a_codeml_m8.tsv", check.names = FALSE)
s10b <- read.delim("results/Supplementary_Table_S10b_hyphy_selection.tsv", check.names = FALSE)
s10c <- read.delim("results/Supplementary_Table_S10c_domain_partitioned_omega.tsv", check.names = FALSE)
s10d <- read.delim("results/Supplementary_Table_S10d_tkp_module_asymmetry.tsv", check.names = FALSE)
s10e <- read.delim("results/Supplementary_Table_S10e_af_surface_mapping.tsv", check.names = FALSE)

med_iqr <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.50, 0.75), type = 7, na.rm = TRUE)
  sprintf("%.3f (IQR = %.3f-%.3f)", q[2], q[1], q[3])
}

# ------------------------------------------------------------------------------
# [1] Domain-Partitioned Selection (omega) Across Structural Classes I-V (S10c)
# ------------------------------------------------------------------------------
c1_lrr <- as.numeric(s10c$omega_M0[s10c$structural_class == "I" & s10c$domain_name == "LRR"])
c2_lrr <- as.numeric(s10c$omega_M0[s10c$structural_class %in% c("I/II", "II") & s10c$domain_name == "LRR"])
c3_kin <- as.numeric(s10c$omega_M0[s10c$structural_class == "III" & s10c$domain_type == "kinase_module"])
c4_ank <- as.numeric(s10c$omega_M0[s10c$structural_class == "IV" & s10c$domain_name != "whole_gene"])
c5_tmd <- as.numeric(s10c$omega_M0[s10c$structural_class == "V" & s10c$domain_name != "whole_gene"])

kw_res <- kruskal.test(list(
  ClassI_LRR = c1_lrr,
  ClassII_LRR = c2_lrr,
  ClassIII_Kinase = c3_kin,
  ClassIV_Domains = c4_ank,
  ClassV_Domains = c5_tmd
))

nlr_lrr_w <- c(c1_lrr, c2_lrr)
non_nlr_dom_w <- c(c3_kin, c4_ank, c5_tmd)
mwu_dom_exact <- suppressWarnings(wilcox.test(nlr_lrr_w, non_nlr_dom_w, exact = TRUE))
mwu_dom_asymp <- wilcox.test(nlr_lrr_w, non_nlr_dom_w, exact = FALSE, correct = TRUE)

# Intra-NLR domain stratification (CC/ID vs NB-ARC vs LRR across n = 20 NLR orthogroups)
nlr_cc_w <- as.numeric(s10c$omega_M0[s10c$structural_class %in% c("I", "I/II", "II") & s10c$domain_name == "N-terminal_CC_or_ID"])
nlr_nb_w <- as.numeric(s10c$omega_M0[s10c$structural_class %in% c("I", "I/II", "II") & s10c$domain_name == "NB-ARC"])
kw_intra_nlr <- kruskal.test(list(CC_ID = nlr_cc_w, NB_ARC = nlr_nb_w, LRR = nlr_lrr_w))
mwu_lrr_nbarc <- wilcox.test(nlr_lrr_w, nlr_nb_w, exact = FALSE, correct = TRUE)

cat(sprintf("[1] Domain-Partitioned dN/dS (omega) Across Structural Classes I-V (S10c, n = %d partitions):\n", nrow(s10c)))
cat(sprintf("    Class I   CC-NLR LRR        (n = %2d): median = %s\n", length(c1_lrr), med_iqr(c1_lrr)))
cat(sprintf("    Class II  NLR-ID LRR        (n = %2d): median = %s\n", length(c2_lrr), med_iqr(c2_lrr)))
cat(sprintf("    Class III TKP Kinase Modules(n = %2d): median = %s\n", length(c3_kin), med_iqr(c3_kin)))
cat(sprintf("    Class IV  ANK-TM Domains    (n = %2d): median = %s\n", length(c4_ank), med_iqr(c4_ank)))
cat(sprintf("    Class V   APR TMD/NBD       (n = %2d): median = %s\n", length(c5_tmd), med_iqr(c5_tmd)))
cat(sprintf("    Kruskal-Wallis test across 5 classes: H = %.2f (df = %d), P = %.4e\n",
            kw_res$statistic, kw_res$parameter, kw_res$p.value))
cat(sprintf("    NLR LRR (n = %d) vs non-NLR functional domains (n = %d): MWU U = %.1f, P = %.4e\n",
            length(nlr_lrr_w), length(non_nlr_dom_w), mwu_dom_asymp$statistic, mwu_dom_asymp$p.value))
cat(sprintf("    Intra-NLR domain hierarchy: CC/ID [%.3f] vs NB-ARC [%.3f] vs LRR [%.3f] (KW H = %.2f, P = %.4e; LRR vs NB-ARC MWU P = %.4e)\n\n",
            median(nlr_cc_w), median(nlr_nb_w), median(nlr_lrr_w), kw_intra_nlr$statistic, kw_intra_nlr$p.value, mwu_lrr_nbarc$p.value))

# ------------------------------------------------------------------------------
# [2] Site-Level Diversifying & Purifying Selection (PAML codeml S10a & HyPhy S10b)
# ------------------------------------------------------------------------------
is_nlr_a <- s10a$structural_class %in% c("I", "I/II", "II")
is_nlr_b <- s10b$structural_class %in% c("I", "I/II", "II")

nlr_beb <- as.numeric(s10a$n_BEB_positive_sites_p95[is_nlr_a])
non_beb <- as.numeric(s10a$n_BEB_positive_sites_p95[!is_nlr_a])
mwu_beb <- wilcox.test(nlr_beb, non_beb, exact = FALSE, correct = TRUE)

nlr_meme <- as.numeric(s10b$n_meme_episodic_p05[is_nlr_b])
non_meme <- as.numeric(s10b$n_meme_episodic_p05[!is_nlr_b])
mwu_meme <- wilcox.test(nlr_meme, non_meme, exact = FALSE, correct = TRUE)

nlr_lrt <- as.numeric(s10a$LRT_M7_M8[is_nlr_a])
non_lrt <- as.numeric(s10a$LRT_M7_M8[!is_nlr_a])
mwu_lrt <- wilcox.test(nlr_lrt, non_lrt, exact = FALSE, correct = TRUE)

non_pur_agb <- as.numeric(s10b$pct_fubar_purifying[!is_nlr_b])
nlr_pur_agb <- as.numeric(s10b$pct_fubar_purifying[is_nlr_b])

# Domain-level FUBAR purifying % (S10c)
nlr_lrr_pur <- as.numeric(s10c$pct_fubar_purifying[s10c$structural_class %in% c("I", "I/II", "II") & s10c$domain_name == "LRR"])
non_dom_pur <- as.numeric(s10c$pct_fubar_purifying[s10c$structural_class %in% c("III", "IV", "V") & s10c$domain_name != "whole_gene"])
mwu_dom_pur <- wilcox.test(non_dom_pur, nlr_lrr_pur, exact = FALSE, correct = TRUE)

cat(sprintf("[2] Site-Level Diversifying (PAML M8 BEB, S10a; HyPhy MEME, S10b) & Purifying Selection (HyPhy FUBAR, n = %d OGs):\n", nrow(s10a)))
cat(sprintf("    PAML M8 BEB (PP >= 0.95) sites: NLR (n = %d) = %.1f +/- %.1f vs non-NLR (n = %d) = %.1f +/- %.1f (MWU U = %.1f, P = %.4e)\n",
            length(nlr_beb), mean(nlr_beb), sd(nlr_beb)/sqrt(length(nlr_beb)),
            length(non_beb), mean(non_beb), sd(non_beb)/sqrt(length(non_beb)),
            mwu_beb$statistic, mwu_beb$p.value))
cat(sprintf("    HyPhy MEME (P < 0.05) sites:     NLR (n = %d) = %.1f +/- %.1f vs non-NLR (n = %d) = %.1f +/- %.1f (MWU U = %.1f, P = %.4e)\n",
            length(nlr_meme), mean(nlr_meme), sd(nlr_meme)/sqrt(length(nlr_meme)),
            length(non_meme), mean(non_meme), sd(non_meme)/sqrt(length(non_meme)),
            mwu_meme$statistic, mwu_meme$p.value))
cat(sprintf("    FUBAR purifying codons (alpha > beta): non-NLR whole-gene = %.1f +/- %.1f%% vs NLR LRR = %.1f +/- %.1f%% (MWU P = %.4e)\n\n",
            mean(non_pur_agb), sd(non_pur_agb)/sqrt(length(non_pur_agb)),
            mean(nlr_lrr_pur), sd(nlr_lrr_pur)/sqrt(length(nlr_lrr_pur)),
            mwu_dom_pur$p.value))

# ------------------------------------------------------------------------------
# [3] Class III Module Asymmetry (S10d)
# ------------------------------------------------------------------------------
cat(sprintf("[3] Class III Tandem Kinase Module Asymmetry (S10d, n = %d pairs):\n", nrow(s10d)))
for (i in seq_len(nrow(s10d))) {
  r <- s10d[i, ]
  cat(sprintf("    %s (%s, %s): %s (w = %s) vs %s (w = %s) | delta_w = %s, site MWU P = %s\n",
              r$reference_query, r$orthogroup, r$pair_category,
              r$module_1, r$module_1_omega_M0, r$module_2, r$module_2_omega_M0,
              r$delta_omega_m2_minus_m1, r$site_MWU_P))
}
cat("\n")

# ------------------------------------------------------------------------------
# [4] AlphaFold 3D Surface Mapping & LRR Concave Beta-Sheet Enrichment (S10e)
# ------------------------------------------------------------------------------
c1_s10e_raw <- s10e[s10e$structural_class == "I", ]
c1_s10e <- c1_s10e_raw[!duplicated(paste(c1_s10e_raw$orthogroup, c1_s10e_raw$uniprot_id, sep = "|")), ]
tot_pos_c1 <- sum(as.integer(c1_s10e$n_selected_sites_BEB95_or_MEME05))
tot_lrr_pos_c1 <- sum(as.integer(c1_s10e$selected_in_LRR))
tot_conc_pos_c1 <- sum(as.integer(c1_s10e$selected_on_concave_face))
tot_conc_res_c1 <- sum(as.integer(c1_s10e$concave_face_residues))
tot_res_c1 <- sum(as.integer(c1_s10e$model_residues))

fisher_ors <- numeric(nrow(c1_s10e))
fisher_pvals <- numeric(nrow(c1_s10e))
for (i in seq_len(nrow(c1_s10e))) {
  r <- c1_s10e[i, ]
  a <- as.integer(r$selected_on_concave_face)
  b <- as.integer(r$n_selected_sites_BEB95_or_MEME05) - a
  c_cell <- max(1L, as.integer(r$concave_face_residues) - a)
  d_cell <- max(1L, (as.integer(r$model_residues) - as.integer(r$concave_face_residues)) - b)
  mat <- matrix(c(a, c_cell, b, d_cell), nrow = 2, byrow = FALSE)
  ft <- fisher.test(mat, alternative = "greater")
  fisher_ors[i] <- (as.numeric(a) * as.numeric(d_cell)) / (as.numeric(b) * as.numeric(c_cell))
  fisher_pvals[i] <- ft$p.value
}

a_pool <- tot_conc_pos_c1
b_pool <- tot_pos_c1 - tot_conc_pos_c1
c_pool <- tot_conc_res_c1 - a_pool
d_pool <- (tot_res_c1 - tot_conc_res_c1) - b_pool
ft_pool <- fisher.test(matrix(c(a_pool, c_pool, b_pool, d_pool), nrow = 2), alternative = "greater")
or_pool_mle <- (as.numeric(a_pool) * as.numeric(d_pool)) / (as.numeric(b_pool) * as.numeric(c_pool))

cat(sprintf("[4] AlphaFold 3D Surface Mapping & LRR Concave Beta-Sheet Enrichment (S10e, n = %d structures):\n", nrow(s10e)))
for (i in seq_len(nrow(c1_s10e))) {
  r <- c1_s10e[i, ]
  cat(sprintf("    %s (Class I, %s): %d/%d (%.1f%%) pos sites on concave beta-sheet (%d/%d residues) | OR = %.2f, Fisher P = %.4e\n",
              r$gene, r$uniprot_id, as.integer(r$selected_on_concave_face), as.integer(r$n_selected_sites_BEB95_or_MEME05),
              100 * as.integer(r$selected_on_concave_face) / as.integer(r$n_selected_sites_BEB95_or_MEME05),
              as.integer(r$concave_face_residues), as.integer(r$model_residues), fisher_ors[i], fisher_pvals[i]))
}
cat(sprintf("    Pooled Class I CC-NLR: %d/%d (%.1f%%) in LRR; %d/%d (%.1f%%) on concave face vs %d/%d (%.2f%%) background\n",
            tot_lrr_pos_c1, tot_pos_c1, 100 * tot_lrr_pos_c1 / tot_pos_c1,
            tot_conc_pos_c1, tot_pos_c1, 100 * tot_conc_pos_c1 / tot_pos_c1,
            tot_conc_res_c1, tot_res_c1, 100 * tot_conc_res_c1 / tot_res_c1))
cat(sprintf("    Pooled Fisher's exact test: MLE OR = %.2f (conditional OR = %.2f), P = %.4e\n\n",
            or_pool_mle, ft_pool$estimate, ft_pool$p.value))

# ------------------------------------------------------------------------------
# Save R Statistical Audit Summary Table
# ------------------------------------------------------------------------------
audit_tbl <- data.frame(
  Test_ID = c(
    "S3.4_01_KW_Domain_Omega_5Classes",
    "S3.4_02_MWU_NLR_LRR_vs_NonNLR_Domains",
    "S3.4_03_KW_IntraNLR_Domains_CC_NBARC_LRR",
    "S3.4_04_MWU_LRR_vs_NBARC",
    "S3.4_05_MWU_PAML_M8_BEB_Sites_NLR_vs_NonNLR",
    "S3.4_06_MWU_HyPhy_MEME_Sites_NLR_vs_NonNLR",
    "S3.4_07_MWU_FUBAR_Purifying_NLR_LRR_vs_NonNLR",
    "S3.4_08_Yr15_Module_Asymmetry_FUBAR_MWU",
    "S3.4_09_Sr60_Module_Asymmetry_FUBAR_MWU",
    "S3.4_10_Sr62_Module_Asymmetry_FUBAR_MWU",
    "S3.4_11_Fisher_Concave_LRR_Sr35",
    "S3.4_12_Fisher_Concave_LRR_Yr27",
    "S3.4_13_Fisher_Concave_LRR_Sr33",
    "S3.4_14_Fisher_Concave_LRR_Sr22",
    "S3.4_15_Fisher_Concave_LRR_Pooled_ClassI"
  ),
  Description = c(
    "Kruskal-Wallis test of domain-partitioned codeml M0 dN/dS across 5 classes (39 domain partitions across 31 OGs: NLR LRR n=20, Class III kinase n=12, Class IV n=4, Class V n=3)",
    "Mann-Whitney U test: Classes I-II LRR domain partitions (n=20 across 20 OGs) vs Classes III-V core functional domain partitions (n=19 across 11 OGs, excluding 3 accessory domains)",
    "Kruskal-Wallis test across NLR domain partitions: CC/ID (n=20) vs NB-ARC (n=20) vs LRR (n=20)",
    "Mann-Whitney U test: NLR LRR (n=20) vs NLR NB-ARC (n=20)",
    "Full multi-assembly alignment test (n=4-26 CDS): Mann-Whitney U test of PAML M8 BEB (PP>=0.95) sites: NLR (n=20) vs non-NLR (n=11)",
    "Primary full-orthogroup site-level test: Mann-Whitney U test of HyPhy MEME (P<0.05) diversifying sites: NLR (n=20) vs non-NLR (n=11)",
    "Mann-Whitney U test: FUBAR purifying codons (%) in non-NLR domain partitions (n=22, including 3 Class III accessory domains) vs NLR LRR (n=20)",
    "Mann-Whitney U test on site-level FUBAR beta-alpha: Yr15 K1_active vs K2_pseudokinase",
    "Mann-Whitney U test on site-level FUBAR beta-alpha: Sr60 K2_active vs K1_pseudokinase",
    "Mann-Whitney U test on site-level FUBAR beta-alpha: Sr62 K1_active vs K2_pseudokinase",
    "One-sided Fisher's exact test for solvent-exposed LRR beta-strand (DSSP E/B, RSA>0.25) enrichment in Sr35",
    "One-sided Fisher's exact test for solvent-exposed LRR beta-strand (DSSP E/B, RSA>0.25) enrichment in Yr27",
    "One-sided Fisher's exact test for solvent-exposed LRR beta-strand (DSSP E/B, RSA>0.25) enrichment in Sr33",
    "One-sided Fisher's exact test for solvent-exposed LRR beta-strand (DSSP E/B, RSA>0.25) enrichment in Sr22",
    "Pooled one-sided Fisher's exact test for solvent-exposed LRR beta-strand (DSSP E/B, RSA>0.25) enrichment across 11 Class I CC-NLRs"
  ),
  Group1_Summary = c(
    sprintf("Class I LRR med=%.3f (IQR %.3f-%.3f); Class II LRR med=%.3f", median(c1_lrr), quantile(c1_lrr,0.25), quantile(c1_lrr,0.75), median(c2_lrr)),
    sprintf("NLR LRR (n=%d): med=%.3f", length(nlr_lrr_w), median(nlr_lrr_w)),
    sprintf("LRR med=%.3f; CC/ID med=%.3f", median(nlr_lrr_w), median(nlr_cc_w)),
    sprintf("LRR (n=%d): med=%.3f", length(nlr_lrr_w), median(nlr_lrr_w)),
    sprintf("NLR (n=%d): %.1f +/- %.1f sites", length(nlr_beb), mean(nlr_beb), sd(nlr_beb)/sqrt(length(nlr_beb))),
    sprintf("NLR (n=%d): %.1f +/- %.1f sites", length(nlr_meme), mean(nlr_meme), sd(nlr_meme)/sqrt(length(nlr_meme))),
    sprintf("Non-NLR domains (n=%d): %.1f +/- %.1f%%", length(non_dom_pur), mean(non_dom_pur), sd(non_dom_pur)/sqrt(length(non_dom_pur))),
    sprintf("Yr15 K2_pseudokinase w=%.4f (n=%s codons)", as.numeric(s10d$module_2_omega_M0[s10d$reference_query=="Yr15"]), s10d$module_2_codons[s10d$reference_query=="Yr15"]),
    sprintf("Sr60 K1_pseudokinase w=%.4f (n=%s codons)", as.numeric(s10d$module_2_omega_M0[s10d$reference_query=="Sr60"]), s10d$module_2_codons[s10d$reference_query=="Sr60"]),
    sprintf("Sr62 K1_active w=%.4f (n=%s codons)", as.numeric(s10d$module_1_omega_M0[s10d$reference_query=="Sr62"]), s10d$module_1_codons[s10d$reference_query=="Sr62"]),
    sprintf("%s concave pos: %d/%d (%.1f%%)", "Sr35", as.integer(c1_s10e$selected_on_concave_face[match("Sr35", c1_s10e$gene)]), as.integer(c1_s10e$n_selected_sites_BEB95_or_MEME05[match("Sr35", c1_s10e$gene)]), 100*as.integer(c1_s10e$selected_on_concave_face[match("Sr35", c1_s10e$gene)])/as.integer(c1_s10e$n_selected_sites_BEB95_or_MEME05[match("Sr35", c1_s10e$gene)])),
    sprintf("%s concave pos: %d/%d (%.1f%%)", "Yr27", as.integer(c1_s10e$selected_on_concave_face[match("Yr27", c1_s10e$gene)]), as.integer(c1_s10e$n_selected_sites_BEB95_or_MEME05[match("Yr27", c1_s10e$gene)]), 100*as.integer(c1_s10e$selected_on_concave_face[match("Yr27", c1_s10e$gene)])/as.integer(c1_s10e$n_selected_sites_BEB95_or_MEME05[match("Yr27", c1_s10e$gene)])),
    sprintf("%s concave pos: %d/%d (%.1f%%)", "Sr33", as.integer(c1_s10e$selected_on_concave_face[match("Sr33", c1_s10e$gene)]), as.integer(c1_s10e$n_selected_sites_BEB95_or_MEME05[match("Sr33", c1_s10e$gene)]), 100*as.integer(c1_s10e$selected_on_concave_face[match("Sr33", c1_s10e$gene)])/as.integer(c1_s10e$n_selected_sites_BEB95_or_MEME05[match("Sr33", c1_s10e$gene)])),
    sprintf("%s concave pos: %d/%d (%.1f%%)", "Sr22", as.integer(c1_s10e$selected_on_concave_face[match("Sr22", c1_s10e$gene)]), as.integer(c1_s10e$n_selected_sites_BEB95_or_MEME05[match("Sr22", c1_s10e$gene)]), 100*as.integer(c1_s10e$selected_on_concave_face[match("Sr22", c1_s10e$gene)])/as.integer(c1_s10e$n_selected_sites_BEB95_or_MEME05[match("Sr22", c1_s10e$gene)])),
    sprintf("Pooled Class I concave pos: %d/%d (%.1f%%)", tot_conc_pos_c1, tot_pos_c1, 100*tot_conc_pos_c1/tot_pos_c1)
  ),
  Group2_Summary = c(
    sprintf("Class III med=%.3f; Class IV med=%.3f; Class V med=%.3f", median(c3_kin), median(c4_ank), median(c5_tmd)),
    sprintf("Non-NLR functional domains (n=%d): med=%.3f", length(non_nlr_dom_w), median(non_nlr_dom_w)),
    sprintf("NB-ARC med=%.3f", median(nlr_nb_w)),
    sprintf("NB-ARC (n=%d): med=%.3f", length(nlr_nb_w), median(nlr_nb_w)),
    sprintf("Non-NLR (n=%d): %.1f +/- %.1f sites", length(non_beb), mean(non_beb), sd(non_beb)/sqrt(length(non_beb))),
    sprintf("Non-NLR (n=%d): %.1f +/- %.1f sites", length(non_meme), mean(non_meme), sd(non_meme)/sqrt(length(non_meme))),
    sprintf("NLR LRR (n=%d): %.1f +/- %.1f%%", length(nlr_lrr_pur), mean(nlr_lrr_pur), sd(nlr_lrr_pur)/sqrt(length(nlr_lrr_pur))),
    sprintf("Yr15 K1_active w=%.4f (n=%s codons)", as.numeric(s10d$module_1_omega_M0[s10d$reference_query=="Yr15"]), s10d$module_1_codons[s10d$reference_query=="Yr15"]),
    sprintf("Sr60 K2_active w=%.4f (n=%s codons)", as.numeric(s10d$module_1_omega_M0[s10d$reference_query=="Sr60"]), s10d$module_1_codons[s10d$reference_query=="Sr60"]),
    sprintf("Sr62 K2_pseudokinase w=%.4f (n=%s codons)", as.numeric(s10d$module_2_omega_M0[s10d$reference_query=="Sr62"]), s10d$module_2_codons[s10d$reference_query=="Sr62"]),
    sprintf("Background concave: %d/%d (%.1f%%)", as.integer(c1_s10e$concave_face_residues[match("Sr35", c1_s10e$gene)]), as.integer(c1_s10e$model_residues[match("Sr35", c1_s10e$gene)]), 100*as.integer(c1_s10e$concave_face_residues[match("Sr35", c1_s10e$gene)])/as.integer(c1_s10e$model_residues[match("Sr35", c1_s10e$gene)])),
    sprintf("Background concave: %d/%d (%.1f%%)", as.integer(c1_s10e$concave_face_residues[match("Yr27", c1_s10e$gene)]), as.integer(c1_s10e$model_residues[match("Yr27", c1_s10e$gene)]), 100*as.integer(c1_s10e$concave_face_residues[match("Yr27", c1_s10e$gene)])/as.integer(c1_s10e$model_residues[match("Yr27", c1_s10e$gene)])),
    sprintf("Background concave: %d/%d (%.1f%%)", as.integer(c1_s10e$concave_face_residues[match("Sr33", c1_s10e$gene)]), as.integer(c1_s10e$model_residues[match("Sr33", c1_s10e$gene)]), 100*as.integer(c1_s10e$concave_face_residues[match("Sr33", c1_s10e$gene)])/as.integer(c1_s10e$model_residues[match("Sr33", c1_s10e$gene)])),
    sprintf("Background concave: %d/%d (%.1f%%)", as.integer(c1_s10e$concave_face_residues[match("Sr22", c1_s10e$gene)]), as.integer(c1_s10e$model_residues[match("Sr22", c1_s10e$gene)]), 100*as.integer(c1_s10e$concave_face_residues[match("Sr22", c1_s10e$gene)])/as.integer(c1_s10e$model_residues[match("Sr22", c1_s10e$gene)])),
    sprintf("Background concave: %d/%d (%.2f%%)", tot_conc_res_c1, tot_res_c1, 100*tot_conc_res_c1/tot_res_c1)
  ),
  Statistic = c(
    sprintf("H = %.2f (df = 4)", kw_res$statistic),
    sprintf("U = %.1f", mwu_dom_asymp$statistic),
    sprintf("H = %.2f (df = 2)", kw_intra_nlr$statistic),
    sprintf("U = %.1f", mwu_lrr_nbarc$statistic),
    sprintf("U = %.1f", mwu_beb$statistic),
    sprintf("U = %.1f", mwu_meme$statistic),
    sprintf("U = %.1f", mwu_dom_pur$statistic),
    sprintf("U = %s", s10d$site_MWU_U[s10d$reference_query=="Yr15"]),
    sprintf("U = %s", s10d$site_MWU_U[s10d$reference_query=="Sr60"]),
    sprintf("U = %s", s10d$site_MWU_U[s10d$reference_query=="Sr62"]),
    sprintf("OR = %.2f", fisher_ors[match("Sr35", c1_s10e$gene)]),
    sprintf("OR = %.2f", fisher_ors[match("Yr27", c1_s10e$gene)]),
    sprintf("OR = %.2f", fisher_ors[match("Sr33", c1_s10e$gene)]),
    sprintf("OR = %.2f", fisher_ors[match("Sr22", c1_s10e$gene)]),
    sprintf("Pooled OR = %.2f", or_pool_mle)
  ),
  P_value = c(
    sprintf("%.4e", kw_res$p.value),
    sprintf("%.4e", mwu_dom_asymp$p.value),
    sprintf("%.4e", kw_intra_nlr$p.value),
    sprintf("%.4e", mwu_lrr_nbarc$p.value),
    sprintf("%.4e", mwu_beb$p.value),
    sprintf("%.4e", mwu_meme$p.value),
    sprintf("%.4e", mwu_dom_pur$p.value),
    sprintf("%.4e", as.numeric(s10d$site_MWU_P[s10d$reference_query=="Yr15"])),
    sprintf("%.4e", as.numeric(s10d$site_MWU_P[s10d$reference_query=="Sr60"])),
    sprintf("%.4e", as.numeric(s10d$site_MWU_P[s10d$reference_query=="Sr62"])),
    sprintf("%.4e", fisher_pvals[match("Sr35", c1_s10e$gene)]),
    sprintf("%.4e", fisher_pvals[match("Yr27", c1_s10e$gene)]),
    sprintf("%.4e", fisher_pvals[match("Sr33", c1_s10e$gene)]),
    sprintf("%.4e", fisher_pvals[match("Sr22", c1_s10e$gene)]),
    sprintf("%.4e", ft_pool$p.value)
  )
)

write.table(audit_tbl, "results/Supplementary_Table_S10f_R_statistical_audit_3.4.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
dir.create("tables", showWarnings = FALSE)
write.table(audit_tbl, "tables/Supplementary_Table_S10f_R_Statistical_Audit_3.4.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------------------------
# Generate Publication-Grade 4-Panel Figure 4 in Base R
# ------------------------------------------------------------------------------
draw_fig4 <- function() {
  layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
  par(oma = c(0.8, 0.8, 2.5, 0.8))

  cols <- c("#D55E00", "#E69F00", "#0072B2", "#009E73", "#CC79A7")
  cols_alpha <- adjustcolor(cols, alpha.f = 0.65)
  class_names <- c("Class I\n(CC-NLR LRR)", "Class II\n(NLR-ID LRR)",
                   "Class III\n(TKP Kinase)", "Class IV\n(ANK-TM)", "Class V\n(APR TMD)")

  # ----------------------------------------------------------------------------
  # Panel a: Domain-Partitioned dN/dS (omega) Across Structural Classes I-V
  # ----------------------------------------------------------------------------
  par(mar = c(4.5, 4.8, 3.5, 1.5))
  omega_list <- list(c1_lrr, c2_lrr, c3_kin, c4_ank, c5_tmd)
  boxplot(omega_list, names = NA, col = cols_alpha, border = "grey15",
          ylim = c(0, 1.25), ylab = expression(bold("Domain-Partitioned " * omega * " (" * italic(d)[N] / italic(d)[S] * ")")),
          las = 1, boxwex = 0.50, outline = FALSE, cex.lab = 1.05, bty = "l")
  abline(h = 1.0, lty = 2, col = "grey45", lwd = 1.2)
  grid(nx = NA, ny = NULL, col = "grey88", lty = 3)
  axis(1, at = 1:5, labels = class_names, tick = TRUE, padj = 0.5, cex.axis = 0.88, font = 2)

  set.seed(42)
  for (i in 1:5) {
    pts <- omega_list[[i]]
    jit <- runif(length(pts), -0.14, 0.14)
    points(rep(i, length(pts)) + jit, pts, pch = 21, bg = cols[i], col = "grey10", cex = 1.15, lwd = 0.8)
    text(i, max(pts) + 0.05, sprintf("med=%.3f", median(pts)), cex = 0.80, font = 2, col = "grey15")
  }
  text(3.1, 1.18, sprintf("Kruskal-Wallis (n=39 domains) H = %.2f, P = %.2e\nNLR LRR (n=20 dom.) vs Non-NLR core (n=19 dom.): U = %.0f, P = %.2e",
                          kw_res$statistic, kw_res$p.value, mwu_dom_asymp$statistic, mwu_dom_asymp$p.value),
       cex = 0.85, font = 3, col = "grey10")
  mtext("a  Domain-Partitioned Selection Ratios Across Classes I-V", side = 3, line = 1.0, adj = 0, font = 2, cex = 0.98)

  # ----------------------------------------------------------------------------
  # Panel b: Pervasive Purifying Selection Constraint (FUBAR) & PAML M8 BEB Sites
  # ----------------------------------------------------------------------------
  par(mar = c(4.5, 4.8, 3.5, 1.5))
  pur_list <- list(
    as.numeric(s10c$pct_fubar_purifying[s10c$structural_class == "I" & s10c$domain_name == "LRR"]),
    as.numeric(s10c$pct_fubar_purifying[s10c$structural_class %in% c("I/II", "II") & s10c$domain_name == "LRR"]),
    as.numeric(s10c$pct_fubar_purifying[s10c$structural_class == "III" & s10c$domain_type == "kinase_module"]),
    as.numeric(s10c$pct_fubar_purifying[s10c$structural_class == "IV" & s10c$domain_name != "whole_gene"]),
    as.numeric(s10c$pct_fubar_purifying[s10c$structural_class == "V" & s10c$domain_name != "whole_gene"])
  )
  pur_means <- sapply(pur_list, mean)
  pur_sems  <- sapply(pur_list, function(x) if (length(x) > 1) sd(x)/sqrt(length(x)) else 0)

  bp <- barplot(pur_means, names.arg = NA, col = cols_alpha, border = "grey15",
                ylim = c(0, 116), ylab = expression(bold("FUBAR Purifying Codons (%) (" * hat(alpha) > hat(beta) * ")")),
                las = 1, width = 0.65, space = 0.35, cex.lab = 1.05)
  grid(nx = NA, ny = NULL, col = "grey88", lty = 3)
  axis(1, at = bp, labels = class_names, tick = TRUE, padj = 0.5, cex.axis = 0.88, font = 2)

  set.seed(42)
  for (i in 1:5) {
    if (pur_sems[i] > 0) {
      arrows(bp[i], pur_means[i] - pur_sems[i], bp[i], pur_means[i] + pur_sems[i],
             angle = 90, code = 3, length = 0.05, lwd = 1.3, col = "grey15")
    }
    pts <- pur_list[[i]]
    jit <- runif(length(pts), -0.12, 0.12)
    points(rep(bp[i], length(pts)) + jit, pts, pch = 21, bg = cols[i], col = "grey10", cex = 1.0)
    y_lbl <- max(pur_means[i] + pur_sems[i], max(pts)) + 3.8
    text(bp[i], y_lbl, sprintf("%.1f%%", pur_means[i]), cex = 0.85, font = 2, col = "grey15")
  }
  text(mean(bp), 110, sprintf("Non-NLR whole-gene purifying: %.1f +/- %.1f%% codons | NLR M8 BEB: %.1f +/- %.1f sites/OG",
                              mean(non_pur_agb), sd(non_pur_agb)/sqrt(length(non_pur_agb)),
                              mean(nlr_beb), sd(nlr_beb)/sqrt(length(nlr_beb))),
       cex = 0.82, font = 3, col = "grey10")
  mtext("b  Pervasive Purifying Selection Constraint vs Positive Selection", side = 3, line = 1.0, adj = 0, font = 2, cex = 0.96)

  # ----------------------------------------------------------------------------
  # Panel c: Class III Tandem Kinase Paired Active vs Pseudokinase Decoy Asymmetry
  # ----------------------------------------------------------------------------
  par(mar = c(4.5, 4.8, 3.5, 1.5))
  mod1_w <- as.numeric(s10d$module_1_omega_M0)
  mod2_w <- as.numeric(s10d$module_2_omega_M0)
  plot(c(0.40, 2.60), c(0.08, 0.98), type = "n", xaxt = "n", las = 1, bty = "l",
       xlab = "", ylab = expression(bold("Kinase Module " * omega * " (" * italic(d)[N] / italic(d)[S] * ")")),
       cex.lab = 1.05)
  grid(nx = NA, ny = NULL, col = "grey88", lty = 3)
  axis(1, at = c(1, 2), labels = c("Primary Module\n(Active / Conserved)", "Secondary Module\n(Pseudokinase / Accessory)"),
       tick = TRUE, padj = 0.5, cex.axis = 0.92, font = 2)

  # Custom y-offsets for labels to prevent collision
  left_y_offsets <- c(Yr15 = 0.025, Lr9 = -0.025, Sr60 = 0.020, Sr43 = -0.022, Sr62 = 0.0)
  right_y_offsets <- c(Yr15 = 0.022, Lr9 = -0.022, Sr60 = 0.0, Sr43 = 0.0, Sr62 = 0.0)

  for (i in seq_len(nrow(s10d))) {
    gene <- s10d$reference_query[i]
    col_line <- if (gene == "Sr62") "#0072B2" else "#D55E00"
    segments(1, mod1_w[i], 2, mod2_w[i], col = col_line, lwd = 2.2)
    points(1, mod1_w[i], pch = 21, bg = "#0072B2", col = "grey15", cex = 1.35)
    points(2, mod2_w[i], pch = 21, bg = "#D55E00", col = "grey15", cex = 1.35)
    
    y1 <- mod1_w[i] + left_y_offsets[gene]
    y2 <- mod2_w[i] + right_y_offsets[gene]
    text(1 - 0.08, y1, gene, adj = 1, cex = 0.85, font = 2, col = "grey15")
    text(2 + 0.08, y2, sprintf("%s: %.3f", gene, mod2_w[i]), adj = 0, cex = 0.82, font = 2, col = "grey15")
  }

  sr62_row <- s10d[s10d$reference_query == "Sr62", ][1, ]
  text(1.5, 0.94, sprintf("Module asymmetry is locus-specific:\nYr15, Lr9, Sr43 (+modest Sr60): Pseudokinase/Accessory w elevated\nSr62: Inverse polarity (K1_active w=%.3f > K2_pseudo w=%.3f)",
                           as.numeric(sr62_row$module_1_omega_M0), as.numeric(sr62_row$module_2_omega_M0)),
       cex = 0.82, font = 3, col = "grey10")
  mtext("c  Class III Tandem Kinase Module Selection Asymmetry", side = 3, line = 1.0, adj = 0, font = 2, cex = 0.96)

  # ----------------------------------------------------------------------------
  # Panel d: AlphaFold 3D Surface Localization of Positively Selected Sites (S10e)
  # ----------------------------------------------------------------------------
  par(mar = c(4.5, 6.2, 3.5, 1.5))
  ord <- rev(seq_len(nrow(s10e)))
  s10e_o <- s10e[ord, ]
  conc_pct <- as.numeric(s10e_o$selected_on_concave_face) / pmax(1, as.numeric(s10e_o$n_selected_sites_BEB95_or_MEME05)) * 100
  conc_pct[is.na(conc_pct)] <- 0
  other_pct <- 100 - conc_pct
  mat_d <- rbind(conc_pct, other_pct)
  gene_lbls <- sprintf("%s (Class %s)", s10e_o$gene, s10e_o$structural_class)

  bp_d <- barplot(mat_d, horiz = TRUE, names.arg = gene_lbls, las = 1,
                  col = c("#D55E00", "#56B4E9"), border = "grey15",
                  xlim = c(0, 215), ylim = c(0, 27.5),
                  xlab = expression(bold("Positively Selected Sites Localization (%)")),
                  cex.names = 0.78, cex.lab = 1.05, font.axis = 2)
  grid(nx = NULL, ny = NA, col = "grey88", lty = 3)

  for (i in seq_len(nrow(s10e_o))) {
    r <- s10e_o[i, ]
    n_sel <- as.integer(r$n_selected_sites_BEB95_or_MEME05)
    n_cnc <- as.integer(r$selected_on_concave_face)
    if (is.na(n_cnc)) n_cnc <- 0
    lbl <- if (r$structural_class == "I") {
      p_fmt <- if (as.numeric(r$concave_P_vs_whole_protein) < 1e-4) sprintf("%.1e", as.numeric(r$concave_P_vs_whole_protein)) else sprintf("%.3f", as.numeric(r$concave_P_vs_whole_protein))
      sprintf("%d/%d concave (%.1f%%) | OR = %.1f, P = %s",
              n_cnc, n_sel, conc_pct[i], as.numeric(r$concave_OR_vs_whole_protein), p_fmt)
    } else if (r$structural_class == "II") {
      sprintf("0/%d concave | %d BED/ID loop sites (%s%% exp.)",
              n_sel, n_sel, r$pct_exposed_rsa25)
    } else if (r$structural_class == "III") {
      sprintf("0/%d concave | %d inter-kinase linker/loop sites",
              n_sel, n_sel)
    } else {
      sprintf("0/%d concave | %d pore-lining/gating residue",
              n_sel, n_sel)
    }
    text(103, bp_d[i], lbl, adj = 0, cex = 0.72, font = 2, col = "grey15")
  }
  legend(0, 27.2, legend = c("LRR Concave Beta-Sheet (Solvent-Exposed)", "Other Surface / Loop / Pore-Lining Sites"),
         fill = c("#D55E00", "#56B4E9"), border = "grey15", cex = 0.80, horiz = TRUE, bty = "n")
  mtext(sprintf("d  AlphaFold 3D Surface Localization (Pooled OR = %.2f, P = %.2e)", or_pool_mle, ft_pool$p.value),
        side = 3, line = 1.0, adj = 0, font = 2, cex = 0.96)

  mtext("Figure 4 | Contrasting Purifying and Diversifying Selection Regimes Across 3D Protein Architectures",
        outer = TRUE, side = 3, line = 0.6, font = 2, cex = 1.15)
}

pdf("figures/Fig4_selection_structures_R.pdf", width = 14, height = 10.5)
draw_fig4()
invisible(dev.off())

png("figures/Fig4_selection_structures_R.png", width = 14, height = 10.5, units = "in", res = 300)
draw_fig4()
invisible(dev.off())

cat("Saved R statistical audit table -> results/Supplementary_Table_S10f_R_statistical_audit_3.4.tsv\n")
cat("Saved R publication figures -> figures/Fig4_selection_structures_R.pdf and .png\n")
cat("======================================================================\n")
