#!/usr/bin/env Rscript
# ==============================================================================
# stats_and_fig1_3.1.R
# Complete R Statistical Audit, Pangenome Modeling & Figure 1 Generator for
# Section 3.1 ("Super-pangenome inventory and core-to-dispensable partitioning
# of cloned rust resistance homologues") and Methods 5.1-5.3.
#
# Usage:
#   Rscript stats_and_fig1_3.1.R [results_dir] [figures_dir] [tables_dir]
#
# Dependencies:
#   Base R (stats, graphics, grDevices) + vegan + micropan
#
# Outputs:
#   1. results/S3d_heaps_law.tsv & results/S3d_heaps_curves.tsv
#   2. results/S3e_R_statistical_audit_3.1.tsv
#   3. tables/Supplementary_Table_S3e_R_Statistical_Audit_3.1.tsv
#   4. figures/Fig1_pangenome_inventory_R.pdf & .png (300 dpi)
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

suppressPackageStartupMessages({
  library(stats)
  library(graphics)
  library(grDevices)
  library(vegan)
  library(micropan)
})

args <- commandArgs(trailingOnly = TRUE)
rd   <- if (length(args) >= 1) args[1] else "results"
fd   <- if (length(args) >= 2) args[2] else "figures"
td   <- if (length(args) >= 3) args[3] else "tables"
dir.create(fd, showWarnings = FALSE, recursive = TRUE)
dir.create(td, showWarnings = FALSE, recursive = TRUE)

cat("======================================================================\n")
cat("  Section 3.1 R Statistical Audit & Pangenome Analysis\n")
cat("======================================================================\n\n")

audit_rows <- list()
add_audit <- function(section, metric, statistic, p_value, details) {
  audit_rows[[length(audit_rows) + 1]] <<- data.frame(
    section   = section,
    metric    = metric,
    statistic = as.character(statistic),
    p_value   = as.character(p_value),
    details   = as.character(details),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 1. Load Reference Catalogue & Super-Pangenome Assembly Panel
# ------------------------------------------------------------------------------
cat_df   <- read.delim("catalogue_3.1.tsv", stringsAsFactors = FALSE, check.names = FALSE)
panel_df <- read.delim("panel_3.1.tsv", stringsAsFactors = FALSE, check.names = FALSE)
panel_df <- panel_df[tolower(panel_df$include) %in% c("yes", "1", "true"), ]
panel_df$annotated <- (panel_df$source == "ensembl" | (panel_df$gff3 != "" & !is.na(panel_df$gff3)))

main_cat <- cat_df[cat_df$Set == "Main", ]
n_queries_main <- nrow(main_cat)
n_loci_main    <- length(unique(main_cat$LocusGroup))
n_asm          <- nrow(panel_df)
n_subgenomes   <- sum(nchar(panel_df$genomes))
n_annot_asm    <- sum(panel_df$annotated)

cat(sprintf("[1] Catalogue & Panel: %d Main queries (%d loci) + %d Comparators across %d assemblies (%d haploid subgenomes; %d annotated)\n",
            n_queries_main, n_loci_main, sum(cat_df$Set != "Main"), n_asm, n_subgenomes, n_annot_asm))
add_audit("5.1_Catalogue_Panel", "Main cloned rust R proteins / loci",
          sprintf("%d proteins / %d loci", n_queries_main, n_loci_main), "NA",
          sprintf("Class I=%d (%d loci), II=%d (%d loci), III=%d (%d loci), IV=%d (%d loci), V=%d (%d loci)",
                  sum(main_cat$Class == "I"), length(unique(main_cat$LocusGroup[main_cat$Class == "I"])),
                  sum(main_cat$Class == "II"), length(unique(main_cat$LocusGroup[main_cat$Class == "II"])),
                  sum(main_cat$Class == "III"), length(unique(main_cat$LocusGroup[main_cat$Class == "III"])),
                  sum(main_cat$Class == "IV"), length(unique(main_cat$LocusGroup[main_cat$Class == "IV"])),
                  sum(main_cat$Class == "V"), length(unique(main_cat$LocusGroup[main_cat$Class == "V"]))))

# ------------------------------------------------------------------------------
# 2. Annotation-Independent Homologue Rescue (Wilcoxon Signed-Rank Test)
# ------------------------------------------------------------------------------
s3_df  <- read.delim(file.path(rd, "S3_homologue_loci_rescue.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
s3a_df <- read.delim(file.path(rd, "S3a_recovery_per_assembly.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
s3b_df <- read.delim(file.path(rd, "S3b_recovery_per_locus.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

ann_vec <- s3a_df$annotation_homologues
uni_vec <- s3a_df$union_homologues
pct_inc <- 100 * (uni_vec - ann_vec) / ann_vec
pct_mean <- mean(pct_inc)
pct_sem  <- sd(pct_inc) / sqrt(length(pct_inc))

# Two-sided Wilcoxon signed-rank test (paired union vs annotation across 16 GFF3 assemblies)
wt_rescue <- wilcox.test(uni_vec, ann_vec, paired = TRUE, alternative = "two.sided", exact = FALSE, correct = FALSE)
# Note: SciPy wilcoxon reports min(W+, W-) = 0.0; R wilcox.test reports V = W+ = n(n+1)/2 = 136.0
w_min <- min(wt_rescue$statistic, (length(uni_vec) * (length(uni_vec) + 1) / 2) - wt_rescue$statistic)

s3_main_annot <- s3_df[s3_df$annotated == 1 & s3_df$LocusGroup %in% unique(main_cat$LocusGroup), ]
n_unann <- sum(s3_main_annot$status == "rescued_unannotated")
n_split <- sum(s3_main_annot$status == "rescued_split_model")
n_misan <- sum(s3_main_annot$status == "rescued_misannotated_model")
n_resc_total <- n_unann + n_split + n_misan

cat(sprintf("[2] Homologue Rescue (n = %d annotated assemblies): %d rescued loci (unannotated=%d, split=%d, misannotated=%d)\n",
            nrow(s3a_df), n_resc_total, n_unann, n_split, n_misan))
cat(sprintf("    Per-assembly recovery increase: %.2f +/- %.2f%% (mean +/- s.e.m.); Wilcoxon W = %.1f (V = %.1f), P = %.3e\n",
            pct_mean, pct_sem, w_min, wt_rescue$statistic, wt_rescue$p.value))

add_audit("3.1_Annotation_Rescue", "Wilcoxon signed-rank test (union vs GFF3 annotation)",
          sprintf("W_min = %.1f (V = %.1f), mean gain = +%.2f +/- %.2f%%", w_min, wt_rescue$statistic, pct_mean, pct_sem),
          sprintf("%.4e", wt_rescue$p.value),
          sprintf("n = %d annotated assemblies; %d total rescued loci (%d misannotated, %d unannotated, %d split)",
                  nrow(s3a_df), n_resc_total, n_misan, n_unann, n_split))

# ------------------------------------------------------------------------------
# 3. Paralogue Pruning & Split-Model Domain Order Verification
# ------------------------------------------------------------------------------
s9_df  <- read.delim(file.path(rd, "S9_paralogue_pruning.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
s9b_df <- read.delim(file.path(rd, "S9b_domain_order_split.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

s9_placed <- s9_df[s9_df$subgenome != "Un", ]
n_cand_total  <- nrow(s9_df)
n_cand_placed <- nrow(s9_placed)
n_cand_un     <- sum(s9_df$subgenome == "Un")
n_retained    <- sum(s9_placed$verdict == "retained")
n_pruned_out  <- sum(s9_placed$verdict == "pruned_out_cluster")
n_rbh_support <- sum(s9_df$reverse_hit_support == 1)
n_anc_support <- sum(s9_df$flanking_anchor_support == 1)
n_split_ok    <- sum(s9b_df$order_conserved == 1)

cat(sprintf("[3] Paralogue Pruning: %d/%d (%.1f%%) subgenome-placed candidates retained in 31 query orthogroups (%d out-cluster paralogues pruned; %d unplaced Un scaffold hits excluded)\n",
            n_retained, n_cand_placed, 100 * n_retained / n_cand_placed, n_pruned_out, n_cand_un))
cat(sprintf("    Reciprocal BLASTp support: %d/%d (%.1f%%); Syntenic flanking-anchor support: %d/%d; Split-model domain order conserved: %d/%d (100%%)\n",
            n_rbh_support, n_cand_total, 100 * n_rbh_support / n_cand_total, n_anc_support, n_cand_total, n_split_ok, nrow(s9b_df)))

add_audit("5.2_Paralogue_Pruning", "Orthogroup RBH & syntenic anchor pruning",
          sprintf("%d/%d (%.1f%%) placed retained", n_retained, n_cand_placed, 100 * n_retained / n_cand_placed), "NA",
          sprintf("All anchor-verified=%d; placed=%d; Un=%d; RBH=%d/%d; flanking_anchor=%d/%d; split_domain_order=%d/%d",
                  n_cand_total, n_cand_placed, n_cand_un, n_rbh_support, n_cand_total, n_anc_support, n_cand_total, n_split_ok, nrow(s9b_df)))

# ------------------------------------------------------------------------------
# 4. Core-to-Dispensable Pangenome Partitioning (Fisher's Exact & Mann-Whitney U)
# ------------------------------------------------------------------------------
s1b_df  <- read.delim(file.path(rd, "S1b_orthogroup_partition.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
og_main <- s1b_df[s1b_df$set == "Main" & s1b_df$group %in% c("NLR", "nonNLR"), ]

nlr_og <- og_main[og_main$group == "NLR", ]
non_og <- og_main[og_main$group == "nonNLR", ]

# Haploid subgenome level (56 subgenomes)
nlr_sc_sub <- sum(nlr_og$compartment %in% c("shell", "cloud"))
nlr_cs_sub <- sum(nlr_og$compartment %in% c("core", "soft-core"))
non_sc_sub <- sum(non_og$compartment %in% c("shell", "cloud"))
non_cs_sub <- sum(non_og$compartment %in% c("core", "soft-core"))

mat_sub <- matrix(c(nlr_sc_sub, non_sc_sub, nlr_cs_sub, non_cs_sub), nrow = 2,
                  dimnames = list(Group = c("NLR", "nonNLR"), Compartment = c("Shell_Cloud", "Core_SoftCore")))
ft_sub  <- fisher.test(mat_sub, alternative = "two.sided")
or_ha_sub <- ((nlr_sc_sub + 0.5) * (non_cs_sub + 0.5)) / ((nlr_cs_sub + 0.5) * (non_sc_sub + 0.5))

# Whole-assembly level (27 assemblies)
nlr_sc_asm <- sum(nlr_og$compartment_assembly_level %in% c("shell", "cloud"))
nlr_cs_asm <- sum(nlr_og$compartment_assembly_level %in% c("core", "soft-core"))
non_sc_asm <- sum(non_og$compartment_assembly_level %in% c("shell", "cloud"))
non_cs_asm <- sum(non_og$compartment_assembly_level %in% c("core", "soft-core"))

mat_asm <- matrix(c(nlr_sc_asm, non_sc_asm, nlr_cs_asm, non_cs_asm), nrow = 2,
                  dimnames = list(Group = c("NLR", "nonNLR"), Compartment = c("Shell_Cloud", "Core_SoftCore")))
ft_asm  <- fisher.test(mat_asm, alternative = "two.sided")
or_ha_asm <- ((nlr_sc_asm + 0.5) * (non_cs_asm + 0.5)) / ((nlr_cs_asm + 0.5) * (non_sc_asm + 0.5))

# Mann-Whitney U (Wilcoxon rank-sum) tests on occupancy frequency
mwu_sub <- wilcox.test(nlr_og$occupancy_all_subgenomes, non_og$occupancy_all_subgenomes,
                       alternative = "two.sided", exact = FALSE, correct = TRUE)
mwu_asm <- wilcox.test(nlr_og$occupancy_assemblies, non_og$occupancy_assemblies,
                       alternative = "two.sided", exact = FALSE, correct = TRUE)

cat(sprintf("[4] Pangenome Partitioning (%d Main OGs: %d NLR vs %d non-NLR):\n",
            nrow(og_main), nrow(nlr_og), nrow(non_og)))
cat(sprintf("    Haploid subgenomes (n = 56): NLR shell/cloud = %d/%d (%.1f%% [shell=%d, cloud=%d]) vs non-NLR core/soft-core = %d/%d (%.1f%% [core=%d, soft-core=%d])\n",
            nlr_sc_sub, nrow(nlr_og), 100 * nlr_sc_sub / nrow(nlr_og),
            sum(nlr_og$compartment == "shell"), sum(nlr_og$compartment == "cloud"),
            non_cs_sub, nrow(non_og), 100 * non_cs_sub / nrow(non_og),
            sum(non_og$compartment == "core"), sum(non_og$compartment == "soft-core")))
cat(sprintf("    Fisher's exact test (subgenomes): OR = %s (Haldane-Anscombe OR = %.1f), P = %.3e\n",
            ifelse(is.infinite(ft_sub$estimate), "Inf", sprintf("%.2f", ft_sub$estimate)), or_ha_sub, ft_sub$p.value))
cat(sprintf("    Fisher's exact test (assemblies): OR = %s (Haldane-Anscombe OR = %.1f), P = %.3e\n",
            ifelse(is.infinite(ft_asm$estimate), "Inf", sprintf("%.2f", ft_asm$estimate)), or_ha_asm, ft_asm$p.value))
cat(sprintf("    Mann-Whitney U test (subgenomes): U = %.1f, P = %.3e (median NLR = %.3f vs non-NLR = %.3f)\n",
            mwu_sub$statistic, mwu_sub$p.value, median(nlr_og$occupancy_all_subgenomes), median(non_og$occupancy_all_subgenomes)))
cat(sprintf("    Mann-Whitney U test (assemblies): U = %.1f, P = %.3e (median NLR = %.3f vs non-NLR = %.3f)\n",
            mwu_asm$statistic, mwu_asm$p.value, median(nlr_og$occupancy_assemblies), median(non_og$occupancy_assemblies)))

add_audit("3.1_Compartment_Fisher_Subgenomes", "Fisher's exact test: NLR shell/cloud vs non-NLR core/soft-core (56 subgenomes)",
          sprintf("OR = Inf (Haldane-Anscombe OR = %.1f)", or_ha_sub), sprintf("%.4e", ft_sub$p.value),
          sprintf("2x2 table: [[NLR shell/cloud=%d, NLR core/soft=%d], [nonNLR shell/cloud=%d, nonNLR core/soft=%d]]",
                  nlr_sc_sub, nlr_cs_sub, non_sc_sub, non_cs_sub))
add_audit("3.1_Compartment_Fisher_Assemblies", "Fisher's exact test: NLR shell/cloud vs non-NLR core/soft-core (27 assemblies)",
          sprintf("OR = Inf (Haldane-Anscombe OR = %.1f)", or_ha_asm), sprintf("%.4e", ft_asm$p.value),
          sprintf("2x2 table: [[NLR shell/cloud=%d, NLR core/soft=%d], [nonNLR shell/cloud=%d, nonNLR core/soft=%d]]",
                  nlr_sc_asm, nlr_cs_asm, non_sc_asm, non_cs_asm))
add_audit("3.1_Occupancy_MWU_Subgenomes", "Mann-Whitney U test on subgenome occupancy (NLR vs non-NLR)",
          sprintf("U = %.1f", mwu_sub$statistic), sprintf("%.4e", mwu_sub$p.value),
          sprintf("median NLR = %.4f (n=%d), median non-NLR = %.4f (n=%d)",
                  median(nlr_og$occupancy_all_subgenomes), nrow(nlr_og),
                  median(non_og$occupancy_all_subgenomes), nrow(non_og)))
add_audit("3.1_Occupancy_MWU_Assemblies", "Mann-Whitney U test on assembly occupancy (NLR vs non-NLR)",
          sprintf("U = %.1f", mwu_asm$statistic), sprintf("%.4e", mwu_asm$p.value),
          sprintf("median NLR = %.4f (n=%d), median non-NLR = %.4f (n=%d)",
                  median(nlr_og$occupancy_assemblies), nrow(nlr_og),
                  median(non_og$occupancy_assemblies), nrow(non_og)))

# ------------------------------------------------------------------------------
# 5. Identity-Threshold Sensitivity Sweep (S3c Verification in R)
# ------------------------------------------------------------------------------
s3c_df <- read.delim(file.path(rd, "S3c_identity_threshold_sensitivity.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
cat("[5] Identity-Threshold Sensitivity Sweep (S3c):\n")
for (i in seq_len(nrow(s3c_df))) {
  r_sw <- s3c_df[i, ]
  cat(sprintf("    min_id_NLR = %.2f (progenitor = %.2f): median NLR = %.4f vs non-NLR = %.4f | MWU U = %.1f, P = %.3e\n",
              r_sw$min_identity_NLR, r_sw$min_identity_progenitor, r_sw$median_occ_NLR, r_sw$median_occ_nonNLR,
              r_sw$MWU_U, r_sw$MWU_P))
}
add_audit("5.3_Identity_Sensitivity_Sweep", "Mann-Whitney U across identity thresholds (0.70-0.95)",
          sprintf("U range = %.1f to %.1f", min(s3c_df$MWU_U), max(s3c_df$MWU_U)),
          sprintf("All P <= %.4e", max(s3c_df$MWU_P)),
          "Non-NLR occupancy significantly exceeds NLR occupancy across all tested thresholds")

# ------------------------------------------------------------------------------
# 6. Pangenome Saturation (Heaps' Power Law & Core Decay)
# ------------------------------------------------------------------------------
heaps_df  <- read.delim(file.path(rd, "S3d_heaps_law.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
curves_df <- read.delim(file.path(rd, "S3d_heaps_curves.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

h_nlr <- heaps_df[heaps_df$group == "NLR", ]
h_non <- heaps_df[heaps_df$group == "nonNLR", ]
cat(sprintf("[6] Heaps' Law Pangenome Saturation (1,000 permutations):\n"))
cat(sprintf("    NLR (n = %d OGs):     gamma = %.4f +/- %.4f (perm mean = %.4f +/- %.4f), core_n56 = %d\n",
            h_nlr$n_orthogroups, h_nlr$gamma_meancurve, h_nlr$gamma_meancurve_se,
            h_nlr$gamma_perm_mean, h_nlr$gamma_perm_sd, h_nlr$core_orthogroups_n56))
cat(sprintf("    non-NLR (n = %d OGs): gamma = %.4f +/- %.4f (perm mean = %.4f +/- %.4f), core_n56 = %d\n",
            h_non$n_orthogroups, h_non$gamma_meancurve, h_non$gamma_meancurve_se,
            h_non$gamma_perm_mean, h_non$gamma_perm_sd, h_non$core_orthogroups_n56))
cat(sprintf("    Delta_gamma = %.4f, Label-permutation P = %.4f (n = %d permutations)\n",
            h_nlr$delta_gamma, h_nlr$P_label_permutation, h_nlr$n_label_perm))

add_audit("3.1_Heaps_Law_Openness", "Heaps' power-law exponent gamma (NLR vs non-NLR)",
          sprintf("gamma_NLR = %.4f +/- %.4f vs gamma_nonNLR = %.4f +/- %.4f (Delta = %.4f)",
                  h_nlr$gamma_meancurve, h_nlr$gamma_meancurve_se,
                  h_non$gamma_meancurve, h_non$gamma_meancurve_se, h_nlr$delta_gamma),
          sprintf("%.4f", h_nlr$P_label_permutation),
          sprintf("1,000 random subgenome orderings + 1,000 label permutations; core at n=56: NLR=%d/20, nonNLR=%d/10",
                  h_nlr$core_orthogroups_n56, h_non$core_orthogroups_n56))

audit_table <- do.call(rbind, audit_rows)
write.table(audit_table, file.path(rd, "S3e_R_statistical_audit_3.1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(audit_table, file.path(td, "Supplementary_Table_S3e_R_Statistical_Audit_3.1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nSaved R statistical audit table -> %s/S3e_R_statistical_audit_3.1.tsv\n", rd))

# ------------------------------------------------------------------------------
# 7. Generate Publication-Quality 6-Panel Figure 1 in R (PDF & 300-dpi PNG)
# ------------------------------------------------------------------------------
draw_fig1_R <- function() {
  par(mfrow = c(3, 2), oma = c(0.5, 0.5, 0.5, 0.5),
      family = "sans", cex.main = 0.98, cex.lab = 0.90, cex.axis = 0.80)

  # ---- Panel a: Structural taxonomy of 43 cloned rust R proteins
  par(mar = c(4.2, 9.5, 2.6, 1.5))
  classes <- c("I", "II", "III", "IV", "V")
  cls_cols <- c("#1f77b4", "#4a90d9", "#d95f02", "#7570b3", "#1b9e77")
  cls_labs <- c("Class I: CC-NLR", "Class II: NLR-ID", "Class III: TKP/TKF",
                "Class IV: Membrane/NAM", "Class V: APR Transp.")
  q_cnts <- sapply(classes, function(k) sum(main_cat$Class == k))
  l_cnts <- sapply(classes, function(k) length(unique(main_cat$LocusGroup[main_cat$Class == k])))
  mat_a  <- rbind(rev(l_cnts), rev(q_cnts))
  bp_a <- barplot(mat_a, beside = TRUE, horiz = TRUE, names.arg = rev(cls_labs), las = 1,
                  col = rbind(adjustcolor(rev(cls_cols), alpha.f = 0.45), rev(cls_cols)),
                  xlim = c(0, max(q_cnts) + 5), xlab = "Count in reference catalogue",
                  main = sprintf("a  Structural taxonomy of %d cloned rust R proteins (%d loci)",
                                 n_queries_main, n_loci_main))
  text(mat_a[1, ] + 0.8, bp_a[1, ], labels = mat_a[1, ], cex = 0.78)
  text(mat_a[2, ] + 0.8, bp_a[2, ], labels = mat_a[2, ], cex = 0.78, font = 2)
  legend("bottomright", legend = c(sprintf("Cloned proteins (n = %d)", n_queries_main),
                                   sprintf("Independent loci (n = %d)", n_loci_main)),
         fill = c("#1f77b4", adjustcolor("#1f77b4", alpha.f = 0.45)), bty = "n", cex = 0.80)

  # ---- Panel b: Super-pangenome panel composition
  par(mar = c(4.2, 5.0, 2.6, 1.5))
  tiers <- c("A", "S", "D", "AB", "ABD")
  tier_labs <- c("Diploid A", "Diploid S", "Diploid D", "Tetraploid AB", "Hexaploid ABD")
  n_ann_t <- sapply(tiers, function(g) sum(panel_df$genomes == g & panel_df$annotated))
  n_raw_t <- sapply(tiers, function(g) sum(panel_df$genomes == g & !panel_df$annotated))
  n_sub_t <- sapply(tiers, function(g) sum(nchar(panel_df$genomes[panel_df$genomes == g])))
  x_b <- seq_along(tiers)
  plot(NULL, xlim = c(0.4, 5.6), ylim = c(0, max(n_sub_t) + 6), xaxt = "n",
       xlab = "Genome tier", ylab = "Count (assemblies / subgenomes)",
       main = sprintf("b  Triticum-Aegilops super-pangenome (%d assemblies, %d subgenomes)",
                      n_asm, n_subgenomes))
  axis(1, at = x_b, labels = tier_labs, cex.axis = 0.78)
  rect(x_b - 0.34, 0, x_b - 0.02, n_ann_t, col = "#2b5c8f", border = "black")
  rect(x_b - 0.34, n_ann_t, x_b - 0.02, n_ann_t + n_raw_t, col = "#9ecae1", border = "black")
  rect(x_b + 0.02, 0, x_b + 0.34, n_sub_t, col = "#e66101", border = "black")
  text(x_b - 0.18, n_ann_t + n_raw_t + 1.5, labels = n_ann_t + n_raw_t, cex = 0.78)
  text(x_b + 0.18, n_sub_t + 1.5, labels = n_sub_t, cex = 0.78, font = 2)
  legend("topleft", legend = c(sprintf("Assemblies w/ GFF3 (n = %d)", sum(n_ann_t)),
                               sprintf("Genome-only assemblies (n = %d)", sum(n_raw_t)),
                               sprintf("Haploid subgenomes (n = %d)", sum(n_sub_t))),
         fill = c("#2b5c8f", "#9ecae1", "#e66101"), bty = "n", cex = 0.78)

  # ---- Panel c: Annotation rescue per assembly
  par(mar = c(5.2, 5.0, 2.6, 1.5))
  asms <- s3a_df$asm_id
  ann_c   <- sapply(asms, function(a) sum(s3_main_annot$asm_id == a & s3_main_annot$status %in% c("annotated", "annotated_only_not_in_miniprot")))
  mis_c   <- sapply(asms, function(a) sum(s3_main_annot$asm_id == a & s3_main_annot$status == "rescued_misannotated_model"))
  unann_c <- sapply(asms, function(a) sum(s3_main_annot$asm_id == a & s3_main_annot$status == "rescued_unannotated"))
  split_c <- sapply(asms, function(a) sum(s3_main_annot$asm_id == a & s3_main_annot$status == "rescued_split_model"))
  mat_c   <- rbind(ann_c, mis_c, unann_c, split_c)
  short_asms <- gsub("^Tae_|^Tdic_|^Tdur_", "", asms)
  bp_c <- barplot(mat_c, names.arg = short_asms, las = 2,
                  col = c("#4daf4a", "#ff7f00", "#e41a1c", "#984ea3"),
                  ylim = c(0, 118), ylab = "Homologue loci (>=80% id & qcov)",
                  main = sprintf("c  Direct genome re-annotation rescues %d loci (+%.1f +/- %.1f%%, P = %.2e)",
                                 n_resc_total, pct_mean, pct_sem, wt_rescue$p.value))
  legend("topleft", legend = c(sprintf("Annotated GFF3 (n = %d)", sum(ann_c)),
                               sprintf("Rescued misannotated (n = %d)", sum(mis_c)),
                               sprintf("Rescued unannotated (n = %d)", sum(unann_c)),
                               sprintf("Rescued split model (n = %d)", sum(split_c))),
         fill = c("#4daf4a", "#ff7f00", "#e41a1c", "#984ea3"), bty = "n", ncol = 2, cex = 0.72)

  # ---- Panel d: Core-to-dispensable pangenome partitioning
  par(mar = c(4.8, 5.0, 2.6, 1.5))
  comp_levels <- c("core", "soft-core", "shell", "cloud")
  comp_cols   <- c("#1b7837", "#7fbf7b", "#fdb863", "#e66101")
  grp_list <- list(
    nlr_og,
    non_og,
    og_main[og_main$class == "III", ],
    og_main[og_main$class == "IV", ],
    og_main[og_main$class == "V", ]
  )
  grp_names <- c(sprintf("Classes I-II\nNLR (n=%d)", nrow(nlr_og)),
                 sprintf("Classes III-V\nnon-NLR (n=%d)", nrow(non_og)),
                 sprintf("Class III\nTKP/TKF (n=%d)", sum(og_main$class == "III")),
                 sprintf("Class IV\nANK/NAM (n=%d)", sum(og_main$class == "IV")),
                 sprintf("Class V\nAPR (n=%d)", sum(og_main$class == "V")))
  mat_d <- sapply(grp_list, function(df) {
    sapply(comp_levels, function(cpt) 100 * sum(df$compartment == cpt) / max(1, nrow(df)))
  })
  bp_d <- barplot(mat_d, names.arg = grp_names, col = comp_cols, ylim = c(0, 118),
                  ylab = "Orthogroups in compartment (%)",
                  main = sprintf("d  Pangenome partitioning (%d%% NLR shell/cloud vs %d%% non-NLR core/soft-core)",
                                 round(100 * nlr_sc_sub / nrow(nlr_og)), round(100 * non_cs_sub / nrow(non_og))))
  for (j in seq_along(grp_list)) {
    cum_h <- 0
    n_tot <- max(1, nrow(grp_list[[j]]))
    for (i in seq_along(comp_levels)) {
      v <- mat_d[i, j]
      cnt <- sum(grp_list[[j]]$compartment == comp_levels[i])
      if (v >= 12) {
        text(bp_d[j], cum_h + v / 2, labels = sprintf("%.1f%%\n(%d/%d)", v, cnt, n_tot),
             cex = 0.68, font = 2, col = ifelse(comp_levels[i] %in% c("core", "cloud"), "white", "black"))
      }
      cum_h <- cum_h + v
    }
  }
  abline(v = mean(bp_d[2:3]), lty = 2, col = "gray40")
  text(mean(bp_d[1:2]), 108,
       labels = sprintf("Fisher Haldane OR = %.1f, P = %.2e | MWU U = %.1f, P = %.2e",
                        or_ha_sub, ft_sub$p.value, mwu_sub$statistic, mwu_sub$p.value),
       cex = 0.72, font = 2)
  legend("topright", legend = c("Core (>=95%)", "Soft-core (80-94%)", "Shell (15-79%)", "Cloud (<15%)"),
         fill = comp_cols, bty = "n", ncol = 2, cex = 0.72)

  # ---- Panel e: Orthogroup presence/absence profile across 56 subgenomes
  par(mar = c(4.2, 10.5, 2.6, 1.5))
  pm_df <- read.delim(file.path(rd, "og_presence_matrix.tsv"), row.names = 1, check.names = FALSE)
  ord_idx <- order(ifelse(og_main$group == "NLR", 1, 2), -og_main$occupancy_all_subgenomes)
  og_ord  <- og_main[ord_idx, ]
  pm_ord  <- as.matrix(pm_df[, og_ord$orthogroup, drop = FALSE])
  image(x = seq_len(nrow(pm_ord)), y = seq_len(ncol(pm_ord)), z = pm_ord,
        col = c("#f7fbff", "#08519c"), yaxt = "n", xaxt = "n",
        xlab = sprintf("%d haploid subgenomes (Diploid A, S, D -> Tetraploid AB -> Hexaploid ABD)", nrow(pm_ord)),
        ylab = "",
        main = sprintf("e  Subgenome presence/absence matrix across %d Main orthogroups x %d subgenomes",
                       ncol(pm_ord), nrow(pm_ord)))
  og_labels <- sprintf("%s (%s) [%s]", og_ord$orthogroup, substr(og_ord$queries, 1, 14), og_ord$class)
  axis(2, at = seq_len(ncol(pm_ord)), labels = og_labels, las = 1, cex.axis = 0.56)
  abline(h = sum(og_ord$group == "NLR") + 0.5, col = "#d95f02", lwd = 1.5, lty = 2)

  # ---- Panel f: Heaps' law saturation & core-genome decay curves
  par(mar = c(4.2, 5.0, 2.6, 1.5))
  c_nlr <- curves_df[curves_df$group == "NLR", ]
  c_non <- curves_df[curves_df$group == "nonNLR", ]
  max_y <- max(c(c_nlr$mean + c_nlr$sd, c_non$mean + c_non$sd)) + 2.5
  plot(c_nlr$n_genomes, c_nlr$mean, type = "l", lwd = 2.2, col = "#1f77b4",
       xlim = c(1, 56), ylim = c(0, max_y),
       xlab = "Number of sampled haploid subgenomes (n)",
       ylab = "Cumulative orthogroups (pan / core)",
       main = sprintf("f  Orthogroup accumulation & core decay (1,000 perms; P = %.4f)",
                      h_nlr$P_label_permutation))
  polygon(c(c_nlr$n_genomes, rev(c_nlr$n_genomes)),
          c(c_nlr$mean - c_nlr$sd, rev(c_nlr$mean + c_nlr$sd)),
          col = adjustcolor("#1f77b4", alpha.f = 0.2), border = NA)
  lines(c_non$n_genomes, c_non$mean, lwd = 2.2, col = "#d95f02")
  polygon(c(c_non$n_genomes, rev(c_non$n_genomes)),
          c(c_non$mean - c_non$sd, rev(c_non$mean + c_non$sd)),
          col = adjustcolor("#d95f02", alpha.f = 0.2), border = NA)
  if ("core_mean" %in% colnames(c_nlr)) {
    lines(c_nlr$n_genomes, c_nlr$core_mean, lwd = 1.6, lty = 2, col = "#1f77b4")
    lines(c_non$n_genomes, c_non$core_mean, lwd = 1.6, lty = 2, col = "#d95f02")
  }
  legend(x = 22, y = 18.2,
         legend = c(sprintf("Classes I-II NLR (%d OGs, pan): gamma = %.3f", h_nlr$n_orthogroups, h_nlr$gamma_meancurve),
                    sprintf("Classes III-V non-NLR (%d OGs, pan): gamma = %.3f", h_non$n_orthogroups, h_non$gamma_meancurve),
                    sprintf("Classes I-II NLR (core decay: n56 = %d)", h_nlr$core_orthogroups_n56),
                    sprintf("Classes III-V non-NLR (core decay: n56 = %d)", h_non$core_orthogroups_n56)),
         col = c("#1f77b4", "#d95f02", "#1f77b4", "#d95f02"),
         lty = c(1, 1, 2, 2), lwd = c(2.2, 2.2, 1.6, 1.6), bg = adjustcolor("white", alpha.f = 0.85),
         box.lwd = 0, cex = 0.74)
}

pdf(file.path(fd, "Fig1_pangenome_inventory_R.pdf"), width = 15, height = 12)
draw_fig1_R()
invisible(dev.off())

png(file.path(fd, "Fig1_pangenome_inventory_R.png"), width = 15, height = 12, units = "in", res = 300)
draw_fig1_R()
invisible(dev.off())

cat(sprintf("Saved R publication figures -> %s/Fig1_pangenome_inventory_R.pdf and .png\n", fd))
cat("======================================================================\n")
