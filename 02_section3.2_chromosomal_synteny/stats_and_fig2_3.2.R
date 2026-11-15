#!/usr/bin/env Rscript
# ==============================================================================
# stats_and_fig2_3.2.R
# Standalone R statistical audit & publication-ready figure generator for
# Results Section 3.2 ("Subgenome asymmetry, chromosomal distribution, and
# microsynteny retention") and Methods 5.4.
#
# Inputs (produced by rust32.py):
#   results/S4_chromosomal_positions.tsv
#   results/S4b_locusgroup_positions.tsv
#   results/S5_hexaploid_triads.tsv
#   results/S5b_orthogroup_triad_summary.tsv
#   results/S5c_ensembl_compara_triads.tsv
#   results/S5d_subgenome_cnv.tsv
#   results/S6a_sri_per_locus.tsv
#   results/S6b_sri_pairwise.tsv
#   results/S6c_microsynteny_ribbons.tsv
#   results/S7a_flanking_te_windows.tsv
#   results/S7b_te_summary_by_class.tsv
#
# Outputs:
#   results/S7c_R_statistical_audit_3.2.tsv
#   tables/Supplementary_Table_S7c_R_Statistical_Audit_3.2.tsv
#   figures/Fig2_chromosomal_distribution_R.pdf & .png
#   figures/Supplementary_Fig1_subgenome_CNV_R.pdf & .png
#   figures/Supplementary_Fig2_TE_vs_SRI_R.pdf & .png
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

options(stringsAsFactors = FALSE, warn = -1)

read_tsv_clean <- function(path) {
  if (!file.exists(path)) stop(paste("Missing required input file:", path))
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "",
                    check.names = FALSE, comment.char = "")
}

pos_all  <- read_tsv_clean("results/S4_chromosomal_positions.tsv")
pos_lg   <- read_tsv_clean("results/S4b_locusgroup_positions.tsv")
triad_lg <- read_tsv_clean("results/S5_hexaploid_triads.tsv")
triad_og <- read_tsv_clean("results/S5b_orthogroup_triad_summary.tsv")
sg_cnv   <- read_tsv_clean("results/S5d_subgenome_cnv.tsv")
sri_lg   <- read_tsv_clean("results/S6a_sri_per_locus.tsv")
sri_pw   <- read_tsv_clean("results/S6b_sri_pairwise.tsv")
ribbons  <- read_tsv_clean("results/S6c_microsynteny_ribbons.tsv")
te_lg    <- read_tsv_clean("results/S7a_flanking_te_windows.tsv")

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

sem <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(0.0)
  stats::sd(x) / sqrt(length(x))
}

iqr_bounds <- function(x) {
  stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, type = 7)
}

haldane_or <- function(a, b, c, d) {
  ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
}

audit_rows <- list()
add_audit <- function(section, metric, group_comparison, estimate, test_statistic, p_value, note = "") {
  audit_rows[[length(audit_rows) + 1]] <<- data.frame(
    section = section,
    metric = metric,
    group_comparison = group_comparison,
    estimate = estimate,
    test_statistic = test_statistic,
    p_value = p_value,
    note = note,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 1. Chromosomal Position (d_norm from centromere 0.0 to telomere 1.0)
# ==============================================================================
d_nlr_all <- as.numeric(pos_all$d_norm[pos_all$Group == "NLR"])
d_non_all <- as.numeric(pos_all$d_norm[pos_all$Group == "nonNLR"])
q_nlr_all <- iqr_bounds(d_nlr_all)
q_non_all <- iqr_bounds(d_non_all)
pct_dist_nlr_all <- 100.0 * mean(d_nlr_all >= 0.85)
pct_dist_non_all <- 100.0 * mean(d_non_all >= 0.85)

ks_all  <- suppressWarnings(stats::ks.test(d_nlr_all, d_non_all, alternative = "two.sided"))
mwu_all <- suppressWarnings(stats::wilcox.test(d_nlr_all, d_non_all, alternative = "two.sided", exact = FALSE, correct = FALSE))

add_audit(
  "3.2.1_Chromosomal_Position", "Arm-normalized position d_norm (all placed loci)",
  sprintf("Classes I-II NLR (n=%d) vs Classes III-V non-NLR (n=%d)", length(d_nlr_all), length(d_non_all)),
  sprintf("NLR median = %.3f (IQR %.3f-%.3f, %.1f%% distal >=0.85) vs non-NLR median = %.3f (IQR %.3f-%.3f, %.1f%% distal >=0.85)",
          stats::median(d_nlr_all), q_nlr_all[1], q_nlr_all[2], pct_dist_nlr_all,
          stats::median(d_non_all), q_non_all[1], q_non_all[2], pct_dist_non_all),
  sprintf("KS D = %.4f; MWU W = %.1f", ks_all$statistic, mwu_all$statistic),
  sprintf("KS P = %.3e; MWU P = %.3e", ks_all$p.value, mwu_all$p.value),
  "Two-sample Kolmogorov-Smirnov and Mann-Whitney U tests across all 12,018 placed candidate loci"
)

d_nlr_lg <- as.numeric(pos_lg$median_d_norm[pos_lg$Group == "NLR" & pos_lg$median_d_norm != ""])
d_non_lg <- as.numeric(pos_lg$median_d_norm[pos_lg$Group == "nonNLR" & pos_lg$median_d_norm != ""])
q_nlr_lg <- iqr_bounds(d_nlr_lg)
q_non_lg <- iqr_bounds(d_non_lg)
n_dist_nlr_lg <- sum(d_nlr_lg >= 0.85)
n_dist_non_lg <- sum(d_non_lg >= 0.85)

ks_lg  <- suppressWarnings(stats::ks.test(d_nlr_lg, d_non_lg, alternative = "two.sided"))
mwu_lg <- suppressWarnings(stats::wilcox.test(d_nlr_lg, d_non_lg, alternative = "two.sided", exact = FALSE, correct = FALSE))

add_audit(
  "3.2.1_Chromosomal_Position", "Per-LocusGroup median d_norm (36 independent loci)",
  sprintf("Classes I-II NLR (n=%d loci) vs Classes III-V non-NLR (n=%d loci)", length(d_nlr_lg), length(d_non_lg)),
  sprintf("NLR median = %.3f (IQR %.3f-%.3f; %d/%d [%.1f%%] distal) vs non-NLR median = %.3f (IQR %.3f-%.3f; %d/%d [%.1f%%] distal)",
          stats::median(d_nlr_lg), q_nlr_lg[1], q_nlr_lg[2], n_dist_nlr_lg, length(d_nlr_lg), 100.0 * n_dist_nlr_lg / length(d_nlr_lg),
          stats::median(d_non_lg), q_non_lg[1], q_non_lg[2], n_dist_non_lg, length(d_non_lg), 100.0 * n_dist_non_lg / length(d_non_lg)),
  sprintf("KS D = %.4f; MWU W = %.1f", ks_lg$statistic, mwu_lg$statistic),
  sprintf("KS P = %.3e; MWU P = %.3e", ks_lg$p.value, mwu_lg$p.value),
  "Independent LocusGroup-level comparison"
)

# ==============================================================================
# 2. Homoeologous A:B:D Triad Balance (13 Hexaploids) & Subgenome CNV Asymmetry
# ==============================================================================
nlr_og <- triad_og[triad_og$Group == "NLR", ]
non_og <- triad_og[triad_og$Group == "nonNLR", ]

n_nlr_con <- sum(as.integer(nlr_og$conserved_1_1_1_triad_regime) == 1)
n_non_con <- sum(as.integer(non_og$conserved_1_1_1_triad_regime) == 1)
mat_con <- matrix(c(n_non_con, nrow(non_og) - n_non_con,
                    n_nlr_con, nrow(nlr_og) - n_nlr_con), nrow = 2, byrow = TRUE)
fe_con <- stats::fisher.test(mat_con, alternative = "two.sided")
ha_con <- haldane_or(n_non_con, nrow(non_og) - n_non_con, n_nlr_con, nrow(nlr_og) - n_nlr_con)

add_audit(
  "3.2.2_Hexaploid_Triads", "Conserved A:B:D triad retention (>=60% [>=8/13] hexaploid assemblies)",
  sprintf("Classes III-V non-NLR OGs (n=%d) vs Classes I-II NLR OGs (n=%d)", nrow(non_og), nrow(nlr_og)),
  sprintf("non-NLR = %d/%d (%.1f%%) vs NLR = %d/%d (%.1f%%)",
          n_non_con, nrow(non_og), 100.0 * n_non_con / nrow(non_og),
          n_nlr_con, nrow(nlr_og), 100.0 * n_nlr_con / nrow(nlr_og)),
  sprintf("Exact OR = %s (Haldane-Anscombe OR = %.1f)",
          ifelse(is.infinite(fe_con$estimate), "Inf", sprintf("%.2f", fe_con$estimate)), ha_con),
  sprintf("%.3e", fe_con$p.value),
  "Two-sided Fisher's exact test on S9-pruned focal orthogroup triads across 13 hexaploid cultivars"
)

n_nlr_ever3 <- sum(as.integer(nlr_og$n_all_three_ABD) >= 1)
n_non_ever3 <- sum(as.integer(non_og$n_all_three_ABD) >= 1)
mat_ever3 <- matrix(c(n_non_ever3, nrow(non_og) - n_non_ever3,
                      n_nlr_ever3, nrow(nlr_og) - n_nlr_ever3), nrow = 2, byrow = TRUE)
fe_ever3 <- stats::fisher.test(mat_ever3, alternative = "two.sided")
ha_ever3 <- haldane_or(n_non_ever3, nrow(non_og) - n_non_ever3, n_nlr_ever3, nrow(nlr_og) - n_nlr_ever3)

add_audit(
  "3.2.2_Hexaploid_Triads", "Ever-present A:B:D triad retention (>=1/13 hexaploid assembly)",
  sprintf("Classes III-V non-NLR OGs (n=%d) vs Classes I-II NLR OGs (n=%d)", nrow(non_og), nrow(nlr_og)),
  sprintf("non-NLR = %d/%d (%.1f%%) vs NLR = %d/%d (%.1f%%)",
          n_non_ever3, nrow(non_og), 100.0 * n_non_ever3 / nrow(non_og),
          n_nlr_ever3, nrow(nlr_og), 100.0 * n_nlr_ever3 / nrow(nlr_og)),
  sprintf("Exact OR = %.2f (Haldane-Anscombe OR = %.1f)", fe_ever3$estimate, ha_ever3),
  sprintf("%.3e", fe_ever3$p.value),
  "Two-sided Fisher's exact test for A+B+D co-occurrence in at least one hexaploid assembly"
)

nlr_cells <- triad_lg[triad_lg$Group == "NLR" & triad_lg$triad_state != "absent", ]
non_cells <- triad_lg[triad_lg$Group == "nonNLR" & triad_lg$triad_state != "absent", ]
c_nlr_111 <- sum(nlr_cells$triad_state == "balanced_1:1:1")
c_non_111 <- sum(non_cells$triad_state == "balanced_1:1:1")
mat_cell  <- matrix(c(c_non_111, nrow(non_cells) - c_non_111,
                      c_nlr_111, nrow(nlr_cells) - c_nlr_111), nrow = 2, byrow = TRUE)
fe_cell <- stats::fisher.test(mat_cell, alternative = "two.sided")

add_audit(
  "3.2.2_Hexaploid_Triads", "Strict 1:1:1 A:B:D triads across present locus x hexaploid cells",
  sprintf("non-NLR cells (n=%d) vs NLR cells (n=%d)", nrow(non_cells), nrow(nlr_cells)),
  sprintf("non-NLR = %d/%d (%.1f%%) vs NLR = %d/%d (%.1f%%)",
          c_non_111, nrow(non_cells), 100.0 * c_non_111 / nrow(non_cells),
          c_nlr_111, nrow(nlr_cells), 100.0 * c_nlr_111 / nrow(nlr_cells)),
  sprintf("Exact OR = %.2f", fe_cell$estimate),
  sprintf("%.3e", fe_cell$p.value),
  "Two-sided Fisher's exact test across present canonical chromosome group triads"
)

cn_A  <- as.numeric(sg_cnv$intact_copies[sg_cnv$Group == "NLR" & sg_cnv$tier == "A"])
cn_BS <- as.numeric(sg_cnv$intact_copies[sg_cnv$Group == "NLR" & sg_cnv$tier == "B/S"])
cn_D  <- as.numeric(sg_cnv$intact_copies[sg_cnv$Group == "NLR" & sg_cnv$tier == "D"])
kw_res <- stats::kruskal.test(list(A = cn_A, BS = cn_BS, D = cn_D))

add_audit(
  "3.2.2_Subgenome_CNV", "NLR intact copy-number expansion across A, B/S, and D subgenomes",
  sprintf("A (n=%d) vs B/S (n=%d) vs D (n=%d) orthogroup x subgenome cells", length(cn_A), length(cn_BS), length(cn_D)),
  sprintf("A = %.2f +/- %.2f (SD), B/S = %.2f +/- %.2f (SD), D = %.2f +/- %.2f (SD)",
          mean(cn_A), stats::sd(cn_A), mean(cn_BS), stats::sd(cn_BS), mean(cn_D), stats::sd(cn_D)),
  sprintf("Kruskal-Wallis H = %.2f (df = 2)", kw_res$statistic),
  sprintf("%.3e", kw_res$p.value),
  "Kruskal-Wallis test across 56 haploid subgenomes (19 A, 21 B/S, 16 D) x 20 NLR orthogroups"
)

# ==============================================================================
# 3. Microsynteny Retention Index (SRI) across +/-500 kb Windows
# ==============================================================================
sri_nlr_loci <- as.numeric(sri_lg$mean_SRI[sri_lg$Group == "NLR" & sri_lg$mean_SRI != ""])
sri_non_loci <- as.numeric(sri_lg$mean_SRI[sri_lg$Group == "nonNLR" & sri_lg$mean_SRI != ""])
mwu_sri_loci <- suppressWarnings(stats::wilcox.test(sri_nlr_loci, sri_non_loci, alternative = "two.sided", exact = FALSE, correct = FALSE))

add_audit(
  "3.2.3_Microsynteny_SRI", "Synteny Retention Index (SRI) per LocusGroup (+/-500 kb)",
  sprintf("Classes I-II NLR (n=%d loci) vs Classes III-V non-NLR (n=%d loci)", length(sri_nlr_loci), length(sri_non_loci)),
  sprintf("NLR = %.3f +/- %.3f (median %.3f) vs non-NLR = %.3f +/- %.3f (median %.3f)",
          mean(sri_nlr_loci), sem(sri_nlr_loci), stats::median(sri_nlr_loci),
          mean(sri_non_loci), sem(sri_non_loci), stats::median(sri_non_loci)),
  sprintf("MWU U = %.1f", mwu_sri_loci$statistic),
  sprintf("%.3e", mwu_sri_loci$p.value),
  sprintf("Based on %d pairwise window comparisons across 16 annotated assemblies", nrow(sri_pw))
)

og_sri_df <- stats::aggregate(as.numeric(sri_lg$mean_SRI) ~ orthogroup + Group, data = sri_lg, FUN = mean)
colnames(og_sri_df)[3] <- "og_mean_SRI"
sri_nlr_og <- og_sri_df$og_mean_SRI[og_sri_df$Group == "NLR"]
sri_non_og <- og_sri_df$og_mean_SRI[og_sri_df$Group == "nonNLR"]
mwu_sri_og <- suppressWarnings(stats::wilcox.test(sri_nlr_og, sri_non_og, alternative = "two.sided", exact = FALSE, correct = FALSE))

add_audit(
  "3.2.3_Microsynteny_SRI", "Synteny Retention Index (SRI) per Main Orthogroup (+/-500 kb)",
  sprintf("Classes I-II NLR (n=%d OGs) vs Classes III-V non-NLR (n=%d OGs)", length(sri_nlr_og), length(sri_non_og)),
  sprintf("NLR = %.3f +/- %.3f (median %.3f) vs non-NLR = %.3f +/- %.3f (median %.3f)",
          mean(sri_nlr_og), sem(sri_nlr_og), stats::median(sri_nlr_og),
          mean(sri_non_og), sem(sri_non_og), stats::median(sri_non_og)),
  sprintf("MWU U = %.1f", mwu_sri_og$statistic),
  sprintf("%.3e", mwu_sri_og$p.value),
  "Aggregated across the 30 Main orthogroups (20 NLR vs 10 non-NLR)"
)

# ==============================================================================
# 4. Flanking TE Coverage & Local Rearrangements (SVs/Indels)
# ==============================================================================
te_nlr <- as.numeric(te_lg$total_te_pct[te_lg$Group == "NLR"])
te_non <- as.numeric(te_lg$total_te_pct[te_lg$Group == "nonNLR"])
mwu_te <- suppressWarnings(stats::wilcox.test(te_nlr, te_non, alternative = "two.sided", exact = FALSE, correct = FALSE))

sv_nlr <- as.numeric(te_lg$mean_sv_events[te_lg$Group == "NLR"])
sv_non <- as.numeric(te_lg$mean_sv_events[te_lg$Group == "nonNLR"])
mwu_sv <- suppressWarnings(stats::wilcox.test(sv_nlr, sv_non, alternative = "two.sided", exact = FALSE, correct = FALSE))

all_te  <- as.numeric(te_lg$total_te_pct)
all_sri <- as.numeric(te_lg$mean_SRI)
all_sv  <- as.numeric(te_lg$mean_sv_events)

sp_sv_sri <- suppressWarnings(stats::cor.test(all_sv, all_sri, method = "spearman", exact = FALSE))
sp_te_sv  <- suppressWarnings(stats::cor.test(all_te, all_sv,  method = "spearman", exact = FALSE))
sp_te_sri <- suppressWarnings(stats::cor.test(all_te, all_sri, method = "spearman", exact = FALSE))

add_audit(
  "3.2.4_Flanking_TE_and_SV", "Flanking +/-200 kb CS Ensembl dna_sm TE coverage (%)",
  sprintf("Classes I-II NLR (n=%d) vs Classes III-V non-NLR (n=%d)", length(te_nlr), length(te_non)),
  sprintf("NLR = %.2f +/- %.2f%% vs non-NLR = %.2f +/- %.2f%%",
          mean(te_nlr), sem(te_nlr), mean(te_non), sem(te_non)),
  sprintf("MWU U = %.1f", mwu_te$statistic),
  sprintf("%.3e", mwu_te$p.value),
  "Reference Chinese Spring Ensembl dna_sm soft-masked repeat coverage in +/-200 kb windows"
)

add_audit(
  "3.2.4_Flanking_TE_and_SV", "Local gene-order rearrangements/indels per +/-500 kb window",
  sprintf("Classes I-II NLR (n=%d) vs Classes III-V non-NLR (n=%d)", length(sv_nlr), length(sv_non)),
  sprintf("NLR = %.2f +/- %.2f vs non-NLR = %.2f +/- %.2f events",
          mean(sv_nlr), sem(sv_nlr), mean(sv_non), sem(sv_non)),
  sprintf("MWU U = %.1f", mwu_sv$statistic),
  sprintf("%.3e", mwu_sv$p.value),
  "Mean local inversions + gene-order indels per +/-500 kb window"
)

add_audit(
  "3.2.4_Flanking_TE_and_SV", "Spearman rank correlations (SV vs SRI; TE vs SV; TE vs SRI)",
  sprintf("All %d cloned R-gene loci", nrow(te_lg)),
  sprintf("rho(SV, SRI) = %.3f; rho(TE, SV) = %.3f; rho(TE, SRI) = %.3f",
          sp_sv_sri$estimate, sp_te_sv$estimate, sp_te_sri$estimate),
  sprintf("S = %.1f / %.1f / %.1f", sp_sv_sri$statistic, sp_te_sv$statistic, sp_te_sri$statistic),
  sprintf("P(SV,SRI) = %.3e; P(TE,SV) = %.3e; P(TE,SRI) = %.3e",
          sp_sv_sri$p.value, sp_te_sv$p.value, sp_te_sri$p.value),
  "Two-sided Spearman rank correlations across n = 36 loci"
)

audit_df <- do.call(rbind, audit_rows)
utils::write.table(audit_df, file = "results/S7c_R_statistical_audit_3.2.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(audit_df, file = "tables/Supplementary_Table_S7c_R_Statistical_Audit_3.2.tsv",
                   sep = "\t", quote = FALSE, row.names = FALSE)

# ==============================================================================
# Print Console Report
# ==============================================================================
cat("======================================================================\n")
cat("  Section 3.2 R Statistical Audit & Chromosomal/Microsynteny Analysis\n")
cat("======================================================================\n\n")

cat(sprintf("[1] Chromosomal Position (d_norm from centromere 0.0 to telomere 1.0):\n"))
cat(sprintf("    All placed loci (n = %d NLR vs %d non-NLR across 36 loci):\n", length(d_nlr_all), length(d_non_all)))
cat(sprintf("      Classes I-II NLR:     median = %.3f (IQR = %.3f-%.3f), %.1f%% distal (d_norm >= 0.85)\n",
            stats::median(d_nlr_all), q_nlr_all[1], q_nlr_all[2], pct_dist_nlr_all))
cat(sprintf("      Classes III-V non-NLR: median = %.3f (IQR = %.3f-%.3f), %.1f%% distal (d_norm >= 0.85)\n",
            stats::median(d_non_all), q_non_all[1], q_non_all[2], pct_dist_non_all))
cat(sprintf("      Two-sample KS test: D = %.4f, P = %.3e | Mann-Whitney U = %.1f, P = %.3e\n",
            ks_all$statistic, ks_all$p.value, mwu_all$statistic, mwu_all$p.value))
cat(sprintf("    Independent LocusGroups (n = %d NLR vs %d non-NLR loci):\n", length(d_nlr_lg), length(d_non_lg)))
cat(sprintf("      Classes I-II NLR:     median = %.3f (IQR = %.3f-%.3f), %d/%d (%.1f%%) distal\n",
            stats::median(d_nlr_lg), q_nlr_lg[1], q_nlr_lg[2], n_dist_nlr_lg, length(d_nlr_lg), 100.0 * n_dist_nlr_lg / length(d_nlr_lg)))
cat(sprintf("      Classes III-V non-NLR: median = %.3f (IQR = %.3f-%.3f), %d/%d (%.1f%%) distal\n",
            stats::median(d_non_lg), q_non_lg[1], q_non_lg[2], n_dist_non_lg, length(d_non_lg), 100.0 * n_dist_non_lg / length(d_non_lg)))
cat(sprintf("      Two-sample KS test: D = %.4f, P = %.3e | Mann-Whitney U = %.1f, P = %.3e\n\n",
            ks_lg$statistic, ks_lg$p.value, mwu_lg$statistic, mwu_lg$p.value))

cat(sprintf("[2] Homoeologous A:B:D Triad Balance (n = 13 hexaploids) & Subgenome CNV Asymmetry:\n"))
cat(sprintf("    Conserved A:B:D triads (>=60%% [>=8/13] hexaploids): non-NLR = %d/%d (%.1f%%) vs NLR = %d/%d (%.1f%%)\n",
            n_non_con, nrow(non_og), 100.0 * n_non_con / nrow(non_og),
            n_nlr_con, nrow(nlr_og), 100.0 * n_nlr_con / nrow(nlr_og)))
cat(sprintf("      Fisher's exact test: OR = %s (Haldane-Anscombe OR = %.1f), P = %.3e\n",
            ifelse(is.infinite(fe_con$estimate), "Inf", sprintf("%.2f", fe_con$estimate)), ha_con, fe_con$p.value))
cat(sprintf("    Ever-present A:B:D triads (>=1/13 hexaploid): non-NLR = %d/%d (%.1f%%) vs NLR = %d/%d (%.1f%%)\n",
            n_non_ever3, nrow(non_og), 100.0 * n_non_ever3 / nrow(non_og),
            n_nlr_ever3, nrow(nlr_og), 100.0 * n_nlr_ever3 / nrow(nlr_og)))
cat(sprintf("      Fisher's exact test: OR = %.2f (Haldane-Anscombe OR = %.1f), P = %.3e\n",
            fe_ever3$estimate, ha_ever3, fe_ever3$p.value))
cat(sprintf("    Strict 1:1:1 A:B:D triads across present locus x cultivar cells: non-NLR = %d/%d (%.1f%%) vs NLR = %d/%d (%.1f%%)\n",
            c_non_111, nrow(non_cells), 100.0 * c_non_111 / nrow(non_cells),
            c_nlr_111, nrow(nlr_cells), 100.0 * c_nlr_111 / nrow(nlr_cells)))
cat(sprintf("      Fisher's exact test: OR = %.2f, P = %.3e\n", fe_cell$estimate, fe_cell$p.value))
cat(sprintf("    NLR intact copies per subgenome tier: A = %.2f +/- %.2f, B/S = %.2f +/- %.2f, D = %.2f +/- %.2f\n",
            mean(cn_A), stats::sd(cn_A), mean(cn_BS), stats::sd(cn_BS), mean(cn_D), stats::sd(cn_D)))
cat(sprintf("      Kruskal-Wallis H = %.2f (df = 2), P = %.3e\n\n", kw_res$statistic, kw_res$p.value))

cat(sprintf("[3] Microsynteny Retention Index (SRI) across +/-500 kb Windows (n = %d pairwise comparisons):\n", nrow(sri_pw)))
cat(sprintf("    Per LocusGroup (n = 26 NLR vs 10 non-NLR): NLR = %.3f +/- %.3f vs non-NLR = %.3f +/- %.3f | MWU U = %.1f, P = %.3e\n",
            mean(sri_nlr_loci), sem(sri_nlr_loci), mean(sri_non_loci), sem(sri_non_loci), mwu_sri_loci$statistic, mwu_sri_loci$p.value))
cat(sprintf("    Per Orthogroup (n = 20 NLR vs 10 non-NLR): NLR = %.3f +/- %.3f vs non-NLR = %.3f +/- %.3f | MWU U = %.1f, P = %.3e\n",
            mean(sri_nlr_og), sem(sri_nlr_og), mean(sri_non_og), sem(sri_non_og), mwu_sri_og$statistic, mwu_sri_og$p.value))
for (k in c("I", "II", "III", "IV", "V")) {
  vk <- as.numeric(sri_lg$mean_SRI[sri_lg$Class == k & sri_lg$mean_SRI != ""])
  cat(sprintf("      Class %-3s (n = %2d loci): mean SRI = %.3f +/- %.3f\n", k, length(vk), mean(vk), sem(vk)))
}
cat("\n")

cat(sprintf("[4] Flanking +/-200 kb RepeatMasker TE Coverage & +/-500 kb Local Structural Rearrangements:\n"))
cat(sprintf("    Flanking TE coverage: NLR = %.2f +/- %.2f%% vs non-NLR = %.2f +/- %.2f%% (MWU U = %.1f, P = %.3e)\n",
            mean(te_nlr), sem(te_nlr), mean(te_non), sem(te_non), mwu_te$statistic, mwu_te$p.value))
cat(sprintf("    Local rearrangements/indels: NLR = %.2f +/- %.2f vs non-NLR = %.2f +/- %.2f (MWU U = %.1f, P = %.3e)\n",
            mean(sv_nlr), sem(sv_nlr), mean(sv_non), sem(sv_non), mwu_sv$statistic, mwu_sv$p.value))
cat(sprintf("    Spearman rho(SV, SRI) = %+.3f (P = %.3e) | rho(TE, SV) = %+.3f (P = %.3e) | rho(TE, SRI) = %+.3f (P = %.3e)\n\n",
            sp_sv_sri$estimate, sp_sv_sri$p.value, sp_te_sv$estimate, sp_te_sv$p.value, sp_te_sri$estimate, sp_te_sri$p.value))

# ==============================================================================
# 5. Publication-Grade R Figures (Main Fig. 2, Suppl. Fig. 1, Suppl. Fig. 2)
# ==============================================================================
class_cols <- c("I" = "#1f77b4", "II" = "#4a90d9", "III" = "#d95f02", "IV" = "#7570b3", "V" = "#1b9e77")

draw_fig2 <- function() {
  # Layout: 3 rows x 3 columns where panels a, b, d, e, f span appropriate cells and panel c has ternary + barplot
  layout_mat <- matrix(c(
    1, 1, 2, 2,
    3, 4, 5, 5,
    6, 6, 7, 7
  ), nrow = 3, byrow = TRUE)
  graphics::layout(layout_mat, widths = c(1.15, 0.90, 1.02, 1.02), heights = c(1.0, 1.05, 1.18))
  graphics::par(oma = c(0.6, 0.6, 0.6, 0.6), cex = 0.82)

  # ---- Panel a: Chromosomal distribution histogram ----
  graphics::par(mar = c(4.2, 4.2, 2.6, 1.2))
  brks <- seq(0, 1, length.out = 21)
  h_nlr <- graphics::hist(d_nlr_all, breaks = brks, plot = FALSE)
  h_non <- graphics::hist(d_non_all, breaks = brks, plot = FALSE)
  ymax_a <- max(c(h_nlr$density, h_non$density)) * 1.18
  graphics::plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, ymax_a),
                 xlab = "Arm-normalized position from centromere (0.0) to telomere (1.0)",
                 ylab = "Normalized locus density", bty = "l",
                 main = sprintf("a  Chromosomal distribution (KS D = %.2f, P = %.2e)", ks_all$statistic, ks_all$p.value),
                 adj = 0, cex.main = 1.05)
  for (i in seq_along(h_nlr$density)) {
    graphics::rect(brks[i], 0, brks[i + 1], h_nlr$density[i],
                   col = grDevices::adjustcolor("#1f77b4", alpha.f = 0.55), border = "black", lwd = 0.5)
    graphics::rect(brks[i], 0, brks[i + 1], h_non$density[i],
                   col = grDevices::adjustcolor("#d95f02", alpha.f = 0.55), border = "black", lwd = 0.5)
  }
  graphics::abline(v = 0.85, col = "#b2182b", lty = 2, lwd = 1.5)
  graphics::legend("topleft", bty = "n", cex = 0.86,
                   legend = c(sprintf("Classes I-II NLR (n = %d, med = %.2f)", length(d_nlr_all), stats::median(d_nlr_all)),
                              sprintf("Classes III-V non-NLR (n = %d, med = %.2f)", length(d_non_all), stats::median(d_non_all)),
                              "Distal 15% subtelomeric boundary (0.85)"),
                   fill = c(grDevices::adjustcolor("#1f77b4", 0.55), grDevices::adjustcolor("#d95f02", 0.55), NA),
                   border = c("black", "black", NA),
                   lty = c(NA, NA, 2), lwd = c(NA, NA, 1.5), col = c(NA, NA, "#b2182b"))

  # ---- Panel b: Boxplots across Classes I-V ----
  graphics::par(mar = c(4.2, 4.2, 2.6, 1.2))
  classes <- c("I", "II", "III", "IV", "V")
  d_by_cls <- lapply(classes, function(k) as.numeric(pos_all$d_norm[pos_all$Class == k]))
  names(d_by_cls) <- sprintf("Class %s\n(n=%d)", classes, sapply(d_by_cls, length))
  graphics::boxplot(d_by_cls, col = grDevices::adjustcolor(class_cols[classes], alpha.f = 0.78),
                    outline = FALSE, ylim = c(0, 1.05), ylab = "Arm-normalized position d_norm",
                    bty = "l", las = 1, boxwex = 0.52,
                    main = "b  Relative chromosomal position across Classes I-V", adj = 0, cex.main = 1.05)
  graphics::abline(h = 0.85, col = "#b2182b", lty = 2, lwd = 1.4)

  # ---- Panel c (Left): Equilateral A:B:D Ternary Plot of 30 Main Orthogroups ----
  graphics::par(mar = c(3.2, 1.2, 2.6, 0.4))
  h_tri <- sqrt(3.0) / 2.0
  graphics::plot(0, 0, type = "n", xlim = c(-0.14, 1.14), ylim = c(-0.26, h_tri + 0.14),
                 asp = 1, axes = FALSE, xlab = "", ylab = "",
                 main = "c  Homoeologous A:B:D triad balance", adj = 0, cex.main = 1.05)
  # Grid lines
  for (frac in c(0.25, 0.50, 0.75)) {
    graphics::segments(0.5 * frac, h_tri * frac, 1.0 - 0.5 * frac, h_tri * frac, col = "#dddddd", lty = 3, lwd = 0.8)
    graphics::segments(frac, 0, 0.5 * frac, h_tri * frac, col = "#dddddd", lty = 3, lwd = 0.8)
    graphics::segments(1.0 - frac, 0, 1.0 - 0.5 * frac, h_tri * frac, col = "#dddddd", lty = 3, lwd = 0.8)
  }
  graphics::polygon(c(0, 1, 0.5), c(0, 0, h_tri), border = "black", lwd = 1.4)
  # Balanced 1:1:1 zone around centroid (0.5, h_tri/3)
  theta <- seq(0, 2 * pi, length.out = 60)
  graphics::polygon(0.5 + 0.13 * cos(theta), (h_tri / 3.0) + 0.13 * sin(theta),
                    col = grDevices::adjustcolor("#1b9e77", alpha.f = 0.18), border = "#1b9e77", lty = 2, lwd = 0.8)
  graphics::points(0.5, h_tri / 3.0, pch = 3, col = "#1b9e77", cex = 1.3, lwd = 1.4)

  set.seed(20261108)
  for (i in seq_len(nrow(triad_og))) {
    ma <- as.numeric(triad_og$mean_copies_A[i])
    mb <- as.numeric(triad_og$mean_copies_B[i])
    md <- as.numeric(triad_og$mean_copies_D[i])
    tot <- ma + mb + md
    if (tot <= 0) next
    pb <- mb / tot
    pd <- md / tot
    tx <- 0.5 * (2.0 * pb + pd) + stats::runif(1, -0.015, 0.015)
    ty <- h_tri * pd + stats::runif(1, -0.015, 0.015)
    is_nlr <- triad_og$Group[i] == "NLR"
    bg_col <- grDevices::adjustcolor(ifelse(is_nlr, "#1f77b4", "#d95f02"), alpha.f = 0.85)
    pch_i  <- ifelse(is_nlr, 21, 23)
    cex_i  <- 1.05 + 0.07 * as.numeric(triad_og$n_all_three_ABD[i])
    graphics::points(tx, ty, pch = pch_i, bg = bg_col, col = "black", cex = cex_i, lwd = 0.7)
  }
  graphics::text(-0.02, -0.04, "A (100%)", font = 2, cex = 0.88, adj = c(0, 1))
  graphics::text(1.02, -0.04, "B (100%)", font = 2, cex = 0.88, adj = c(1, 1))
  graphics::text(0.50, h_tri + 0.05, "D (100%)", font = 2, cex = 0.88, adj = c(0.5, 0))
  graphics::legend(0.5, -0.10, xjust = 0.5, yjust = 1, bty = "n", cex = 0.80,
                   legend = c("NLR OGs (0/20 >=60% A:B:D)", "non-NLR OGs (9/10 >=60% A:B:D)"),
                   pch = c(21, 23), pt.bg = c("#1f77b4", "#d95f02"), col = "black", pt.cex = 1.1)

  # ---- Panel c (Right): Stacked Bar Chart of Hexaploid Triad Regimes ----
  graphics::par(mar = c(4.2, 4.3, 2.6, 1.0))
  states <- c("balanced_1:1:1", "present_A+B+D_multicopy", "subgenome_asymmetric", "absent")
  st_cols <- c("#1b9e77", "#66a61e", "#e7298a", "#cccccc")
  st_labs <- c("Balanced 1:1:1", "A+B+D multicopy", "Asymmetric PAV", "Absent")
  mat_c <- sapply(c("NLR", "nonNLR"), function(g) {
    sub_t <- triad_lg[triad_lg$Group == g, ]
    sapply(states, function(st) 100.0 * sum(sub_t$triad_state == st) / nrow(sub_t))
  })
  bp_c <- graphics::barplot(mat_c, col = st_cols, border = "black", ylim = c(0, 118),
                            names.arg = c("NLR\n(26x13)", "non-NLR\n(10x13)"),
                            ylab = "Hexaploid cultivar loci (%)", las = 1, width = 0.55, space = 0.45)
  for (col_i in 1:2) {
    cum_h <- 0
    for (row_i in seq_len(nrow(mat_c))) {
      v <- mat_c[row_i, col_i]
      if (v >= 6.0) {
        graphics::text(bp_c[col_i], cum_h + v / 2.0, sprintf("%.1f%%", v), font = 2, cex = 0.78)
      }
      cum_h <- cum_h + v
    }
  }
  graphics::legend("topright", bty = "n", cex = 0.72, fill = st_cols, legend = st_labs)

  # ---- Panel d: Synteny Retention Index (SRI) across Classes I-V ----
  graphics::par(mar = c(4.2, 4.2, 2.6, 1.2))
  sri_cls <- lapply(classes, function(k) as.numeric(sri_lg$mean_SRI[sri_lg$Class == k & sri_lg$mean_SRI != ""]))
  means_d <- sapply(sri_cls, mean)
  sems_d  <- sapply(sri_cls, sem)
  bp_d <- graphics::barplot(means_d, col = grDevices::adjustcolor(class_cols[classes], alpha.f = 0.80),
                            border = "black", ylim = c(0, 1.05), las = 1,
                            names.arg = sprintf("Class %s\n(n=%d)", classes, sapply(sri_cls, length)),
                            ylab = "Synteny Retention Index (SRI, +/-500 kb)",
                            main = sprintf("d  Microsynteny retention across +/-500 kb (MWU U = %.1f, P = %.2e)",
                                           mwu_sri_loci$statistic, mwu_sri_loci$p.value),
                            adj = 0, cex.main = 1.05)
  graphics::arrows(bp_d, means_d - sems_d, bp_d, means_d + sems_d, angle = 90, code = 3, length = 0.05, lwd = 1.2)
  set.seed(20261101)
  for (i in seq_along(sri_cls)) {
    pts <- sri_cls[[i]]
    if (length(pts) > 0) {
      jx <- bp_d[i] + stats::runif(length(pts), -0.14, 0.14)
      graphics::points(jx, pts, pch = 19, col = grDevices::adjustcolor("black", alpha.f = 0.75), cex = 0.85)
    }
  }

  # ---- Panel e: Flanking TE & Intergenic Coverage ----
  graphics::par(mar = c(4.2, 4.2, 2.6, 1.2))
  sm_means <- sapply(classes, function(k) mean(as.numeric(te_lg$total_te_pct[te_lg$Class == k])))
  sm_sems  <- sapply(classes, function(k) sem(as.numeric(te_lg$total_te_pct[te_lg$Class == k])))
  ig_means <- sapply(classes, function(k) {
    sub <- te_lg[te_lg$Class == k, ]
    mean(pmax(0, as.numeric(sub$intergenic_noncoding_pct) - as.numeric(sub$total_te_pct)))
  })
  mat_e <- rbind(sm_means, ig_means)
  bp_e <- graphics::barplot(mat_e, col = c("#d73027", "#fddbc7"), border = "black",
                            ylim = c(0, 118), las = 1,
                            names.arg = c("Class I\nCC-NLR", "Class II\nNLR-ID", "Class III\nTKP/TKF", "Class IV\nANK-TM/NAM", "Class V\nTransporter"),
                            ylab = "Flanking window sequence coverage (%)",
                            main = "e  Flanking RepeatMasker TE & intergenic sequence coverage",
                            adj = 0, cex.main = 1.05)
  graphics::arrows(bp_e, sm_means - sm_sems, bp_e, sm_means + sm_sems, angle = 90, code = 3, length = 0.05, lwd = 1.1)
  for (i in seq_along(bp_e)) {
    graphics::text(bp_e[i], sm_means[i] / 2.0, sprintf("%.1f%%", sm_means[i]), col = "white", font = 2, cex = 0.82)
    graphics::text(bp_e[i], sm_means[i] + ig_means[i] + 4.0, sprintf("%.1f%%", sm_means[i] + ig_means[i]), col = "black", font = 2, cex = 0.80)
  }
  graphics::legend("topright", bty = "n", cex = 0.82, fill = c("#d73027", "#fddbc7"),
                   legend = c("RepeatMasker soft-masked TE (%)", "Unmasked intergenic non-coding (%)"))

  # ---- Panel f: Representative Microsynteny Ribbons ----
  graphics::par(mar = c(4.2, 6.8, 2.6, 1.2))
  rep_lgs  <- c("Sr35", "Lr13/Yr27", "Yr15", "Lr14a", "Lr34")
  rep_asms <- c("Tu_G1812", "Aet_AL878", "Tdic_Zavitan", "Tae_CSv21", "Tae_Kariega")
  tracks <- list()
  for (lg in rep_lgs) {
    sub_r <- ribbons[ribbons$LocusGroup == lg, ]
    for (a in rep_asms) {
      ar <- sub_r[sub_r$asm_id == a & abs(as.numeric(sub_r$rel_start_kb)) <= 300, ]
      if (nrow(ar) > 0) {
        short_a <- gsub("^(Tae_|Tdic_|Tu_|Aet_)", "", a)
        short_lg <- strsplit(lg, "/")[[1]][1]
        tracks[[length(tracks) + 1]] <- list(label = sprintf("%s (%s)", short_lg, short_a), rows = ar)
      }
    }
  }
  n_tr <- length(tracks)
  graphics::plot(0, 0, type = "n", xlim = c(-310, 310), ylim = c(n_tr + 0.8, 0.2),
                 yaxt = "n", bty = "l",
                 xlab = "Relative distance from R-gene midpoint (kb)", ylab = "",
                 main = "f  Microsynteny conservation across diploid, tetraploid, and hexaploid assemblies",
                 adj = 0, cex.main = 1.05)
  graphics::abline(v = 0, col = "#e41a1c", lty = 3, lwd = 1.0)
  for (idx in seq_along(tracks)) {
    graphics::segments(-300, idx, 300, idx, col = "#dddddd", lwd = 1.2)
    ar <- tracks[[idx]]$rows
    mids <- 0.5 * (as.numeric(ar$rel_start_kb) + as.numeric(ar$rel_end_kb))
    is_rg <- ar$feature_type == "R_gene"
    if (any(!is_rg)) {
      cols_f <- ifelse(as.character(ar$collinear_to_ref[!is_rg]) == "1", "#2b8cbe", "#bababa")
      graphics::points(mids[!is_rg], rep(idx, sum(!is_rg)), pch = 15, col = cols_f, cex = 0.85)
    }
    if (any(is_rg)) {
      graphics::points(mids[is_rg], rep(idx, sum(is_rg)), pch = 18, col = "#e41a1c", cex = 1.35)
    }
  }
  graphics::axis(2, at = seq_len(n_tr), labels = sapply(tracks, `[[`, "label"), las = 1, cex.axis = 0.72)
  graphics::legend("bottomright", bty = "n", cex = 0.78,
                   legend = c("Cloned R-gene locus (0 kb)", "Collinear flanking orthologue", "Non-collinear / PAV gene"),
                   pch = c(18, 15, 15), col = c("#e41a1c", "#2b8cbe", "#bababa"), pt.cex = c(1.2, 0.9, 0.9))
}

grDevices::pdf("figures/Fig2_chromosomal_distribution_R.pdf", width = 15.5, height = 11.8)
draw_fig2()
invisible(grDevices::dev.off())

grDevices::png("figures/Fig2_chromosomal_distribution_R.png", width = 15.5, height = 11.8, units = "in", res = 300)
draw_fig2()
invisible(grDevices::dev.off())

# ---- Supplementary Fig. 1 (Subgenome CNV) ----
draw_sfig1 <- function() {
  graphics::par(mar = c(4.5, 4.5, 3.0, 1.2), cex = 0.9)
  tiers <- c("A", "B/S", "D")
  nlr_m <- sapply(tiers, function(t) mean(as.numeric(sg_cnv$intact_copies[sg_cnv$Group == "NLR" & sg_cnv$tier == t])))
  nlr_s <- sapply(tiers, function(t) sem(as.numeric(sg_cnv$intact_copies[sg_cnv$Group == "NLR" & sg_cnv$tier == t])))
  non_m <- sapply(tiers, function(t) mean(as.numeric(sg_cnv$intact_copies[sg_cnv$Group == "nonNLR" & sg_cnv$tier == t])))
  non_s <- sapply(tiers, function(t) sem(as.numeric(sg_cnv$intact_copies[sg_cnv$Group == "nonNLR" & sg_cnv$tier == t])))
  mat_s1 <- rbind(nlr_m, non_m)
  bp_s1 <- graphics::barplot(mat_s1, beside = TRUE, col = c("#1f77b4", "#d95f02"), border = "black",
                             ylim = c(0, max(nlr_m + nlr_s) * 1.22), las = 1,
                             names.arg = c("A subgenome (n = 19)", "B/S subgenome (n = 21)", "D subgenome (n = 16)"),
                             ylab = "Intact copies per haploid subgenome (mean +/- s.e.m.)",
                             main = "Supplementary Fig. 1 | Subgenome copy-number distribution across A, B/S, and D lineages",
                             adj = 0, cex.main = 0.95)
  graphics::arrows(bp_s1[1, ], nlr_m - nlr_s, bp_s1[1, ], nlr_m + nlr_s, angle = 90, code = 3, length = 0.05, lwd = 1.2)
  graphics::arrows(bp_s1[2, ], non_m - non_s, bp_s1[2, ], non_m + non_s, angle = 90, code = 3, length = 0.05, lwd = 1.2)
  n_og_nlr <- length(unique(sg_cnv$orthogroup[sg_cnv$Group == "NLR"]))
  n_og_non <- length(unique(sg_cnv$orthogroup[sg_cnv$Group == "nonNLR"]))
  graphics::legend("topright", bty = "n", fill = c("#1f77b4", "#d95f02"),
                   legend = c(sprintf("Classes I-II NLR orthogroups (n = %d)", n_og_nlr),
                              sprintf("Classes III-V non-NLR orthogroups (n = %d)", n_og_non)))
}

grDevices::pdf("figures/Supplementary_Fig1_subgenome_CNV_R.pdf", width = 7.2, height = 4.8)
draw_sfig1()
invisible(grDevices::dev.off())

grDevices::png("figures/Supplementary_Fig1_subgenome_CNV_R.png", width = 7.2, height = 4.8, units = "in", res = 300)
draw_sfig1()
invisible(grDevices::dev.off())

# ---- Supplementary Fig. 2 (SV & TE vs SRI) ----
draw_sfig2 <- function() {
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.0, 1.2), cex = 0.88)
  cols_pt <- ifelse(te_lg$Group == "NLR", "#1f77b4", "#d95f02")
  graphics::plot(all_sv, all_sri, pch = 21, bg = grDevices::adjustcolor(cols_pt, alpha.f = 0.85),
                 col = "black", cex = 1.25, bty = "l",
                 xlab = "Mean local SV / gene rearrangement events (+/-500 kb)",
                 ylab = "Synteny Retention Index (SRI)",
                 main = sprintf("a  Local rearrangements vs. SRI (rho = %.2f, P = %.2e)",
                                sp_sv_sri$estimate, sp_sv_sri$p.value),
                 adj = 0, cex.main = 0.98)
  graphics::abline(stats::lm(all_sri ~ all_sv), col = "black", lty = 2, lwd = 1.3)

  graphics::plot(all_te, all_sv, pch = 21, bg = grDevices::adjustcolor(cols_pt, alpha.f = 0.85),
                 col = "black", cex = 1.25, bty = "l",
                 xlab = "Flanking +/-200 kb TE coverage (%)",
                 ylab = "Mean local SV / gene rearrangement events",
                 main = sprintf("b  TE density vs. rearrangements (rho = %.2f, P = %.2e)",
                                sp_te_sv$estimate, sp_te_sv$p.value),
                 adj = 0, cex.main = 0.98)
  graphics::abline(stats::lm(all_sv ~ all_te), col = "black", lty = 2, lwd = 1.3)
  graphics::legend("topleft", bty = "n", pch = 21, pt.bg = c("#1f77b4", "#d95f02"), col = "black",
                   legend = c("Classes I-II NLR", "Classes III-V non-NLR"))
}

grDevices::pdf("figures/Supplementary_Fig2_TE_vs_SRI_R.pdf", width = 11.8, height = 4.6)
draw_sfig2()
invisible(grDevices::dev.off())

grDevices::png("figures/Supplementary_Fig2_TE_vs_SRI_R.png", width = 11.8, height = 4.6, units = "in", res = 300)
draw_sfig2()
invisible(grDevices::dev.off())

cat("Saved R statistical audit table -> results/S7c_R_statistical_audit_3.2.tsv\n")
cat("Saved R publication figures -> figures/Fig2_chromosomal_distribution_R.pdf/.png,\n")
cat("   figures/Supplementary_Fig1_subgenome_CNV_R.pdf/.png, and\n")
cat("   figures/Supplementary_Fig2_TE_vs_SRI_R.pdf/.png\n")
cat("======================================================================\n")
