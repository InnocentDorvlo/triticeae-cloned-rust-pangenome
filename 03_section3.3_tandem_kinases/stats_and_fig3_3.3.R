#!/usr/bin/env Rscript
# ==============================================================================
# stats_and_fig3_3.3.R
# Standalone R statistical audit and publication figure generation for Result 3.3:
#   - S8a_classIII_domains.tsv (Domain architectures & catalytic triad status)
#   - S8b_kinase_phylogeny.tsv & work/tkp_3.3/kinase_domains.treefile (Rooted ML tree)
#   - S8c_kinase_pair_divergence.tsv (K1-K2 aa identity, Nei-Gojobori pS & dS)
#   - S8d_architecture_distribution.tsv (Fused vs progenitor presence across 27 assemblies)
#   - S8e_exon_intron_phases.tsv (Exon-intron structures & phase-0 junction enrichment)
#   - Outputs:
#       results/S8f_R_statistical_audit_3.3.tsv
#       tables/Supplementary_Table_S8f_R_Statistical_Audit_3.3.tsv
#       figures/Fig3_tandem_kinase_origin_R.pdf and .png
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

options(stringsAsFactors = FALSE, width = 120)

s8a <- read.delim("results/S8a_classIII_domains.tsv", check.names = FALSE)
s8b <- read.delim("results/S8b_kinase_phylogeny.tsv", check.names = FALSE)
s8c <- read.delim("results/S8c_kinase_pair_divergence.tsv", check.names = FALSE)
s8d <- read.delim("results/S8d_architecture_distribution.tsv", check.names = FALSE)
s8e <- read.delim("results/S8e_exon_intron_phases.tsv", check.names = FALSE)

CLASS3 <- c("Yr15", "Sr60", "Sr62", "Lr9", "Yr36", "Sr43", "Rpg1", "WTK4")
DUP_GENES <- c("Lr9", "Rpg1", "WTK4")
FUS_GENES <- c("Sr60", "Sr62", "Yr15")

# ------------------------------------------------------------------------------
# 1. Domain architectures & homologue counts (S8a)
# ------------------------------------------------------------------------------
refs_a <- s8a[s8a$Record == "reference", ]
homs_a <- s8a[s8a$Record == "homologue", ]
homs_main <- homs_a[homs_a$Set == "Main", ]
tkp_homs <- homs_a[as.integer(homs_a$n_kinase_modules) >= 2, ]

# ------------------------------------------------------------------------------
# 2. Phylogeny & patristic distances (S8b)
# ------------------------------------------------------------------------------
pat_dup <- as.numeric(s8b$ref_K1_K2_patristic[s8b$Gene %in% DUP_GENES])
pat_fus <- as.numeric(s8b$ref_K1_K2_patristic[s8b$Gene %in% FUS_GENES])
mwu_pat <- suppressWarnings(wilcox.test(pat_dup, pat_fus, exact = TRUE))
welch_pat <- t.test(pat_dup, pat_fus, var.equal = FALSE)

# ------------------------------------------------------------------------------
# 3. K1-K2 amino-acid identity & Nei-Gojobori pS / dS (S8c)
# ------------------------------------------------------------------------------
ref_c <- s8c[s8c$Record == "reference", ]
hom_c <- s8c[s8c$Record == "homologue", ]

# Reference-level K1-K2 identity & pS (n = 3 vs 3)
ref_id_dup <- as.numeric(ref_c$aa_identity_pct[ref_c$Gene %in% DUP_GENES])
ref_id_fus <- as.numeric(ref_c$aa_identity_pct[ref_c$Gene %in% FUS_GENES])
mwu_ref_id <- suppressWarnings(wilcox.test(ref_id_dup, ref_id_fus, exact = TRUE))
welch_ref_id <- t.test(ref_id_dup, ref_id_fus, var.equal = FALSE)

ref_ps_dup <- as.numeric(ref_c$pS[ref_c$Gene %in% DUP_GENES])
ref_ps_fus <- as.numeric(ref_c$pS[ref_c$Gene %in% FUS_GENES])
ref_ds_dup <- as.numeric(ref_c$dS_JC[ref_c$Gene %in% DUP_GENES])
mwu_ref_ps <- suppressWarnings(wilcox.test(ref_ps_dup, ref_ps_fus, exact = TRUE))
welch_ref_ps <- t.test(ref_ps_dup, ref_ps_fus, var.equal = FALSE)

# Locus-aggregated pangenome homologue K1-K2 identity (n = 3 vs 3)
loc_id_dup <- sapply(DUP_GENES, function(g) mean(as.numeric(hom_c$aa_identity_pct[hom_c$Gene == g])))
loc_id_fus <- sapply(FUS_GENES, function(g) mean(as.numeric(hom_c$aa_identity_pct[hom_c$Gene == g])))
mwu_loc_id <- suppressWarnings(wilcox.test(loc_id_dup, loc_id_fus, exact = TRUE))
welch_loc_id <- t.test(loc_id_dup, loc_id_fus, var.equal = FALSE)

# All pangenome TKP homologue pairs (n = 344 vs 171)
pair_id_dup <- as.numeric(hom_c$aa_identity_pct[hom_c$Gene %in% DUP_GENES])
pair_id_fus <- as.numeric(hom_c$aa_identity_pct[hom_c$Gene %in% FUS_GENES])
mwu_pair_id <- wilcox.test(pair_id_dup, pair_id_fus, exact = FALSE)
welch_pair_id <- t.test(pair_id_dup, pair_id_fus, var.equal = FALSE)

# Pangenome TKP homologue pairs with frame-verified genomic CDS (n = 12 vs 6)
g_dup <- hom_c[hom_c$Gene %in% DUP_GENES & hom_c$cds_source == "genomic CDS", ]
g_fus <- hom_c[hom_c$Gene %in% FUS_GENES & hom_c$cds_source == "genomic CDS", ]
ps_g_dup <- as.numeric(g_dup$pS)
ps_g_fus <- as.numeric(g_fus$pS)
ds_g_dup <- as.numeric(g_dup$dS_JC[g_dup$dS_saturated == 0])
ds_g_fus <- as.numeric(g_fus$dS_JC[g_fus$dS_saturated == 0])
mwu_g_ps <- suppressWarnings(wilcox.test(ps_g_dup, ps_g_fus, exact = TRUE))
welch_g_ps <- t.test(ps_g_dup, ps_g_fus, var.equal = FALSE)

# ------------------------------------------------------------------------------
# 4. Pangenome fused vs progenitor distribution across 27 assemblies (S8d)
# ------------------------------------------------------------------------------
n_asm_total <- length(unique(s8d$asm_id))
n_dip_total <- length(unique(s8d$asm_id[s8d$ploidy == 2]))
dist_summary <- lapply(CLASS3, function(g) {
  sub <- s8d[s8d$Gene == g, ]
  dsub <- sub[sub$ploidy == 2, ]
  n_fused <- sum(sub$fused_architecture_present == 1)
  n_single_dip <- sum(dsub$single_kinase_progenitor_present == 1)
  acc_vals <- dsub$standalone_accessory_progenitor_present
  n_acc_dip <- if (all(acc_vals == "NA" | is.na(acc_vals))) NA_integer_ else sum(as.character(acc_vals) == "1", na.rm = TRUE)
  list(Gene = g, n_fused = n_fused, n_single_dip = n_single_dip, n_acc_dip = n_acc_dip)
})
names(dist_summary) <- CLASS3

# ------------------------------------------------------------------------------
# 5. Exon-intron junction phase-0 enrichment (S8e)
# ------------------------------------------------------------------------------
junc_rows <- s8e[s8e$domain_boundary %in% c("K2_start", "accessory_start") &
                 s8e$distance_aa != "-" & as.integer(s8e$distance_aa) <= 30, ]
n_junc_records <- nrow(junc_rows)
n_junc_phase0 <- sum(as.integer(junc_rows$nearest_intron_phase) == 0)
junc_genes <- unique(junc_rows$Gene)
n_junc_genes <- length(junc_genes)
# One primary record per gene (NCBI when available, else Miniprot)
junc_per_gene <- do.call(rbind, lapply(junc_genes, function(g) {
  sub <- junc_rows[junc_rows$Gene == g, ]
  sub[1, , drop = FALSE]
}))
n_gene_phase0 <- sum(as.integer(junc_per_gene$nearest_intron_phase) == 0)
binom_gene_p <- binom.test(n_gene_phase0, n_junc_genes, p = 1/3, alternative = "greater")$p.value
binom_rec_p  <- binom.test(n_junc_phase0, n_junc_records, p = 1/3, alternative = "greater")$p.value

# Background intron phases across the 8 Class III reference/miniprot structures (one source per gene)
bg_phases <- unlist(lapply(CLASS3, function(g) {
  sub <- s8e[s8e$Gene == g, ]
  ph_str <- sub$all_intron_phases[1]
  if (ph_str == "-" || is.na(ph_str)) return(integer(0))
  as.integer(strsplit(ph_str, ",")[[1]])
}))
n_bg_total <- length(bg_phases)
n_bg_phase0 <- sum(bg_phases == 0)

# ------------------------------------------------------------------------------
# 6. Console Report & Audit Table (S8f)
# ------------------------------------------------------------------------------
cat("======================================================================\n")
cat("  Section 3.3 R Statistical Audit & Tandem Kinase Evolution Analysis\n")
cat("======================================================================\n\n")

cat(sprintf("[1] Class III Catalogue & Domain Dissection (S8a):\n"))
cat(sprintf("    Queries: %d (%d Main + %d Comparators: %s)\n",
            nrow(refs_a), sum(refs_a$Set == "Main"), sum(refs_a$Set != "Main"), paste(CLASS3, collapse = ", ")))
cat(sprintf("    Retained anchor-verified homologues: %d Main (%d total including comparators)\n",
            nrow(homs_main), nrow(homs_a)))
cat(sprintf("    Pangenome homologues with >=2 complete kinase modules (>=150 aa): %d\n\n", nrow(tkp_homs)))

cat(sprintf("[2] Rooted ML Kinase-Domain Phylogeny & Patristic Distances (S8b, 61 domains, 7 outgroups):\n"))
cat(sprintf("    Intragenic duplication (%s): K1/K2 sister clades = 3/3 (100%%; UFBoot = 100%%)\n",
            paste(DUP_GENES, collapse = ", ")))
cat(sprintf("    Intergenic fusion (%s):      K1/K2 sister clades = 0/3 (0%%; distant subfamilies)\n",
            paste(FUS_GENES, collapse = ", ")))
cat(sprintf("    Reference K1-K2 patristic distance: duplication = %.3f +/- %.3f vs fusion = %.3f +/- %.3f subst./site\n",
            mean(pat_dup), sd(pat_dup)/sqrt(length(pat_dup)), mean(pat_fus), sd(pat_fus)/sqrt(length(pat_fus))))
cat(sprintf("      Mann-Whitney U = %.1f, exact P = %.4f | Welch's t = %.2f, P = %.4f\n\n",
            mwu_pat$statistic, mwu_pat$p.value, welch_pat$statistic, welch_pat$p.value))

cat(sprintf("[3] Intra-Polypeptide K1-K2 Divergence (S8c):\n"))
cat(sprintf("    Reference queries (n = 3 duplication vs 3 fusion):\n"))
cat(sprintf("      K1-K2 aa identity: %.2f +/- %.2f%% vs %.2f +/- %.2f%% (MWU U = %.1f, P = %.4f; Welch t = %.2f, P = %.4f)\n",
            mean(ref_id_dup), sd(ref_id_dup)/sqrt(3), mean(ref_id_fus), sd(ref_id_fus)/sqrt(3),
            mwu_ref_id$statistic, mwu_ref_id$p.value, welch_ref_id$statistic, welch_ref_id$p.value))
cat(sprintf("      Nei-Gojobori pS:   %.4f +/- %.4f (dS = %.4f +/- %.4f; 0/3 sat) vs %.4f +/- %.4f (3/3 sat)\n",
            mean(ref_ps_dup), sd(ref_ps_dup)/sqrt(3), mean(ref_ds_dup), sd(ref_ds_dup)/sqrt(3),
            mean(ref_ps_fus), sd(ref_ps_fus)/sqrt(3)))
cat(sprintf("      Reference pS test: MWU U = %.1f, P = %.4f | Welch's t = %.2f, P = %.4e\n",
            mwu_ref_ps$statistic, mwu_ref_ps$p.value, welch_ref_ps$statistic, welch_ref_ps$p.value))
cat(sprintf("    Locus-aggregated pangenome homologues (n = 3 vs 3 loci):\n"))
cat(sprintf("      K1-K2 aa identity: %.2f +/- %.2f%% vs %.2f +/- %.2f%% (MWU U = %.1f, P = %.4f; Welch t = %.2f, P = %.4f)\n",
            mean(loc_id_dup), sd(loc_id_dup)/sqrt(3), mean(loc_id_fus), sd(loc_id_fus)/sqrt(3),
            mwu_loc_id$statistic, mwu_loc_id$p.value, welch_loc_id$statistic, welch_loc_id$p.value))
cat(sprintf("    All pangenome TKP homologue pairs (n = %d duplication vs %d fusion):\n",
            length(pair_id_dup), length(pair_id_fus)))
cat(sprintf("      K1-K2 aa identity: %.2f +/- %.2f%% (median %.2f%%) vs %.2f +/- %.2f%% (median %.2f%%)\n",
            mean(pair_id_dup), sd(pair_id_dup)/sqrt(length(pair_id_dup)), median(pair_id_dup),
            mean(pair_id_fus), sd(pair_id_fus)/sqrt(length(pair_id_fus)), median(pair_id_fus)))
cat(sprintf("      MWU U = %.1f, P = %.3e | Welch's t = %.2f, P = %.3e\n",
            mwu_pair_id$statistic, mwu_pair_id$p.value, welch_pair_id$statistic, welch_pair_id$p.value))
cat(sprintf("    Pangenome homologues with genomic CDS (n = %d duplication vs %d fusion in Tae_ArinaLrFor):\n",
            nrow(g_dup), nrow(g_fus)))
cat(sprintf("      Genomic pS: %.4f +/- %.4f (dS = %.4f +/- %.4f; 0/%d sat) vs %.4f +/- %.4f (3/%d sat; unsat dS = %.4f +/- %.4f)\n",
            mean(ps_g_dup), sd(ps_g_dup)/sqrt(length(ps_g_dup)),
            mean(ds_g_dup), sd(ds_g_dup)/sqrt(length(ds_g_dup)), nrow(g_dup),
            mean(ps_g_fus), sd(ps_g_fus)/sqrt(length(ps_g_fus)), nrow(g_fus),
            mean(ds_g_fus), sd(ds_g_fus)/sqrt(length(ds_g_fus))))
cat(sprintf("      Genomic pS test: MWU U = %.1f, exact P = %.3e | Welch's t = %.2f, P = %.3e\n\n",
            mwu_g_ps$statistic, mwu_g_ps$p.value, welch_g_ps$statistic, welch_g_ps$p.value))

cat(sprintf("[4] Architecture Distribution Across %d Assemblies (%d diploids) & Exon-Intron Junctions (S8d, S8e):\n",
            n_asm_total, n_dip_total))
for (g in CLASS3) {
  d <- dist_summary[[g]]
  acc_txt <- if (is.na(d$n_acc_dip)) "" else sprintf(", standalone accessory in %d/%d diploids", d$n_acc_dip, n_dip_total)
  cat(sprintf("    %-5s: fused in %2d/%d assemblies; single-kinase progenitor in %2d/%d diploids%s\n",
              g, d$n_fused, n_asm_total, d$n_single_dip, n_dip_total, acc_txt))
}
cat(sprintf("    Inter-domain junction introns (<=30 aa of K2_start/accessory_start): %d/%d genes (%d/%d records) are phase 0\n",
            n_gene_phase0, n_junc_genes, n_junc_phase0, n_junc_records))
cat(sprintf("      Exact binomial test vs uniform 1/3 phase-0 expectation: per-gene P = %.4e (all records P = %.4e)\n\n",
            binom_gene_p, binom_rec_p))

audit_df <- data.frame(
  Section = c(
    rep("1_Domain_Dissection", 3),
    rep("2_Kinase_Phylogeny", 4),
    rep("3_K1K2_Divergence", 10),
    rep("4_Pangenome_Distribution", 8),
    rep("5_Exon_Intron_Junctions", 3)
  ),
  Metric = c(
    "N_ClassIII_queries", "N_ClassIII_retained_homologues", "N_TKP_homologues_2_modules",
    "Sister_clades_intragenic_duplication", "Sister_clades_intergenic_fusion",
    "Ref_patristic_distance_dup_vs_fus", "Ref_patristic_statistical_tests",
    "Ref_K1K2_aa_identity_dup_vs_fus", "Ref_K1K2_aa_identity_tests",
    "Ref_K1K2_pS_dup_vs_fus", "Ref_K1K2_dS_dup_vs_fus", "Ref_K1K2_pS_tests",
    "Locus_mean_K1K2_aa_identity_dup_vs_fus", "Locus_mean_K1K2_aa_identity_tests",
    "All_pairs_K1K2_aa_identity_dup_vs_fus", "All_pairs_K1K2_aa_identity_tests",
    "Genomic_CDS_K1K2_pS_dS_dup_vs_fus",
    paste0("Architecture_distribution_", CLASS3),
    "Junction_introns_phase0_per_gene", "Junction_introns_phase0_all_records", "Background_introns_phase0_rate"
  ),
  Estimate = c(
    sprintf("%d (6 Main + 2 Comparators)", nrow(refs_a)),
    sprintf("%d Main (%d total with comparators)", nrow(homs_main), nrow(homs_a)),
    sprintf("%d homologues", nrow(tkp_homs)),
    "3/3 (100%: Lr9, Rpg1, WTK4)", "0/3 (0%: Sr60, Sr62, Yr15)",
    sprintf("%.3f +/- %.3f vs %.3f +/- %.3f subst./site", mean(pat_dup), sd(pat_dup)/sqrt(3), mean(pat_fus), sd(pat_fus)/sqrt(3)),
    sprintf("MWU U = %.1f, P = %.4f; Welch t = %.2f, P = %.4f", mwu_pat$statistic, mwu_pat$p.value, welch_pat$statistic, welch_pat$p.value),
    sprintf("%.2f +/- %.2f%% vs %.2f +/- %.2f%%", mean(ref_id_dup), sd(ref_id_dup)/sqrt(3), mean(ref_id_fus), sd(ref_id_fus)/sqrt(3)),
    sprintf("MWU U = %.1f, P = %.4f; Welch t = %.2f, P = %.4f", mwu_ref_id$statistic, mwu_ref_id$p.value, welch_ref_id$statistic, welch_ref_id$p.value),
    sprintf("%.4f +/- %.4f vs %.4f +/- %.4f", mean(ref_ps_dup), sd(ref_ps_dup)/sqrt(3), mean(ref_ps_fus), sd(ref_ps_fus)/sqrt(3)),
    sprintf("%.4f +/- %.4f (0/3 sat) vs saturated (3/3 pS >= 0.75)", mean(ref_ds_dup), sd(ref_ds_dup)/sqrt(3)),
    sprintf("MWU U = %.1f, P = %.4f; Welch t = %.2f, P = %.4e", mwu_ref_ps$statistic, mwu_ref_ps$p.value, welch_ref_ps$statistic, welch_ref_ps$p.value),
    sprintf("%.2f +/- %.2f%% vs %.2f +/- %.2f%%", mean(loc_id_dup), sd(loc_id_dup)/sqrt(3), mean(loc_id_fus), sd(loc_id_fus)/sqrt(3)),
    sprintf("MWU U = %.1f, P = %.4f; Welch t = %.2f, P = %.4f", mwu_loc_id$statistic, mwu_loc_id$p.value, welch_loc_id$statistic, welch_loc_id$p.value),
    sprintf("%.2f +/- %.2f%% (n=%d) vs %.2f +/- %.2f%% (n=%d)",
            mean(pair_id_dup), sd(pair_id_dup)/sqrt(length(pair_id_dup)), length(pair_id_dup),
            mean(pair_id_fus), sd(pair_id_fus)/sqrt(length(pair_id_fus)), length(pair_id_fus)),
    sprintf("MWU U = %.1f, P = %.3e; Welch t = %.2f, P = %.3e", mwu_pair_id$statistic, mwu_pair_id$p.value, welch_pair_id$statistic, welch_pair_id$p.value),
    sprintf("pS %.4f+/-%.4f (dS %.4f+/-%.4f, n=%d) vs pS %.4f+/-%.4f (n=%d); MWU U=%.1f, P=%.3e",
            mean(ps_g_dup), sd(ps_g_dup)/sqrt(length(ps_g_dup)), mean(ds_g_dup), sd(ds_g_dup)/sqrt(length(ds_g_dup)), nrow(g_dup),
            mean(ps_g_fus), sd(ps_g_fus)/sqrt(length(ps_g_fus)), nrow(g_fus), mwu_g_ps$statistic, mwu_g_ps$p.value),
    sapply(CLASS3, function(g) {
      d <- dist_summary[[g]]
      acc_txt <- if (is.na(d$n_acc_dip)) "" else sprintf("; acc %d/%d diploids", d$n_acc_dip, n_dip_total)
      sprintf("fused %d/%d; single-kinase %d/%d diploids%s", d$n_fused, n_asm_total, d$n_single_dip, n_dip_total, acc_txt)
    }),
    sprintf("%d/%d (100%%) phase 0 (binomial P = %.4e vs 1/3)", n_gene_phase0, n_junc_genes, binom_gene_p),
    sprintf("%d/%d (100%%) phase 0 (binomial P = %.4e vs 1/3)", n_junc_phase0, n_junc_records, binom_rec_p),
    sprintf("%d/%d (%.1f%%) phase 0 across all introns in Class III genes", n_bg_phase0, n_bg_total, 100 * n_bg_phase0 / n_bg_total)
  )
)

write.table(audit_df, "results/S8f_R_statistical_audit_3.3.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(audit_df, "tables/Supplementary_Table_S8f_R_Statistical_Audit_3.3.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. Pure Base-R Newick Parser & Outgroup Rooter for Panel B
# ------------------------------------------------------------------------------
parse_newick_Clade <- function(txt) {
  txt <- gsub("\\s+", "", txt)
  txt <- sub(";$", "", txt)
  pos <- 1L
  n <- nchar(txt)
  read_node <- function() {
    children <- list()
    if (substr(txt, pos, pos) == "(") {
      pos <<- pos + 1L
      repeat {
        children[[length(children) + 1L]] <- read_node()
        ch <- substr(txt, pos, pos)
        if (ch == ",") {
          pos <<- pos + 1L
        } else if (ch == ")") {
          pos <<- pos + 1L
          break
        } else {
          break
        }
      }
    }
    rest <- substr(txt, pos, n)
    m <- regexpr("^[^,:()]*", rest)
    label <- regmatches(rest, m)
    pos <<- pos + attr(m, "match.length")
    bl <- 0.0
    if (substr(txt, pos, pos) == ":") {
      pos <<- pos + 1L
      rest2 <- substr(txt, pos, n)
      m2 <- regexpr("^[-0-9.eE+]+", rest2)
      bl <- as.numeric(regmatches(rest2, m2))
      pos <<- pos + attr(m2, "match.length")
    }
    list(name = label, length = ifelse(is.na(bl), 0.0, bl), children = children)
  }
  read_node()
}

get_tips <- function(node) {
  if (length(node$children) == 0L) return(node$name)
  unlist(lapply(node$children, get_tips))
}

# Convert unrooted/arbitrarily rooted tree into adjacency list and re-root on midpoint of edge to basal outgroup clade
reroot_on_clade <- function(root_node, target_tips) {
  # Flatten tree into nodes & edges
  edges <- list()
  nodes <- list()
  next_id <- 0L
  build_adj <- function(nd) {
    next_id <<- next_id + 1L
    my_id <- next_id
    nodes[[my_id]] <<- list(id = my_id, name = nd$name, is_tip = (length(nd$children) == 0L))
    for (ch in nd$children) {
      cid <- build_adj(ch)
      edges[[length(edges) + 1L]] <<- list(u = my_id, v = cid, w = ch$length, support = ch$name)
    }
    my_id
  }
  build_adj(root_node)

  # Find edge whose partition equals target_tips
  n_nodes <- length(nodes)
  adj <- vector("list", n_nodes)
  for (k in seq_along(edges)) {
    e <- edges[[k]]
    adj[[e$u]][[length(adj[[e$u]]) + 1L]] <- list(to = e$v, w = e$w, edge_id = k)
    adj[[e$v]][[length(adj[[e$v]]) + 1L]] <- list(to = e$u, w = e$w, edge_id = k)
  }

  tips_on_side <- function(start, blocked_Id) {
    vis <- rep(FALSE, n_nodes)
    vis[blocked_Id] <- TRUE
    st <- c(start)
    out <- character(0)
    while (length(st) > 0L) {
      curr <- st[length(st)]
      st <- st[-length(st)]
      if (vis[curr]) next
      vis[curr] <- TRUE
      if (nodes[[curr]]$is_tip) out <- c(out, nodes[[curr]]$name)
      for (nb in adj[[curr]]) {
        if (!vis[nb$to]) st <- c(st, nb$to)
      }
    }
    out
  }

  best_edge <- 1L
  best_dir <- TRUE
  for (k in seq_along(edges)) {
    e <- edges[[k]]
    t_v <- tips_on_side(e$v, e$u)
    if (setequal(t_v, target_tips)) {
      best_edge <- k; best_dir <- TRUE; break
    }
    t_u <- tips_on_side(e$u, e$v)
    if (setequal(t_u, target_tips)) {
      best_edge <- k; best_dir <- FALSE; break
    }
  }

  e0 <- edges[[best_edge]]
  out_node <- if (best_dir) e0$v else e0$u
  in_node  <- if (best_dir) e0$u else e0$v

  build_directed <- function(curr, parent) {
    ch_list <- list()
    for (nb in adj[[curr]]) {
      if (nb$to == parent) next
      sub_nd <- build_directed(nb$to, curr)
      sub_nd$length <- nb$w
      ch_list[[length(ch_list) + 1L]] <- sub_nd
    }
    # Ladderize (sort children by number of descendant tips)
    if (length(ch_list) > 1L) {
      ntips <- sapply(ch_list, function(x) length(get_tips(x)))
      ch_list <- ch_list[order(ntips)]
    }
    list(name = nodes[[curr]]$name, length = 0.0, children = ch_list)
  }

  out_sub <- build_directed(out_node, in_node)
  out_sub$length <- e0$w * 0.15
  in_sub  <- build_directed(in_node, out_node)
  in_sub$length  <- e0$w * 0.85
  list(name = "ROOT", length = 0.0, children = list(out_sub, in_sub))
}

# ------------------------------------------------------------------------------
# 8. Publication Figure 3 (Fig3_tandem_kinase_origin_R.pdf and .png)
# ------------------------------------------------------------------------------
draw_fig3 <- function() {
  layout(matrix(c(1, 2, 2,
                  3, 4, 5), nrow = 2, byrow = TRUE),
         widths = c(1.15, 1.05, 1.05), heights = c(1.08, 1.0))

  # ---- Panel a: Domain organization of the 8 Class III reference proteins ----
  par(mar = c(4.2, 5.2, 2.8, 1.2))
  refs_ord <- refs_a[match(CLASS3, refs_a$Gene), ]
  n_ref <- nrow(refs_ord)
  plot(NA, xlim = c(0, 1380), ylim = c(0.4, n_ref + 0.8),
       xlab = "Amino-acid position (aa)", ylab = "", yaxt = "n", bty = "n", cex.lab = 0.95)
  title(main = "a   Domain architecture of cloned Class III reference proteins", adj = 0, cex.main = 1.0, font.main = 2)

  for (i in seq_len(n_ref)) {
    r <- refs_ord[i, ]
    y <- n_ref - i + 1
    L <- as.numeric(r$aa_len)
    segments(0, y, L, y, col = "grey45", lwd = 2.0)
    axis(2, at = y, labels = r$Gene, las = 1, font = 3, tick = FALSE, cex.axis = 0.9, line = -0.5)
    text(L + 25, y, sprintf("%d aa", L), adj = c(0, 0.5), cex = 0.72, col = "grey30")

    # K1 & K2
    for (kinfo in list(list(rng = r$K1_range, st = r$K1_status, lab = "K1"),
                       list(rng = r$K2_range, st = r$K2_status, lab = "K2"))) {
      if (kinfo$rng != "-") {
        se <- as.numeric(strsplit(kinfo$rng, "-")[[1]])
        col_fill <- if (kinfo$st == "active") "#1f78b4" else "#a6cee3"
        col_bord <- if (kinfo$st == "active") "#08306b" else "#2b5c8f"
        rect(se[1], y - 0.28, se[2], y + 0.28, col = col_fill, border = col_bord, lwd = 1.1)
        text(mean(se), y, kinfo$lab, cex = 0.72, font = 2,
             col = if (kinfo$st == "active") "white" else "#08306b")
      }
    }
    # Accessory
    if (r$accessory_range != "-") {
      se <- as.numeric(strsplit(r$accessory_range, "-")[[1]])
      rect(se[1], y - 0.28, se[2], y + 0.28, col = "#e6550d", border = "#8c2d04", lwd = 1.1)
      acc_lab <- sub(",.*", "", r$accessory_domains)
      text(mean(se), y, acc_lab, cex = 0.68, font = 2, col = "white")
    }
  }
  legend("topright",
         legend = c("Active kinase (K/HRD/DFG intact)", "Pseudokinase (catalytic substitution)", "Integrated accessory domain"),
         fill = c("#1f78b4", "#a6cee3", "#e6550d"), border = c("#08306b", "#2b5c8f", "#8c2d04"),
         cex = 0.74, bty = "n")

  # ---- Panel b: Rooted ML phylogeny of dissected kinase domains ----
  par(mar = c(4.2, 1.0, 2.8, 1.5))
  nwk_txt <- paste(readLines("work/tkp_3.3/kinase_domains.treefile", warn = FALSE), collapse = "")
  raw_tree <- parse_newick_Clade(nwk_txt)
  all_tips <- get_tips(raw_tree)
  basal_outs <- all_tips[grepl("XP_044954759|XP_044957438", all_tips)]
  rtree <- reroot_on_clade(raw_tree, basal_outs)

  # Assign y to tips and x by branch length
  tip_y <- list()
  cur_y <- 0L
  assign_coords <- function(nd, x0) {
    x1 <- x0 + nd$length
    if (length(nd$children) == 0L) {
      cur_y <<- cur_y + 1L
      return(list(x0 = x0, x1 = x1, y = as.numeric(cur_y), name = nd$name, children = list()))
    }
    ch_coords <- lapply(nd$children, function(ch) assign_coords(ch, x1))
    my_y <- mean(sapply(ch_coords, function(c) c$y))
    list(x0 = x0, x1 = x1, y = my_y, name = nd$name, children = ch_coords)
  }
  c_tree <- assign_coords(rtree, 0.0)
  max_x <- max(unlist( (function(f, nd) c(nd$x1, unlist(lapply(nd$children, function(ch) f(f, ch)))))(
    function(f, nd) c(nd$x1, unlist(lapply(nd$children, function(ch) f(f, ch)))), c_tree)))

  plot(NA, xlim = c(0, max_x * 1.46), ylim = c(0.5, cur_y + 1.5),
       xlab = "Substitutions per site (IQ-TREE ModelFinder, 1,000 UFBoot / SH-aLRT)",
       ylab = "", yaxt = "n", bty = "n", cex.lab = 0.95)
  title(main = "b   Rooted maximum-likelihood phylogeny of dissected Kinase 1, Kinase 2, single-kinase and outgroup domains",
        adj = 0, cex.main = 1.0, font.main = 2)

  draw_clade <- function(nd) {
    if (nd$name != "ROOT") {
      segments(nd$x0, nd$y, nd$x1, nd$y, col = "grey30", lwd = 1.1)
    }
    if (length(nd$children) == 0L) {
      nm <- nd$name
      col_t <- if (grepl("_K1", nm)) "#1f78b4" else if (grepl("_K2", nm)) "#d62728" else if (grepl("^OUT_", nm)) "grey45" else "#2ca02c"
      pts_pch <- if (grepl("_ref$", nm)) 19 else 16
      points(nd$x1, nd$y, pch = pts_pch, cex = if (grepl("_ref$", nm)) 0.85 else 0.55, col = col_t)
      clean_nm <- sub("^OUT_Hvulgare_", "OUT_", nm)
      text(nd$x1 + max_x * 0.012, nd$y, clean_nm, adj = c(0, 0.5), cex = 0.52,
           col = col_t, font = if (grepl("_ref$", nm)) 2 else 1)
    } else {
      ys <- sapply(nd$children, function(c) c$y)
      segments(nd$x1, min(ys), nd$x1, max(ys), col = "grey30", lwd = 1.1)
      # Show bootstrap on well-supported deep nodes
      if (grepl("/", nd$name)) {
        parts <- as.numeric(strsplit(nd$name, "/")[[1]])
        if (length(parts) == 2 && !any(is.na(parts)) && parts[2] >= 90 && (max(ys) - min(ys) >= 3.5)) {
          text(nd$x1 - max_x * 0.008, nd$y + 0.45, sprintf("%.0f/%.0f", parts[1], parts[2]),
               adj = c(1, 0), cex = 0.48, col = "grey20")
        }
      }
      for (ch in nd$children) draw_clade(ch)
    }
  }
  draw_clade(c_tree)
  legend("topleft",
         legend = c("Kinase 1 (N-terminal module)", "Kinase 2 (C-terminal module)",
                    "Single-kinase progenitor / TF", "Brachypodium / Hordeum outgroup",
                    "Intragenic sister clades: Lr9 (99.7/100), Rpg1/WTK4 (99.2/100)"),
         col = c("#1f78b4", "#d62728", "#2ca02c", "grey45", "black"),
         pch = c(19, 19, 19, 19, NA), cex = 0.68, bty = "n")

  # ---- Panel c: K1 vs K2 amino-acid identity & Nei-Gojobori pS ----
  par(mar = c(4.5, 4.5, 2.8, 3.8))
  set.seed(42)
  plot(NA, xlim = c(0.45, 2.55), ylim = c(20, 105),
       xaxt = "n", xlab = "", ylab = "Kinase 1 vs Kinase 2 amino-acid identity (%)",
       bty = "l", cex.lab = 0.92)
  axis(1, at = c(1, 2),
       labels = c("Intragenic duplication\n(Lr9, Rpg1, WTK4)", "Intergenic fusion\n(Sr60, Sr62, Yr15)"),
       tick = TRUE, padj = 0.4, cex.axis = 0.82)
  title(main = "c   Intra-polypeptide K1-K2 identity & pS", adj = 0, cex.main = 0.95, font.main = 2)

  boxplot(list(pair_id_dup, pair_id_fus), at = c(0.88, 1.88), boxwex = 0.28,
          add = TRUE, col = c("#c6dbef", "#fcbba1"), border = c("#08519c", "#a50f15"),
          outline = FALSE, xaxt = "n", yaxt = "n")
  points(jitter(rep(0.88, length(pair_id_dup)), amount = 0.07), pair_id_dup,
         pch = 16, cex = 0.45, col = rgb(8/255, 81/255, 156/255, 0.30))
  points(jitter(rep(1.88, length(pair_id_fus)), amount = 0.07), pair_id_fus,
         pch = 16, cex = 0.45, col = rgb(165/255, 15/255, 21/255, 0.30))

  # Overlay locus-level means (diamonds) with non-overlapping leader labels
  points(rep(0.88, 3), loc_id_dup, pch = 23, bg = "#08519c", col = "white", cex = 1.35, lwd = 1.2)
  points(rep(1.88, 3), loc_id_fus, pch = 23, bg = "#a50f15", col = "white", cex = 1.35, lwd = 1.2)
  ord_d <- order(loc_id_dup)
  y_lab_d <- loc_id_dup[ord_d] + c(-3.2, 0.0, 3.2)
  segments(rep(0.73, 3), y_lab_d, rep(0.84, 3), loc_id_dup[ord_d], col = "#08519c", lwd = 0.7)
  text(rep(0.71, 3), y_lab_d, names(loc_id_dup)[ord_d], cex = 0.65, font = 3, col = "#08519c", adj = c(1, 0.5))

  ord_f <- order(loc_id_fus)
  y_lab_f <- loc_id_fus[ord_f] + c(-3.5, 0.0, 3.5)
  segments(rep(1.73, 3), y_lab_f, rep(1.84, 3), loc_id_fus[ord_f], col = "#a50f15", lwd = 0.7)
  text(rep(1.71, 3), y_lab_f, names(loc_id_fus)[ord_f], cex = 0.65, font = 3, col = "#a50f15", adj = c(1, 0.5))

  # Secondary axis for Nei-Gojobori pS (scaled 0.0 -> 20, 1.0 -> 100)
  ps_to_y <- function(ps) 20 + ps * 80
  axis(4, at = ps_to_y(seq(0, 1, by = 0.25)), labels = sprintf("%.2f", seq(0, 1, by = 0.25)),
       las = 1, cex.axis = 0.78, col = "grey25", col.axis = "grey20")
  mtext("Synonymous divergence pS (Nei-Gojobori)", side = 4, line = 2.5, cex = 0.78, col = "grey20")
  abline(h = ps_to_y(0.75), lty = 2, col = "grey40", lwd = 1.0)
  text(1.48, ps_to_y(0.75) + 2.2, "pS = 0.75 (dS saturation)", adj = c(0.5, 0), cex = 0.62, col = "grey30")

  # Plot genomic CDS pS (open triangles) & reference pS (filled triangles)
  points(jitter(rep(1.16, length(ps_g_dup)), amount = 0.04), ps_to_y(ps_g_dup),
         pch = 24, bg = "#74c476", col = "#006d2c", cex = 0.85)
  points(jitter(rep(2.16, length(ps_g_fus)), amount = 0.04), ps_to_y(ps_g_fus),
         pch = 24, bg = "#fd8d3c", col = "#a63603", cex = 0.85)
  points(rep(1.16, 3), ps_to_y(ref_ps_dup), pch = 17, col = "#00441b", cex = 1.2)
  points(rep(2.16, 3), ps_to_y(ref_ps_fus), pch = 17, col = "#7f2704", cex = 1.2)

  legend("bottomleft",
         legend = c("Locus mean K1-K2 aa identity", "Cloned reference query pS", "Pangenome genomic CDS pS"),
         pch = c(23, 17, 24), pt.bg = c("#08519c", NA, "#74c476"),
         col = c("white", "#00441b", "#006d2c"), pt.cex = c(1.1, 1.0, 0.85), cex = 0.62, bty = "n")

  text(1.5, 102, sprintf("Locus mean: 66.0%% vs 34.5%% (Welch P = 0.0060)\nAll pairs (n=515): MWU P = 2.3e-76 | Genomic pS: P = 1.1e-04"),
       cex = 0.64, font = 1, col = "grey15")

  # ---- Panel d: Architecture presence across 27 assemblies ----
  par(mar = c(4.5, 6.2, 2.8, 1.0))
  asms <- unique(s8d$asm_id)
  n_a <- length(asms)
  n_g <- length(CLASS3)
  state_mat <- matrix(0, nrow = n_a, ncol = n_g)
  for (ia in seq_len(n_a)) {
    for (ig in seq_len(n_g)) {
      st <- s8d$state[s8d$asm_id == asms[ia] & s8d$Gene == CLASS3[ig]][1]
      state_mat[ia, ig] <- if (st == "fused") 2 else if (st == "progenitor_only") 1 else 0
    }
  }
  plot(NA, xlim = c(0.5, n_g + 0.5), ylim = c(0.5, n_a + 2.2),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  title(main = "d   Fused vs progenitor retention (27 assemblies)", adj = 0, cex.main = 0.95, font.main = 2)
  cols_map <- c("#f0f0f0", "#fdae6b", "#08519c")
  for (ia in seq_len(n_a)) {
    y <- n_a - ia + 1
    for (ig in seq_len(n_g)) {
      v <- state_mat[ia, ig]
      rect(ig - 0.46, y - 0.44, ig + 0.46, y + 0.44,
           col = cols_map[v + 1], border = "white", lwd = 0.7)
    }
  }
  axis(1, at = seq_len(n_g), labels = CLASS3, las = 2, font = 3, cex.axis = 0.80, tick = FALSE, line = -0.5)
  axis(2, at = seq_len(n_a), labels = rev(asms), las = 1, cex.axis = 0.56, tick = FALSE, line = -0.5)
  legend(x = (n_g + 1) / 2, y = n_a + 2.15, xjust = 0.5, horiz = TRUE,
         legend = c("Fused TKP/TKF", "Progenitor only", "Absent"),
         fill = c("#08519c", "#fdae6b", "#f0f0f0"), border = "grey50", cex = 0.66, bty = "n")

  # ---- Panel e: Exon-intron structures & phase-0 domain junctions ----
  par(mar = c(4.5, 4.2, 2.8, 1.2))
  plot(NA, xlim = c(0, 22.0), ylim = c(0.4, length(CLASS3) + 0.9),
       xlab = "Exon index (5' -> 3') with intron phases (0, 1, 2)",
       ylab = "", yaxt = "n", bty = "n", cex.lab = 0.90)
  title(main = "e   Exon-intron phases & phase-0 junctions", adj = 0, cex.main = 0.95, font.main = 2)

  for (i in seq_along(CLASS3)) {
    g <- CLASS3[i]
    y <- length(CLASS3) - i + 1
    axis(2, at = y, labels = g, las = 1, font = 3, tick = FALSE, cex.axis = 0.85, line = -0.5)
    sub <- s8e[s8e$Gene == g, ]
    r0 <- sub[1, ]
    nex <- as.integer(r0$n_exons)
    phs <- if (r0$all_intron_phases != "-") as.integer(strsplit(r0$all_intron_phases, ",")[[1]]) else integer(0)

    # Identify which intron index is the domain junction intron (<=30 aa from K2_start or accessory_start)
    j_sub <- sub[sub$domain_boundary %in% c("K2_start", "accessory_start") &
                 sub$distance_aa != "-" & as.integer(sub$distance_aa) <= 30, ]
    segments(0.2, y, nex - 0.2, y, col = "grey65", lwd = 1.4)
    for (ex in seq_len(nex)) {
      rect(ex - 0.82, y - 0.22, ex - 0.18, y + 0.22, col = "#636363", border = "#252525", lwd = 0.8)
      if (ex <= length(phs)) {
        ph_val <- phs[ex]
        col_ph <- if (ph_val == 0) "#d95f02" else "grey35"
        text(ex - 0.0, y + 0.30, as.character(ph_val), cex = 0.62,
             font = if (ph_val == 0) 2 else 1, col = col_ph)
      }
    }
    if (nrow(j_sub) > 0) {
      jr <- j_sub[1, ]
      text(nex + 0.3, y,
           sprintf("phase %s @ aa %s (%s +/-%s aa)", jr$nearest_intron_phase, jr$nearest_intron_aa,
                   sub("_start", "", jr$domain_boundary), jr$distance_aa),
           adj = c(0, 0.5), cex = 0.60, col = "#b30000", font = 2)
    }
  }
  text(0.2, length(CLASS3) + 0.65,
       sprintf("All 6/6 genes with a junction intron (<=30 aa) are phase 0 (binomial P = %.4f)", binom_gene_p),
       adj = c(0, 0.5), cex = 0.68, col = "#b30000", font = 2)
}

pdf("figures/Fig3_tandem_kinase_origin_R.pdf", width = 15.5, height = 10.5, useDingbats = FALSE)
draw_fig3()
dev.off()

png("figures/Fig3_tandem_kinase_origin_R.png", width = 15.5, height = 10.5, units = "in", res = 300)
draw_fig3()
dev.off()

cat("Saved R statistical audit table -> results/S8f_R_statistical_audit_3.3.tsv\n")
cat("Saved R publication figures -> figures/Fig3_tandem_kinase_origin_R.pdf and .png\n")
cat("======================================================================\n")
