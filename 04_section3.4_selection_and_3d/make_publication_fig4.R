#!/usr/bin/env Rscript
# ==============================================================================
# make_publication_fig4.R
# Publication-Grade Figure 4 Generator embedding PyMOL 3D AlphaFold Renders
# and Statistical Selection Metrics using ggplot2, png, grid, and patchwork
#
# Generates:
#   figures/Fig4_selection_structures_publication.pdf
#   figures/Fig4_selection_structures_publication.png (300 dpi)
#   (also updates Fig4_selection_structures.png and Fig4_selection_structures_R.png)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(png)
  library(grid)
  library(patchwork)
  library(cowplot)
})

cat("=======================================================================\n")
cat("   Building Publication-Grade Figure 4 (R + PyMOL 3D AlphaFold)\n")
cat("=======================================================================\n\n")

.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", .args[grep("^--file=", .args)])
base_dir <- if (length(.file_arg) > 0) {
  dirname(dirname(normalizePath(.file_arg[1])))
} else if (file.exists("config_3.1.json")) {
  getwd()
} else if (file.exists("../config_3.1.json")) {
  normalizePath("..")
} else {
  getwd()
}
res_dir  <- file.path(base_dir, "results")
fig_dir  <- file.path(base_dir, "figures")
struct_dir <- file.path(fig_dir, "structures")

# ------------------------------------------------------------------------------
# Panel a: Domain-Partitioned Codon Selection (omega)
# ------------------------------------------------------------------------------
cat("Processing Panel a: Domain-partitioned omega distributions...\n")
s10c_df <- read.delim(file.path(res_dir, "Supplementary_Table_S10c_domain_partitioned_omega.tsv"), stringsAsFactors = FALSE)
sec34_vals <- read.delim(file.path(res_dir, "section3.4_values.tsv"), stringsAsFactors = FALSE)
get_v <- function(k) { col <- setdiff(names(sec34_vals), c("value", "Value", "note", "Notes"))[1]; sec34_vals$value[sec34_vals[[col]] == k][1] }

s10c_clean <- s10c_df %>%
  filter(omega_M0 < 10 & domain_name != "whole_gene") %>%
  mutate(
    Domain_Label = case_when(
      domain_name == "N-terminal_CC_or_ID" ~ "CC / ID\n(Signaling)\nn = 20",
      domain_name == "NB-ARC" ~ "NB-ARC\n(Nucleotide)\nn = 20",
      domain_name == "LRR" ~ "LRR\n(Recognition)\nn = 20",
      structural_class == "III" & domain_type == "kinase_module" ~ "Class III\n(Kinase)\nn = 12",
      structural_class == "IV" ~ "Class IV\n(Ank/TM)\nn = 4",
      structural_class == "V" ~ "Class V\n(Transp.)\nn = 3",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Domain_Label)) %>%
  mutate(
    Domain_Label = factor(Domain_Label, levels = c(
      "CC / ID\n(Signaling)\nn = 20", "NB-ARC\n(Nucleotide)\nn = 20",
      "LRR\n(Recognition)\nn = 20", "Class III\n(Kinase)\nn = 12",
      "Class IV\n(Ank/TM)\nn = 4", "Class V\n(Transp.)\nn = 3"
    ))
  )

kw_p_num <- as.numeric(get_v("kw_pval"))
kw_exp <- floor(log10(kw_p_num))
kw_mant <- kw_p_num / (10^kw_exp)

lrr_vals <- s10c_clean$omega_M0[grepl("^LRR", s10c_clean$Domain_Label)]
nb_vals  <- s10c_clean$omega_M0[grepl("^NB-ARC", s10c_clean$Domain_Label)]
wt_lrr_nb <- wilcox.test(lrr_vals, nb_vals, exact = FALSE, correct = TRUE)

p_a <- ggplot(s10c_clean, aes(x = Domain_Label, y = omega_M0, fill = Domain_Label)) +
  geom_violin(alpha = 0.40, color = NA, width = 0.82, scale = "width") +
  geom_boxplot(outlier.shape = NA, alpha = 0.78, width = 0.34, color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.7, shape = 21, color = "black", stroke = 0.4, alpha = 0.75) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "#b2182b", linewidth = 0.8) +
  annotate("segment", x = 2, xend = 3, y = 2.34, yend = 2.34, linewidth = 0.55, color = "black") +
  annotate("segment", x = 2, xend = 2, y = 2.28, yend = 2.34, linewidth = 0.55, color = "black") +
  annotate("segment", x = 3, xend = 3, y = 2.28, yend = 2.34, linewidth = 0.55, color = "black") +
  annotate("text", x = 2.5, y = 2.49,
           label = sprintf("LRR vs NB-ARC:\nU = %.1f, P = %.2e",
                           as.numeric(wt_lrr_nb$statistic), wt_lrr_nb$p.value),
           size = 2.45, fontface = "bold", color = "#b2182b", lineheight = 0.9) +
  scale_fill_manual(values = c(
    "CC / ID\n(Signaling)\nn = 20"   = "#4a90d9",
    "NB-ARC\n(Nucleotide)\nn = 20"   = "#9ecae1",
    "LRR\n(Recognition)\nn = 20"     = "#d73027",
    "Class III\n(Kinase)\nn = 12"    = "#d95f02",
    "Class IV\n(Ank/TM)\nn = 4"      = "#7570b3",
    "Class V\n(Transp.)\nn = 3"      = "#1b7837"
  )) +
  scale_y_continuous(limits = c(0, 2.90), breaks = seq(0, 2.5, 0.5)) +
  annotate("text", x = 3.5, y = 2.75,
           label = sprintf("list('5-class KW (39 dom. across 31 OGs):'~italic(H) == %s, italic(P) == %.2f %%*%% 10^%d)",
                           get_v("kw_H_stat"), kw_mant, kw_exp),
           parse = TRUE, size = 2.55, fontface = "bold", color = "black") +
  labs(
    title = expression(bold("a  Domain-partitioned ") * bolditalic(omega) * bold(" (n = 79 dom. across 31 OGs)")),
    x = NULL,
    y = expression(bold("Codeml M0 rate ratio (") * bolditalic(omega) * bold(")"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    legend.position = "none",
    axis.text.x = element_text(size = 7.4, face = "bold", color = "black"),
    axis.text.y = element_text(size = 9, color = "black")
  )

# ------------------------------------------------------------------------------
# Panel b: Paired Tandem Kinase Module Asymmetry (Dumbbell Asymmetry Plot)
# ------------------------------------------------------------------------------
cat("Processing Panel b: Paired tandem kinase asymmetry...\n")
s10d_df <- read.delim(file.path(res_dir, "Supplementary_Table_S10d_tkp_module_asymmetry.tsv"), stringsAsFactors = FALSE)

dumbbell_df <- s10d_df %>%
  arrange(delta_omega_m2_minus_m1) %>%
  mutate(
    gene_ord = factor(reference_query, levels = reference_query),
    lbl_text = sprintf("Delta*omega == '%+.2f (%s)'", delta_omega_m2_minus_m1,
                       ifelse(site_MWU_P < 0.001, "P < 0.001 ***",
                              ifelse(site_MWU_P < 0.01, "P < 0.01 **",
                                     ifelse(site_MWU_P < 0.05, "P < 0.05 *",
                                            sprintf("P = %.2f ns", site_MWU_P))))),
    sig_col = ifelse(site_MWU_P < 0.05, "#b30000", "#525252"),
    m1_short = case_when(
      module_1 == "K1_active" ~ "K1 act",
      module_1 == "K2_active" ~ "K2 act",
      module_1 == "K1_pseudokinase" ~ "pK1",
      TRUE ~ module_1
    ),
    m2_short = case_when(
      module_2 == "K2_pseudokinase" ~ "pK2",
      module_2 == "K1_pseudokinase" ~ "pK1",
      module_2 == "Accessory_DUF3475" ~ "DUF3475",
      TRUE ~ module_2
    ),
    m1_x = ifelse(abs(delta_omega_m2_minus_m1) < 0.18,
                  module_1_omega_M0 + ifelse(delta_omega_m2_minus_m1 >= 0, -0.02, 0.02),
                  module_1_omega_M0),
    m2_x = ifelse(abs(delta_omega_m2_minus_m1) < 0.18,
                  module_2_omega_M0 + ifelse(delta_omega_m2_minus_m1 >= 0, 0.02, -0.02),
                  module_2_omega_M0),
    m1_hjust = ifelse(abs(delta_omega_m2_minus_m1) < 0.18,
                      ifelse(delta_omega_m2_minus_m1 >= 0, 1, 0), 0.5),
    m2_hjust = ifelse(abs(delta_omega_m2_minus_m1) < 0.18,
                      ifelse(delta_omega_m2_minus_m1 >= 0, 0, 1), 0.5)
  )

wilc_tkp <- wilcox.test(s10d_df$module_1_omega_M0, s10d_df$module_2_omega_M0, paired = TRUE)

p_b <- ggplot(dumbbell_df) +
  # Connecting arrow
  geom_segment(aes(x = module_1_omega_M0, xend = module_2_omega_M0, y = gene_ord, yend = gene_ord),
               arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
               color = "#525252", linewidth = 1.15) +
  # Points
  geom_point(aes(x = module_1_omega_M0, y = gene_ord, fill = "Module 1 (Catalytic/Primary)"),
             size = 3.8, shape = 21, color = "black", stroke = 0.5) +
  geom_point(aes(x = module_2_omega_M0, y = gene_ord, fill = "Module 2 (Decoy/Accessory)"),
             size = 3.8, shape = 21, color = "black", stroke = 0.5) +
  # Domain labels above points
  geom_text(aes(x = m1_x, y = as.numeric(gene_ord) + 0.24, label = m1_short, hjust = m1_hjust),
            size = 2.4, fontface = "bold", color = "#1b7837") +
  geom_text(aes(x = m2_x, y = as.numeric(gene_ord) + 0.24, label = m2_short, hjust = m2_hjust),
            size = 2.4, fontface = "bold", color = "#d95f02") +
  # Statistics on the right
  geom_text(aes(x = pmax(module_1_omega_M0, module_2_omega_M0) + 0.04, y = gene_ord, label = lbl_text, color = sig_col),
            parse = TRUE, hjust = 0, size = 2.40, fontface = "bold") +
  scale_color_identity(guide = "none") +
  scale_fill_manual(
    values = c("Module 1 (Catalytic/Primary)" = "#1b7837", "Module 2 (Decoy/Accessory)" = "#d95f02"),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(0.02, 1.56),
    breaks = seq(0.2, 1.0, 0.2),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_discrete(expand = expansion(mult = c(0.14, 0.20))) +
  labs(
    title = expression(bold("b  Paired Class III kinase module asymmetry (") * bolditalic("n") * bold(" = 5 loci)")),
    subtitle = "M1 (green) -> M2 (orange) | FUBAR MWU (* <0.05, *** <0.001, ns)",
    x = expression(bold("Domain evolutionary rate ratio (") * bolditalic(omega) * bold(")")),
    y = "Cloned Class III locus"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10.2),
    plot.subtitle = element_text(size = 7.8, color = "#444444"),
    legend.position = "top",
    legend.text = element_text(size = 8.0, face = "bold"),
    legend.margin = margin(t = 1, b = 2),
    axis.title.x = element_text(size = 8.8, face = "bold", color = "black", margin = margin(t = 4)),
    axis.title.y = element_text(size = 9, face = "bold", color = "black", margin = margin(r = 5)),
    axis.text.y = element_text(size = 8.8, face = "bold.italic", color = "black"),
    axis.text.x = element_text(size = 8, color = "black")
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3.5)))

# ------------------------------------------------------------------------------
# Panel c: LRR Structural Domain Distribution & Concave Face Enrichment
# ------------------------------------------------------------------------------
cat("Processing Panel c: Structural domain distribution & concave enrichment...\n")
s10e_df <- read.delim(file.path(res_dir, "Supplementary_Table_S10e_af_surface_mapping.tsv"), stringsAsFactors = FALSE)
c1 <- s10e_df %>%
  filter(structural_class == "I") %>%
  group_by(orthogroup, uniprot_id) %>%
  summarise(gene = paste(gene, collapse = "/"), across(-gene, first), .groups = "drop") %>%
  arrange(concave_OR_vs_whole_protein)

plot_df <- list()
for (i in 1:nrow(c1)) {
  tot   <- c1$n_selected_sites_BEB95_or_MEME05[i]
  conc  <- c1$selected_on_concave_face[i]
  other <- c1$selected_in_LRR[i] - conc
  non   <- tot - c1$selected_in_LRR[i]
  g     <- c1$gene[i]
  
  plot_df[[length(plot_df) + 1]] <- data.frame(
    gene = g, category = "Concave LRR beta-strand", count = conc, pct = 100 * conc / tot,
    stringsAsFactors = FALSE
  )
  plot_df[[length(plot_df) + 1]] <- data.frame(
    gene = g, category = "Other LRR (loops/convex)", count = other, pct = 100 * other / tot,
    stringsAsFactors = FALSE
  )
  plot_df[[length(plot_df) + 1]] <- data.frame(
    gene = g, category = "Scaffold (CC / NB-ARC)", count = non, pct = 100 * non / tot,
    stringsAsFactors = FALSE
  )
}
df_long <- do.call(rbind, plot_df)
df_long$gene <- factor(df_long$gene, levels = c1$gene)
df_long$category <- factor(df_long$category, levels = c("Scaffold (CC / NB-ARC)", "Other LRR (loops/convex)", "Concave LRR beta-strand"))

bg_pct <- as.numeric(get_v("classI_concave_bg_pct"))
fisher_p_val <- as.numeric(get_v("classI_concave_fisher_pval"))

anno_df <- c1 %>% mutate(
  gene = factor(gene, levels = gene),
  p_sig = ifelse(concave_P_vs_whole_protein < 0.001, " ***",
                 ifelse(concave_P_vs_whole_protein < 0.01, " **",
                        ifelse(concave_P_vs_whole_protein < 0.05, " *", " ns"))),
  p_fmt = ifelse(concave_P_vs_whole_protein < 0.001,
                 sprintf("%.2e", concave_P_vs_whole_protein),
                 sprintf("%.3f", concave_P_vs_whole_protein)),
  lbl_txt = sprintf("%d/%d | OR=%.1f (P=%s)%s",
                    selected_on_concave_face, n_selected_sites_BEB95_or_MEME05,
                    concave_OR_vs_whole_protein, p_fmt, p_sig),
  Sig = ifelse(concave_P_vs_whole_protein < 0.05, "Significant", "Not significant")
)

p_c <- ggplot(df_long, aes(x = pct, y = gene, fill = category)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.28) +
  geom_vline(xintercept = bg_pct, linetype = "dotted", color = "#b2182b", linewidth = 0.75) +
  geom_text(data = anno_df, aes(x = 102, y = gene, label = lbl_txt, color = Sig),
            inherit.aes = FALSE, hjust = 0, size = 2.35, fontface = "bold") +
  scale_fill_manual(
    values = c(
      "Concave LRR beta-strand"  = "#d73027",
      "Other LRR (loops/convex)" = "#fdae61",
      "Scaffold (CC / NB-ARC)"   = "#d9d9d9"
    ),
    labels = c(
      "Scaffold (CC / NB-ARC)"   = "Scaffold (CC/NB-ARC)",
      "Other LRR (loops/convex)" = "Other LRR loops",
      "Concave LRR beta-strand"  = "Concave LRR beta"
    ),
    name = NULL
  ) +
  scale_color_manual(
    values = c("Significant" = "#b2182b", "Not significant" = "#525252"),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 178),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  annotate("text", x = bg_pct + 1.2, y = 1.2,
           label = sprintf("%.1f%% proteome bg", bg_pct),
           hjust = 0, size = 2.15, fontface = "italic", color = "#b2182b") +
  labs(
    title = sprintf("c  3D partitioning of positive selection (OR = %s, P = %.2e)",
                    get_v("classI_concave_odds_ratio"), fisher_p_val),
    subtitle = "Concave beta (red) vs LRR loops (orange) vs Scaffold (gray)",
    x = "% of positive selection sites in receptor",
    y = "Cloned Class I NLR"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 9.3),
    plot.subtitle = element_text(size = 7.8, color = "#444444"),
    legend.position = "top",
    legend.text = element_text(size = 8.0),
    legend.title = element_text(size = 8.2, face = "bold"),
    legend.margin = margin(t = 1, b = 2),
    axis.title.x = element_text(size = 8.8, face = "bold", color = "black", margin = margin(t = 4)),
    axis.title.y = element_text(size = 9, face = "bold", color = "black", margin = margin(r = 5)),
    axis.text.y = element_text(size = 8.5, face = "bold.italic", color = "black"),
    axis.text.x = element_text(size = 8, color = "black")
  ) +
  guides(fill = guide_legend(nrow = 1))

# ------------------------------------------------------------------------------
# Panels d, e, f: Embedded Ray-Traced 3D AlphaFold Models (PyMOL)
# ------------------------------------------------------------------------------
cat("Loading and embedding PyMOL 3D AlphaFold renders...\n")

make_struct_panel <- function(png_path, title_text, subtitle_text) {
  if (file.exists(png_path)) {
    img <- readPNG(png_path)
    # Trim transparent alpha margins so protein sits tight to title
    if (dim(img)[3] == 4) {
      alpha <- img[, , 4]
      rows <- which(rowSums(alpha > 0.05) > 0)
      cols <- which(colSums(alpha > 0.05) > 0)
      if (length(rows) > 0 && length(cols) > 0) {
        r1 <- max(1, min(rows) - 8)
        r2 <- min(nrow(img), max(rows) + 8)
        c1 <- max(1, min(cols) - 8)
        c2 <- min(ncol(img), max(cols) + 8)
        img <- img[r1:r2, c1:c2, ]
      }
    }
    g <- rasterGrob(img, interpolate = TRUE, width = unit(0.96, "npc"), height = unit(0.96, "npc"))
    ggplot() +
      annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      theme_void() +
      labs(title = title_text, subtitle = subtitle_text) +
      theme(
        plot.title = element_text(face = "bold", size = 10.5, hjust = 0, margin = margin(b = 2)),
        plot.subtitle = element_text(size = 8.2, hjust = 0, color = "#525252", margin = margin(b = 4)),
        plot.margin = margin(t = 4, r = 6, b = 4, l = 6)
      )
  } else {
    ggplot() + theme_void() + labs(title = paste0(title_text, " (Structure pending)"))
  }
}

sr35_n <- as.integer(s10e_df$n_selected_sites_BEB95_or_MEME05[s10e_df$gene == "Sr35"][1])
yr27_n <- as.integer(s10e_df$n_selected_sites_BEB95_or_MEME05[s10e_df$gene == "Yr27"][1])

p_d <- make_struct_panel(
  file.path(struct_dir, "Sr35_selection_surface.png"),
  "d  Sr35 AlphaFold model with exposed LRR beta-strand sites",
  sprintf("%d sites under diversifying selection (red spheres) on solvent-exposed LRR beta-strands", sr35_n)
)

p_e <- make_struct_panel(
  file.path(struct_dir, "Yr15_wtk_asymmetry.png"),
  "e  Yr15 AlphaFold model showing WTK module asymmetry",
  "Active Kinase 1 (green, catalytic triad) vs Pseudokinase 2 decoy (orange, degenerated)"
)

p_f <- make_struct_panel(
  file.path(struct_dir, "Yr27_selection_surface.png"),
  "f  Yr27 AlphaFold model showing dense LRR recognition cluster",
  sprintf("%d positive selection sites (red spheres) mapping to solvent-exposed surface", yr27_n)
)

# ------------------------------------------------------------------------------
# Compose Full Figure 4 with Patchwork
# ------------------------------------------------------------------------------
cat("Composing full Figure 4 with patchwork...\n")

# Top row: Statistical metrics (a: domain omega, b: kinase asymmetry, c: surface OR)
# Bottom row: 3D AlphaFold ray-traced structures (d: Sr35, e: Yr15, f: Yr27)
top_row <- (p_a | p_b | p_c) + plot_layout(widths = c(1.3, 1.0, 1.2))
bot_row <- (p_d | p_e | p_f) + plot_layout(widths = c(1.0, 1.0, 1.0))

fig4_composed <- top_row / bot_row + plot_layout(heights = c(1.0, 1.05))

out_pdf <- file.path(fig_dir, "Fig4_selection_structures_publication.pdf")
out_png <- file.path(fig_dir, "Fig4_selection_structures_publication.png")

cat(sprintf("Saving PDF -> %s ...\n", out_pdf))
pdf(out_pdf, width = 16, height = 13.5)
print(fig4_composed)
invisible(dev.off())

cat(sprintf("Saving 300-DPI PNG -> %s ...\n", out_png))
png(out_png, width = 16, height = 13.5, units = "in", res = 300)
print(fig4_composed)
invisible(dev.off())

file.copy(out_png, file.path(fig_dir, "Fig4_selection_structures.png"), overwrite = TRUE)
file.copy(out_png, file.path(fig_dir, "Fig4_selection_structures_R.png"), overwrite = TRUE)
file.copy(out_pdf, file.path(fig_dir, "Fig4_selection_structures.pdf"), overwrite = TRUE)
file.copy(out_pdf, file.path(fig_dir, "Fig4_selection_structures_R.pdf"), overwrite = TRUE)

cat("\n>>> Figure 4 successfully generated and verified! <<<\n")
cat("=======================================================================\n")
