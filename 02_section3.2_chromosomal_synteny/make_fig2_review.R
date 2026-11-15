#!/usr/bin/env Rscript
# ==============================================================================
# make_fig2_review.R
# Publication-Grade Figure 2 Generator using ggplot2, scales, and patchwork
#
# Generates:
#   figures/Fig2_chromosomal_distribution_review.pdf
#   figures/Fig2_chromosomal_distribution_review.png (300 dpi)
#   (also updates Fig2_chromosomal_distribution.png and Fig2_chromosomal_distribution_R.png)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(patchwork)
  library(grid)
  library(cowplot)
  library(gggenes)
})

cat("=======================================================================\n")
cat("   Building Publication-Grade Figure 2 (R + ggplot2 + patchwork)\n")
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
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

class_cols <- c(
  "Class I"   = "#1f77b4",
  "Class II"  = "#4a90d9",
  "Class III" = "#d95f02",
  "Class IV"  = "#7570b3",
  "Class V"   = "#1b9e77"
)

# ------------------------------------------------------------------------------
# Panel a: Centromere-to-telomere chromosomal density + LocusGroup rug
# ------------------------------------------------------------------------------
cat("Processing Panel a: Centromere-to-telomere density...\n")
s4_df <- read.delim(file.path(res_dir, "S4_chromosomal_positions.tsv"), stringsAsFactors = FALSE)
s4_clean <- s4_df %>% filter(!is.na(d_norm))

s4b_df <- read.delim(file.path(res_dir, "S4b_locusgroup_positions.tsv"), stringsAsFactors = FALSE) %>%
  filter(!is.na(median_d_norm))

p_a <- ggplot(s4_clean, aes(x = d_norm, fill = Group, color = Group)) +
  geom_density(alpha = 0.42, linewidth = 0.95) +
  geom_rug(data = s4b_df, aes(x = median_d_norm, color = Group),
           inherit.aes = FALSE, sides = "b", length = unit(0.035, "npc"),
           linewidth = 0.7, alpha = 0.85) +
  geom_vline(xintercept = 0.85, linetype = "dashed", color = "#b2182b", linewidth = 0.95) +
  scale_fill_manual(values = c("NLR" = "#1f77b4", "nonNLR" = "#d95f02"),
                    labels = c("NLR (Classes I-II, n = 9,776)", "non-NLR (Classes III-V, n = 2,242)")) +
  scale_color_manual(values = c("NLR" = "#1f77b4", "nonNLR" = "#d95f02"),
                     labels = c("NLR (Classes I-II, n = 9,776)", "non-NLR (Classes III-V, n = 2,242)")) +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.2), expand = expansion(mult = c(0.01, 0.01))) +
  annotate("text", x = 0.83, y = 3.55,
           label = "Distal 15% boundary (d_norm = 0.85)\n77.6% vs 36.4% loci (KS D = 0.451, P = 3.95e-323)\nLocus medians (rug): 80.8% vs 40.0% distal",
           hjust = 1, size = 2.8, fontface = "italic", color = "#b2182b") +
  labs(
    title = "a  Chromosomal position from centromere (0.0) to telomere (1.0)",
    x = "Arm-normalized distance from centromere (0.0) to telomere (1.0)",
    y = "Normalized locus density",
    fill = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    legend.position = c(0.32, 0.86),
    legend.background = element_rect(fill = alpha("white", 0.88), color = NA),
    legend.text = element_text(size = 8.5, face = "bold"),
    axis.text = element_text(size = 9, color = "black")
  )

# ------------------------------------------------------------------------------
# Panel b: Violin + Boxplot across Classes I-V with LocusGroup median diamonds
# ------------------------------------------------------------------------------
cat("Processing Panel b: Violin + Boxplots across Classes I-V...\n")
s4_clean <- s4_clean %>%
  mutate(
    Class_Label = paste0("Class ", Class),
    Class_Label = factor(Class_Label, levels = c("Class I", "Class II", "Class III", "Class IV", "Class V"))
  )

s4b_clean <- s4b_df %>%
  mutate(
    Class_Label = paste0("Class ", Class),
    Class_Label = factor(Class_Label, levels = c("Class I", "Class II", "Class III", "Class IV", "Class V"))
  )

class_counts <- s4_clean %>%
  count(Class_Label) %>%
  mutate(label = paste0("n=", comma(n)))

p_b <- ggplot(s4_clean, aes(x = Class_Label, y = d_norm, fill = Class_Label)) +
  geom_violin(alpha = 0.45, color = NA, width = 0.85, scale = "width") +
  geom_boxplot(outlier.shape = NA, alpha = 0.82, width = 0.28, color = "black", linewidth = 0.45) +
  geom_point(data = s4b_clean, aes(x = Class_Label, y = median_d_norm),
             inherit.aes = FALSE, shape = 23, size = 2.1, fill = "white", color = "black", stroke = 0.6,
             position = position_jitter(width = 0.10, height = 0, seed = 42)) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "#b2182b", linewidth = 0.9) +
  geom_text(data = class_counts, aes(x = Class_Label, y = -0.05, label = label),
            inherit.aes = FALSE, size = 2.7, color = "#333333", fontface = "bold") +
  scale_fill_manual(values = class_cols) +
  scale_y_continuous(limits = c(-0.09, 1.05), breaks = seq(0, 1.0, 0.2)) +
  labs(
    title = "b  Relative chromosomal position d_norm across Classes I-V",
    subtitle = "Violins: all 12,018 placed loci; white diamonds: 36 LocusGroup medians",
    x = NULL,
    y = "Arm-normalized position d_norm"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, color = "#444444"),
    legend.position = "none",
    axis.text.x = element_text(size = 9.5, face = "bold", color = "black"),
    axis.text.y = element_text(size = 9, color = "black")
  )

# ------------------------------------------------------------------------------
# Panel c: 30-Orthogroup A:B:D Copy-Number Heatmap + 13-Cultivar Triad-State
#          Stacked Bar (Walkowiak 2020 / Jiao 2025 style) + Retention Summary
# ------------------------------------------------------------------------------
cat("Processing Panel c: A:B:D copy-number heatmap & 13-cultivar triad states...\n")
s5b_df <- read.delim(file.path(res_dir, "S5b_orthogroup_triad_summary.tsv"), stringsAsFactors = FALSE)
s5_df  <- read.delim(file.path(res_dir, "S5_hexaploid_triads.tsv"), stringsAsFactors = FALSE) %>%
  distinct(orthogroup, asm_id, .keep_all = TRUE)

# Clean short labels for the 30 Main orthogroups
og_order_df <- s5b_df %>%
  mutate(
    short_label = case_when(
      grepl("^Sr21,", queries) ~ "Sr21/Yr5/Yr7*",
      grepl("^Sr33,", queries) ~ "Sr33/Sr50*",
      grepl("^Sr13,", queries) ~ "Sr13/Sr27*",
      grepl("^Sr22,", queries) ~ "Sr22/Sr22b",
      TRUE ~ queries
    ),
    total_ABD = mean_copies_A + mean_copies_B + mean_copies_D,
    Group_ord = ifelse(Group == "nonNLR", 1, 2)
  ) %>%
  arrange(Group_ord, desc(n_strict_1_1_1), desc(n_all_three_ABD), desc(total_ABD))

# Bottom-to-top factor order so non-NLRs appear at the top of the heatmap
og_levels <- rev(og_order_df$orthogroup)
label_map <- setNames(og_order_df$short_label, og_order_df$orthogroup)

hm_df <- og_order_df %>%
  select(orthogroup, short_label, Group, A = mean_copies_A, B = mean_copies_B, D = mean_copies_D) %>%
  pivot_longer(cols = c("A", "B", "D"), names_to = "Subgenome", values_to = "MeanCopies") %>%
  mutate(
    orthogroup = factor(orthogroup, levels = og_levels),
    Subgenome  = factor(Subgenome, levels = c("A", "B", "D")),
    txt_lbl    = ifelse(MeanCopies == 0, "0", sprintf("%.1f", MeanCopies)),
    txt_col    = ifelse(MeanCopies >= 3.5, "white", "black")
  )

# Separator y-position between NLR (20) and non-NLR (10)
sep_y <- 20.5

p_c1a <- ggplot(hm_df, aes(x = Subgenome, y = orthogroup, fill = log1p(MeanCopies))) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = txt_lbl, color = txt_col), size = 2.15, fontface = "bold") +
  geom_hline(yintercept = sep_y, color = "black", linewidth = 0.75, linetype = "dashed") +
  scale_color_identity() +
  scale_fill_gradientn(
    colors = c("#f7f7f7", "#deebf7", "#6baed6", "#2171b5", "#08306b"),
    values = rescale(c(0, log1p(0.5), log1p(1.5), log1p(3.5), log1p(7.2))),
    name   = "Mean copies\n(log1p)"
  ) +
  scale_y_discrete(labels = label_map) +
  scale_x_discrete(position = "top") +
  labs(
    title = "c  Subgenome A:B:D copy number",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7.0, face = "bold", color = "black"),
    axis.text.x = element_text(size = 8.5, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.key.height = unit(2.8, "mm"),
    legend.key.width  = unit(5.5, "mm"),
    legend.title = element_text(size = 7.2, face = "bold"),
    legend.text  = element_text(size = 6.8)
  )

# 13-cultivar triad-state stacked bar aligned row-for-row with p_c1a
state_df <- s5_df %>%
  count(orthogroup, triad_state) %>%
  mutate(
    orthogroup  = factor(orthogroup, levels = og_levels),
    triad_state = factor(
      triad_state,
      levels = c("balanced_1:1:1", "present_A+B+D_multicopy", "subgenome_asymmetric", "absent"),
      labels = c("1:1:1 triad", "A+B+D multicopy", "Asymmetric PAV", "Absent")
    )
  )

p_c1b <- ggplot(state_df, aes(x = n / 13 * 100, y = orthogroup, fill = triad_state)) +
  geom_col(width = 0.82, color = "white", linewidth = 0.2) +
  geom_hline(yintercept = sep_y, color = "black", linewidth = 0.75, linetype = "dashed") +
  scale_fill_manual(values = c(
    "1:1:1 triad"     = "#1b9e77",
    "A+B+D multicopy" = "#66c2a5",
    "Asymmetric PAV"  = "#fc8d62",
    "Absent"          = "#d9d9d9"
  )) +
  scale_x_continuous(limits = c(0, 100.5), breaks = c(0, 50, 100), expand = c(0, 0)) +
  labs(
    title = "13-cultivar triad state (%)",
    x = "% of 13 hexaploids", y = NULL, fill = NULL
  ) +
  guides(fill = guide_legend(nrow = 2)) +
  theme_classic(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 9.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 7.8, color = "black"),
    axis.title.x = element_text(size = 8.0),
    legend.position = "bottom",
    legend.key.size = unit(3.0, "mm"),
    legend.text = element_text(size = 7.0)
  )

# Orthogroup-level triad summaries synced with S7c_R_statistical_audit_3.2.tsv
triad_df <- data.frame(
  Group = c("NLR (Classes I-II)", "non-NLR (Classes III-V)"),
  Conserved = c(0.0, 90.0),
  Strict_111 = c(3.0, 18.6)
) %>%
  pivot_longer(cols = c("Conserved", "Strict_111"), names_to = "Metric", values_to = "Pct") %>%
  mutate(Metric = factor(Metric, levels = c("Conserved", "Strict_111"),
                         labels = c("A:B:D conserved (>=8/13)", "Strict 1:1:1 cells (7/231 vs 11/59)")),
         Group = factor(Group, levels = c("NLR (Classes I-II)", "non-NLR (Classes III-V)")))

p_c2 <- ggplot(triad_df, aes(x = Group, y = Pct, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", Pct)), position = position_dodge(width = 0.72),
            vjust = -0.35, size = 2.65, fontface = "bold") +
  scale_fill_manual(values = c("A:B:D conserved (>=8/13)" = "#1b9e77",
                               "Strict 1:1:1 cells (7/231 vs 11/59)" = "#d95f02")) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.03))) +
  labs(title = "Triad retention (P <= 1.05e-04)", x = NULL, y = "Proportion (%)", fill = NULL) +
  guides(fill = guide_legend(nrow = 2)) +
  theme_classic(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 9.5),
    axis.text.x = element_text(size = 7.8, face = "bold", color = "black"),
    axis.text.y = element_text(size = 8.0, color = "black"),
    legend.position = "bottom",
    legend.key.size = unit(3.0, "mm"),
    legend.text = element_text(size = 7.0)
  )

p_c <- wrap_elements(full = (p_c1a + p_c1b + p_c2) + plot_layout(ncol = 3, widths = c(1.05, 0.95, 1.0)))

# ------------------------------------------------------------------------------
# Panel d: Syntenic Retention Index (SRI) across Classes I-V
# ------------------------------------------------------------------------------
cat("Processing Panel d: SRI across Classes I-V...\n")
s6a_df <- read.delim(file.path(res_dir, "S6a_sri_per_locus.tsv"), stringsAsFactors = FALSE)
s6a_clean <- s6a_df %>%
  mutate(
    Class_Label = paste0("Class ", Class),
    Class_Label = factor(Class_Label, levels = c("Class I", "Class II", "Class III", "Class IV", "Class V"))
  )

sri_counts <- s6a_clean %>%
  count(Class_Label) %>%
  mutate(label = paste0("n=", n))

sri_extremes <- s6a_clean %>%
  filter(LocusGroup %in% c("Lr21", "Lr34"))

p_d <- ggplot(s6a_clean, aes(x = Class_Label, y = mean_SRI, fill = Class_Label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.55, color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.16, size = 2.1, shape = 21, color = "black", stroke = 0.5) +
  geom_text(data = sri_counts, aes(x = Class_Label, y = -0.05, label = label),
            inherit.aes = FALSE, size = 2.8, fontface = "bold", color = "#333333") +
  geom_text(data = sri_extremes, aes(x = Class_Label, y = mean_SRI, label = LocusGroup),
            vjust = -1.0, hjust = -0.15, size = 2.7, fontface = "bold.italic", color = "black") +
  scale_fill_manual(values = class_cols) +
  scale_y_continuous(limits = c(-0.08, 1.08), breaks = seq(0, 1.0, 0.2)) +
  annotate("text", x = 3.0, y = 1.03,
           label = sprintf("NLR (mean = %.2f) vs non-NLR (mean = %.2f)",
                           mean(s6a_clean$mean_SRI[s6a_clean$Group == "NLR"], na.rm = TRUE),
                           mean(s6a_clean$mean_SRI[s6a_clean$Group != "NLR"], na.rm = TRUE)),
           size = 2.95, fontface = "bold", color = "#b30000") +
  annotate("text", x = 3.0, y = 0.94,
           label = "Mann-Whitney U = 33.0, P = 6.11e-04",
           size = 2.95, fontface = "bold", color = "#b30000") +
  labs(
    title = "d  Microsynteny retention across +/-500 kb windows (n = 36 loci)",
    x = NULL,
    y = "Synteny Retention Index (SRI)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "none",
    axis.text.x = element_text(size = 9.5, face = "bold", color = "black"),
    axis.text.y = element_text(size = 9, color = "black")
  )

# ------------------------------------------------------------------------------
# Panel e: Flanking TE vs SRI (null) paired with Local SV Events vs SRI (driver)
# ------------------------------------------------------------------------------
cat("Processing Panel e: Flanking TE vs SV mechanistic drivers...\n")
s7a_df <- read.delim(file.path(res_dir, "S7a_flanking_te_windows.tsv"), stringsAsFactors = FALSE) %>%
  mutate(Group_Label = ifelse(Group == "NLR", "NLR (n = 26)", "non-NLR (n = 10)"))

p_e1 <- ggplot(s7a_df, aes(x = total_te_pct, y = mean_SRI, color = Group_Label, shape = Group_Label)) +
  geom_smooth(method = "lm", se = TRUE, color = "#666666", fill = "#dddddd", linewidth = 0.7, linetype = "dashed", inherit.aes = FALSE,
              aes(x = total_te_pct, y = mean_SRI)) +
  geom_point(size = 2.3, alpha = 0.85) +
  scale_color_manual(values = c("NLR (n = 26)" = "#1f77b4", "non-NLR (n = 10)" = "#d95f02")) +
  scale_shape_manual(values = c("NLR (n = 26)" = 16, "non-NLR (n = 10)" = 17)) +
  annotate("text", x = 42, y = 0.96, hjust = 0, size = 2.7, fontface = "bold",
           label = "Spearman rho = 0.107, P = 0.533\nTE%: 70.9% vs 74.1% (P = 0.417)") +
  scale_y_continuous(limits = c(0, 1.04), breaks = seq(0, 1.0, 0.25)) +
  labs(
    title = "e  Flanking TE (%) vs SRI",
    x = "Flanking +/-200 kb TE (%)",
    y = "Synteny Retention Index (SRI)",
    color = NULL, shape = NULL
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    axis.text = element_text(size = 8.5, color = "black")
  )

p_e2 <- ggplot(s7a_df, aes(x = mean_sv_events, y = mean_SRI, color = Group_Label, shape = Group_Label)) +
  geom_smooth(method = "lm", se = TRUE, color = "#b2182b", fill = "#f4a582", alpha = 0.25, linewidth = 0.8, inherit.aes = FALSE,
              aes(x = mean_sv_events, y = mean_SRI)) +
  geom_point(size = 2.3, alpha = 0.85) +
  scale_color_manual(values = c("NLR (n = 26)" = "#1f77b4", "non-NLR (n = 10)" = "#d95f02")) +
  scale_shape_manual(values = c("NLR (n = 26)" = 16, "non-NLR (n = 10)" = 17)) +
  annotate("text", x = 31, y = 0.96, hjust = 1, size = 2.7, fontface = "bold", color = "#b2182b",
           label = "Spearman rho = -0.332, P = 0.048\nSV events: 17.6 vs 12.2") +
  scale_y_continuous(limits = c(0, 1.04), breaks = seq(0, 1.0, 0.25)) +
  labs(
    title = "Local SV rearrangements vs SRI",
    x = "Mean SV events / +/-500 kb",
    y = NULL,
    color = NULL, shape = NULL
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    axis.text = element_text(size = 8.5, color = "black")
  )

p_e <- (p_e1 + p_e2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# Panel f: All-36-Locus Synteny Retention (SRI +/- SEM) & SV Burden Ranking
# ------------------------------------------------------------------------------
cat("Processing Panel f: All-36-locus SRI ranking & SV burden...\n")
s6a_rank <- s6a_df %>%
  mutate(
    Group_Label = ifelse(Group == "NLR", "NLR (Classes I-II, n = 26)", "non-NLR (Classes III-V, n = 10)"),
    Locus_ord   = reorder(LocusGroup, mean_SRI),
    sri_low     = pmax(0, mean_SRI - sem_SRI),
    sri_high    = pmin(1, mean_SRI + sem_SRI)
  )

p_f <- ggplot(s6a_rank, aes(x = mean_SRI, y = Locus_ord, color = Group_Label)) +
  geom_vline(xintercept = mean(s6a_rank$mean_SRI[s6a_rank$Group == "NLR"], na.rm = TRUE),
             linetype = "dashed", color = "#1f77b4", alpha = 0.7, linewidth = 0.6) +
  geom_vline(xintercept = mean(s6a_rank$mean_SRI[s6a_rank$Group != "NLR"], na.rm = TRUE),
             linetype = "dashed", color = "#d95f02", alpha = 0.7, linewidth = 0.6) +
  geom_segment(aes(x = 0, xend = mean_SRI, y = Locus_ord, yend = Locus_ord),
               linewidth = 0.5, alpha = 0.45) +
  geom_errorbarh(aes(xmin = sri_low, xmax = sri_high), height = 0.3, linewidth = 0.45, alpha = 0.85) +
  geom_point(aes(size = mean_sv_events), alpha = 0.92) +
  scale_color_manual(values = c("NLR (Classes I-II, n = 26)" = "#1f77b4",
                                "non-NLR (Classes III-V, n = 10)" = "#d95f02")) +
  scale_size_continuous(range = c(1.6, 4.4), breaks = c(5, 15, 25), name = "Mean SV events") +
  scale_x_continuous(limits = c(0, 0.90), breaks = seq(0, 0.8, 0.2), expand = expansion(mult = c(0.01, 0.03))) +
  labs(
    title = "f  Per-locus microsynteny retention across all 36 cloned R-gene loci",
    subtitle = "Mean SRI +/- s.e.m. (dashed lines: NLR 0.20 vs non-NLR 0.55); size = SV events / +/-500 kb",
    x = "Synteny Retention Index (SRI)",
    y = "Cloned R-gene locus",
    color = NULL
  ) +
  guides(
    color = guide_legend(order = 1, nrow = 1),
    size  = guide_legend(order = 2, nrow = 1)
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.0, color = "#444444"),
    axis.title.x = element_text(size = 9.5, face = "bold", color = "black", margin = margin(t = 4)),
    axis.title.y = element_text(size = 9.5, face = "bold", color = "black", margin = margin(r = 5)),
    axis.text.y = element_text(size = 7.0, face = "bold", color = "black"),
    axis.text.x = element_text(size = 8.8, color = "black"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 7.8),
    legend.title = element_text(size = 7.8, face = "bold"),
    legend.spacing.y = unit(0.5, "mm")
  )

# ------------------------------------------------------------------------------
# Compose Full Figure 2 with Patchwork
# ------------------------------------------------------------------------------
cat("Composing full 6-panel Figure 2 with patchwork...\n")

fig2_composed <- (p_a | p_b) /
                 (p_c | p_d) /
                 (p_e | p_f) +
                 plot_layout(heights = c(0.95, 1.15, 1.0))

out_pdf <- file.path(fig_dir, "Fig2_chromosomal_distribution_review.pdf")
out_png <- file.path(fig_dir, "Fig2_chromosomal_distribution_review.png")

cat(sprintf("Saving PDF -> %s ...\n", out_pdf))
pdf(out_pdf, width = 16.5, height = 16)
print(fig2_composed)
invisible(dev.off())

cat(sprintf("Saving 300-DPI PNG -> %s ...\n", out_png))
png(out_png, width = 16.5, height = 16, units = "in", res = 300)
print(fig2_composed)
invisible(dev.off())

file.copy(out_png, file.path(fig_dir, "Fig2_chromosomal_distribution.png"), overwrite = TRUE)
file.copy(out_png, file.path(fig_dir, "Fig2_chromosomal_distribution_R.png"), overwrite = TRUE)
file.copy(out_pdf, file.path(fig_dir, "Fig2_chromosomal_distribution.pdf"), overwrite = TRUE)
file.copy(out_pdf, file.path(fig_dir, "Fig2_chromosomal_distribution_R.pdf"), overwrite = TRUE)

cat("\n>>> Figure 2 successfully generated and verified! <<<\n")
cat("=======================================================================\n")
