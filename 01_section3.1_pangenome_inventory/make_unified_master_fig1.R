#!/usr/bin/env Rscript
# ==============================================================================
# Script: make_unified_master_fig1.R
# Purpose: Generate Unified Master Figure 1:
#   Left Column (Methodological & Evolutionary Foundation):
#     - Panel a: Triticeae super-pangenome phylogeny (42 subgenomes, ggtree)
#     - Panel b: Miniprot raw-genome annotation rescue (+76.4% gain, P = 4.36e-4)
#     - Panel c: Permutation Heaps' law saturation curves (gamma = 0.134 vs 0.007)
#     - Panel d: Pangenome compartment breakdown (NLR vs non-NLR, Fisher OR = 259.7)
#   Right Column (Nature Genetics-style Discovery Panels):
#     - Panel e: Landmark cloned rust resistance PAV matrix across subgenomes
#     - Panel f: Core-to-Dispensable Waterfall ranking plot with gene callout arrows
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(ape)
  library(ggtree)
  library(patchwork)
  library(grid)
  library(ggtext)
})

cat("=======================================================================\n")
cat("   Building Unified Master Figure 1 (R + ggtree + patchwork)\n")
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

# Color palettes
tier_colors <- c(
  "Diploid A (2x)"     = "#e41a1c",
  "Diploid S (2x)"     = "#ff7f00",
  "Diploid D (2x)"     = "#4daf4a",
  "Tetraploid AB (4x)" = "#984ea3",
  "Hexaploid ABD (6x)" = "#377eb8"
)

compartment_colors <- c(
  "Core"      = "#1b7837",
  "Soft-core" = "#7fbf7b",
  "Shell"     = "#fdb863",
  "Cloud"     = "#e66101"
)

rescue_colors <- c(
  "Annotated (GFF3)"       = "#525252",
  "Rescued (misannotated)" = "#fe9929",
  "Rescued (unannotated)"  = "#d95f02",
  "Rescued (split model)"  = "#7570b3"
)

pav_colors <- c(
  "Diploid A (2x)"     = "#e41a1c",
  "Diploid S (2x)"     = "#ff7f00",
  "Diploid D (2x)"     = "#4daf4a",
  "Tetraploid AB (4x)" = "#984ea3",
  "Hexaploid ABD (6x)" = "#377eb8",
  "Absent"             = "#f0f0f0"
)

# ------------------------------------------------------------------------------
# 1. Load Metadata & Presence Matrix
# ------------------------------------------------------------------------------
pm_df   <- read.delim(file.path(res_dir, "og_presence_matrix.tsv"), row.names = 1, check.names = FALSE)
og_df   <- read.delim(file.path(res_dir, "S1b_orthogroup_partition.tsv"), stringsAsFactors = FALSE)
tips_df <- read.delim(file.path(res_dir, "S2b_species_tree_tips.tsv"), stringsAsFactors = FALSE)
s2_df   <- read.delim(file.path(res_dir, "S2_assembly_panel_qc.tsv"), stringsAsFactors = FALSE)
cat_df  <- read.delim(file.path(base_dir, "catalogue_3.1.tsv"), stringsAsFactors = FALSE)

subg_meta <- data.frame(subgenome = rownames(pm_df), stringsAsFactors = FALSE) %>%
  mutate(
    asm_id = gsub("__[ABDSUn]+$", "", subgenome),
    sub = gsub("^.*__", "", subgenome)
  ) %>%
  left_join(s2_df %>% select(asm_id, ploidy, species) %>% distinct(), by = "asm_id") %>%
  mutate(
    Tier = case_when(
      ploidy == 2 & sub == "A" ~ "Diploid A (2x)",
      ploidy == 2 & sub == "S" ~ "Diploid S (2x)",
      ploidy == 2 & sub == "D" ~ "Diploid D (2x)",
      ploidy == 4 ~ "Tetraploid AB (4x)",
      ploidy == 6 ~ "Hexaploid ABD (6x)",
      TRUE ~ "Hexaploid ABD (6x)"
    ),
    Tier = factor(Tier, levels = names(tier_colors))
  )

subg_order <- subg_meta %>%
  arrange(Tier, asm_id, sub) %>%
  pull(subgenome)

# ------------------------------------------------------------------------------
# Panel a: Pangenome species tree (ggtree)
# ------------------------------------------------------------------------------
cat("Step 1/6: Building Panel a (Species phylogeny)...\n")
tree_file <- file.path(res_dir, "species_tree.contree")
tr <- read.tree(tree_file)
if ("Aet_AL878__D" %in% tr$tip.label) {
  tr <- root(tr, outgroup = "Aet_AL878__D", resolve.root = TRUE)
} else {
  tr <- midpoint(tr)
}

tip_meta <- data.frame(tip = tr$tip.label, stringsAsFactors = FALSE) %>%
  left_join(tips_df %>% select(tip, asm_id, subgenome, ploidy, species), by = "tip") %>%
  mutate(
    Ploidy_Tier = case_when(
      ploidy == 2 ~ "Diploid (2x)",
      ploidy == 4 ~ "Tetraploid (4x)",
      ploidy == 6 ~ "Hexaploid (6x)",
      TRUE ~ "Other"
    ),
    Display_Label = gsub("__", " [", paste0(tip, "]"))
  )

p_a_base <- ggtree(tr, layout = "rectangular", branch.length = "none") %<+% tip_meta
max_x <- max(p_a_base$data$x, na.rm = TRUE)
clade_ranges <- p_a_base$data %>%
  filter(isTip) %>%
  mutate(sub_clade = gsub("^.*__", "", label)) %>%
  group_by(sub_clade) %>%
  summarise(ymin = min(y) - 0.4, ymax = max(y) + 0.4, ymid = mean(y), .groups = "drop") %>%
  mutate(
    clade_label = paste0(sub_clade, " subgenome"),
    clade_col = case_when(
      sub_clade == "A" ~ "#e41a1c",
      sub_clade == "B" ~ "#984ea3",
      sub_clade == "D" ~ "#4daf4a",
      TRUE ~ "#666666"
    )
  )

p_a <- p_a_base +
  geom_rect(data = clade_ranges,
            aes(xmin = max_x * 0.55, xmax = max_x * 1.44, ymin = ymin, ymax = ymax, fill = clade_col),
            inherit.aes = FALSE, alpha = 0.08) +
  scale_fill_identity() +
  geom_tippoint(aes(color = Ploidy_Tier), size = 2.0) +
  geom_tiplab(aes(label = Display_Label), size = 2.1, offset = max_x * 0.025, align = FALSE) +
  geom_segment(data = clade_ranges,
               aes(x = max_x * 1.46, xend = max_x * 1.46, y = ymin + 0.15, yend = ymax - 0.15, color = clade_col),
               inherit.aes = FALSE, linewidth = 1.4, show.legend = FALSE) +
  geom_text(data = clade_ranges,
            aes(x = max_x * 1.50, y = ymid, label = sub_clade, color = clade_col),
            inherit.aes = FALSE, size = 3.8, fontface = "bold", hjust = 0, show.legend = FALSE) +
  scale_color_manual(values = c(
    "Diploid (2x)"    = "#2ca02c",
    "Tetraploid (4x)" = "#1f77b4",
    "Hexaploid (6x)"  = "#ff7f0e",
    "#e41a1c" = "#e41a1c", "#984ea3" = "#984ea3", "#4daf4a" = "#4daf4a"
  ), breaks = c("Diploid (2x)", "Tetraploid (4x)", "Hexaploid (6x)")) +
  labs(
    title = "a  Triticeae super-pangenome phylogeny (42 subgenomes)",
    subtitle = "Cladogram of 250 BUSCO proteins (16 GFF3 assemblies; +11 genome-only in f)",
    color = NULL
  ) +
  theme_tree(legend.position = "top") +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8.0, face = "italic", color = "#444444"),
    legend.text = element_text(size = 8.0, face = "bold"),
    legend.key.size = unit(0.3, "cm")
  ) +
  xlim(0, max_x * 1.60)

# ------------------------------------------------------------------------------
# Panel b: Annotation rescue per assembly
# ------------------------------------------------------------------------------
cat("Step 2/6: Building Panel b (Annotation rescue)...\n")
s3_df  <- read.delim(file.path(res_dir, "S3_homologue_loci_rescue.tsv"), stringsAsFactors = FALSE)
s3a_df <- read.delim(file.path(res_dir, "S3a_recovery_per_assembly.tsv"), stringsAsFactors = FALSE)

rescue_summary <- s3_df %>%
  mutate(
    Status = case_when(
      status %in% c("annotated", "annotated_only_not_in_miniprot") ~ "Annotated (GFF3)",
      status == "rescued_misannotated_model" ~ "Rescued (misannotated)",
      status == "rescued_unannotated" ~ "Rescued (unannotated)",
      status == "rescued_split_model" ~ "Rescued (split model)",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Status != "Other") %>%
  count(asm_id, Status) %>%
  group_by(asm_id) %>%
  mutate(Total = sum(n)) %>%
  ungroup() %>%
  left_join(s3a_df %>% select(asm_id, pct_increase), by = "asm_id") %>%
  mutate(
    asm_id = reorder(asm_id, Total),
    Status = factor(Status, levels = c("Annotated (GFF3)", "Rescued (misannotated)", 
                                       "Rescued (unannotated)", "Rescued (split model)"))
  )

rescue_totals <- rescue_summary %>%
  select(asm_id, Total, pct_increase) %>%
  distinct()

p_b <- ggplot(rescue_summary, aes(x = n, y = asm_id, fill = Status)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7, color = "black", linewidth = 0.2) +
  geom_text(data = rescue_totals, aes(x = Total + 1.5, y = asm_id, label = sprintf("+%.0f%%", pct_increase)),
            inherit.aes = FALSE, hjust = 0, size = 2.2, fontface = "bold", color = "#b30000") +
  scale_fill_manual(values = rescue_colors) +
  scale_x_continuous(limits = c(0, 108), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
  annotate("text", x = 106, y = 2.0,
           label = "Direct spliced alignment rescues +76.4% loci",
           size = 2.35, fontface = "bold", color = "#b30000", hjust = 1) +
  annotate("text", x = 106, y = 1.0,
           label = "list('Wilcoxon signed-rank', italic(P) == 4.36 %*% 10^-4)",
           parse = TRUE, size = 2.25, fontface = "bold", color = "#b30000", hjust = 1) +
  labs(
    title = "b  Homologue recovery: GFF3 vs Miniprot",
    x = "Number of verified homologue loci",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 7.0, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 7.0),
    legend.key.size = unit(0.3, "cm")
  ) +
  guides(fill = guide_legend(nrow = 2))

# ------------------------------------------------------------------------------
# Panel c: Permutation Heaps' law curves
# ------------------------------------------------------------------------------
cat("Step 3/6: Building Panel c (Heaps' law curves)...\n")
curves_df <- read.delim(file.path(res_dir, "S3d_heaps_curves.tsv"), stringsAsFactors = FALSE) %>%
  mutate(
    Group_Label = ifelse(group == "NLR", "Classes I-II NLR", "Classes III-V non-NLR")
  )

p_c <- ggplot(curves_df, aes(x = n_genomes, color = Group_Label, fill = Group_Label)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
  geom_line(aes(y = mean), linewidth = 1.1) +
  geom_line(aes(y = core_mean), linewidth = 0.9, linetype = "dashed") +
  scale_color_manual(
    values = c("Classes I-II NLR" = "#1f77b4", "Classes III-V non-NLR" = "#d95f02"),
    labels = c(
      "Classes I-II NLR"      = expression("Classes I-II NLR (" * gamma == 0.134 %+-% 0.011 * ")"),
      "Classes III-V non-NLR" = expression("Classes III-V non-NLR (" * gamma == 0.007 %+-% 0.001 * ")")
    )
  ) +
  scale_fill_manual(
    values = c("Classes I-II NLR" = "#1f77b4", "Classes III-V non-NLR" = "#d95f02"),
    labels = c(
      "Classes I-II NLR"      = expression("Classes I-II NLR (" * gamma == 0.134 %+-% 0.011 * ")"),
      "Classes III-V non-NLR" = expression("Classes III-V non-NLR (" * gamma == 0.007 %+-% 0.001 * ")")
    )
  ) +
  scale_x_continuous(breaks = seq(0, 56, 10), limits = c(1, 56)) +
  scale_y_continuous(limits = c(0, 24), breaks = seq(0, 24, 4)) +
  annotate("text", x = 28, y = 22, 
           label = "Open NLR (solid) vs Saturating Non-NLR (dashed)\nPermutation difference P = 0.0010", 
           size = 2.6, fontface = "italic", color = "black") +
  labs(
    title = "c  Heaps' law pan-genome accumulation & core decay",
    x = "Sampled haploid subgenomes (n)",
    y = "Cumulative orthogroups count",
    color = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 8, color = "black"),
    legend.position = c(0.48, 0.58),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
    legend.text = element_text(size = 7.5)
  )

# ------------------------------------------------------------------------------
# Panel d: Pangenome compartment breakdown
# ------------------------------------------------------------------------------
cat("Step 4/6: Building Panel d (Compartment breakdown)...\n")
part_calc <- function(subset_df, group_label) {
  n_tot <- nrow(subset_df)
  data.frame(
    Group = group_label,
    Compartment = c("Core", "Soft-core", "Shell", "Cloud"),
    Pct = c(
      100 * sum(subset_df$compartment == "core") / n_tot,
      100 * sum(subset_df$compartment == "soft-core") / n_tot,
      100 * sum(subset_df$compartment == "shell") / n_tot,
      100 * sum(subset_df$compartment == "cloud") / n_tot
    ),
    Count = c(
      sum(subset_df$compartment == "core"),
      sum(subset_df$compartment == "soft-core"),
      sum(subset_df$compartment == "shell"),
      sum(subset_df$compartment == "cloud")
    )
  )
}

og_main  <- og_df %>% filter(set == "Main")
og_nlr   <- og_main %>% filter(group == "NLR")
og_non   <- og_main %>% filter(group != "NLR")
lbl_nlr  <- sprintf("NLR (I-II)\n(n = %d)", nrow(og_nlr))
lbl_non  <- sprintf("non-NLR (III-V)\n(n = %d)", nrow(og_non))

d_nlr    <- part_calc(og_nlr, lbl_nlr)
d_nonnlr <- part_calc(og_non, lbl_non)

d_combined <- bind_rows(d_nlr, d_nonnlr) %>%
  mutate(
    Group = factor(Group, levels = c(lbl_nlr, lbl_non)),
    Compartment = factor(Compartment, levels = rev(c("Core", "Soft-core", "Shell", "Cloud")))
  )

core_nlr <- sum(og_nlr$compartment %in% c("core", "soft-core"))
disp_nlr <- nrow(og_nlr) - core_nlr
core_non <- sum(og_non$compartment %in% c("core", "soft-core"))
disp_non <- nrow(og_non) - core_non
ft_part  <- fisher.test(matrix(c(core_non, disp_non, core_nlr, disp_nlr), nrow = 2))
or_part  <- if (is.infinite(ft_part$estimate)) "Inf" else sprintf("%.1f", ft_part$estimate)
p_exp    <- floor(log10(ft_part$p.value))
p_mant   <- ft_part$p.value / (10^p_exp)
part_lbl <- sprintf("list('Fisher OR' == %.1f, italic(P) == %.2f %%*%% 10^%d)",
                    ifelse(is.infinite(ft_part$estimate), 0, ft_part$estimate), p_mant, p_exp)
if (core_nlr == 0 || disp_non == 0) {
  # Haldane-corrected empirical Fisher label
  or_h <- ((core_non + 0.5) * (disp_nlr + 0.5)) / ((disp_non + 0.5) * (core_nlr + 0.5))
  part_lbl <- sprintf("list('Fisher OR' == %.1f~'(Haldane)', italic(P) == %.2f %%*%% 10^%d)", or_h, p_mant, p_exp)
}

p_d_bar <- ggplot(d_combined, aes(x = Group, y = Pct, fill = Compartment)) +
  geom_col(width = 0.55, color = "black", linewidth = 0.2) +
  geom_text(
    data = d_combined %>% filter(Pct >= 9),
    aes(label = sprintf("%.0f%%\n(%d)", Pct, Count)),
    position = position_stack(vjust = 0.5),
    size = 2.3, fontface = "bold", lineheight = 0.9,
    color = ifelse(d_combined %>% filter(Pct >= 9) %>% pull(Compartment) %in% c("Core", "Cloud"), "white", "black")
  ) +
  scale_fill_manual(values = compartment_colors) +
  scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 25)) +
  annotate("text", x = 1.5, y = 110,
           label = part_lbl, parse = TRUE,
           size = 2.25, fontface = "bold", color = "#b30000") +
  labs(
    title = "d  Compartment & occupancy (n = 30 Main OGs)",
    x = NULL,
    y = "OGs (%)",
    fill = NULL
  ) +
  theme_classic(base_size = 9.0) +
  theme(
    plot.title = element_text(face = "bold", size = 9.5),
    axis.text.x = element_text(size = 7.2, face = "bold", color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 7.0),
    legend.key.size = unit(0.26, "cm")
  ) +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE))

og_occ <- og_main %>%
  mutate(
    Group_Short = ifelse(group == "NLR", "NLR\n(n=20)", "non-NLR\n(n=10)"),
    Group_Short = factor(Group_Short, levels = c("NLR\n(n=20)", "non-NLR\n(n=10)")),
    occ_pct = occupancy_all_subgenomes * 100
  )

wt_occ <- wilcox.test(occ_pct ~ Group_Short, data = og_occ)
occ_p_exp <- floor(log10(wt_occ$p.value))
occ_p_mant <- wt_occ$p.value / (10^occ_p_exp)

p_d_occ <- ggplot(og_occ, aes(x = Group_Short, y = occ_pct, fill = Group_Short)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.48, color = "black", linewidth = 0.35) +
  geom_jitter(width = 0.14, size = 1.8, shape = 21, color = "black", stroke = 0.35, alpha = 0.85) +
  geom_hline(yintercept = c(15, 80, 95), linetype = "dashed", color = "#888888", linewidth = 0.3) +
  scale_fill_manual(values = c("NLR\n(n=20)" = "#fdb863", "non-NLR\n(n=10)" = "#7fbf7b")) +
  scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 25)) +
  annotate("text", x = 1.5, y = 110,
           label = sprintf("list('MWU', italic(P) == %.2f %%*%% 10^%d)", occ_p_mant, occ_p_exp),
           parse = TRUE, size = 2.25, fontface = "bold", color = "#b30000") +
  labs(
    title = "",
    x = NULL,
    y = "Subgenome occupancy (%)"
  ) +
  theme_classic(base_size = 9.0) +
  theme(
    plot.title = element_text(face = "bold", size = 9.5),
    axis.text.x = element_text(size = 7.2, face = "bold", color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.position = "none"
  )

p_d <- p_d_bar + p_d_occ + plot_layout(widths = c(1.15, 0.85))

# ------------------------------------------------------------------------------
# Panel e: Per-locus Miniprot recovery gain across Classes I-V (S3b)
# ------------------------------------------------------------------------------
cat("Step 5/6: Building Panel e (Per-locus Miniprot recovery gain)...\n")
s3b_df <- read.delim(file.path(res_dir, "S3b_recovery_per_locus.tsv"), stringsAsFactors = FALSE) %>%
  filter(!is.na(pct_increase) & is.finite(pct_increase)) %>%
  mutate(
    Class_Label = paste0("Class ", Class),
    Class_Label = factor(Class_Label, levels = c("Class I", "Class II", "Class III", "Class IV", "Class V"))
  )

s3b_labels <- s3b_df %>%
  filter(pct_increase >= 180 | LocusGroup %in% c("Yr15", "Lr34", "Lr67"))

p_e <- ggplot(s3b_df, aes(x = Class_Label, y = pct_increase + 1, fill = Class_Label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.72, width = 0.52, color = "black", linewidth = 0.35) +
  geom_jitter(width = 0.16, size = 2.0, shape = 21, color = "black", stroke = 0.4, alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = s3b_labels,
    aes(label = sprintf("%s (+%.0f%%)", LocusGroup, pct_increase)),
    size = 2.3, fontface = "bold.italic", color = "#b30000",
    box.padding = 0.25, point.padding = 0.15, min.segment.length = 0, seed = 42
  ) +
  scale_fill_manual(values = c(
    "Class I"   = "#386cb0",
    "Class II"  = "#7fc97f",
    "Class III" = "#f0027f",
    "Class IV"  = "#beaed4",
    "Class V"   = "#fdc086"
  )) +
  scale_y_log10(
    breaks = c(1, 11, 101, 1001),
    labels = c("0%", "+10%", "+100%", "+1,000%")
  ) +
  labs(
    title = "e  Per-locus Miniprot recovery gain across structural classes (n = 36 loci)",
    x = NULL,
    y = "Recovery gain (log10)"
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    plot.title.position = "plot",
    axis.text.x = element_text(size = 8.5, face = "bold", color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.position = "none"
  )

# ------------------------------------------------------------------------------
# Panel f: Core-to-Dispensable Waterfall (56 subgenomes across 27 assemblies)
# ------------------------------------------------------------------------------
cat("Step 6/6: Building Panel f (Waterfall plot)...\n")

short_og_name <- function(q) {
  dplyr::case_when(
    grepl("Sr21", q) ~ "Sr21/Yr5/Yr7",
    grepl("Sr33", q) ~ "Sr33/Sr50",
    grepl("Sr13", q) ~ "Sr13/Sr27",
    grepl("Rpg1", q) ~ "Rpg1/WTK4 [comp]",
    grepl("Yr27", q) ~ "Yr27/Lr13",
    grepl("Sr22", q) ~ "Sr22",
    q == "YrNAM_Yr10" ~ "YrNAM/Yr10",
    q == "Yr87-Lr85" ~ "Yr87/Lr85",
    TRUE ~ q
  )
}

og_waterfall <- og_df %>%
  arrange(desc(occupancy_all_subgenomes), desc(n_members)) %>%
  mutate(
    rank = row_number(),
    pct_occupancy = occupancy_all_subgenomes * 100,
    short_gene = short_og_name(queries),
    is_landmark = short_gene %in% c("Lr34", "Sr60", "Yr15", "Sr35", "Yr27/Lr13", "Lr10", "Lr21", "Sr26"),
    y_label = ifelse(
      is_landmark,
      sprintf("<b style='color:#b30000;'>%d.<i> %s</i></b>", rank, short_gene),
      sprintf("<span style='color:#333333;'>%d.<i> %s</i></span>", rank, short_gene)
    ),
    Comp_Label = case_when(
      pct_occupancy >= 95 ~ "Core",
      pct_occupancy >= 80 ~ "Soft-core",
      pct_occupancy >= 15 ~ "Shell",
      TRUE                ~ "Cloud"
    ),
    Comp_Label = factor(Comp_Label, levels = c("Core", "Soft-core", "Shell", "Cloud"))
  )

s5d_df <- read.delim(file.path(res_dir, "S5d_subgenome_cnv.tsv"), stringsAsFactors = FALSE) %>%
  mutate(subgenome_id = paste0(asm_id, "__", subgenome)) %>%
  select(orthogroup, subgenome = subgenome_id, intact_copies)

# Map ordered subgenomes to exact integer columns 1..56 and centered 27 assembly ticks
col_map <- subg_meta %>%
  filter(subgenome %in% subg_order) %>%
  mutate(subgenome = factor(subgenome, levels = subg_order)) %>%
  arrange(subgenome) %>%
  mutate(x_idx = row_number()) %>%
  left_join(s2_df %>% select(asm_id, annotation_gff3) %>% distinct(), by = "asm_id")

short_asm_name <- function(a, gff) {
  nm <- sub("^T(ae|sp)_", "", a)
  ifelse(gff == "no", paste0(nm, "*"), nm)
}

asm_ticks <- col_map %>%
  group_by(Tier, asm_id, annotation_gff3) %>%
  summarise(
    x_start = min(x_idx),
    x_end = max(x_idx),
    x_mid = mean(x_idx),
    .groups = "drop"
  ) %>%
  arrange(x_start) %>%
  mutate(
    lbl = short_asm_name(asm_id, annotation_gff3),
    lbl_md = ifelse(
      annotation_gff3 == "no",
      sprintf("<b style='color:#b30000;'>%s</b>", lbl),
      sprintf("<span style='color:#222222;'>%s</span>", lbl)
    )
  )

tier_bounds <- col_map %>%
  group_by(Tier) %>%
  summarise(x_end = max(x_idx), .groups = "drop") %>%
  arrange(x_end) %>%
  filter(x_end < nrow(col_map)) %>%
  pull(x_end) + 0.5

asm_bounds <- asm_ticks %>%
  filter(x_end < nrow(col_map)) %>%
  pull(x_end) + 0.5
asm_bounds <- setdiff(asm_bounds, tier_bounds)

waterfall_long <- pm_df[subg_order, og_waterfall$orthogroup, drop = FALSE] %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(subgenome = rownames(pm_df[subg_order, ])) %>%
  pivot_longer(-subgenome, names_to = "orthogroup", values_to = "presence") %>%
  left_join(col_map %>% select(subgenome, x_idx, Tier), by = "subgenome") %>%
  left_join(og_waterfall %>% select(orthogroup, rank, pct_occupancy), by = "orthogroup") %>%
  left_join(s5d_df, by = c("orthogroup", "subgenome")) %>%
  mutate(
    cn_val = ifelse(is.na(intact_copies), ifelse(presence > 0, 1L, 0L),
                    ifelse(presence > 0 & intact_copies == 0, 1L, intact_copies)),
    cn_bin = case_when(
      presence == 0 ~ "0",
      cn_val == 1   ~ "1",
      cn_val <= 4   ~ "2-4",
      cn_val <= 9   ~ "5-9",
      TRUE          ~ ">=10"
    ),
    cn_bin = factor(cn_bin, levels = c("0", "1", "2-4", "5-9", ">=10"))
  )

cn_cols <- c("0" = "#f0f0f0", "1" = "#bdd7e7", "2-4" = "#6baed6", "5-9" = "#2171b5", ">=10" = "#08306b")

n_core <- sum(og_waterfall$pct_occupancy >= 95)
n_sc   <- sum(og_waterfall$pct_occupancy >= 80)
n_sh   <- sum(og_waterfall$pct_occupancy >= 15)
n_tot_og <- nrow(og_waterfall)

# Top ploidy strip locked strictly to x = 0.5..56.5
p_f_strip <- ggplot(col_map, aes(x = x_idx, y = 1, fill = Tier)) +
  geom_tile(width = 1, height = 1, color = "white", linewidth = 0.12) +
  geom_vline(xintercept = tier_bounds, color = "white", linewidth = 0.9) +
  scale_fill_manual(values = pav_colors[names(pav_colors) != "Absent"], name = "Ploidy Tier:") +
  scale_x_continuous(limits = c(0.5, nrow(col_map) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  labs(title = sprintf("f  Copy-number-resolved waterfall across 27 assemblies (%d OGs x 56 subgenomes)", n_tot_og)) +
  theme_void(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 2)),
    plot.title.position = "plot",
    legend.position = "top",
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.text = element_text(size = 7.0),
    legend.key.size = unit(0.26, "cm"),
    legend.margin = margin(t = 0, b = 2),
    plot.margin = margin(t = 2, r = 1, b = 1, l = 2)
  ) +
  guides(fill = guide_legend(nrow = 1))

# Main 31 x 56 copy-number heatmap locked strictly to x = 0.5..56.5, y = 31.5..0.5
p_f_main <- ggplot(waterfall_long, aes(x = x_idx, y = rank, fill = cn_bin)) +
  geom_tile(width = 1, height = 1, color = "white", linewidth = 0.08) +
  geom_vline(xintercept = asm_bounds, color = "grey78", linewidth = 0.22) +
  geom_vline(xintercept = tier_bounds, color = "#222222", linewidth = 0.55) +
  geom_hline(yintercept = c(n_core + 0.5, n_sc + 0.5, n_sh + 0.5),
             linetype = c("dashed", "dotted", "dashed"), color = "#111111", linewidth = 0.45) +
  scale_fill_manual(values = cn_cols, name = "Intact Copies:") +
  scale_x_continuous(
    limits = c(0.5, nrow(col_map) + 0.5),
    expand = c(0, 0),
    breaks = asm_ticks$x_mid,
    labels = asm_ticks$lbl_md
  ) +
  scale_y_reverse(
    limits = c(n_tot_og + 0.5, 0.5),
    expand = c(0, 0),
    breaks = og_waterfall$rank,
    labels = og_waterfall$y_label
  ) +
  labs(
    x = "27 Triticeae assemblies sorted by ploidy tier (* = 11 genome-only assemblies absent from a-b)",
    y = NULL
  ) +
  theme_classic(base_size = 9.0) +
  theme(
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    axis.ticks.length = unit(0.08, "cm"),
    axis.text.y = ggtext::element_markdown(size = 6.2, hjust = 1, margin = margin(r = 2)),
    axis.text.x = ggtext::element_markdown(size = 5.8, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_text(size = 8.2, face = "bold", color = "black", margin = margin(t = 3)),
    legend.position = "bottom",
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.text = element_text(size = 7.2),
    legend.key.size = unit(0.28, "cm"),
    legend.margin = margin(t = 1, b = 0),
    plot.margin = margin(t = 0, r = 2, b = 2, l = 2)
  ) +
  guides(fill = guide_legend(nrow = 1))

# Aligned right-hand occupancy (%) + compartment bracket track (0..100%)
comp_brackets <- data.frame(
  Comp_Label = factor(c("Core", "Soft-core", "Shell", "Cloud"), levels = c("Core", "Soft-core", "Shell", "Cloud")),
  ymin = c(1 - 0.35, n_core + 1 - 0.35, n_sc + 1 - 0.35, n_sh + 1 - 0.35),
  ymax = c(n_core + 0.35, n_sc + 0.35, n_sh + 0.35, n_tot_og + 0.35),
  ymid = c((1 + n_core) / 2, (n_core + 1 + n_sc) / 2, (n_sc + 1 + n_sh) / 2, (n_sh + 1 + n_tot_og) / 2),
  txt = c(
    sprintf("Core (%d)", n_core),
    sprintf("Soft-core (%d)", n_sc - n_core),
    sprintf("Shell (%d)", n_sh - n_sc),
    sprintf("Cloud (%d)", n_tot_og - n_sh)
  )
)

p_f_occ <- ggplot(og_waterfall, aes(y = rank, x = pct_occupancy, fill = Comp_Label)) +
  geom_vline(xintercept = c(15, 80, 95), linetype = "dashed", color = "#999999", linewidth = 0.3) +
  geom_rect(aes(xmin = 0, xmax = pct_occupancy, ymin = rank - 0.36, ymax = rank + 0.36),
            color = "black", linewidth = 0.12) +
  geom_hline(yintercept = c(n_core + 0.5, n_sc + 0.5, n_sh + 0.5),
             linetype = c("dashed", "dotted", "dashed"), color = "#111111", linewidth = 0.45) +
  geom_segment(data = comp_brackets,
               aes(x = 104, xend = 104, y = ymin, yend = ymax, color = Comp_Label),
               inherit.aes = FALSE, linewidth = 1.6, show.legend = FALSE) +
  geom_text(data = comp_brackets,
            aes(x = 107, y = ymid, label = txt, color = Comp_Label),
            inherit.aes = FALSE, hjust = 0, size = 2.35, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = compartment_colors, guide = "none") +
  scale_color_manual(values = compartment_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 160), breaks = c(0, 50, 100), labels = c("0", "50", "100%"), expand = c(0, 0)) +
  scale_y_reverse(limits = c(n_tot_og + 0.5, 0.5), expand = c(0, 0)) +
  labs(x = "Occupancy (%)", y = NULL) +
  theme_classic(base_size = 9.0) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.4, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.35, color = "black"),
    axis.text.x = element_text(size = 7.0, color = "black"),
    axis.title.x = element_text(size = 8.0, face = "bold", color = "black", margin = margin(t = 3)),
    plot.margin = margin(t = 0, r = 4, b = 2, l = 1)
  )

# ------------------------------------------------------------------------------
# Composing the Unified Master Layout (100% pure patchwork)
# ------------------------------------------------------------------------------
cat("Composing unified layout with patchwork...\n")

top_left <- (p_a | p_b) + plot_layout(widths = c(1.22, 0.88))
bot_left <- (p_c | p_d)
left_col <- (top_left / bot_left) + plot_layout(heights = c(1.26, 0.84))

right_design <- "
  AA
  BC
  DE
"
right_col <- wrap_plots(
  A = p_e + theme(plot.title.position = "plot"),
  B = p_f_strip,
  C = plot_spacer(),
  D = p_f_main,
  E = p_f_occ,
  design = right_design,
  widths = c(1.0, 0.38),
  heights = c(0.40, 0.045, 1.0)
)

master_fig1 <- (left_col | right_col) +
  plot_layout(widths = c(1.12, 1.06))

out_pdf <- file.path(fig_dir, "Fig1_unified_master_publication.pdf")
out_png <- file.path(fig_dir, "Fig1_unified_master_publication.png")

cat(sprintf("Saving PDF -> %s ...\n", out_pdf))
pdf(out_pdf, width = 18, height = 12)
print(master_fig1)
invisible(dev.off())

cat(sprintf("Saving 300-DPI PNG -> %s ...\n", out_png))
png(out_png, width = 18, height = 12, units = "in", res = 300)
print(master_fig1)
invisible(dev.off())

# Also update the primary Fig1_pangenome_inventory.png
file.copy(out_png, file.path(fig_dir, "Fig1_pangenome_inventory.png"), overwrite = TRUE)
file.copy(out_pdf, file.path(fig_dir, "Fig1_pangenome_inventory.pdf"), overwrite = TRUE)

cat("\n>>> Unified Master Figure 1 successfully generated! <<<\n")
cat("=======================================================================\n")
