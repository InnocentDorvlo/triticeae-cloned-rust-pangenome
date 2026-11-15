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
  library(cowplot)
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
    y = "Homologue recovery gain (log10)"
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 8.5, face = "bold", color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.position = "none"
  )

# ------------------------------------------------------------------------------
# Panel f: Core-to-Dispensable Waterfall Plot (Image 1c style)
# ------------------------------------------------------------------------------
cat("Step 6/6: Building Panel f (Waterfall plot)...\n")
og_waterfall <- og_df %>%
  arrange(desc(occupancy_all_subgenomes), desc(n_members)) %>%
  mutate(
    rank = row_number(),
    pct_occupancy = occupancy_all_subgenomes * 100,
    cloned_marker = case_when(
      grepl("Lr34", queries) ~ "Lr34 OG (100% Core)",
      grepl("Yr15", queries) ~ "Yr15 progenitor OG (94.6%)",
      grepl("Sr60", queries) ~ "Sr60 progenitor OG (94.6%)",
      grepl("Sr35", queries) ~ "Sr35 (51.8% Shell)",
      grepl("Yr27", queries) ~ "Yr27 (26.8% Shell)",
      grepl("Lr10", queries) ~ "Lr10 (17.9% Shell)",
      grepl("Lr21", queries) ~ "Lr21 (12.5% Cloud)",
      grepl("Sr26", queries) ~ "Sr26 (7.1% Cloud)",
      TRUE ~ NA_character_
    )
  )

s5d_df <- read.delim(file.path(res_dir, "S5d_subgenome_cnv.tsv"), stringsAsFactors = FALSE) %>%
  mutate(subgenome_id = paste0(asm_id, "__", subgenome)) %>%
  select(orthogroup, subgenome = subgenome_id, intact_copies)

waterfall_long <- pm_df[subg_order, og_waterfall$orthogroup, drop = FALSE] %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(subgenome = rownames(pm_df[subg_order, ])) %>%
  pivot_longer(-subgenome, names_to = "orthogroup", values_to = "presence") %>%
  left_join(subg_meta, by = "subgenome") %>%
  left_join(og_waterfall %>% select(orthogroup, rank, pct_occupancy, cloned_marker), by = "orthogroup") %>%
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
    cn_bin = factor(cn_bin, levels = c("0", "1", "2-4", "5-9", ">=10")),
    subgenome = factor(subgenome, levels = subg_order)
  )

cn_cols <- c("0" = "#f0f0f0", "1" = "#bdd7e7", "2-4" = "#6baed6", "5-9" = "#2171b5", ">=10" = "#08306b")

n_core <- sum(og_waterfall$pct_occupancy >= 95)
n_sc   <- sum(og_waterfall$pct_occupancy >= 80)
n_sh   <- sum(og_waterfall$pct_occupancy >= 15)
n_tot_og <- nrow(og_waterfall)

tier_strip <- subg_meta %>%
  mutate(subgenome = factor(subgenome, levels = subg_order))

p_f_strip <- ggplot(tier_strip, aes(x = subgenome, y = 1, fill = Tier)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_manual(values = pav_colors[names(pav_colors) != "Absent"], name = "Ploidy Tier") +
  coord_cartesian(xlim = c(-26, 76), clip = "off") +
  labs(title = sprintf("f  Copy-number-resolved Core-to-Dispensable waterfall (%d OGs x 56 subgenomes)", n_tot_og)) +
  theme_void(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    legend.position = "top",
    legend.text = element_text(size = 7.2),
    legend.key.size = unit(0.28, "cm")
  ) +
  guides(fill = guide_legend(nrow = 1))

callouts <- og_waterfall %>% filter(!is.na(cloned_marker))

p_f_main <- ggplot() +
  geom_tile(data = waterfall_long, aes(x = subgenome, y = rank, fill = cn_bin), linewidth = 0.1) +
  scale_fill_manual(values = cn_cols, name = "Intact Copies") +
  scale_y_reverse(breaks = c(1, 10, 20, n_tot_og), expand = expansion(mult = c(0.01, 0.02))) +
  geom_hline(yintercept = n_core + 0.5, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = n_sc + 0.5, linetype = "dotted", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = n_sh + 0.5, linetype = "dashed", color = "black", linewidth = 0.6) +
  annotate("text", x = 58, y = max(1, n_core / 2), label = "} Core (>=95%)", hjust = 0, size = 2.6, fontface = "bold", color = "#1b7837") +
  annotate("text", x = 58, y = (n_core + n_sc + 1) / 2, label = "} Soft-core (80-94%)", hjust = 0, size = 2.5, fontface = "bold", color = "#7fbf7b") +
  annotate("text", x = 58, y = (n_sc + n_sh + 1) / 2, label = "} Shell (15-79%)", hjust = 0, size = 2.6, fontface = "bold", color = "#fdb863") +
  annotate("text", x = 58, y = (n_sh + n_tot_og + 1) / 2, label = "} Cloud (<15%)", hjust = 0, size = 2.6, fontface = "bold", color = "#e66101") +
  geom_segment(data = callouts, aes(x = -0.5, xend = -4.0, y = rank, yend = rank),
               arrow = arrow(length = unit(0.10, "cm"), ends = "first", type = "closed"),
               color = "#b30000", linewidth = 0.5) +
  geom_text(data = callouts, aes(x = -4.5, y = rank, label = cloned_marker),
            hjust = 1, size = 2.2, fontface = "bold.italic", color = "#b30000") +
  coord_cartesian(xlim = c(-26, 76), clip = "off") +
  labs(
    x = "56 haploid subgenomes sorted by evolutionary tier",
    y = "Orthogroup Rank (Ordered by Subgenome Occupancy)"
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    axis.title.y = element_text(size = 9, face = "bold", color = "black", margin = margin(r = 3)),
    axis.title.x = element_text(size = 9, face = "bold", color = "black", margin = margin(t = 4)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(0.3, "cm")
  ) +
  guides(fill = guide_legend(nrow = 1))

p_f <- p_f_strip / p_f_main + plot_layout(heights = c(0.08, 1.0))

# ------------------------------------------------------------------------------
# Composing the Unified Master Layout
# ------------------------------------------------------------------------------
cat("Composing unified layout with patchwork...\n")

top_left <- (p_a | p_b) + plot_layout(widths = c(1.22, 0.88))
bot_left <- (p_c | p_d)
left_col <- top_left / bot_left + plot_layout(heights = c(1.28, 0.82))
right_col <- cowplot::plot_grid(p_e, p_f, ncol = 1, rel_heights = c(0.45, 1.0))

master_fig1 <- (left_col | right_col) + 
  plot_layout(widths = c(1.15, 1.0))

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
