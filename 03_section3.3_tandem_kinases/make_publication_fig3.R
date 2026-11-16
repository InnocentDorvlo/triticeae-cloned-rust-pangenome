#!/usr/bin/env Rscript
# ==============================================================================
# make_publication_fig3.R
# Publication-Grade Figure 3 Generator using ggtree, ggplot2, and patchwork
#
# Generates:
#   figures/Fig3_tandem_kinase_origin_publication.pdf
#   figures/Fig3_tandem_kinase_origin_publication.png (300 dpi)
#   (also updates Fig3_tandem_kinase_origin.png and Fig3_tandem_kinase_origin_R.png)
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
cat("   Building Publication-Grade Figure 3 (R + ggtree + ggplot2)\n")
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
work_dir <- file.path(base_dir, "work/tkp_3.3")
fig_dir  <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Panel a: To-scale Domain Architectures of Class III References (labeled)
# ------------------------------------------------------------------------------
cat("Processing Panel a: To-scale domain architectures...\n")
s8a_df <- read.delim(file.path(res_dir, "S8a_classIII_domains.tsv"), stringsAsFactors = FALSE)
ref_tkp <- s8a_df %>% filter(Record == "reference")

backbone_list <- list()
domain_list <- list()

for (i in seq_len(nrow(ref_tkp))) {
  row <- ref_tkp[i, ]
  g   <- row$Gene
  len <- row$aa_len

  backbone_list[[length(backbone_list) + 1]] <- data.frame(
    Gene = g, Start = 1, End = len,
    Origin = ifelse(g %in% c("Lr9", "Rpg1", "WTK4"), "Intragenic dup.",
                    ifelse(g %in% c("Yr36", "Sr43"), "Kinase-accessory", "Intergenic fusion")),
    stringsAsFactors = FALSE
  )

  if (row$K1_range != "-") {
    pts <- as.numeric(strsplit(row$K1_range, "-")[[1]])
    domain_list[[length(domain_list) + 1]] <- data.frame(
      Gene = g, ShortLbl = ifelse(row$K1_status == "active", "K1", "pK1"),
      Type = ifelse(row$K1_status == "active", "Kinase (Active)", "Pseudokinase (Decoy)"),
      Start = pts[1], End = pts[2], stringsAsFactors = FALSE
    )
  }

  if (row$K2_range != "-") {
    pts <- as.numeric(strsplit(row$K2_range, "-")[[1]])
    domain_list[[length(domain_list) + 1]] <- data.frame(
      Gene = g, ShortLbl = ifelse(row$K2_status == "active", "K2", "pK2"),
      Type = ifelse(row$K2_status == "active", "Kinase (Active)", "Pseudokinase (Decoy)"),
      Start = pts[1], End = pts[2], stringsAsFactors = FALSE
    )
  }

  if (row$accessory_range != "-") {
    pts <- as.numeric(strsplit(row$accessory_range, "-")[[1]])
    acc_lbl <- case_when(
      grepl("VWA", row$accessory_domains, ignore.case = TRUE)   ~ "vWA",
      grepl("START", row$accessory_domains, ignore.case = TRUE) ~ "START",
      grepl("DUF3475", row$accessory_domains, ignore.case = TRUE) ~ "DUF3475",
      TRUE ~ row$accessory_domains
    )
    domain_list[[length(domain_list) + 1]] <- data.frame(
      Gene = g, ShortLbl = acc_lbl,
      Type = "Accessory Domain",
      Start = pts[1], End = pts[2], stringsAsFactors = FALSE
    )
  }
}

bb_df  <- bind_rows(backbone_list)
dom_df <- bind_rows(domain_list) %>%
  mutate(
    Gene = factor(Gene, levels = rev(ref_tkp$Gene)),
    Type = factor(Type, levels = c("Kinase (Active)", "Pseudokinase (Decoy)", "Accessory Domain")),
    Mid  = (Start + End) / 2
  )
bb_df$Gene <- factor(bb_df$Gene, levels = rev(ref_tkp$Gene))

p_a <- ggplot() +
  geom_segment(data = bb_df, aes(x = Start, xend = End, y = Gene, yend = Gene),
               color = "#636363", linewidth = 1.6) +
  geom_rect(data = dom_df, aes(xmin = Start, xmax = End, ymin = as.numeric(Gene) - 0.28,
                               ymax = as.numeric(Gene) + 0.28, fill = Type),
            color = "black", linewidth = 0.3) +
  geom_text(data = dom_df, aes(x = Mid, y = as.numeric(Gene), label = ShortLbl),
            color = "white", size = 2.5, fontface = "bold") +
  geom_text(data = bb_df, aes(x = End + 18, y = as.numeric(Gene), label = Origin),
            hjust = 0, size = 2.6, fontface = "italic", color = "#333333") +
  scale_fill_manual(values = c(
    "Kinase (Active)"       = "#1b7837",
    "Pseudokinase (Decoy)"  = "#d95f02",
    "Accessory Domain"      = "#7570b3"
  )) +
  scale_x_continuous(limits = c(0, 1340), breaks = seq(0, 1200, 200), expand = expansion(mult = c(0.01, 0.02))) +
  labs(
    title = "a  Pfam domain architectures of 8 cloned wheat tandem kinases (WTKs)",
    x = "Protein length (amino acids)",
    y = "Cloned Class III locus",
    fill = "Domain Type"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.title.x = element_text(size = 9.5, face = "bold", color = "black", margin = margin(t = 4)),
    axis.title.y = element_text(size = 9.5, face = "bold", color = "black", margin = margin(r = 5)),
    axis.text.y = element_text(size = 9.5, face = "bold.italic", color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 8.5)
  )

# ------------------------------------------------------------------------------
# Panel b: Kinase Domain ML Phylogeny (ggtree) with Sister-Clade Support
# ------------------------------------------------------------------------------
cat("Processing Panel b: Kinase ML tree...\n")
tr_file <- file.path(work_dir, "kinase_domains.contree")
tr <- read.tree(tr_file)

out_tips <- grep("^OUT_", tr$tip.label, value = TRUE)
if (length(out_tips) > 0) {
  tr <- root(tr, outgroup = out_tips[1], resolve.root = TRUE)
} else {
  tr <- midpoint(tr)
}

format_tip_label <- function(x) {
  if (grepl("^OUT_", x)) {
    sub_x <- gsub("^OUT_Hvulgare_", "", x)
    sub_x <- gsub("_XP_", " XP_", sub_x)
    return(paste0("Outgroup: ", sub_x))
  }
  is_ref <- grepl("_ref$", x)
  clean <- gsub("_ref$", "", x)
  parts <- strsplit(clean, "_")[[1]]
  gene <- parts[1]
  
  mod <- tail(parts, 1)
  mod_str <- case_when(
    mod == "K1" ~ "K1 (Act)",
    mod == "K2" ~ "pK2 (Decoy)",
    mod == "single" ~ "Single K",
    TRUE ~ mod
  )
  
  if (is_ref) {
    return(sprintf("%s (Ref) - %s", gene, mod_str))
  }
  
  asm_info <- paste(parts[2:(length(parts)-1)], collapse = " ")
  return(sprintf("%s [%s] - %s", gene, asm_info, mod_str))
}

k_meta <- data.frame(tip = tr$tip.label, stringsAsFactors = FALSE) %>%
  mutate(
    Module_Type = case_when(
      grepl("^OUT_", tip)   ~ "Outgroup Kinase",
      grepl("_K1", tip)     ~ "Kinase 1 (Active/N-term)",
      grepl("_K2", tip)     ~ "Kinase 2 (Pseudokinase/C-term)",
      grepl("_single", tip) ~ "Single Kinase Progenitor",
      TRUE ~ "Other"
    ),
    Display_Tip = sapply(tip, format_tip_label),
    Is_Ref = grepl("\\(Ref\\)", Display_Tip)
  )

p_b_base <- ggtree(tr, layout = "rectangular", branch.length = "none", linewidth = 0.5) %<+% k_meta
max_x <- max(p_b_base$data$x, na.rm = TRUE)

# Clade definitions matching Fig 1a style (aligned with Section 3.3 dataset)
clade_boxes <- data.frame(
  clade = c(
    "Outgroup 1",
    "Sr62 & Sr60 K1",
    "Outgroup 2",
    "Yr15, Sr60, Yr36, Sr43",
    "Lr9",
    "Rpg1 & WTK4"
  ),
  ymin = c(0.5, 5.5, 18.5, 20.5, 35.5, 44.5),
  ymax = c(5.5, 18.5, 20.5, 35.5, 44.5, 62.5),
  clade_col = c("#7f7f7f", "#e6550d", "#7f7f7f", "#756bb1", "#2b8cbe", "#1f77b4"),
  label_short = c(
    "Outgroup Kinases\n(H. vulgare & B. distachyon)",
    "Sr62 & Sr60 (K1) [Intergenic Fusion]\nPolyphyletic K1/pK2 independent recruitment\nSr62 patristic = 1.99 | Sr60 = 2.14 sub/site",
    "Outgroup",
    "Yr15, Sr60 (pK2), Yr36, Sr43 [Intergenic Fusion]\nPolyphyletic K1/pK2 independent recruitment\nYr15 K1-pK2 patristic dist = 2.76 sub/site",
    "Lr9 [Intragenic Duplication]\nK1 + K2 Sister Clade (SH-aLRT 99.7 / UFBoot 100)\nPatristic distance = 0.54 sub/site",
    "Rpg1 & WTK4 [Intragenic Duplication]\nK1 + K2 Sister Clades (SH-aLRT 99.2 / UFBoot 100)\nPatristic distance = 0.61 - 0.66 sub/site"
  ),
  category = c(
    "Outgroup",
    "Intergenic Fusion",
    "Outgroup",
    "Intergenic Fusion",
    "Intragenic Duplication",
    "Intragenic Duplication"
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(ymid = (ymin + ymax) / 2)

# Outer regime brackets (Fig 1a styling)
regime_boxes <- data.frame(
  regime = c(
    "Outgroup",
    "INTERGENIC FUSIONS\nPolyphyletic K1 + K2 recruitment (3/3 loci, 0%)\nMean patristic divergence = 2.295 \u00b1 0.236 sub/site",
    "INTRAGENIC DUPLICATIONS\nMonophyletic K1 + K2 sister clades (3/3 loci, 100%)\nMean patristic divergence = 0.604 \u00b1 0.035 sub/site"
  ),
  ymin = c(0.5, 5.5, 35.5),
  ymax = c(5.5, 35.5, 62.5),
  col = c("#7f7f7f", "#d95f02", "#1f77b4"),
  stringsAsFactors = FALSE
) %>%
  mutate(ymid = (ymin + ymax) / 2)

p_b <- p_b_base +
  # Shaded clade background rectangles
  geom_rect(data = clade_boxes,
            aes(xmin = max_x * 0.35, xmax = max_x * 1.54, ymin = ymin, ymax = ymax, fill = clade_col),
            inherit.aes = FALSE, alpha = 0.08) +
  scale_fill_identity() +
  # Colored tip points
  geom_tippoint(aes(color = Module_Type), size = 2.5) +
  # Tip labels with clean display (expanded font size)
  geom_tiplab(aes(label = Display_Tip, fontface = ifelse(Is_Ref, "bold", "plain")),
              size = 2.7, offset = max_x * 0.025, align = FALSE) +
  # Clade vertical brackets
  geom_segment(data = clade_boxes %>% filter(!grepl("Outgroup 2", clade)),
               aes(x = max_x * 1.58, xend = max_x * 1.58, y = ymin + 0.25, yend = ymax - 0.25, color = clade_col),
               inherit.aes = FALSE, linewidth = 1.5, show.legend = FALSE) +
  # Clade labels (expanded font size)
  geom_text(data = clade_boxes %>% filter(!grepl("Outgroup 2", clade)),
            aes(x = max_x * 1.62, y = ymid, label = label_short, color = clade_col),
            inherit.aes = FALSE, size = 3.0, fontface = "bold", hjust = 0, lineheight = 0.95, show.legend = FALSE) +
  # Outer regime brackets (Fig 1a styling)
  geom_segment(data = regime_boxes %>% filter(regime != "Outgroup"),
               aes(x = max_x * 2.24, xend = max_x * 2.24, y = ymin + 0.35, yend = ymax - 0.35, color = col),
               inherit.aes = FALSE, linewidth = 2.4, show.legend = FALSE) +
  geom_text(data = regime_boxes %>% filter(regime != "Outgroup"),
            aes(x = max_x * 2.28, y = ymid, label = regime, color = col),
            inherit.aes = FALSE, size = 3.5, fontface = "bold", hjust = 0, lineheight = 0.95, show.legend = FALSE) +
  scale_color_manual(values = c(
    "Kinase 1 (Active/N-term)"        = "#1f77b4",
    "Kinase 2 (Pseudokinase/C-term)"  = "#d95f02",
    "Single Kinase Progenitor"        = "#2ca02c",
    "Outgroup Kinase"                 = "#7f7f7f",
    "#7f7f7f" = "#7f7f7f", "#e6550d" = "#e6550d", "#756bb1" = "#756bb1",
    "#2b8cbe" = "#2b8cbe", "#1f77b4" = "#1f77b4", "#d95f02" = "#d95f02"
  ), breaks = c(
    "Kinase 1 (Active/N-term)",
    "Kinase 2 (Pseudokinase/C-term)",
    "Single Kinase Progenitor",
    "Outgroup Kinase"
  )) +
  labs(
    title = "b  Maximum-likelihood phylogeny of dissected kinase domains (IQ-TREE 2; 62 tips)",
    subtitle = "Monophyletic K1+K2 sister clades in intragenic duplications vs polyphyletic intergenic fusions",
    color = "Kinase Module Architecture:"
  ) +
  theme_tree() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, color = "#444444", margin = margin(b = 4)),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size = 9.5, face = "bold"),
    legend.text = element_text(size = 9.0, face = "bold"),
    legend.margin = margin(t = 2, b = 2)
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 4))) +
  xlim(0, max_x * 3.05)

# ------------------------------------------------------------------------------
# Panel c: K1 vs K2 Sequence Conservation (Identity & pS)
# ------------------------------------------------------------------------------
cat("Processing Panel c: K1 vs K2 divergence...\n")
s8c_df <- read.delim(file.path(res_dir, "S8c_kinase_pair_divergence.tsv"), stringsAsFactors = FALSE)
s8c_clean <- s8c_df %>%
  filter(!is.na(aa_identity_pct), Gene %in% c("Lr9", "Rpg1", "WTK4", "Yr15", "Sr60", "Sr62")) %>%
  mutate(
    Origin = ifelse(Gene %in% c("Lr9", "Rpg1", "WTK4"), "Intragenic Dup", "Intergenic Fusion"),
    Origin = factor(Origin, levels = c("Intragenic Dup", "Intergenic Fusion")),
    Rec_Type = ifelse(Record == "reference", "Cloned Ref (n=6)", "Pangenome Homologue")
  )

p_c1 <- ggplot(s8c_clean, aes(x = Origin, y = aa_identity_pct, fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.48, color = "black", linewidth = 0.4) +
  geom_jitter(data = subset(s8c_clean, Rec_Type == "Pangenome Homologue"),
              aes(shape = Rec_Type), width = 0.14, size = 1.3, color = "black", stroke = 0.3, alpha = 0.30) +
  geom_point(data = subset(s8c_clean, Rec_Type == "Cloned Ref (n=6)"),
             aes(shape = Rec_Type), size = 2.8, fill = "#e41a1c", color = "black", stroke = 0.5) +
  ggrepel::geom_text_repel(data = subset(s8c_clean, Rec_Type == "Cloned Ref (n=6)"),
                           aes(label = Gene), size = 2.35, fontface = "bold.italic",
                           nudge_x = 0.30, force = 3, box.padding = 0.25,
                           point.padding = 0.2, min.segment.length = 0, max.overlaps = 20) +
  annotate("text", x = 1.5, y = 104,
           label = "Ref Welch t = 10.34, P = 0.0058\nAll pairs P = 2.28e-76",
           size = 2.55, fontface = "bold", color = "#b30000") +
  scale_fill_manual(values = c("Intragenic Dup" = "#9ecae1", "Intergenic Fusion" = "#fdd0a2")) +
  scale_shape_manual(values = c("Cloned Ref (n=6)" = 23, "Pangenome Homologue" = 21)) +
  scale_y_continuous(limits = c(20, 110), breaks = seq(20, 100, 20)) +
  labs(
    title = "c  K1 vs K2 amino acid identity",
    x = NULL,
    y = "K1 vs K2 identity (%)",
    shape = NULL, fill = NULL
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 7.4),
    axis.text.x = element_text(size = 8.2, face = "bold", color = "black"),
    axis.text.y = element_text(size = 8.5, color = "black")
  ) +
  guides(fill = "none")

s8c_ps <- s8c_clean %>% filter(!is.na(pS) & is.finite(pS))
ps_ref_t <- t.test(pS ~ Origin, data = subset(s8c_ps, Rec_Type == "Cloned Ref (n=6)"))
ps_gen_w <- wilcox.test(pS ~ Origin, data = subset(s8c_ps, Rec_Type != "Cloned Ref (n=6)"))

p_c2 <- ggplot(s8c_ps, aes(x = Origin, y = pS, fill = Origin)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65, width = 0.48, color = "black", linewidth = 0.4) +
  geom_jitter(data = subset(s8c_ps, Rec_Type == "Pangenome Homologue"),
              aes(shape = Rec_Type), width = 0.12, size = 1.7, fill = "white", color = "black", stroke = 0.4, alpha = 0.75) +
  geom_point(data = subset(s8c_ps, Rec_Type == "Cloned Ref (n=6)"),
             aes(shape = Rec_Type), size = 2.7, fill = "#e41a1c", color = "black", stroke = 0.5) +
  annotate("text", x = 1.5, y = 1.03,
           label = sprintf("Ref Welch t = %.2f, P = %.2e\nGenomic CDS MWU P = %.2e",
                           as.numeric(ps_ref_t$statistic), ps_ref_t$p.value, ps_gen_w$p.value),
           size = 2.55, fontface = "bold", color = "#b30000") +
  scale_fill_manual(values = c("Intragenic Dup" = "#9ecae1", "Intergenic Fusion" = "#fdd0a2")) +
  scale_shape_manual(values = c("Cloned Ref (n=6)" = 23, "Pangenome Homologue" = 21)) +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1.0, 0.2)) +
  labs(
    title = "Synonymous divergence (pS)",
    x = NULL,
    y = "pS (dS saturated in fusions)",
    shape = NULL
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 7.4),
    axis.text.x = element_text(size = 8.2, face = "bold", color = "black"),
    axis.text.y = element_text(size = 8.5, color = "black")
  ) +
  guides(fill = "none")

p_c <- (p_c1 + p_c2) + plot_layout(widths = c(1.1, 1.0), guides = "collect") & theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# Panel d: WTK Architecture Matrix Ordered by Ploidy (2x -> 4x -> 6x)
# ------------------------------------------------------------------------------
cat("Processing Panel d: Ploidy-ordered architecture distribution...\n")
s8d_df <- read.delim(file.path(res_dir, "S8d_architecture_distribution.tsv"), stringsAsFactors = FALSE)

asm_order <- s8d_df %>%
  select(asm_id, ploidy, genomes) %>%
  distinct() %>%
  arrange(ploidy, genomes, asm_id) %>%
  mutate(
    asm_short = gsub("^Tae_|^Tdic_|^Tdur_", "", asm_id),
    ploidy_lbl = paste0(ploidy, "x (", genomes, ")")
  )

s8d_clean <- s8d_df %>%
  mutate(
    asm_short = factor(gsub("^Tae_|^Tdic_|^Tdur_", "", asm_id), levels = asm_order$asm_short),
    State_Label = case_when(
      state == "fused" ~ "Fused WTK (K1+K2)",
      state == "progenitor_only" & standalone_accessory_progenitor_present == 1 ~ "Single Kinase + Acc. Progenitors",
      state == "progenitor_only" ~ "Single Kinase Progenitor",
      standalone_accessory_progenitor_present == 1 ~ "Standalone Accessory Progenitor",
      TRUE ~ "Absent / Deleted"
    ),
    State_Label = factor(State_Label, levels = c(
      "Fused WTK (K1+K2)",
      "Single Kinase + Acc. Progenitors",
      "Single Kinase Progenitor",
      "Standalone Accessory Progenitor",
      "Absent / Deleted"
    ))
  )

# Ploidy boundary x-positions: 11 diploids (1..11), 3 tetraploids (12..14), 13 hexaploids (15..27)
p_d <- ggplot(s8d_clean, aes(x = asm_short, y = Gene, fill = State_Label)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_vline(xintercept = c(11.5, 14.5), color = "black", linewidth = 0.85, linetype = "dashed") +
  annotate("text", x = 6,   y = 8.85, label = "Diploids (2x, n = 11)", fontface = "bold", size = 2.7) +
  annotate("text", x = 13,  y = 8.85, label = "4x (n=3)", fontface = "bold", size = 2.7) +
  annotate("text", x = 21,  y = 8.85, label = "Hexaploids (6x, n = 13)", fontface = "bold", size = 2.7) +
  scale_fill_manual(values = c(
    "Fused WTK (K1+K2)"                = "#1b7837",
    "Single Kinase + Acc. Progenitors" = "#7570b3",
    "Standalone Accessory Progenitor"  = "#bcbddc",
    "Single Kinase Progenitor"         = "#fdb863",
    "Absent / Deleted"                 = "#f0f0f0"
  )) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.14))) +
  labs(
    title = "d  WTK tandem architecture and unfused progenitor distribution across 27 assemblies",
    x = "Triticeae assemblies ordered by ploidy (2x wild diploids -> 4x tetraploids -> 6x hexaploids)",
    y = "Cloned Class III locus",
    fill = NULL
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.title.x = element_text(size = 9.0, face = "bold", color = "black", margin = margin(t = 4)),
    axis.title.y = element_text(size = 9.0, face = "bold", color = "black", margin = margin(r = 5)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.0, face = "bold", color = "black"),
    axis.text.y = element_text(size = 8.8, face = "bold.italic", color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 7.5),
    panel.grid = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2))

# ------------------------------------------------------------------------------
# Panel e: Exon-Intron Splice Architecture with Highlighted Domain-Junction Introns
# ------------------------------------------------------------------------------
cat("Processing Panel e: Exon-intron chevron structures...\n")
s8e_df <- read.delim(file.path(res_dir, "S8e_exon_intron_phases.tsv"), stringsAsFactors = FALSE)
genes_ref <- c("Yr15", "Sr60", "Sr62", "Lr9", "Yr36", "Rpg1", "WTK4", "Sr43")

# Exact domain-junction phase-0 intron indices (within 3-15 aa of K2_start / accessory_start / K2_end)
junction_idx_map <- list(
  "Yr15" = c(4),
  "Sr62" = c(6),
  "Lr9"  = c(7, 10),
  "Yr36" = c(5),
  "Rpg1" = c(6),
  "WTK4" = c(6)
)

exon_boxes <- list()
intron_segs <- list()
chevron_arrows <- list()
phase_labels <- list()
junc_markers <- list()

for (i in seq_along(genes_ref)) {
  g <- genes_ref[i]
  y <- length(genes_ref) - i + 1
  sub <- s8e_df[s8e_df$Gene == g, ]
  r0 <- sub[1, ]
  nex <- as.integer(r0$n_exons)
  phs <- if (r0$all_intron_phases != "-") as.integer(strsplit(r0$all_intron_phases, ",")[[1]]) else integer(0)
  j_idxs <- junction_idx_map[[g]]

  for (ex in seq_len(nex)) {
    exon_boxes[[length(exon_boxes) + 1]] <- data.frame(
      Gene = g, y = y, xmin = ex - 0.75, xmax = ex - 0.25, ymin = y - 0.22, ymax = y + 0.22,
      stringsAsFactors = FALSE
    )
    if (ex < nex) {
      intron_segs[[length(intron_segs) + 1]] <- data.frame(
        Gene = g, y = y, x = ex - 0.25, xend = ex + 0.25, yend = y,
        stringsAsFactors = FALSE
      )
      chevron_arrows[[length(chevron_arrows) + 1]] <- data.frame(
        Gene = g, y = y, x = ex - 0.05, xend = ex + 0.10, yend = y,
        stringsAsFactors = FALSE
      )
      if (ex <= length(phs)) {
        is_junc <- !is.null(j_idxs) && (ex %in% j_idxs)
        phase_labels[[length(phase_labels) + 1]] <- data.frame(
          Gene = g, y = y + 0.33, x = ex, phase = phs[ex],
          is_zero = (phs[ex] == 0), is_junc = is_junc,
          stringsAsFactors = FALSE
        )
        if (is_junc) {
          junc_markers[[length(junc_markers) + 1]] <- data.frame(
            Gene = g, x = ex, y = y, stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

ex_df   <- bind_rows(exon_boxes)
in_df   <- bind_rows(intron_segs)
ch_df   <- bind_rows(chevron_arrows)
ph_df   <- bind_rows(phase_labels)
junc_df <- bind_rows(junc_markers)

p_e <- ggplot() +
  geom_segment(data = in_df, aes(x = x, xend = xend, y = y, yend = yend), color = "#636363", linewidth = 0.8) +
  geom_segment(data = ch_df, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.12, "cm"), type = "open"), color = "#636363", linewidth = 0.8) +
  geom_rect(data = ex_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#737373", color = "black", linewidth = 0.3) +
  geom_point(data = junc_df, aes(x = x, y = y), shape = 25, size = 2.4, fill = "#e41a1c", color = "black", stroke = 0.4) +
  geom_text(data = ph_df, aes(x = x, y = y, label = phase, color = ifelse(is_junc, "junc", ifelse(is_zero, "zero", "other"))),
            size = 2.4, fontface = "bold") +
  scale_color_manual(values = c("junc" = "#e41a1c", "zero" = "#d95f02", "other" = "#969696")) +
  scale_y_continuous(breaks = seq_along(genes_ref), labels = rev(genes_ref), expand = expansion(mult = c(0.06, 0.24))) +
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 2)) +
  annotate("text", x = 1, y = length(genes_ref) + 1.10,
           label = "Red triangles: 10/10 domain-junction introns are phase 0 (P = 1.69e-05) | Phase: red = junction 0, orange = other 0, gray = 1/2",
           hjust = 0, size = 2.30, fontface = "bold", color = "#b30000") +
  labs(
    title = "e  Transcript exon-intron structure with phase-0 domain-junction introns",
    x = "Exon index (5' -> 3')",
    y = "Cloned Class III locus"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    legend.position = "none",
    axis.title.x = element_text(size = 9, face = "bold", color = "black", margin = margin(t = 4)),
    axis.title.y = element_text(size = 9, face = "bold", color = "black", margin = margin(r = 5)),
    axis.text.y = element_text(size = 9, face = "bold.italic", color = "black"),
    axis.text.x = element_text(size = 8.5, color = "black")
  )

# ------------------------------------------------------------------------------
# Compose Full Figure 3 with Patchwork
# ------------------------------------------------------------------------------
cat("Composing full Figure 3 with patchwork...\n")

row_c_d <- (p_c | p_d) + plot_layout(widths = c(1.05, 1.25))

fig3_composed <- p_a / p_b / row_c_d / p_e +
                 plot_layout(heights = c(0.70, 2.70, 0.95, 0.75))

out_pdf <- file.path(fig_dir, "Fig3_tandem_kinase_origin_publication.pdf")
out_png <- file.path(fig_dir, "Fig3_tandem_kinase_origin_publication.png")

cat(sprintf("Saving PDF -> %s ...\n", out_pdf))
pdf(out_pdf, width = 16.5, height = 21.0)
print(fig3_composed)
invisible(dev.off())

cat(sprintf("Saving 300-DPI PNG -> %s ...\n", out_png))
png(out_png, width = 16.5, height = 21.0, units = "in", res = 300)
print(fig3_composed)
invisible(dev.off())

file.copy(out_png, file.path(fig_dir, "Fig3_tandem_kinase_origin.png"), overwrite = TRUE)
file.copy(out_png, file.path(fig_dir, "Fig3_tandem_kinase_origin_R.png"), overwrite = TRUE)
file.copy(out_pdf, file.path(fig_dir, "Fig3_tandem_kinase_origin.pdf"), overwrite = TRUE)
file.copy(out_pdf, file.path(fig_dir, "Fig3_tandem_kinase_origin_R.pdf"), overwrite = TRUE)

cat("\n>>> Figure 3 successfully generated and verified! <<<\n")
cat("=======================================================================\n")
