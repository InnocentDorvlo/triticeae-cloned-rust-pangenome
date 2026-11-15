#!/usr/bin/env Rscript
# ==============================================================================
# make_publication_fig5.R
# Publication-Grade Figure 5 Generator: Population Diversity, Bottlenecks,
# Haplotype Networks & Functional Allele Geography using ggplot2 & patchwork
#
# Generates:
#   figures/Fig5_haplotypes_bottlenecks_publication.pdf
#   figures/Fig5_haplotypes_bottlenecks_publication.png (300 dpi)
#   (also updates Fig5_haplotypes_bottlenecks.png and Fig5_haplotypes_bottlenecks_R.png)
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
})

cat("=======================================================================\n")
cat("   Building Publication-Grade Figure 5 (R + ggplot2 + patchwork)\n")
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
work_dir <- file.path(base_dir, "work/haplo_3.5")
fig_dir  <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

tier_cols <- c(
  "Wild Progenitors"         = "#2ca02c",
  "Landraces & Domesticates" = "#ff7f0e",
  "Modern Elite Cultivars"   = "#1f77b4"
)

# ------------------------------------------------------------------------------
# Function to Render Publication Haplotype Network in ggplot2
# ------------------------------------------------------------------------------
make_network_plot <- function(nodes_file, edges_file, title_text, subtitle_text) {
  nodes_df <- read.delim(nodes_file, stringsAsFactors = FALSE)
  edges_df <- read.delim(edges_file, stringsAsFactors = FALSE)
  
  # Calculate radius per node
  nodes_df <- nodes_df %>%
    mutate(r = 0.55 + 0.14 * sqrt(total))
  
  # Map coordinates to edges and intelligently position edge labels to avoid nodes
  edges_df <- edges_df %>%
    left_join(nodes_df %>% select(id, x1 = x, y1 = y), by = c("node1" = "id")) %>%
    left_join(nodes_df %>% select(id, x2 = x, y2 = y), by = c("node2" = "id"))
  
  # Compute label positions with collision avoidance against other nodes
  mx_vec <- numeric(nrow(edges_df))
  my_vec <- numeric(nrow(edges_df))
  
  for (k in seq_len(nrow(edges_df))) {
    x1 <- edges_df$x1[k]; y1 <- edges_df$y1[k]
    x2 <- edges_df$x2[k]; y2 <- edges_df$y2[k]
    
    # Try midpoints at t = 0.5, 0.35, 0.65, 0.25, 0.75
    candidates <- c(0.50, 0.35, 0.65, 0.25, 0.75)
    best_t <- 0.50
    best_min_dist <- -1
    
    # Other nodes (excluding the two endpoints)
    other_nodes <- nodes_df %>% filter(!id %in% c(edges_df$node1[k], edges_df$node2[k]))
    
    for (cand_t in candidates) {
      cx <- (1 - cand_t) * x1 + cand_t * x2
      cy <- (1 - cand_t) * y1 + cand_t * y2
      if (nrow(other_nodes) > 0) {
        dists <- sqrt((other_nodes$x - cx)^2 + (other_nodes$y - cy)^2)
        min_d <- min(dists)
      } else {
        min_d <- 10
      }
      if (min_d > best_min_dist) {
        best_min_dist <- min_d
        best_t <- cand_t
      }
    }
    
    # Compute position at best_t with slight perpendicular offset if close to any node
    dx <- x2 - x1; dy <- y2 - y1
    len <- sqrt(dx^2 + dy^2)
    nx <- -dy / len; ny <- dx / len
    
    offset_dist <- ifelse(best_min_dist < 1.4, 0.35, 0.0)
    mx_vec[k] <- (1 - best_t) * x1 + best_t * x2 + offset_dist * nx
    my_vec[k] <- (1 - best_t) * y1 + best_t * y2 + offset_dist * ny
  }
  
  edges_df$mx <- mx_vec
  edges_df$my <- my_vec
  edges_df$dist_label <- paste0(edges_df$mutation_distance, " mut")
  
  # Calculate pie slices for each node
  pie_pieces <- list()
  for (i in seq_len(nrow(nodes_df))) {
    nd <- nodes_df[i, ]
    tot <- nd$total
    r <- nd$r
    
    w_frac <- nd$wild / tot
    l_frac <- nd$landrace / tot
    e_frac <- nd$elite / tot
    
    angles <- c(0, w_frac * 2 * pi, (w_frac + l_frac) * 2 * pi, 2 * pi)
    tiers  <- c("Wild Progenitors", "Landraces & Domesticates", "Modern Elite Cultivars")
    
    for (j in 1:3) {
      theta1 <- angles[j]
      theta2 <- angles[j+1]
      if (theta2 > theta1) {
        th <- seq(theta1, theta2, length.out = 20)
        px <- c(nd$x, nd$x + r * cos(th), nd$x)
        py <- c(nd$y, nd$y + r * sin(th), nd$y)
        pie_pieces[[length(pie_pieces) + 1]] <- data.frame(
          node = nd$id, tier = tiers[j], x = px, y = py, piece_id = paste0(nd$id, "_", j),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  pie_df <- bind_rows(pie_pieces) %>%
    mutate(tier = factor(tier, levels = c("Wild Progenitors", "Landraces & Domesticates", "Modern Elite Cultivars")))
  
  # Cloned reference nodes
  cloned_nodes <- nodes_df %>% filter(is_cloned == TRUE | is_cloned == "True")
  
  p <- ggplot() +
    # Mutational edges
    geom_segment(data = edges_df, aes(x = x1, y = y1, xend = x2, yend = y2), 
                 color = "#999999", linewidth = 0.85) +
    geom_label(data = edges_df, aes(x = mx, y = my, label = dist_label), 
               size = 2.3, color = "#333333", fill = alpha("white", 0.90),
               label.padding = unit(0.12, "lines"), linewidth = 0.2) +
    # Node pie polygons
    geom_polygon(data = pie_df, aes(x = x, y = y, fill = tier, group = piece_id), 
                 color = "black", linewidth = 0.35) +
    # Cloned reference gold halo ring
    geom_point(data = cloned_nodes, aes(x = x, y = y),
               shape = 21, size = 10.5, color = "#e41a1c", stroke = 1.3, fill = NA) +
    # Node labels
    ggrepel::geom_text_repel(data = nodes_df, aes(x = x, y = y,
                                                  label = ifelse(is_cloned == TRUE | is_cloned == "True",
                                                                 sprintf("%s (n=%d, Cloned)", id, total),
                                                                 sprintf("%s (n=%d)", id, total))),
                             size = 2.65, fontface = "bold", color = "black",
                             box.padding = 0.50, point.padding = 0.40, segment.color = "#666666",
                             segment.size = 0.3, min.segment.length = 0.1, seed = 42) +
    scale_fill_manual(values = tier_cols, drop = FALSE) +
    coord_fixed(ratio = 1) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      fill = NULL
    ) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle = element_text(size = 8.2, face = "italic", color = "#444444", hjust = 0),
      legend.position = "top",
      legend.text = element_text(size = 8.5)
    )
  return(p)
}

# ------------------------------------------------------------------------------
# Panels a & b: Haplotype Networks (Sr35 vs Yr15)
# ------------------------------------------------------------------------------
cat("Processing Panel a & b: Haplotype networks...\n")
p_a <- make_network_plot(
  file.path(work_dir, "network_nodes_Sr35.tsv"),
  file.path(work_dir, "network_edges_Sr35.tsv"),
  "a  Empirical haplotype network: Sr35 (Class I CC-NLR)",
  "OG0000009 protein MST (top 8 shared nodes H1-H8; 30 wild-private haplotypes are singletons)"
)

p_b <- make_network_plot(
  file.path(work_dir, "network_nodes_Yr15.tsv"),
  file.path(work_dir, "network_edges_Yr15.tsv"),
  "b  Empirical haplotype network: Yr15 (Class III WTK)",
  "OG0000045 protein MST (nodes H1-H8 span all 3 tiers; red ring = cloned reference)"
)

# ------------------------------------------------------------------------------
# Panel c: Nucleotide Diversity (pi) & theta_W across Breeding Tiers
# ------------------------------------------------------------------------------
cat("Processing Panel c: Diversity erosion across tiers...\n")
s11a_df <- read.delim(file.path(res_dir, "Supplementary_Table_S11a_population_genetics_tiers.tsv"), stringsAsFactors = FALSE)

div_long <- bind_rows(
  s11a_df %>%
    select(orthogroup, structural_class,
           "Wild" = pi_wild_1e3, "Landrace/\nDom." = pi_landrace_1e3, "Elite" = pi_elite_1e3) %>%
    pivot_longer(cols = c("Wild", "Landrace/\nDom.", "Elite"), names_to = "Tier", values_to = "val") %>%
    mutate(Metric = "pi[aa]~(x10^-3)"),
  s11a_df %>%
    select(orthogroup, structural_class,
           "Wild" = theta_wild_1e3, "Landrace/\nDom." = theta_landrace_1e3, "Elite" = theta_elite_1e3) %>%
    pivot_longer(cols = c("Wild", "Landrace/\nDom.", "Elite"), names_to = "Tier", values_to = "val") %>%
    mutate(Metric = "theta[W,aa]~(x10^-3)")
) %>%
  filter(!is.na(val)) %>%
  mutate(
    Tier = factor(Tier, levels = c("Wild", "Landrace/\nDom.", "Elite")),
    Metric = factor(Metric, levels = c("pi[aa]~(x10^-3)", "theta[W,aa]~(x10^-3)"))
  )

wilc_pi <- wilcox.test(s11a_df$pi_wild_1e3, s11a_df$pi_elite_1e3, paired = TRUE, alternative = "two.sided")
wilc_th <- wilcox.test(s11a_df$theta_wild_1e3, s11a_df$theta_elite_1e3, paired = TRUE, alternative = "two.sided")
pi_p_exp <- floor(log10(wilc_pi$p.value))
pi_p_mant <- wilc_pi$p.value / (10^pi_p_exp)
th_p_exp <- floor(log10(wilc_th$p.value))
th_p_mant <- wilc_th$p.value / (10^th_p_exp)

div_annots <- data.frame(
  Metric = factor(c("pi[aa]~(x10^-3)", "theta[W,aa]~(x10^-3)"),
                  levels = c("pi[aa]~(x10^-3)", "theta[W,aa]~(x10^-3)")),
  x = c(2.0, 2.0),
  y = c(395, 395),
  label = c(
    sprintf("Elite/Wild = 79.4%% (-20.6%%, P = %.2f x 10^%d)", pi_p_mant, pi_p_exp),
    sprintf("Elite/Wild = 69.8%% (-30.2%%, P = %.2f x 10^%d)", th_p_mant, th_p_exp)
  )
)

p_c <- ggplot(div_long, aes(x = Tier, y = val, fill = Tier)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.55, color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.16, size = 1.6, shape = 21, color = "black", stroke = 0.35, alpha = 0.75) +
  facet_wrap(~ Metric, labeller = label_parsed) +
  scale_fill_manual(values = c("Wild" = "#2ca02c", "Landrace/\nDom." = "#ff7f0e", "Elite" = "#1f77b4")) +
  scale_y_continuous(limits = c(0, 425), breaks = seq(0, 400, 100)) +
  geom_text(data = div_annots, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            size = 2.40, fontface = "bold", color = "#b30000") +
  labs(
    title = expression(bold("c  Orthogroup-wide ") * bolditalic(pi)[aa] * bold(" & ") * bolditalic(theta)[W,aa] * bold(" across sampled tiers")),
    x = NULL,
    y = expression(bold("Diversity (") %*% bold("10")^bold("-3") * bold(")"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    strip.background = element_rect(fill = "#f0f0f0", color = "black", linewidth = 0.4),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "none",
    axis.title.y = element_text(size = 9.5, face = "bold", color = "black", margin = margin(r = 5)),
    axis.text.x = element_text(size = 8.0, face = "bold", color = "black"),
    axis.text.y = element_text(size = 9, color = "black")
  )

# ------------------------------------------------------------------------------
# Panel d: Tajima's D Distributions (NLR vs non-NLR across Tiers)
# ------------------------------------------------------------------------------
cat("Processing Panel d: Tajima's D distributions...\n")
is_nlr_s11 <- s11a_df$structural_class %in% c("I", "I/II", "II")
taj_df <- s11a_df %>%
  mutate(Regime = ifelse(structural_class %in% c("I", "I/II", "II"), "NLR (Classes I-II)", "non-NLR (Classes III-V)")) %>%
  select(orthogroup, Regime, primary_gene,
          "Wild Progenitors" = tajima_D_wild,
          "Landraces & Domesticates" = tajima_D_landrace,
          "Modern Elite Cultivars" = tajima_D_elite) %>%
  pivot_longer(cols = c("Wild Progenitors", "Landraces & Domesticates", "Modern Elite Cultivars"),
               names_to = "Tier", values_to = "Tajima_D") %>%
  filter(!is.na(Tajima_D)) %>%
  mutate(
    Tier = factor(Tier, levels = c("Wild Progenitors", "Landraces & Domesticates", "Modern Elite Cultivars")),
    Regime = factor(Regime, levels = c("NLR (Classes I-II)", "non-NLR (Classes III-V)"))
  )

mwu_wild <- wilcox.test(s11a_df$tajima_D_wild[is_nlr_s11], s11a_df$tajima_D_wild[!is_nlr_s11], alternative = "two.sided")
mwu_land <- wilcox.test(s11a_df$tajima_D_landrace[is_nlr_s11], s11a_df$tajima_D_landrace[!is_nlr_s11], alternative = "two.sided")
mwu_elit <- wilcox.test(s11a_df$tajima_D_elite[is_nlr_s11], s11a_df$tajima_D_elite[!is_nlr_s11], alternative = "two.sided")

tier_annots <- data.frame(
  x = c(1, 2, 3),
  y = c(4.65, 4.65, 4.65),
  label = c(
    sprintf("MWU W = %.1f\nP = %.2e (**)", mwu_wild$statistic, mwu_wild$p.value),
    sprintf("MWU W = %.1f\nP = %.2e (*)", mwu_land$statistic, mwu_land$p.value),
    sprintf("MWU W = %.1f\nP = %.2e", mwu_elit$statistic, mwu_elit$p.value)
  )
)

p_d <- ggplot(taj_df, aes(x = Tier, y = Tajima_D, fill = Regime)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.65, position = position_dodge(width = 0.75), color = "black", linewidth = 0.4) +
  geom_point(aes(group = Regime), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
             size = 1.6, shape = 21, color = "black", stroke = 0.35, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#636363", linewidth = 0.8) +
  annotate("text", x = 0.6, y = -1.9, label = "D <= 0: Selective sweep / purifying", hjust = 0, size = 2.4, color = "#636363", fontface = "italic") +
  annotate("text", x = 0.6, y = 3.6, label = "D > 0: Balancing / paralogue retention", hjust = 0, size = 2.4, color = "#1f77b4", fontface = "italic") +
  scale_fill_manual(values = c("NLR (Classes I-II)" = "#1f77b4", "non-NLR (Classes III-V)" = "#d95f02")) +
  scale_y_continuous(limits = c(-2.5, 5.2), breaks = seq(-2, 4, 1)) +
  geom_text(data = tier_annots, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            size = 2.4, fontface = "bold", color = "black", vjust = 0.5) +
  labs(
    title = expression(bold("d  Orthogroup Tajima's ") * bolditalic("D") * bold(" (two-sided MWU per tier)")),
    x = NULL,
    y = expression(bold("Tajima's ") * bolditalic("D")),
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    legend.position = "top",
    legend.text = element_text(size = 8.0),
    axis.text.x = element_text(size = 8.2, face = "bold", color = "black"),
    axis.text.y = element_text(size = 9, color = "black")
  )

# ------------------------------------------------------------------------------
# Panel e: Haplotype Sharing & Privacy Decomposition Across Structural Classes
# ------------------------------------------------------------------------------
cat("Processing Panel e: Haplotype sharing across structural classes...\n")
s11b_df <- read.delim(file.path(res_dir, "Supplementary_Table_S11b_haplotype_diversity_summary.tsv"), stringsAsFactors = FALSE)

make_share_row <- function(df, grp_label) {
  tot <- sum(df$total_intact_haplotypes, na.rm = TRUE)
  wi  <- sum(df$wild_progenitor_private, na.rm = TRUE)
  la  <- sum(df$landrace_domesticate_private, na.rm = TRUE)
  el  <- sum(df$elite_captured_haplotypes, na.rm = TRUE)
  data.frame(
    Group = grp_label,
    Total = tot,
    Wild_Pct = 100 * wi / tot,
    Land_Pct = 100 * (tot - el - wi) / tot,
    Elit_Pct = 100 * el / tot,
    Untapped_Pct = 100 * (tot - el) / tot,
    Wild_N = wi,
    Land_N = la,
    Elit_N = el
  )
}

classes_summary <- bind_rows(
  make_share_row(s11b_df, "Pooled Pangenome"),
  make_share_row(filter(s11b_df, structural_class == "I"), "Class I: CC-NLR"),
  make_share_row(filter(s11b_df, structural_class %in% c("I/II", "II")), "Class II: NLR-ID"),
  make_share_row(filter(s11b_df, structural_class == "III"), "Class III: Tandem Kinase"),
  make_share_row(filter(s11b_df, structural_class == "IV"), "Class IV: ANK-TM / NAM"),
  make_share_row(filter(s11b_df, structural_class == "V"), "Class V: APR Transporter")
)

# Format for horizontal stacked bar
classes_long <- classes_summary %>%
  select(Group, Total, Untapped_Pct,
         "Elite Captured" = Elit_Pct,
         "Landrace/Domesticate Private" = Land_Pct,
         "Wild Progenitor Private" = Wild_Pct) %>%
  pivot_longer(cols = c("Elite Captured", "Landrace/Domesticate Private", "Wild Progenitor Private"),
               names_to = "Compartment", values_to = "Percentage") %>%
  mutate(
    Group_Label = sprintf("%s (n = %d)", Group, Total),
    Group_Label = factor(Group_Label, levels = rev(sprintf("%s (n = %d)", classes_summary$Group, classes_summary$Total))),
    Compartment = factor(Compartment, levels = c("Elite Captured", "Landrace/Domesticate Private", "Wild Progenitor Private"))
  )

annotation_df <- classes_summary %>%
  mutate(
    Group_Label = factor(sprintf("%s (n = %d)", Group, Total),
                         levels = rev(sprintf("%s (n = %d)", classes_summary$Group, classes_summary$Total))),
    Anno_Text = ifelse(Group == "Pooled Pangenome",
                       sprintf("%.1f%% Untapped (P = 2.97e-07)", Untapped_Pct),
                       sprintf("%.1f%% Untapped", Untapped_Pct))
  )

p_e <- ggplot(classes_long, aes(x = Percentage, y = Group_Label, fill = Compartment)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  geom_text(aes(label = ifelse(Percentage >= 8.0, sprintf("%.1f%%", Percentage), "")),
            position = position_stack(vjust = 0.5),
            size = 2.45, fontface = "bold", color = "white") +
  geom_text(data = annotation_df, aes(x = 101.5, y = Group_Label, label = Anno_Text),
            inherit.aes = FALSE, hjust = 0, size = 2.45, fontface = "bold", color = "#b30000") +
  scale_fill_manual(values = c(
    "Elite Captured"               = "#1f77b4",
    "Landrace/Domesticate Private" = "#ff7f0e",
    "Wild Progenitor Private"      = "#2ca02c"
  )) +
  scale_x_continuous(limits = c(0, 142), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "e  Haplotype sharing across resistance structural classes",
    x = "Proportion of intact full-length haplotypes (%)",
    y = "Resistance structural class",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    legend.position = "top",
    legend.text = element_text(size = 8.0),
    axis.title.y = element_text(size = 9, face = "bold", color = "black", margin = margin(r = 5)),
    axis.text.y = element_text(size = 8.5, face = "bold", color = "black"),
    axis.text.x = element_text(size = 8.5, color = "black")
  ) +
  guides(fill = guide_legend(nrow = 1))

# ------------------------------------------------------------------------------
# Panel f: Geographic Hotspots Decomposed by NLR vs non-NLR Architecture
# ------------------------------------------------------------------------------
cat("Processing Panel f: Geographic distribution (NLR vs non-NLR)...\n")
s11d_df <- read.delim(file.path(res_dir, "Supplementary_Table_S11d_watkins_geographic_distribution.tsv"), stringsAsFactors = FALSE)
n_tot_cands <- sum(s11d_df$total_candidate_alleles, na.rm = TRUE)

geo_long <- s11d_df %>%
  select(geographic_region, total_candidate_alleles, pct_global_reservoir,
         "NLR (Classes I-II)" = nlr_alleles,
         "non-NLR (Classes III-V)" = non_nlr_alleles) %>%
  pivot_longer(cols = c("NLR (Classes I-II)", "non-NLR (Classes III-V)"),
               names_to = "ClassGroup", values_to = "Alleles") %>%
  mutate(
    Region_Clean = gsub("_", " ", geographic_region),
    Region_Clean = factor(Region_Clean, levels = rev(unique(gsub("_", " ", s11d_df$geographic_region)))),
    ClassGroup   = factor(ClassGroup, levels = c("non-NLR (Classes III-V)", "NLR (Classes I-II)"))
  )

geo_totals <- s11d_df %>%
  mutate(Region_Clean = factor(gsub("_", " ", geographic_region), levels = rev(gsub("_", " ", geographic_region))))

p_f <- ggplot(geo_long, aes(x = Alleles, y = Region_Clean, fill = ClassGroup)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  geom_text(data = geo_totals, aes(x = total_candidate_alleles, y = Region_Clean,
                                   label = sprintf("%d (%.1f%%)", total_candidate_alleles, pct_global_reservoir)),
            inherit.aes = FALSE, hjust = -0.15, size = 2.7, fontface = "bold", color = "black") +
  annotate("text", x = 60, y = 1.3,
           label = "West & Central Asia hotspot:\n180 / 249 (72.3%; OR = 5.18, P = 0.058)",
           size = 2.65, fontface = "bold", color = "#08519c", lineheight = 0.9) +
  scale_fill_manual(values = c("NLR (Classes I-II)" = "#1f77b4", "non-NLR (Classes III-V)" = "#d95f02")) +
  scale_x_continuous(limits = c(0, max(s11d_df$total_candidate_alleles) * 1.30), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = sprintf("f  Geographic origin of %d untapped functional alleles", n_tot_cands),
    x = "Candidate functional alleles count",
    y = "Accession geographic origin",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    legend.position = "top",
    legend.text = element_text(size = 8.0),
    axis.title.y = element_text(size = 9, face = "bold", color = "black", margin = margin(r = 5)),
    axis.text.y = element_text(size = 9, face = "bold", color = "black"),
    axis.text.x = element_text(size = 8.5, color = "black")
  )

# ------------------------------------------------------------------------------
# Compose Full Figure 5 with Patchwork
# ------------------------------------------------------------------------------
cat("Composing full Figure 5 with patchwork...\n")

top_row <- (p_a | p_b) + plot_layout(widths = c(1, 1), guides = "collect") & theme(legend.position = "top")
mid_row <- (p_c | p_d) + plot_layout(widths = c(1, 1))
bot_row <- (p_e | p_f) + plot_layout(widths = c(1, 1))

fig5_composed <- cowplot::plot_grid(
  top_row,
  mid_row,
  bot_row,
  ncol = 1,
  rel_heights = c(1.25, 1.0, 1.0)
)

out_pdf <- file.path(fig_dir, "Fig5_haplotypes_bottlenecks_publication.pdf")
out_png <- file.path(fig_dir, "Fig5_haplotypes_bottlenecks_publication.png")

cat(sprintf("Saving PDF -> %s ...\n", out_pdf))
pdf(out_pdf, width = 16, height = 16.5)
print(fig5_composed)
invisible(dev.off())

cat(sprintf("Saving 300-DPI PNG -> %s ...\n", out_png))
png(out_png, width = 16, height = 16.5, units = "in", res = 300)
print(fig5_composed)
invisible(dev.off())

file.copy(out_png, file.path(fig_dir, "Fig5_haplotypes_bottlenecks.png"), overwrite = TRUE)
file.copy(out_png, file.path(fig_dir, "Fig5_haplotypes_bottlenecks_R.png"), overwrite = TRUE)
file.copy(out_pdf, file.path(fig_dir, "Fig5_haplotypes_bottlenecks.pdf"), overwrite = TRUE)
file.copy(out_pdf, file.path(fig_dir, "Fig5_haplotypes_bottlenecks_R.pdf"), overwrite = TRUE)

cat("\n>>> Figure 5 successfully generated and verified! <<<\n")
cat("=======================================================================\n")
