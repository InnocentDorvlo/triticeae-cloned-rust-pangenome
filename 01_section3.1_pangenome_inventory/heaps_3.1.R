#!/usr/bin/env Rscript
# heaps_3.1.R - pangenome openness (Heaps' law) and core-genome decay for
#               NLR (Class I-II) vs non-NLR (Class III-V) orthogroups.
# Input : results/og_presence_matrix.tsv, results/S1b_orthogroup_partition.tsv  (written by rust31.py stats)
# Output: results/S3d_heaps_law.tsv, results/S3d_heaps_curves.tsv
# Usage : Rscript heaps_3.1.R [resultsdir] [n_perm] [n_label_perm]
#
# NOTE: Models both cumulative pan-orthogroup discovery (N_pan = kappa * n^gamma)
#       and core-orthogroup decay (N_core) across n = 1..56 haploid subgenomes
#       using 1,000 random addition permutations, vegan::specaccum exact rarefaction,
#       micropan::heaps, and a 1,000-replicate label-permutation test on Delta_gamma.
suppressPackageStartupMessages({ library(micropan); library(vegan) })
.cmd_args <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", .cmd_args[grep("^--file=", .cmd_args)])
.root_dir <- if (length(.file_arg) > 0) dirname(dirname(normalizePath(.file_arg[1]))) else getwd()
args  <- commandArgs(trailingOnly = TRUE)
rd    <- if (length(args) >= 1) args[1] else file.path(.root_dir, "results")
nperm <- if (length(args) >= 2) as.integer(args[2]) else 1000L
nlab  <- if (length(args) >= 3) as.integer(args[3]) else 1000L
set.seed(20261024)

pm  <- read.delim(file.path(rd, "og_presence_matrix.tsv"), row.names = 1, check.names = FALSE)
ogt <- read.delim(file.path(rd, "S1b_orthogroup_partition.tsv"), check.names = FALSE)
ogt <- ogt[ogt$set == "Main" & ogt$group %in% c("NLR", "nonNLR"), ]
pm  <- as.matrix(pm[, ogt$orthogroup, drop = FALSE])

# cumulative pan (union) and core (intersection) orthogroups over random genome orders
cum_curves <- function(m, n) {
  g <- nrow(m)
  pan_out  <- matrix(NA_real_, n, g)
  core_out <- matrix(NA_real_, n, g)
  for (i in seq_len(n)) {
    o <- sample(g)
    seen_pan  <- rep(FALSE, ncol(m))
    seen_core <- rep(TRUE,  ncol(m))
    for (j in seq_len(g)) {
      row_pos       <- (m[o[j], ] > 0)
      seen_pan      <- seen_pan  | row_pos
      seen_core     <- seen_core & row_pos
      pan_out[i, j]  <- sum(seen_pan)
      core_out[i, j] <- sum(seen_core)
    }
  }
  list(pan = pan_out, core = core_out)
}

fit_gamma <- function(y) {
  x <- seq_along(y)
  if (diff(range(y)) == 0) return(c(k = y[1], gamma = 0, gamma_se = 0))
  f <- try(nls(y ~ k * x^g, start = list(k = max(1, y[1]), g = 0.1),
               control = nls.control(maxiter = 200, warnOnly = TRUE)), silent = TRUE)
  if (inherits(f, "try-error")) return(c(k = NA, gamma = NA, gamma_se = NA))
  s <- tryCatch(summary(f)$coefficients, error = function(e) NULL)
  if (is.null(s)) return(c(k = coef(f)[["k"]], gamma = coef(f)[["g"]], gamma_se = NA))
  c(k = s["k", 1], gamma = s["g", 1], gamma_se = s["g", 2])
}

res <- list(); curves <- list()
for (grp in c("NLR", "nonNLR")) {
  m  <- pm[, ogt$group == grp, drop = FALSE]
  cc <- cum_curves(m, nperm)
  sp <- suppressWarnings(vegan::specaccum(m, method = "exact"))
  per <- t(apply(cc$pan, 1, fit_gamma))
  mean_fit <- fit_gamma(colMeans(cc$pan))
  hp <- tryCatch(micropan::heaps(m, n.perm = nperm), error = function(e) c(Intercept = NA, alpha = NA))
  res[[grp]] <- data.frame(
    group              = grp,
    n_orthogroups      = ncol(m),
    n_subgenomes       = nrow(m),
    kappa_meancurve    = round(unname(mean_fit["k"]), 4),
    gamma_meancurve    = unname(mean_fit["gamma"]),
    gamma_meancurve_se = unname(mean_fit["gamma_se"]),
    gamma_perm_mean    = mean(per[, "gamma"], na.rm = TRUE),
    gamma_perm_sd      = sd(per[, "gamma"], na.rm = TRUE),
    micropan_alpha     = unname(hp["alpha"]),
    micropan_intercept = unname(hp["Intercept"]),
    core_orthogroups_n56 = sum(colSums(m > 0) == nrow(m))
  )
  curves[[grp]] <- data.frame(
    group           = grp,
    n_genomes       = seq_len(ncol(cc$pan)),
    mean            = colMeans(cc$pan),
    sd              = apply(cc$pan, 2, sd),
    specaccum_exact = as.numeric(sp$richness),
    specaccum_sd    = as.numeric(sp$sd),
    core_mean       = colMeans(cc$core),
    core_sd         = apply(cc$core, 2, sd)
  )
}
obs <- res$NLR$gamma_meancurve - res$nonNLR$gamma_meancurve

# label permutation test: shuffle NLR/nonNLR labels among orthogroups
lab <- ogt$group; null <- numeric(nlab)
for (i in seq_len(nlab)) {
  l <- sample(lab)
  a <- fit_gamma(colMeans(cum_curves(pm[, l == "NLR", drop = FALSE], 100)$pan))["gamma"]
  b <- fit_gamma(colMeans(cum_curves(pm[, l == "nonNLR", drop = FALSE], 100)$pan))["gamma"]
  null[i] <- a - b
}
p <- (sum(abs(null) >= abs(obs), na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)
out <- do.call(rbind, res)
out$delta_gamma         <- obs
out$P_label_permutation <- p
out$n_label_perm        <- nlab
write.table(out, file.path(rd, "S3d_heaps_law.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, curves), file.path(rd, "S3d_heaps_curves.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(out)
