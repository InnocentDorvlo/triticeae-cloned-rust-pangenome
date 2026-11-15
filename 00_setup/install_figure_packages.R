#!/usr/bin/env Rscript
# ==============================================================================
# Script: install_figure_packages.R
# Purpose: Automatically verify and install all R and Bioconductor packages
#          required for publication-grade Figures 1-5.
# ==============================================================================

cat("=======================================================================\n")
cat("      Wheat Resistance Pangenome - Figure Package Installer\n")
cat("=======================================================================\n\n")

# Set CRAN repository mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 1. CRAN Packages required
cran_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "ggrepel",
  "patchwork",
  "cowplot",
  "scales",
  "ape",
  "igraph",
  "pegas",
  "ggtern",
  "ggnewscale",
  "png",
  "gridExtra",
  "BiocManager"
)

# 2. Bioconductor Packages required
bioc_packages <- c(
  "ComplexHeatmap",
  "ggtree"
)

cat("Step 1: Checking and installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  -> Installing CRAN package: %s ...\n", pkg))
    tryCatch({
      install.packages(pkg, dependencies = TRUE, quiet = FALSE)
    }, error = function(e) {
      cat(sprintf("     [ERROR] Failed to install %s: %s\n", pkg, e$message))
    })
  } else {
    cat(sprintf("  [OK] %s is already installed.\n", pkg))
  }
}

cat("\nStep 2: Checking and installing Bioconductor packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = FALSE)
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  -> Installing Bioconductor package: %s ...\n", pkg))
    tryCatch({
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }, error = function(e) {
      cat(sprintf("     [ERROR] Failed to install %s: %s\n", pkg, e$message))
    })
  } else {
    cat(sprintf("  [OK] %s is already installed.\n", pkg))
  }
}

cat("\n=======================================================================\n")
cat("                    Verification & Status Report\n")
cat("=======================================================================\n")

all_packages <- c(cran_packages, bioc_packages)
all_ok <- TRUE

cat(sprintf("%-20s %-15s %s\n", "Package", "Source", "Status"))
cat(paste(rep("-", 50), collapse = ""), "\n")

for (pkg in all_packages) {
  src <- if (pkg %in% bioc_packages) "Bioconductor" else "CRAN"
  is_installed <- requireNamespace(pkg, quietly = TRUE)
  ver <- if (is_installed) as.character(packageVersion(pkg)) else "NOT INSTALLED"
  status_str <- if (is_installed) paste0("SUCCESS (v", ver, ")") else "FAILED"
  
  if (!is_installed) all_ok <- FALSE
  cat(sprintf("%-20s %-15s %s\n", pkg, src, status_str))
}

cat(paste(rep("-", 50), collapse = ""), "\n")
if (all_ok) {
  cat("\n>>> ALL PACKAGES SUCCESSFULLY INSTALLED AND READY! <<<\n")
} else {
  cat("\n>>> SOME PACKAGES FAILED TO INSTALL. PLEASE REVIEW OUTPUT ABOVE. <<<\n")
}
cat("=======================================================================\n")
