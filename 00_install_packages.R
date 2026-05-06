################################################################################
# 00_install_packages.R
# Helper script for installing packages required by the comparative oncogenomics
# pipeline.
################################################################################

core_cran_pkgs <- c(
  "dplyr",
  "readr",
  "tidyr",
  "stringr",
  "ggplot2",
  "httr",
  "jsonlite",
  "purrr",
  "tibble"
)

core_bioc_pkgs <- c(
  "limma",
  "biomaRt",
  "fgsea",
  "msigdbr",
  "clusterProfiler",
  "BiocParallel"
)

optional_cran_pkgs <- c(
  "igraph",
  "ggrepel",
  "survival",
  "survminer"
)

optional_bioc_pkgs <- c(
  "STRINGdb",
  "edgeR",
  "TCGAbiolinks",
  "GSVA"
)

install_if_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    install.packages(missing)
  }
  invisible(missing)
}

install_if_missing_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
  invisible(missing)
}

install_core_packages <- function() {
  install_if_missing_cran(core_cran_pkgs)
  install_if_missing_bioc(core_bioc_pkgs)
  cat("Core packages installed/available.\n")
  invisible(TRUE)
}

install_optional_packages <- function() {
  install_if_missing_cran(optional_cran_pkgs)
  install_if_missing_bioc(optional_bioc_pkgs)
  cat("Optional network/TCGA/patient-level validation packages installed/available.\n")
  invisible(TRUE)
}

install_all_packages <- function() {
  install_core_packages()
  install_optional_packages()
  cat("All core and optional packages installed/available.\n")
  invisible(TRUE)
}

# Default behaviour when sourcing this file: install the main pipeline dependencies.
install_core_packages()

cat("\nTo install optional network, TCGA and patient-level validation dependencies, run:\n")
cat("  install_optional_packages()\n")
cat("\nTo install everything at once, run:\n")
cat("  install_all_packages()\n")
