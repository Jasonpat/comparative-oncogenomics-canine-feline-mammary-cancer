################################################################################
# CONFIGURATION & GLOBAL SETUP
# Comparative oncogenomics of canine and feline mammary cancer reveals a conserved
# and tractable oncogenic core with human breast cancer concordance
################################################################################

# ==============================================================================
# PROJECT ROOT
# ==============================================================================

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

# If the working directory is accidentally set to the R/ folder,
# move one level up to the actual project root.
if (basename(PROJECT_ROOT) == "R") {
  PROJECT_ROOT <- dirname(PROJECT_ROOT)
}

assign("PROJECT_ROOT", PROJECT_ROOT, envir = .GlobalEnv)

.path <- function(...) {
  normalizePath(file.path(PROJECT_ROOT, ...), winslash = "/", mustWork = FALSE)
}
# ==============================================================================
# PACKAGE GROUPS
# ==============================================================================
# Core packages are required for the main manuscript pipeline.
CORE_PACKAGES <- c(
  "limma",
  "biomaRt",
  "fgsea",
  "msigdbr",
  "clusterProfiler",
  "BiocParallel",
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

# Optional packages are checked only when the corresponding optional module is enabled.
OPTIONAL_PACKAGES_NETWORK <- c("igraph", "STRINGdb")
OPTIONAL_PACKAGES_TCGA <- c("TCGAbiolinks", "edgeR", "ggrepel")
OPTIONAL_PACKAGES_TCGA_SIGNATURE <- c("TCGAbiolinks", "edgeR", "GSVA", "survival", "survminer")

check_required_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
      "\nInstall CRAN packages with install.packages(...), and Bioconductor packages with BiocManager::install(...).",
      "\nA helper installer is provided in 00_install_packages.R.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Check core packages before loading them, so missing dependencies produce a clear message.
check_required_packages(CORE_PACKAGES)

suppressPackageStartupMessages({
  library(limma)
  library(biomaRt)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(BiocParallel)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(httr)
  library(jsonlite)
  library(purrr)
  library(ggplot2)
  library(tibble)
})

# ==============================================================================
# PATH CONFIGURATION
# ==============================================================================
PATH_CONFIG <- list(
  project_root = PROJECT_ROOT,
  data_dir     = .path("data"),
  dog_expr     = .path("data", "GSE119810_CMT_222S_FPKM.csv"),
  cat_up       = .path("data", "CAT_up.txt"),
  cat_down     = .path("data", "CAT_down.txt"),
  tcga_dir     = .path("data", "TCGA_BRCA_STAR_Counts"),
  results_dir  = .path("results"),
  cache_dir    = .path("cache"),
  fig_dir      = .path("figures"),
  create_dirs  = TRUE
)

# ==============================================================================
# PIPELINE EXECUTION CONTROL
# ==============================================================================
EXECUTION_CONFIG <- list(
  # Core analysis steps
  run_step1_dog_gsea = TRUE,
  run_step2_cat_hallmarks = TRUE,
  run_step3_conserved_hallmarks = TRUE,
  run_step4_conserved_targets = TRUE,
  run_step5_drug_targets = TRUE,

  # Step 5B is part of the main workflow because Step 6/Figure 3 are intended to
  # use the integrated drug-level/probe-linked indication layer.
  run_step5b_probe_indications = TRUE,

  run_step6_prioritisation = TRUE,

  # Optional sensitivity analysis comparing primary target-level evidence with
  # the integrated Step 5B evidence layer.
  run_sensitivity_prioritisation = TRUE,

  run_step7_outputs = TRUE,

  # Figure/Table generation
  run_step8_figure1 = TRUE,
  run_step9_figure2_table1 = TRUE,
  run_step10_figure3 = TRUE,

  # Optional validation layers
  run_validation_background_light = FALSE,
  run_validation_permutation = FALSE,
  run_validation_network = FALSE,
  run_validation_tcga_brca = FALSE,

  # Optional TCGA-BRCA patient-level validation layer
  # Runs ssGSEA scoring of the strict conserved UP-core, PAM50 subtype
  # association, and Kaplan-Meier OS visualisation. Cox models are not run.
  run_validation_tcga_signature = FALSE
)

# ==============================================================================
# ANALYSIS PARAMETERS
# ==============================================================================
PARAM_CONFIG <- list(
  # Annotation parameters
  ensembl_version = 104,
  ortholog_type_filter = "ortholog_one2one",
  ortholog_conf_filter = 1,

  # Dog differential-expression / GSEA parameters
  # Main manuscript mode: all_blockaware.
  # Sensitivity mode: complete_pairs_only.
  dog_analysis_mode = "all_blockaware",
  dog_gsea_min_size = 15,
  dog_gsea_max_size = 500,
  dog_hallmark_padj = 0.05,

  # Cat ORA parameters
  cat_ora_padj = 0.05,
  conservation_padj = 0.05,

  # Drug evidence parameters
  drug_min_phase = 1,
  drug_keep_approved_only = FALSE,
  drug_keep_na_phase = FALSE,

  # Step 5B probe-linked/drug-level indication expansion
  probe_min_phase = 1,
  probe_overwrite_target_master = TRUE,

  # Target prioritisation
  top_n_targets = 15,
  use_string_network = TRUE,
  canonical_mitosis_filter = TRUE,

  # Permutation validation
  permutation_n = 10000,

  # Reproducibility-sensitive stochastic modules
  random_seed = 1
)

# ==============================================================================
# REPRODUCIBILITY & METADATA
# ==============================================================================
REPRODUCIBILITY_CONFIG <- list(
  analysis_name = "Comparative oncogenomics of canine and feline mammary cancer reveals a conserved and tractable oncogenic core with human breast cancer concordance",
  analysis_version = "2.0",
  publication_ready = TRUE,
  q1_submission = TRUE,

  # Reference genomes
  reference_genome_human = "Ensembl 104 (GRCh38.p14)",
  reference_genome_dog = "Ensembl 104 (ROS_Cfam_1.0)",
  reference_genome_cat = "Ensembl 104 (Felis_catus_9.0)",

  # Required R version
  min_r_version = "4.0",

  # Required and optional packages
  required_packages = CORE_PACKAGES,
  optional_packages_network = OPTIONAL_PACKAGES_NETWORK,
  optional_packages_tcga = OPTIONAL_PACKAGES_TCGA,
  optional_packages_tcga_signature = OPTIONAL_PACKAGES_TCGA_SIGNATURE,

  # Contact and citation
  authors = "Iason-Spyridon Patergiannakis; Ioannis S. Pappas",
  contact_email = "ipatergiannakis@vet.uth.gr",
  doi_or_preprint = "pending",
  github_url = "pending",
  license = "MIT"
)

# ==============================================================================
# QUALITY CONTROL THRESHOLDS
# ==============================================================================
QC_CONFIG <- list(
  # Input file QC
  check_input_files = TRUE,

  # Mapping QC
  min_dog_genes_mapped = 500,
  min_cat_genes_mapped = 500,
  warn_symbol_collisions_pct = 5,

  # Conservation QC
  min_conserved_hallmarks = 1,
  min_core_genes = 5,

  # Drug evidence QC
  warn_ot_failed_pct = 10,
  min_drug_targets = 5,

  # Output QC
  expect_top_targets = 15,
  min_shortlist_size = 10,

  # Fail-hard vs warn-soft
  fail_on_missing_conserved = TRUE,
  fail_on_no_core_genes = TRUE,
  warn_on_high_ot_failure = TRUE
)

# ==============================================================================
# ANALYSIS COMPONENTS METADATA
# ==============================================================================
COMPONENT_DESCRIPTIONS <- list(
  step1 = "Dog transcriptome analysis and Hallmark GSEA",
  step2 = "Cat gene list parsing, mapping to human orthologs, Hallmark ORA with fixed universe",
  step3 = "Intersection of significant dog and cat hallmarks by direction",
  step4 = "Extraction of conserved core genes from hallmark leading edges",
  step5 = "Druggability assessment via Open Targets and ChEMBL APIs",
  step5b = "Drug-level and high-quality probe-linked indication expansion",
  step6 = "Multi-component target prioritisation and ranking",
  step6b = "Sensitivity analysis: primary target-level evidence versus integrated Step 5B evidence",
  step7 = "Generation of final summary tables and reports",
  step8 = "Figure 1 panels: conserved Hallmark programs",
  step9 = "Figure 2 and Table 1: gene-level overlap support metrics",
  step10 = "Figure 3: target prioritisation and translational landscape",
  tcga_concordance = "Optional TCGA-BRCA tumour-vs-normal concordance validation and Figure 4",
  tcga_signature = "Optional TCGA-BRCA patient-level ssGSEA, PAM50 subtype association, and Kaplan-Meier OS visualisation without Cox models"
)

# ==============================================================================
# HELPERS
# ==============================================================================
`%s*%` <- function(x, n) paste(rep(x, n), collapse = "")

init_directories <- function() {
  if (isTRUE(PATH_CONFIG$create_dirs)) {
    dir.create(PATH_CONFIG$data_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(PATH_CONFIG$results_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(PATH_CONFIG$cache_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(PATH_CONFIG$fig_dir, showWarnings = FALSE, recursive = TRUE)
  }
}

save_result <- function(obj, filename, format = "csv") {
  outfile <- file.path(PATH_CONFIG$results_dir, filename)
  dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

  if (format == "csv") {
    readr::write_csv(obj, outfile)
  } else if (format == "csv_rownames") {
    if (is.data.frame(obj) && !is.null(rownames(obj))) {
      obj <- tibble::rownames_to_column(as.data.frame(obj), "rowname")
    }
    readr::write_csv(obj, outfile)
  } else if (format == "rds") {
    saveRDS(obj, outfile)
  } else if (format == "txt") {
    writeLines(as.character(obj), outfile)
  } else {
    stop("Unsupported format: ", format, call. = FALSE)
  }

  cat(" ✓ Saved:", normalizePath(outfile, winslash = "/", mustWork = FALSE), "\n")
  invisible(outfile)
}

load_cached <- function(name, compute_fn) {
  cache_file <- file.path(PATH_CONFIG$cache_dir, paste0(name, ".rds"))

  if (file.exists(cache_file)) {
    cat("  ↻ Loading from cache:", name, "\n")
    return(readRDS(cache_file))
  }

  cat("  ⧗ Computing:", name, "\n")
  result <- compute_fn()

  if (!is.null(result)) {
    dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
    saveRDS(result, cache_file)
    cat(" ✓ Cached:", name, "\n")
  }

  result
}

check_input_files <- function() {
  if (!isTRUE(QC_CONFIG$check_input_files)) {
    return(invisible(NULL))
  }

  req <- c(PATH_CONFIG$dog_expr, PATH_CONFIG$cat_up, PATH_CONFIG$cat_down)
  missing <- req[!file.exists(req)]

  if (length(missing) > 0) {
    stop("Missing required input file(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  cat("  ✓ All input files found\n")
  invisible(TRUE)
}

check_packages <- function() {
  check_required_packages(REPRODUCIBILITY_CONFIG$required_packages)
}

print_pipeline_header <- function() {
  cat("\n╔════════════════════════════════════════════════════════════╗\n")
  cat("║  CANINE-FELINE COMPARATIVE ONCOGENOMICS PIPELINE          ║\n")
  cat("║  Structured • Modular • Reproducible                     ║\n")
  cat("║  Q1 Publication Ready                                    ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n\n")
}

print_config <- function() {
  cat("EXECUTION CONTROL:\n")
  for (step in names(EXECUTION_CONFIG)) {
    status <- ifelse(isTRUE(EXECUTION_CONFIG[[step]]), "✓", "✗")
    step_clean <- gsub("run_", "", step)
    cat("  ", status, " ", step_clean, "\n", sep = "")
  }

  cat("\nREPRODUCIBILITY:\n")
  cat("  Analysis:", REPRODUCIBILITY_CONFIG$analysis_name, "\n")
  cat("  Version:", REPRODUCIBILITY_CONFIG$analysis_version, "\n")
  cat("  Ensembl:", PARAM_CONFIG$ensembl_version, "\n")
  cat("  Dog analysis mode:", PARAM_CONFIG$dog_analysis_mode, "\n")
  cat("  Step 5B overwrites target_master:", PARAM_CONFIG$probe_overwrite_target_master, "\n")
  cat("  Q1 submission ready:", REPRODUCIBILITY_CONFIG$q1_submission, "\n")

  cat("\nQC THRESHOLDS:\n")
  cat("  Min conserved hallmarks:", QC_CONFIG$min_conserved_hallmarks, "\n")
  cat("  Min core genes:", QC_CONFIG$min_core_genes, "\n")
  cat("  Min drug targets:", QC_CONFIG$min_drug_targets, "\n")
  cat("  Warn if Open Targets failed >", QC_CONFIG$warn_ot_failed_pct, "%\n")

  cat("\nDIRECTORIES:\n")
  cat("  Results:", normalizePath(PATH_CONFIG$results_dir, winslash = "/", mustWork = FALSE), "\n")
  cat("  Cache:", normalizePath(PATH_CONFIG$cache_dir, winslash = "/", mustWork = FALSE), "\n")
  cat("  Figures:", normalizePath(PATH_CONFIG$fig_dir, winslash = "/", mustWork = FALSE), "\n\n")
}

# ==============================================================================
# GLOBAL INITIALIZATION
# ==============================================================================
check_packages()
init_directories()
check_input_files()

if (isTRUE(EXECUTION_CONFIG$run_validation_tcga_brca)) {
  check_required_packages(OPTIONAL_PACKAGES_TCGA)
}

if (isTRUE(EXECUTION_CONFIG$run_validation_tcga_signature)) {
  check_required_packages(OPTIONAL_PACKAGES_TCGA_SIGNATURE)
}

if (isTRUE(EXECUTION_CONFIG$run_validation_network)) {
  check_required_packages(OPTIONAL_PACKAGES_NETWORK)
}

print_pipeline_header()
print_config()
set.seed(PARAM_CONFIG$random_seed)

cat("✓ Configuration loaded successfully\n\n")
