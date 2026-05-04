################################################################################
# 00_run_pipeline.R
# Standalone orchestrator for the canine-feline comparative oncogenomics pipeline.
################################################################################

# ==============================================================================
# PROJECT ROOT DETECTION
# ==============================================================================
.detect_this_file <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  }

  for (i in rev(seq_along(sys.frames()))) {
    ofile <- tryCatch(sys.frames()[[i]]$ofile, error = function(e) NULL)
    if (!is.null(ofile) && nzchar(ofile)) {
      return(normalizePath(ofile, winslash = "/", mustWork = FALSE))
    }
  }

  normalizePath("00_run_pipeline.R", winslash = "/", mustWork = FALSE)
}

THIS_FILE <- .detect_this_file()
PROJECT_ROOT <- dirname(THIS_FILE)
assign("PROJECT_ROOT", PROJECT_ROOT, envir = .GlobalEnv)

source_pipeline_file <- function(filename) {
  candidate <- file.path(PROJECT_ROOT, "R", filename)
  if (!file.exists(candidate)) {
    stop("Cannot find pipeline file: ", candidate, call. = FALSE)
  }
  source(candidate, local = .GlobalEnv)
}

# ==============================================================================
# SOURCE CONFIGURATION & CORE MODULES
# ==============================================================================
source_pipeline_file("00_config.R")
source_pipeline_file("01_module_dog_gsea.R")
source_pipeline_file("02_module_cat_hallmarks.R")
source_pipeline_file("03_module_conserved_hallmarks.R")
source_pipeline_file("04_module_conserved_targets.R")
source_pipeline_file("05_module_drug_targets.R")
source_pipeline_file("05b_module_probe_linked_indications.R")
source_pipeline_file("06_module_target_prioritisation.R")
source_pipeline_file("07_module_outputs.R")
source_pipeline_file("08_figure1_conserved_hallmarks.R")
source_pipeline_file("09_figure2_overlap_support.R")
source_pipeline_file("10_figure3_target_prioritisation.R")

# ==============================================================================
# SOURCE OPTIONAL MODULES WHEN ENABLED
# ==============================================================================
if (isTRUE(EXECUTION_CONFIG$run_sensitivity_prioritisation)) {
  source_pipeline_file("06b_prioritisation_sensitivity_analysis.R")
}

if (isTRUE(EXECUTION_CONFIG$run_validation_background_light) ||
    isTRUE(EXECUTION_CONFIG$run_validation_permutation)) {
  source_pipeline_file("11a_build_background_tractability_universe_light.R")
}

if (isTRUE(EXECUTION_CONFIG$run_validation_permutation)) {
  source_pipeline_file("11b_permutation_validation_light.R")
}

if (isTRUE(EXECUTION_CONFIG$run_validation_network)) {
  check_required_packages(OPTIONAL_PACKAGES_NETWORK)
  source_pipeline_file("13_network_validation_prioritised_targets.R")
}

if (isTRUE(EXECUTION_CONFIG$run_validation_tcga_brca)) {
  check_required_packages(OPTIONAL_PACKAGES_TCGA)
  source_pipeline_file("14_module_tcga_brca_external_concordance.R")
  source_pipeline_file("15_figure4_tcga_brca_concordance.R")
}

# ==============================================================================
# REQUIRED FUNCTION CHECK
# ==============================================================================
required_funs <- c(
  "module_dog_gsea",
  "module_cat_hallmarks",
  "module_conserved_hallmarks",
  "module_conserved_targets",
  "module_drug_targets",
  "module_probe_linked_indications",
  "module_target_prioritisation",
  "module_outputs"
)

missing_funs <- required_funs[!vapply(required_funs, exists, logical(1), mode = "function")]
if (length(missing_funs) > 0) {
  stop("Missing required module function(s): ", paste(missing_funs, collapse = ", "), call. = FALSE)
}

# ==============================================================================
# STEP RUNNERS
# ==============================================================================
run_step <- function(step_name, step_title, fun) {
  config_key <- paste0("run_", step_name)

  if ((config_key %in% names(EXECUTION_CONFIG)) && !isTRUE(EXECUTION_CONFIG[[config_key]])) {
    cat(strrep("═", 70), "\n", sep = "")
    cat(step_title, "\n")
    cat(strrep("═", 70), "\n", sep = "")
    cat("  ⊘ SKIPPED\n\n")
    return(invisible(NULL))
  }

  cat(strrep("═", 70), "\n", sep = "")
  cat(step_title, "\n")
  cat(strrep("═", 70), "\n\n", sep = "")

  tryCatch(
    fun(),
    error = function(e) {
      cat("\n✗ PIPELINE FAILED AT: ", step_title, "\n", sep = "")
      cat("Error message:\n  ", conditionMessage(e), "\n\n", sep = "")
      stop(e)
    }
  )
}

run_optional_step <- function(enabled, step_title, fun) {
  if (!isTRUE(enabled)) return(invisible(NULL))

  cat(strrep("═", 70), "\n", sep = "")
  cat(step_title, "\n")
  cat(strrep("═", 70), "\n\n", sep = "")

  tryCatch(
    fun(),
    error = function(e) {
      cat("\n✗ OPTIONAL STEP FAILED AT: ", step_title, "\n", sep = "")
      cat("Error message:\n  ", conditionMessage(e), "\n\n", sep = "")
      stop(e)
    }
  )
}

# ==============================================================================
# POST-PIPELINE QC
# ==============================================================================
run_qc_checks <- function() {
  cat("\n", strrep("═", 70), "\n", sep = "")
  cat("POST-PIPELINE QUALITY CONTROL CHECKS\n")
  cat(strrep("═", 70), "\n\n", sep = "")

  issues <- list(errors = character(), warnings = character())

  if (exists("cons_up", envir = .GlobalEnv) && exists("cons_down", envir = .GlobalEnv)) {
    n_hallmarks <- length(get("cons_up", envir = .GlobalEnv)) + length(get("cons_down", envir = .GlobalEnv))
    if (n_hallmarks < QC_CONFIG$min_conserved_hallmarks) {
      msg <- paste0("Fewer conserved hallmarks than expected: ", n_hallmarks)
      if (isTRUE(QC_CONFIG$fail_on_missing_conserved)) issues$errors <- c(issues$errors, msg) else issues$warnings <- c(issues$warnings, msg)
    } else {
      cat("  ✓ Conserved hallmarks:", n_hallmarks, "\n")
    }
  } else {
    issues$warnings <- c(issues$warnings, "Conserved hallmark objects not found")
  }

  if (exists("core_up_genes_u", envir = .GlobalEnv)) {
    n_core <- length(get("core_up_genes_u", envir = .GlobalEnv))
    if (n_core < QC_CONFIG$min_core_genes) {
      msg <- paste0("Fewer core UP genes than expected: ", n_core)
      if (isTRUE(QC_CONFIG$fail_on_no_core_genes)) issues$errors <- c(issues$errors, msg) else issues$warnings <- c(issues$warnings, msg)
    } else {
      cat("  ✓ Core UP genes:", n_core, "\n")
    }
  } else {
    issues$warnings <- c(issues$warnings, "core_up_genes_u not found")
  }

  if (exists("target_ok", envir = .GlobalEnv) && exists("target_fail", envir = .GlobalEnv)) {
    target_ok <- get("target_ok", envir = .GlobalEnv)
    target_fail <- get("target_fail", envir = .GlobalEnv)
    total_ot <- length(target_ok) + length(target_fail)
    if (total_ot > 0) {
      fail_pct <- 100 * length(target_fail) / total_ot
      if (fail_pct > QC_CONFIG$warn_ot_failed_pct) {
        issues$warnings <- c(issues$warnings, paste0("High Open Targets failure rate: ", round(fail_pct, 1), "%"))
      } else {
        cat("  ✓ Open Targets query success rate:", round(100 - fail_pct, 1), "%\n")
      }
    }
  }

  if (exists("target_deep_indication_summary", envir = .GlobalEnv)) {
    cat("  ✓ Step 5B integrated indication summary available\n")
  } else if (isTRUE(EXECUTION_CONFIG$run_step5b_probe_indications)) {
    issues$warnings <- c(issues$warnings, "Step 5B was enabled but target_deep_indication_summary was not found")
  }

  if (exists("shortlist_overall", envir = .GlobalEnv)) {
    shortlist_overall <- get("shortlist_overall", envir = .GlobalEnv)
    if (nrow(shortlist_overall) < QC_CONFIG$min_shortlist_size) {
      issues$warnings <- c(issues$warnings, paste0("Small shortlist: ", nrow(shortlist_overall), " targets"))
    } else {
      cat("  ✓ Shortlist size:", nrow(shortlist_overall), "targets\n")
    }
  }

  if (exists("target_master", envir = .GlobalEnv)) {
    target_master <- get("target_master", envir = .GlobalEnv)
    if (nrow(target_master) < QC_CONFIG$min_drug_targets) {
      issues$warnings <- c(issues$warnings, paste0("Few target-master rows: ", nrow(target_master)))
    } else {
      cat("  ✓ Target-master rows:", nrow(target_master), "\n")
    }
  }

  cat("\n")
  if (length(issues$errors) > 0) {
    cat("ERRORS:\n")
    for (err in issues$errors) cat("  ✗ ", err, "\n", sep = "")
  }
  if (length(issues$warnings) > 0) {
    cat("WARNINGS:\n")
    for (warn in issues$warnings) cat("  ⚠ ", warn, "\n", sep = "")
  }
  if (length(issues$errors) == 0 && length(issues$warnings) == 0) {
    cat("✓ All QC checks passed with no issues\n")
  }
  cat("\n")

  if (length(issues$errors) > 0) stop("Fatal QC issue(s) detected.", call. = FALSE)
  invisible(issues)
}

# ==============================================================================
# RUN PIPELINE
# ==============================================================================
pipeline_start <- Sys.time()
cat("Pipeline started at:", as.character(pipeline_start), "\n\n")

run_step("step1_dog_gsea", "STEP 1: DOG ANALYSIS, BLOCK-AWARE LIMMA & HALLMARK GSEA", module_dog_gsea)
run_step("step2_cat_hallmarks", "STEP 2: CAT GENE LISTS, HUMAN MAPPING & HALLMARK ORA", module_cat_hallmarks)
run_step("step3_conserved_hallmarks", "STEP 3: IDENTIFY DIRECTIONALLY CONSERVED HALLMARKS", module_conserved_hallmarks)
run_step("step4_conserved_targets", "STEP 4: EXTRACT STRICT CONSERVED CORE TARGETS + OVERLAP SUPPORT", module_conserved_targets)
run_step("step5_drug_targets", "STEP 5: DRUGGABILITY ASSESSMENT (OPEN TARGETS + CHEMBL)", module_drug_targets)
run_step(
  "step5b_probe_indications",
  "STEP 5B: DRUG-LEVEL AND PROBE-LINKED INDICATION EXPANSION",
  function() {
    module_probe_linked_indications(
      overwrite_target_master = isTRUE(PARAM_CONFIG$probe_overwrite_target_master),
      min_probe_phase = PARAM_CONFIG$probe_min_phase
    )
  }
)
run_step("step6_prioritisation", "STEP 6: CONSENSUS TARGET PRIORITISATION", module_target_prioritisation)

if (isTRUE(EXECUTION_CONFIG$run_sensitivity_prioritisation)) {
  run_optional_step(
    TRUE,
    "STEP 6B: PRIORITISATION SENSITIVITY ANALYSIS",
    module_prioritisation_sensitivity_primary_vs_integrated
  )
}

run_step("step7_outputs", "STEP 7: SUMMARY OUTPUTS AND TRANSLATIONAL GAP TABLES", module_outputs)

if (exists("module_figure1_hallmarks", mode = "function")) {
  run_step("step8_figure1", "STEP 8: FIGURE 1 - CONSERVED HALLMARK PROGRAMS", module_figure1_hallmarks)
}
if (exists("module_figure2_overlap_support", mode = "function")) {
  run_step("step9_figure2_table1", "STEP 9: FIGURE 2 + TABLE 1 - GENE-LEVEL OVERLAP SUPPORT", module_figure2_overlap_support)
}
if (exists("module_figure3_target_prioritisation", mode = "function")) {
  run_step("step10_figure3", "STEP 10: FIGURE 3 - TARGET PRIORITISATION AND TRANSLATIONAL LANDSCAPE", module_figure3_target_prioritisation)
}

if (isTRUE(EXECUTION_CONFIG$run_validation_background_light)) {
  run_optional_step(
    TRUE,
    "VALIDATION: BUILD HALLMARK-CONSTRAINED BACKGROUND TRACTABILITY UNIVERSE",
    module_build_background_tractability_universe_light
  )
}

if (isTRUE(EXECUTION_CONFIG$run_validation_permutation)) {
  run_optional_step(
    TRUE,
    "VALIDATION: PERMUTATION VALIDATION AGAINST BACKGROUND UNIVERSE",
    function() module_permutation_validation_light(B = PARAM_CONFIG$permutation_n, seed = PARAM_CONFIG$random_seed)
  )
}

if (isTRUE(EXECUTION_CONFIG$run_validation_network)) {
  run_optional_step(
    TRUE,
    "VALIDATION: NETWORK VALIDATION OF PRIORITISED TARGETS",
    function() module_network_validation_prioritised_targets(
      top_n = PARAM_CONFIG$top_n_targets,
      n_permutations = PARAM_CONFIG$permutation_n,
      seed = PARAM_CONFIG$random_seed,
      save_pdf = FALSE
    )
  )
}

if (isTRUE(EXECUTION_CONFIG$run_validation_tcga_brca)) {
  run_optional_step(
    TRUE,
    "VALIDATION: TCGA-BRCA EXTERNAL HUMAN CONCORDANCE",
    function() {
      module_tcga_brca_external_concordance_local(
        tcga_dir = PATH_CONFIG$tcga_dir,
        use_cached_query = TRUE
      )
      figure4_tcga_brca_concordance()
    }
  )
}

qc_issues <- run_qc_checks()

# ==============================================================================
# SESSION LOGGING
# ==============================================================================
pipeline_end <- Sys.time()
pipeline_runtime <- difftime(pipeline_end, pipeline_start, units = "mins")

session_info_text <- c(
  "PIPELINE EXECUTION SUMMARY",
  "═════════════════════════════════════════════════════",
  paste("Analysis Name:", REPRODUCIBILITY_CONFIG$analysis_name),
  paste("Pipeline Version:", REPRODUCIBILITY_CONFIG$analysis_version),
  paste("Project Root:", PROJECT_ROOT),
  paste("Start Time:", as.character(pipeline_start)),
  paste("End Time:", as.character(pipeline_end)),
  paste("Runtime (minutes):", round(as.numeric(pipeline_runtime), 2)),
  "",
  "LIVE RESOURCE ACCESS NOTE",
  "═════════════════════════════════════════════════════",
  paste("Live-resource run timestamp:", as.character(pipeline_start)),
  "Open Targets, ChEMBL, STRINGdb, Ensembl BioMart and GDC/TCGA may update over time.",
  "For Open Targets and ChEMBL, retain cache/ plus OT_* query logs for reproducibility.",
  "",
  "ANALYSIS PARAMETERS",
  "═════════════════════════════════════════════════════",
  paste("Dog analysis mode:", PARAM_CONFIG$dog_analysis_mode),
  paste("Ensembl Version:", PARAM_CONFIG$ensembl_version),
  paste("Dog Hallmark FDR:", PARAM_CONFIG$dog_hallmark_padj),
  paste("Cat ORA FDR:", PARAM_CONFIG$cat_ora_padj),
  paste("Drug minimum phase:", PARAM_CONFIG$drug_min_phase),
  paste("Step 5B enabled:", EXECUTION_CONFIG$run_step5b_probe_indications),
  paste("Step 5B overwrote target_master:", PARAM_CONFIG$probe_overwrite_target_master),
  paste("Probe minimum phase:", PARAM_CONFIG$probe_min_phase),
  paste("Top N targets:", PARAM_CONFIG$top_n_targets),
  paste("Random seed:", PARAM_CONFIG$random_seed),
  "",
  "SESSION INFO",
  "═════════════════════════════════════════════════════"
)

writeLines(
  c(session_info_text, capture.output(sessionInfo())),
  file.path(PATH_CONFIG$results_dir, "REPRODUCIBILITY_RECORD.txt")
)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║  ✓ PIPELINE FINISHED                                     ║\n")
cat("║  Check results/ for tables, logs, and figure panels      ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

cat("KEY OUTPUTS:\n")
cat("  • results/DOG_sample_structure_QC.csv\n")
cat("  • results/", paste0("DOG_DE_limma_", PARAM_CONFIG$dog_analysis_mode, "_full.csv"), "\n", sep = "")
cat("  • results/", paste0("DOG_DE_limma_", PARAM_CONFIG$dog_analysis_mode, ".rds"), "\n", sep = "")
cat("  • results/", paste0("DOG_DE_limma_", PARAM_CONFIG$dog_analysis_mode, "_sig_FDR0.05_logFC1.csv"), "\n", sep = "")
cat("  • results/DOG_Hallmark_fgsea_simple_plus.csv\n")
cat("  • results/CAT_Hallmark_ORA_UP_universeFixed.csv\n")
cat("  • results/CAT_Hallmark_ORA_DOWN_universeFixed.csv\n")
cat("  • results/CONSERVED_HALLMARKS_summary.csv\n")
cat("  • results/CONSERVED_CORE_TARGETS_by_hallmark.csv\n")
cat("  • results/DOG_LE_vs_CAT_UP_Fisher_byHallmark.csv\n")
cat("  • results/DOG_LE_vs_CAT_DOWN_Fisher_byHallmark.csv\n")
cat("  • results/Conserved_Core_UP_genes_strict.txt\n")
cat("  • results/target_master_primary_targetlevel.csv\n")
cat("  • results/target_master_integrated_deep_indications.csv\n")
cat("  • results/target_master.csv\n")
cat("  • results/target_prioritisation_scored.csv\n")
cat("  • results/TOP_TARGETS_FORMATTED.csv\n")
cat("  • results/TRANSLATIONAL_GAP_summary.csv\n")
cat("  • results/TRANSLATIONAL_GAP_target_tiers.csv\n")
cat("  • results/TRANSLATIONAL_GAP_clinical_leverage_counts.csv\n")
cat("  • results/REPRODUCIBILITY_RECORD.txt\n\n")

cat("Pipeline finished at:", as.character(pipeline_end), "\n")
cat("Total runtime (minutes):", round(as.numeric(pipeline_runtime), 2), "\n\n")
