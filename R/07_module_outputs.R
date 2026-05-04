################################################################################
# STEP 7: SUMMARY OUTPUTS AND TRANSLATIONAL-EVIDENCE TABLES
################################################################################

module_outputs <- function() {
  
  cat("  [7A] Checking required upstream objects...\n")
  
  required_objs <- c("target_master", "target_scored")
  missing_objs <- required_objs[
    !vapply(required_objs, exists, logical(1), envir = .GlobalEnv)
  ]
  
  if (length(missing_objs) > 0) {
    stop(
      "Missing required object(s) for output step: ",
      paste(missing_objs, collapse = ", "),
      ". Run previous modules first.",
      call. = FALSE
    )
  }
  
  target_master <- get("target_master", envir = .GlobalEnv)
  target_scored <- get("target_scored", envir = .GlobalEnv)
  
  results_dir <- if (exists("PATH_CONFIG", envir = .GlobalEnv)) {
    get("PATH_CONFIG", envir = .GlobalEnv)$results_dir
  } else {
    "results"
  }
  
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  
  save_output <- function(x, filename) {
    out_path <- file.path(results_dir, filename)
    readr::write_csv(x, out_path)
    cat("  ✓ Saved: ", normalizePath(out_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
    invisible(out_path)
  }
  
  # ---------------------------------------------------------------------------
  # 7B. Pipeline metadata
  # ---------------------------------------------------------------------------
  
  cat("  [7B] Writing pipeline metadata...\n")
  
  metadata <- data.frame(
    field = c(
      "analysis_name",
      "analysis_version",
      "dog_analysis_mode",
      "ensembl_version",
      "dog_hallmark_fdr",
      "cat_ora_fdr",
      "drug_min_phase",
      "top_n_targets",
      "translational_tier_system"
    ),
    value = c(
      if (exists("REPRODUCIBILITY_CONFIG", envir = .GlobalEnv)) get("REPRODUCIBILITY_CONFIG", envir = .GlobalEnv)$analysis_name else NA_character_,
      if (exists("REPRODUCIBILITY_CONFIG", envir = .GlobalEnv)) get("REPRODUCIBILITY_CONFIG", envir = .GlobalEnv)$analysis_version else NA_character_,
      if (exists("PARAM_CONFIG", envir = .GlobalEnv)) get("PARAM_CONFIG", envir = .GlobalEnv)$dog_analysis_mode else NA_character_,
      if (exists("PARAM_CONFIG", envir = .GlobalEnv)) get("PARAM_CONFIG", envir = .GlobalEnv)$ensembl_version else NA_character_,
      if (exists("PARAM_CONFIG", envir = .GlobalEnv)) get("PARAM_CONFIG", envir = .GlobalEnv)$dog_hallmark_padj else NA_character_,
      if (exists("PARAM_CONFIG", envir = .GlobalEnv)) get("PARAM_CONFIG", envir = .GlobalEnv)$cat_ora_padj else NA_character_,
      if (exists("PARAM_CONFIG", envir = .GlobalEnv)) get("PARAM_CONFIG", envir = .GlobalEnv)$drug_min_phase else NA_character_,
      if (exists("PARAM_CONFIG", envir = .GlobalEnv)) get("PARAM_CONFIG", envir = .GlobalEnv)$top_n_targets else NA_character_,
      "Tier 1: approved/phase >=2 evidence; Tier 2: early oncology-linked evidence; Tier 3: tractable, no retained drug evidence"
    ),
    stringsAsFactors = FALSE
  )
  
  save_output(metadata, "PIPELINE_METADATA.csv")
  
  # ---------------------------------------------------------------------------
  # 7C. Formatted top-target table
  # ---------------------------------------------------------------------------
  
  cat("  [7C] Creating formatted top-targets table...\n")
  
  top_n <- if (exists("PARAM_CONFIG", envir = .GlobalEnv)) {
    get("PARAM_CONFIG", envir = .GlobalEnv)$top_n_targets
  } else {
    15
  }
  
  target_scored_ordered <- target_scored %>%
    dplyr::arrange(rank, dplyr::desc(consensus_priority_score), target_symbol)
  
  # Inclusive top-N cutoff by rank, consistent with tie-aware prioritisation.
  top_targets_formatted <- target_scored_ordered %>%
    dplyr::filter(rank <= top_n) %>%
    dplyr::mutate(
      approved_or_phase2 = (any_approved %in% TRUE) |
        (!is.na(max_drug_phase) & suppressWarnings(as.numeric(max_drug_phase)) >= 2),
      early_oncology_linked = !(approved_or_phase2 %in% TRUE) &
        (any_onco_drug %in% TRUE),
      translational_tier = dplyr::case_when(
        approved_or_phase2 ~ "Approved/phase >=2 evidence",
        early_oncology_linked ~ "Early oncology-linked evidence",
        (has_tractability_data %in% TRUE) |
          (any_druggable %in% TRUE) |
          (any_sm %in% TRUE) |
          (any_ab %in% TRUE) |
          (any_pr %in% TRUE) |
          (any_oc %in% TRUE) ~ "Tractable, no retained drug evidence",
        TRUE ~ "Biology only"
      )
    ) %>%
    dplyr::select(
      rank,
      target_symbol,
      target_ensembl,
      consensus_priority_score,
      translational_tier,
      conservation_score,
      n_drugs,
      n_indications,
      n_oncology_drug_indications,
      n_breast_drug_indications,
      max_drug_phase,
      any_approved,
      any_onco_drug,
      any_breast_drug,
      has_tractability_data,
      any_sm,
      any_ab,
      any_pr,
      any_oc,
      any_druggable
    )
  
  save_output(top_targets_formatted, "TOP_TARGETS_FORMATTED.csv")
  
  # ---------------------------------------------------------------------------
  # 7D. Translational-evidence tiers
  # ---------------------------------------------------------------------------
  
  cat("  [7D] Building translational-evidence tier summaries...\n")
  
  target_tiers <- target_master %>%
    dplyr::mutate(
      max_drug_phase_num = suppressWarnings(as.numeric(max_drug_phase)),
      approved_or_phase2 = (any_approved %in% TRUE) |
        (!is.na(max_drug_phase_num) & max_drug_phase_num >= 2),
      early_oncology_linked = !(approved_or_phase2 %in% TRUE) &
        (any_onco_drug %in% TRUE),
      tractable_any = (has_tractability_data %in% TRUE) |
        (any_druggable %in% TRUE) |
        (any_sm %in% TRUE) |
        (any_ab %in% TRUE) |
        (any_pr %in% TRUE) |
        (any_oc %in% TRUE),
      translational_tier = dplyr::case_when(
        approved_or_phase2 ~ "Approved/phase >=2 evidence",
        early_oncology_linked ~ "Early oncology-linked evidence",
        tractable_any ~ "Tractable, no retained drug evidence",
        TRUE ~ "Biology only"
      ),
      translational_tier_order = dplyr::case_when(
        translational_tier == "Approved/phase >=2 evidence" ~ 1L,
        translational_tier == "Early oncology-linked evidence" ~ 2L,
        translational_tier == "Tractable, no retained drug evidence" ~ 3L,
        TRUE ~ 4L
      )
    ) %>%
    dplyr::arrange(translational_tier_order, target_symbol)
  
  tier_counts <- target_tiers %>%
    dplyr::count(translational_tier, translational_tier_order, name = "n") %>%
    dplyr::mutate(
      pct = round(100 * n / sum(n), 1)
    ) %>%
    dplyr::arrange(translational_tier_order)
  
  translational_gap_summary <- data.frame(
    n_total = nrow(target_tiers),
    n_approved_or_phase2 = sum(target_tiers$translational_tier == "Approved/phase >=2 evidence"),
    n_early_oncology_linked = sum(target_tiers$translational_tier == "Early oncology-linked evidence"),
    n_tractable_no_retained_drug_evidence = sum(target_tiers$translational_tier == "Tractable, no retained drug evidence"),
    n_biology_only = sum(target_tiers$translational_tier == "Biology only"),
    n_tractable = sum(target_tiers$tractable_any %in% TRUE),
    n_small_molecule_tractable = sum(target_tiers$any_sm %in% TRUE),
    pct_approved_or_phase2 = round(100 * sum(target_tiers$translational_tier == "Approved/phase >=2 evidence") / nrow(target_tiers), 1),
    pct_early_oncology_linked = round(100 * sum(target_tiers$translational_tier == "Early oncology-linked evidence") / nrow(target_tiers), 1),
    pct_tractable_no_retained_drug_evidence = round(100 * sum(target_tiers$translational_tier == "Tractable, no retained drug evidence") / nrow(target_tiers), 1),
    pct_biology_only = round(100 * sum(target_tiers$translational_tier == "Biology only") / nrow(target_tiers), 1),
    pct_tractable = round(100 * sum(target_tiers$tractable_any %in% TRUE) / nrow(target_tiers), 1),
    pct_small_molecule_tractable = round(100 * sum(target_tiers$any_sm %in% TRUE) / nrow(target_tiers), 1),
    stringsAsFactors = FALSE
  )
  
  # Compatibility table for downstream/legacy naming, but with updated tier language.
  clinical_evidence_counts <- tier_counts %>%
    dplyr::select(tier_group = translational_tier, n, pct)
  
  save_output(translational_gap_summary, "TRANSLATIONAL_GAP_summary.csv")
  save_output(target_tiers, "TRANSLATIONAL_GAP_target_tiers.csv")
  save_output(tier_counts, "TRANSLATIONAL_GAP_tier_counts.csv")
  save_output(clinical_evidence_counts, "TRANSLATIONAL_GAP_clinical_leverage_counts.csv")
  
  # ---------------------------------------------------------------------------
  # 7E. Shortlist-specific drug evidence
  # ---------------------------------------------------------------------------
  
  cat("  [7E] Exporting shortlist-specific drug evidence...\n")
  
  drug_evidence_file <- file.path(results_dir, "drug_target_table_with_ATC_and_fallback.csv")
  
  if (file.exists(drug_evidence_file)) {
    drug_evidence <- readr::read_csv(drug_evidence_file, show_col_types = FALSE)
    
    shortlist_symbols <- unique(top_targets_formatted$target_symbol)
    
    shortlist_drug_evidence <- drug_evidence %>%
      dplyr::filter(target_symbol %in% shortlist_symbols) %>%
      dplyr::arrange(target_symbol, dplyr::desc(phase_for_filter), drug_name, disease_name)
    
    save_output(shortlist_drug_evidence, "SHORTLIST_drug_evidence.csv")
  } else {
    warning("drug_target_table_with_ATC_and_fallback.csv not found; skipping SHORTLIST_drug_evidence.csv")
  }
  
  # ---------------------------------------------------------------------------
  # 7F. Human-readable summary
  # ---------------------------------------------------------------------------
  
  cat("  [7F] Writing human-readable pipeline summary...\n")
  
  summary_lines <- c(
    "DOG-CAT COMPARATIVE ONCOGENOMICS PIPELINE SUMMARY",
    "=================================================",
    "",
    paste0("Strict conserved UP targets: ", nrow(target_master)),
    paste0("Tractable targets: ", translational_gap_summary$n_tractable, "/", translational_gap_summary$n_total,
           " (", translational_gap_summary$pct_tractable, "%)"),
    paste0("Small-molecule tractable targets: ", translational_gap_summary$n_small_molecule_tractable, "/", translational_gap_summary$n_total,
           " (", translational_gap_summary$pct_small_molecule_tractable, "%)"),
    "",
    "Translational-evidence tiers:",
    paste0(
      "  - ", tier_counts$translational_tier, ": ", tier_counts$n,
      "/", translational_gap_summary$n_total, " (", tier_counts$pct, "%)"
    ),
    "",
    "Top ranked targets:",
    paste0(
      "  ", top_targets_formatted$rank, ". ", top_targets_formatted$target_symbol,
      " | score=", round(top_targets_formatted$consensus_priority_score, 2),
      " | tier=", top_targets_formatted$translational_tier
    )
  )
  
  writeLines(summary_lines, file.path(results_dir, "PIPELINE_SUMMARY.txt"))
  
  # ---------------------------------------------------------------------------
  # Export objects
  # ---------------------------------------------------------------------------
  
  assign("top_targets_formatted", top_targets_formatted, envir = .GlobalEnv)
  assign("target_tiers", target_tiers, envir = .GlobalEnv)
  assign("tier_counts", tier_counts, envir = .GlobalEnv)
  assign("translational_gap_summary", translational_gap_summary, envir = .GlobalEnv)
  
  cat("  ✓ STEP 7 COMPLETE\n\n")
  
  invisible(list(
    top_targets_formatted = top_targets_formatted,
    target_tiers = target_tiers,
    tier_counts = tier_counts,
    translational_gap_summary = translational_gap_summary
  ))
}
