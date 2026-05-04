################################################################################
# STEP 6: CONSENSUS TARGET PRIORITISATION
# Rank-based feature scoring with deterministic exact top-N selection.
################################################################################

module_target_prioritisation <- function() {
  
  cat("  [6A] Checking required upstream objects...\n")
  
  required_objs <- c("target_master")
  missing_objs <- required_objs[
    !vapply(required_objs, exists, logical(1), envir = .GlobalEnv)
  ]
  
  if (length(missing_objs) > 0) {
    stop(
      "Missing required object(s) for target prioritisation: ",
      paste(missing_objs, collapse = ", "),
      ". Run previous modules first.",
      call. = FALSE
    )
  }
  
  target_master <- get("target_master", envir = .GlobalEnv)
  
  top_n <- if (exists("PARAM_CONFIG", envir = .GlobalEnv)) {
    get("PARAM_CONFIG", envir = .GlobalEnv)$top_n_targets
  } else {
    15
  }
  
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
  # Helper: rank-based points with tie-aware averaging
  # ---------------------------------------------------------------------------
  # Quantitative/count-based features are ranked in descending order.
  # Only non-zero/non-missing values can receive points.
  # Point positions range from top_n to 1.
  # Tied feature values receive the average of the positions they occupy.
  #
  # Important distinction:
  # - Feature-level points are tie-aware.
  # - Final shortlist selection is exact top-N using deterministic tie-breakers,
  #   not rank <= top_n by dense score rank.
  
  rank_points <- function(x, top_n = 15, higher_is_better = TRUE, zero_gets_zero = TRUE) {
    x_num <- suppressWarnings(as.numeric(x))
    
    valid <- !is.na(x_num)
    
    if (zero_gets_zero) {
      valid <- valid & x_num > 0
    }
    
    out <- rep(0, length(x_num))
    
    if (!any(valid)) {
      return(out)
    }
    
    vals <- x_num[valid]
    
    if (!higher_is_better) {
      vals <- -vals
    }
    
    feature_dense_rank <- dplyr::dense_rank(dplyr::desc(vals))
    keep <- feature_dense_rank <= top_n
    
    if (!any(keep)) {
      return(out)
    }
    
    kept_vals <- vals[keep]
    kept_original_idx <- which(valid)[keep]
    
    ord <- order(kept_vals, decreasing = TRUE)
    sorted_vals <- kept_vals[ord]
    sorted_idx <- kept_original_idx[ord]
    
    point_values <- rep(0, length(sorted_vals))
    
    r <- rle(sorted_vals)
    ends <- cumsum(r$lengths)
    starts <- c(1, head(ends, -1) + 1)
    
    for (i in seq_along(r$lengths)) {
      pos <- starts[i]:ends[i]
      raw_points <- pmax(top_n - pos + 1, 1)
      point_values[pos] <- mean(raw_points)
    }
    
    out[sorted_idx] <- point_values
    out
  }
  
  # ---------------------------------------------------------------------------
  # 6B. Build prioritisation table
  # ---------------------------------------------------------------------------
  
  cat("  [6B] Building integrated prioritisation table...\n")
  
  target_tbl <- target_master %>%
    dplyr::mutate(
      conservation_score = dplyr::coalesce(suppressWarnings(as.numeric(conservation_score)), 0),
      n_drugs = dplyr::coalesce(suppressWarnings(as.numeric(n_drugs)), 0),
      n_indications = dplyr::coalesce(suppressWarnings(as.numeric(n_indications)), 0),
      n_oncology_drug_indications = dplyr::coalesce(suppressWarnings(as.numeric(n_oncology_drug_indications)), 0),
      n_breast_drug_indications = dplyr::coalesce(suppressWarnings(as.numeric(n_breast_drug_indications)), 0),
      max_drug_phase = dplyr::coalesce(suppressWarnings(as.numeric(max_drug_phase)), 0),
      any_approved = any_approved %in% TRUE,
      any_onco_drug = any_onco_drug %in% TRUE,
      any_breast_drug = any_breast_drug %in% TRUE,
      has_tractability_data = has_tractability_data %in% TRUE,
      any_sm = any_sm %in% TRUE,
      any_ab = any_ab %in% TRUE,
      any_pr = any_pr %in% TRUE,
      any_oc = any_oc %in% TRUE,
      any_druggable = any_druggable %in% TRUE
    )
  
  # ---------------------------------------------------------------------------
  # 6C. Compute rank-based feature points
  # ---------------------------------------------------------------------------
  
  cat("  [6C] Computing rank-based consensus feature points...\n")
  
  target_scored <- target_tbl %>%
    dplyr::mutate(
      # Quantitative/count-based features.
      pts_conservation = rank_points(conservation_score, top_n = top_n),
      pts_phase = rank_points(max_drug_phase, top_n = top_n),
      pts_n_drugs = rank_points(n_drugs, top_n = top_n),
      pts_n_indications = rank_points(n_indications, top_n = top_n),
      pts_oncology_indications = rank_points(n_oncology_drug_indications, top_n = top_n),
      pts_breast_indications = rank_points(n_breast_drug_indications, top_n = top_n),
      
      # Small-molecule tractability is retained only as a minor binary support flag.
      # It is not ranked because TRUE/FALSE is not an ordinal quantitative variable.
      pts_sm_flag = ifelse(any_sm, 1, 0),
      
      consensus_priority_score =
        pts_conservation +
        pts_phase +
        pts_n_drugs +
        pts_n_indications +
        pts_oncology_indications +
        pts_breast_indications +
        pts_sm_flag
    ) %>%
    dplyr::arrange(
      dplyr::desc(consensus_priority_score),
      dplyr::desc(conservation_score),
      dplyr::desc(max_drug_phase),
      dplyr::desc(n_drugs),
      dplyr::desc(n_indications),
      dplyr::desc(n_oncology_drug_indications),
      dplyr::desc(n_breast_drug_indications),
      dplyr::desc(any_sm),
      target_symbol
    ) %>%
    dplyr::mutate(
      # Keeps tied score groups visible for auditing.
      score_rank = dplyr::dense_rank(dplyr::desc(consensus_priority_score)),
      
      # Unique deterministic order used for exact top-N selection.
      priority_rank = dplyr::row_number(),
      
      # Compatibility column for downstream code.
      rank = priority_rank
    ) %>%
    dplyr::select(
      rank,
      score_rank,
      priority_rank,
      target_symbol,
      target_ensembl,
      consensus_priority_score,
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
      any_druggable,
      pts_conservation,
      pts_phase,
      pts_n_drugs,
      pts_n_indications,
      pts_oncology_indications,
      pts_breast_indications,
      pts_sm_flag
    )
  
  save_output(target_scored, "target_prioritisation_scored.csv")
  
  # ---------------------------------------------------------------------------
  # 6D. Generate exact top-N shortlists
  # ---------------------------------------------------------------------------
  
  cat("  [6D] Generating target shortlists...\n")
  
  shortlist_overall <- target_scored %>%
    dplyr::filter(priority_rank <= top_n) %>%
    dplyr::arrange(priority_rank)
  
  shortlist_onco_or_approved <- target_scored %>%
    dplyr::filter((any_onco_drug %in% TRUE) | (any_approved %in% TRUE)) %>%
    dplyr::arrange(priority_rank)
  
  save_output(shortlist_overall, "shortlist_overall_top15_exact.csv")
  save_output(shortlist_onco_or_approved, "shortlist_onco_or_approved_all.csv")
  
  # Compatibility aliases for existing downstream scripts.
  # These now contain exact top-N output, despite the legacy filename.
  save_output(shortlist_overall, "shortlist_overall_inclusive_top15_cutoff.csv")
  save_output(shortlist_onco_or_approved, "shortlist_onco_or_approved_inclusive_top15_cutoff.csv")
  
  cat("      Targets scored:", nrow(target_scored), "\n")
  cat("      Exact top-", top_n, " shortlist: ", nrow(shortlist_overall), " targets\n", sep = "")
  cat("      Onco/approved evidence set:", nrow(shortlist_onco_or_approved), "targets\n")
  cat("      Top 5 targets:", paste(head(target_scored$target_symbol, 5), collapse = ", "), "\n")
  
  assign("target_scored", target_scored, envir = .GlobalEnv)
  assign("shortlist_overall", shortlist_overall, envir = .GlobalEnv)
  assign("shortlist_onco_or_approved", shortlist_onco_or_approved, envir = .GlobalEnv)
  
  cat("  ✓ STEP 6 COMPLETE\n\n")
  
  invisible(list(
    target_scored = target_scored,
    shortlist_overall = shortlist_overall,
    shortlist_onco_or_approved = shortlist_onco_or_approved
  ))
}
