################################################################################
# STEP 6B: PRIORITISATION SENSITIVITY ANALYSIS
# Previous target-level evidence vs integrated drug/probe-linked evidence
################################################################################

module_prioritisation_sensitivity_primary_vs_integrated <- function(
    previous_target_master_file = file.path(PATH_CONFIG$results_dir, "target_master_primary_targetlevel.csv"),
    current_target_master_file  = file.path(PATH_CONFIG$results_dir, "target_master.csv"),
    top_n = NULL
) {
  
  cat("\n")
  cat("======================================================================\n")
  cat("STEP 6B: PRIORITISATION SENSITIVITY ANALYSIS\n")
  cat("======================================================================\n\n")
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tibble)
    library(stringr)
  })
  
  results_dir <- if (exists("PATH_CONFIG", envir = .GlobalEnv)) {
    get("PATH_CONFIG", envir = .GlobalEnv)$results_dir
  } else {
    "results"
  }
  
  if (is.null(top_n)) {
    top_n <- if (exists("PARAM_CONFIG", envir = .GlobalEnv)) {
      get("PARAM_CONFIG", envir = .GlobalEnv)$top_n_targets
    } else {
      15
    }
  }
  
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  
  save_output <- function(x, filename) {
    out_path <- file.path(results_dir, filename)
    readr::write_csv(x, out_path)
    cat("  ✓ Saved: ", normalizePath(out_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
    invisible(out_path)
  }
  
  as_bool_safe <- function(x) {
    if (is.logical(x)) return(dplyr::coalesce(x, FALSE))
    if (is.numeric(x)) return(dplyr::coalesce(x != 0, FALSE))
    x_chr <- tolower(trimws(as.character(x)))
    out <- x_chr %in% c("true", "t", "1", "yes", "y")
    out[is.na(x_chr)] <- FALSE
    out
  }
  
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
  
  score_target_master <- function(target_master, evidence_label, top_n = 15) {
    
    required_cols <- c(
      "target_symbol",
      "target_ensembl",
      "conservation_score",
      "n_drugs",
      "n_indications",
      "n_oncology_drug_indications",
      "n_breast_drug_indications",
      "max_drug_phase",
      "any_approved",
      "any_onco_drug",
      "any_breast_drug",
      "has_tractability_data",
      "any_sm",
      "any_ab",
      "any_pr",
      "any_oc",
      "any_druggable"
    )
    
    for (col in required_cols) {
      if (!col %in% names(target_master)) {
        if (col %in% c("target_symbol", "target_ensembl")) {
          target_master[[col]] <- NA_character_
        } else if (str_starts(col, "any_") || col == "has_tractability_data") {
          target_master[[col]] <- FALSE
        } else {
          target_master[[col]] <- 0
        }
      }
    }
    
    target_tbl <- target_master %>%
      mutate(
        evidence_label = evidence_label,
        conservation_score = coalesce(suppressWarnings(as.numeric(conservation_score)), 0),
        n_drugs = coalesce(suppressWarnings(as.numeric(n_drugs)), 0),
        n_indications = coalesce(suppressWarnings(as.numeric(n_indications)), 0),
        n_oncology_drug_indications = coalesce(suppressWarnings(as.numeric(n_oncology_drug_indications)), 0),
        n_breast_drug_indications = coalesce(suppressWarnings(as.numeric(n_breast_drug_indications)), 0),
        max_drug_phase = coalesce(suppressWarnings(as.numeric(max_drug_phase)), 0),
        any_approved = as_bool_safe(any_approved),
        any_onco_drug = as_bool_safe(any_onco_drug),
        any_breast_drug = as_bool_safe(any_breast_drug),
        has_tractability_data = as_bool_safe(has_tractability_data),
        any_sm = as_bool_safe(any_sm),
        any_ab = as_bool_safe(any_ab),
        any_pr = as_bool_safe(any_pr),
        any_oc = as_bool_safe(any_oc),
        any_druggable = as_bool_safe(any_druggable)
      )
    
    target_scored <- target_tbl %>%
      mutate(
        pts_conservation = rank_points(conservation_score, top_n = top_n),
        pts_phase = rank_points(max_drug_phase, top_n = top_n),
        pts_n_drugs = rank_points(n_drugs, top_n = top_n),
        pts_n_indications = rank_points(n_indications, top_n = top_n),
        pts_oncology_indications = rank_points(n_oncology_drug_indications, top_n = top_n),
        pts_breast_indications = rank_points(n_breast_drug_indications, top_n = top_n),
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
      arrange(
        desc(consensus_priority_score),
        desc(conservation_score),
        desc(max_drug_phase),
        desc(n_drugs),
        desc(n_indications),
        desc(n_oncology_drug_indications),
        desc(n_breast_drug_indications),
        desc(any_sm),
        target_symbol
      ) %>%
      mutate(
        score_rank = dense_rank(desc(consensus_priority_score)),
        priority_rank = row_number(),
        rank = priority_rank,
        top15 = priority_rank <= top_n
      )
    
    target_scored
  }
  
  # ---------------------------------------------------------------------------
  # Load previous and current target master
  # ---------------------------------------------------------------------------
  
  if (!file.exists(previous_target_master_file)) {
    stop("Previous target_master file not found: ", previous_target_master_file, call. = FALSE)
  }
  
  if (!file.exists(current_target_master_file)) {
    stop("Current target_master file not found: ", current_target_master_file, call. = FALSE)
  }
  
  cat("  [6B-1] Loading target master files...\n")
  cat("      Previous: ", normalizePath(previous_target_master_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("      Current:  ", normalizePath(current_target_master_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
  
  previous_master <- readr::read_csv(previous_target_master_file, show_col_types = FALSE)
  current_master  <- readr::read_csv(current_target_master_file, show_col_types = FALSE)
  
  # ---------------------------------------------------------------------------
  # Score both using the same framework
  # ---------------------------------------------------------------------------
  
  cat("  [6B-2] Scoring previous and integrated evidence tables...\n")
  
  previous_scored <- score_target_master(
    previous_master,
    evidence_label = "Primary target-level evidence",
    top_n = top_n
  )
  
  integrated_scored <- score_target_master(
    current_master,
    evidence_label = "Integrated drug/probe-linked evidence",
    top_n = top_n
  )
  
  save_output(previous_scored, "target_prioritisation_previous_scored.csv")
  save_output(integrated_scored, "target_prioritisation_integrated_rescored_for_sensitivity.csv")
  
  previous_top15 <- previous_scored %>%
    filter(priority_rank <= top_n) %>%
    arrange(priority_rank)
  
  integrated_top15 <- integrated_scored %>%
    filter(priority_rank <= top_n) %>%
    arrange(priority_rank)
  
  save_output(previous_top15, "shortlist_previous_top15_for_sensitivity.csv")
  save_output(integrated_top15, "shortlist_integrated_top15_for_sensitivity.csv")
  
  # ---------------------------------------------------------------------------
  # Build comparison table
  # ---------------------------------------------------------------------------
  
  cat("  [6B-3] Building sensitivity comparison table...\n")
  
  previous_cmp <- previous_scored %>%
    select(
      target_symbol,
      previous_rank = priority_rank,
      previous_score = consensus_priority_score,
      previous_conservation_score = conservation_score,
      previous_max_phase = max_drug_phase,
      previous_n_drug_records = n_drugs,
      previous_n_indication_records = n_indications,
      previous_n_oncology_records = n_oncology_drug_indications,
      previous_n_breast_records = n_breast_drug_indications,
      previous_any_approved = any_approved,
      previous_any_onco = any_onco_drug,
      previous_any_breast = any_breast_drug,
      previous_any_sm = any_sm,
      previous_top15 = top15
    )
  
  integrated_cmp <- integrated_scored %>%
    select(
      target_symbol,
      integrated_rank = priority_rank,
      integrated_score = consensus_priority_score,
      integrated_conservation_score = conservation_score,
      integrated_max_phase = max_drug_phase,
      integrated_n_drug_records = n_drugs,
      integrated_n_indication_records = n_indications,
      integrated_n_oncology_records = n_oncology_drug_indications,
      integrated_n_breast_records = n_breast_drug_indications,
      integrated_any_approved = any_approved,
      integrated_any_onco = any_onco_drug,
      integrated_any_breast = any_breast_drug,
      integrated_any_sm = any_sm,
      integrated_top15 = top15
    )
  
  sensitivity_table <- previous_cmp %>%
    full_join(integrated_cmp, by = "target_symbol") %>%
    mutate(
      rank_change = previous_rank - integrated_rank,
      score_change = integrated_score - previous_score,
      max_phase_change = integrated_max_phase - previous_max_phase,
      n_drug_records_change = integrated_n_drug_records - previous_n_drug_records,
      n_indication_records_change = integrated_n_indication_records - previous_n_indication_records,
      n_oncology_records_change = integrated_n_oncology_records - previous_n_oncology_records,
      n_breast_records_change = integrated_n_breast_records - previous_n_breast_records,
      top15_status = case_when(
        previous_top15 %in% TRUE & integrated_top15 %in% TRUE ~ "Top-15 in both",
        previous_top15 %in% TRUE & !(integrated_top15 %in% TRUE) ~ "Dropped from top-15",
        !(previous_top15 %in% TRUE) & integrated_top15 %in% TRUE ~ "Entered top-15",
        TRUE ~ "Not top-15"
      )
    ) %>%
    arrange(integrated_rank)
  
  save_output(sensitivity_table, "TABLE_S8_prioritisation_sensitivity_previous_vs_integrated.csv")
  
  # ---------------------------------------------------------------------------
  # Summary table
  # ---------------------------------------------------------------------------
  
  previous_top15_symbols <- previous_top15$target_symbol
  integrated_top15_symbols <- integrated_top15$target_symbol
  
  top15_overlap <- intersect(previous_top15_symbols, integrated_top15_symbols)
  
  same_top15_set <- setequal(previous_top15_symbols, integrated_top15_symbols)
  same_top15_order <- identical(previous_top15_symbols, integrated_top15_symbols)
  
  stability_summary <- tibble(
    metric = c(
      "Number of targets scored in previous evidence table",
      "Number of targets scored in integrated evidence table",
      "Top-N threshold",
      "Number of overlapping top-N targets",
      "Same top-N set",
      "Same top-N order",
      "Targets entering top-N after integration",
      "Targets dropping from top-N after integration"
    ),
    value = c(
      as.character(nrow(previous_scored)),
      as.character(nrow(integrated_scored)),
      as.character(top_n),
      as.character(length(top15_overlap)),
      as.character(same_top15_set),
      as.character(same_top15_order),
      paste(setdiff(integrated_top15_symbols, previous_top15_symbols), collapse = ", "),
      paste(setdiff(previous_top15_symbols, integrated_top15_symbols), collapse = ", ")
    )
  )
  
  save_output(stability_summary, "TABLE_S8_top15_stability_summary.csv")
  
  cat("\n")
  cat("  Top-", top_n, " overlap: ", length(top15_overlap), "/", top_n, "\n", sep = "")
  cat("  Same top-", top_n, " set: ", same_top15_set, "\n", sep = "")
  cat("  Same top-", top_n, " order: ", same_top15_order, "\n", sep = "")
  
  cat("\n  Previous top-", top_n, ":\n  ", paste(previous_top15_symbols, collapse = ", "), "\n", sep = "")
  cat("\n  Integrated top-", top_n, ":\n  ", paste(integrated_top15_symbols, collapse = ", "), "\n", sep = "")
  
  assign("previous_prioritisation_scored", previous_scored, envir = .GlobalEnv)
  assign("integrated_prioritisation_sensitivity_scored", integrated_scored, envir = .GlobalEnv)
  assign("prioritisation_sensitivity_table", sensitivity_table, envir = .GlobalEnv)
  assign("prioritisation_top15_stability_summary", stability_summary, envir = .GlobalEnv)
  
  cat("\n  ✓ STEP 6B COMPLETE\n\n")
  
  invisible(list(
    previous_scored = previous_scored,
    integrated_scored = integrated_scored,
    previous_top15 = previous_top15,
    integrated_top15 = integrated_top15,
    sensitivity_table = sensitivity_table,
    stability_summary = stability_summary
  ))
}