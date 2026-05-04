module_permutation_validation_light <- function(
    B = 10000,
    seed = 1
) {
  
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("PERMUTATION VALIDATION AGAINST LIGHT BACKGROUND (FIXED)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(tidyr)
  })
  
  # ---------------------------------------------------------------------------
  # 1. FILE CHECKS
  # ---------------------------------------------------------------------------
  req_files <- c(
    file.path(PATH_CONFIG$results_dir, "comparative_background_tractability_light.csv"),
    file.path(PATH_CONFIG$results_dir, "Conserved_Core_UP_genes_strict.txt")
  )
  
  missing_files <- req_files[!file.exists(req_files)]
  if (length(missing_files) > 0) {
    stop(
      "Missing required file(s) for permutation validation:\n",
      paste(" - ", missing_files, collapse = "\n")
    )
  }
  
  # ---------------------------------------------------------------------------
  # 2. LOAD INPUTS
  # ---------------------------------------------------------------------------
  cat("  [PV-L1] Loading corrected background and strict core...\n")
  
  bg <- read_csv(
    file.path(PATH_CONFIG$results_dir, "comparative_background_tractability_light.csv"),
    show_col_types = FALSE
  )
  
  strict_core <- readLines(
    file.path(PATH_CONFIG$results_dir, "Conserved_Core_UP_genes_strict.txt"),
    warn = FALSE
  )
  strict_core <- unique(trimws(strict_core))
  strict_core <- strict_core[strict_core != ""]
  
  if (length(strict_core) == 0) {
    stop("Strict core is empty.")
  }
  
  # ---------------------------------------------------------------------------
  # 3. USE BACKGROUND AS-IS
  # ---------------------------------------------------------------------------
  cat("  [PV-L2] Preparing light background from corrected CSV...\n")
  
  bg2 <- bg %>%
    transmute(
      target_ensembl = as.character(target_ensembl),
      target_symbol = as.character(target_symbol),
      any_sm = as.logical(any_sm),
      any_ab = as.logical(any_ab),
      any_pr = as.logical(any_pr),
      any_oc = as.logical(any_oc),
      any_druggable = as.logical(any_druggable),
      n_drugs = suppressWarnings(as.numeric(n_drugs)),
      max_drug_phase = suppressWarnings(as.numeric(max_drug_phase)),
      any_approved = as.logical(any_approved)
    ) %>%
    filter(
      !is.na(target_ensembl), target_ensembl != "",
      !is.na(target_symbol), target_symbol != ""
    ) %>%
    distinct(target_ensembl, .keep_all = TRUE) %>%
    mutate(
      n_drugs = coalesce(n_drugs, 0),
      max_drug_phase = coalesce(max_drug_phase, 0),
      any_approved = coalesce(any_approved, FALSE),
      any_sm = coalesce(any_sm, FALSE),
      any_ab = coalesce(any_ab, FALSE),
      any_pr = coalesce(any_pr, FALSE),
      any_oc = coalesce(any_oc, FALSE),
      any_druggable = coalesce(any_druggable, FALSE),
      clinical_evidence = any_approved | (max_drug_phase >= 2),
      sm_tractable = any_sm %in% TRUE,
      sm_unexploited = sm_tractable & !clinical_evidence
    )
  obs_df <- bg2 %>%
    filter(target_symbol %in% strict_core)
  
  k <- nrow(obs_df)
  
  missing_from_bg <- setdiff(strict_core, bg2$target_symbol)
  if (length(missing_from_bg) > 0) {
    cat("      Strict core genes missing from background:", length(missing_from_bg), "\n")
  }
  
  if (k < 5) {
    stop("Too few strict core genes remain in the background universe (k = ", k, ").")
  }
  
  cat("      Background size:", nrow(bg2), "\n")
  cat("      Strict core size:", length(strict_core), "\n")
  cat("      Strict core in background:", k, "\n")
  cat("      Background SM tractable %:", round(100 * mean(bg2$sm_tractable), 1), "\n")
  cat("      Background clinical evidence %:", round(100 * mean(bg2$clinical_evidence), 1), "\n")
  cat("      Background SM unexploited %:", round(100 * mean(bg2$sm_unexploited), 1), "\n")
  
  # ---------------------------------------------------------------------------
  # 4. OBSERVED METRICS
  # ---------------------------------------------------------------------------
  cat("  [PV-L3] Computing observed metrics...\n")
  
  obs <- tibble(
    metric = c("pct_sm_tractable", "pct_clinical_evidence", "pct_sm_unexploited"),
    observed = c(
      mean(obs_df$sm_tractable),
      mean(obs_df$clinical_evidence),
      mean(obs_df$sm_unexploited)
    )
  )
  
  # ---------------------------------------------------------------------------
  # 5. NULL POOL
  # ---------------------------------------------------------------------------
  cat("  [PV-L4] Building null pool excluding strict core genes...\n")
  
  bg_null <- bg2 %>%
    filter(!target_symbol %in% strict_core)
  
  N_null <- nrow(bg_null)
  
  if (N_null < k) {
    stop(
      "Null background after excluding strict core is smaller than k.\n",
      "N_null = ", N_null, ", k = ", k
    )
  }
  
  cat("      Null pool size:", N_null, "\n")
  cat("      Null clinical evidence count:", sum(bg_null$clinical_evidence, na.rm = TRUE), "\n")
  
  # ---------------------------------------------------------------------------
  # 6. PERMUTATION NULL
  # ---------------------------------------------------------------------------
  cat("  [PV-L5] Running permutations...\n")
  
  set.seed(seed)
  
  perm_mat <- replicate(B, {
    idx <- sample.int(N_null, k, replace = FALSE)
    samp <- bg_null[idx, ]
    
    c(
      pct_sm_tractable = mean(samp$sm_tractable),
      pct_clinical_evidence = mean(samp$clinical_evidence),
      pct_sm_unexploited = mean(samp$sm_unexploited)
    )
  })
  
  perm_df <- as.data.frame(t(perm_mat))
  perm_df$iter <- seq_len(nrow(perm_df))
  
  write_csv(
    perm_df,
    file.path(PATH_CONFIG$results_dir, "VALIDATION_permutation_null_distribution_light.csv")
  )
  
  # ---------------------------------------------------------------------------
  # 7. EMPIRICAL P-VALUES
  # ---------------------------------------------------------------------------
  cat("  [PV-L6] Calculating empirical p-values...\n")
  
  p_empirical <- function(null, obs) {
    (sum(null >= obs) + 1) / (length(null) + 1)
  }
  
  pvals <- tibble(
    metric = c("pct_sm_tractable", "pct_clinical_evidence", "pct_sm_unexploited"),
    p_value = c(
      p_empirical(perm_df$pct_sm_tractable,
                  obs$observed[obs$metric == "pct_sm_tractable"]),
      p_empirical(perm_df$pct_clinical_evidence,
                  obs$observed[obs$metric == "pct_clinical_evidence"]),
      p_empirical(perm_df$pct_sm_unexploited,
                  obs$observed[obs$metric == "pct_sm_unexploited"])
    )
  )
  
  null_summary <- bind_rows(
    tibble(metric = "pct_sm_tractable",
           null_mean = mean(perm_df$pct_sm_tractable),
           null_sd = sd(perm_df$pct_sm_tractable)),
    tibble(metric = "pct_clinical_evidence",
           null_mean = mean(perm_df$pct_clinical_evidence),
           null_sd = sd(perm_df$pct_clinical_evidence)),
    tibble(metric = "pct_sm_unexploited",
           null_mean = mean(perm_df$pct_sm_unexploited),
           null_sd = sd(perm_df$pct_sm_unexploited))
  )
  
  obs_vs_null <- obs %>%
    left_join(null_summary, by = "metric") %>%
    left_join(pvals, by = "metric") %>%
    mutate(
      observed_pct = round(100 * observed, 1),
      null_mean_pct = round(100 * null_mean, 1),
      null_sd_pct = round(100 * null_sd, 1),
      z_score = ifelse(null_sd > 0, (observed - null_mean) / null_sd, NA_real_)
    )
  
  write_csv(
    obs_vs_null,
    file.path(PATH_CONFIG$results_dir, "VALIDATION_permutation_observed_vs_null_light.csv")
  )

  metric_definitions <- tibble::tibble(
    metric = c("pct_sm_tractable", "pct_clinical_evidence", "pct_sm_unexploited"),
    definition = c(
      "Proportion of sampled targets with small-molecule tractability in Open Targets.",
      "Proportion of sampled targets with approved-drug evidence or maximum clinical phase >= 2.",
      "Proportion of sampled targets with small-molecule tractability but without established clinical leverage."
    )
  )
  write_csv(
    metric_definitions,
    file.path(PATH_CONFIG$results_dir, "VALIDATION_permutation_metric_definitions.csv")
  )
  
  write_csv(
    pvals,
    file.path(PATH_CONFIG$results_dir, "VALIDATION_permutation_pvalues_light.csv")
  )
  
  # ---------------------------------------------------------------------------
  # 8. DIAGNOSTICS
  # ---------------------------------------------------------------------------
  cat("  [PV-L7] Diagnostics...\n")
  
  diag_df <- tibble(
    metric = c("pct_sm_tractable", "pct_clinical_evidence", "pct_sm_unexploited"),
    n_unique_null_values = c(
      length(unique(perm_df$pct_sm_tractable)),
      length(unique(perm_df$pct_clinical_evidence)),
      length(unique(perm_df$pct_sm_unexploited))
    )
  )
  
  write_csv(
    diag_df,
    file.path(PATH_CONFIG$results_dir, "VALIDATION_permutation_diagnostics_light.csv")
  )
  
  # ---------------------------------------------------------------------------
  # 9. PLOT
  # ---------------------------------------------------------------------------
  cat("  [PV-L8] Building permutation validation panel...\n")
  
  plot_df <- perm_df %>%
    pivot_longer(
      cols = c(pct_sm_tractable, pct_clinical_evidence, pct_sm_unexploited),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c("pct_sm_tractable", "pct_clinical_evidence", "pct_sm_unexploited"),
        labels = c("SM tractable", "Clinical evidence", "SM tractable but\nclinically unexploited")
      )
    )
  
  obs_plot_df <- obs_vs_null %>%
    mutate(
      metric = factor(
        metric,
        levels = c("pct_sm_tractable", "pct_clinical_evidence", "pct_sm_unexploited"),
        labels = c("SM tractable", "Clinical evidence", "SM tractable but\nclinically unexploited")
      ),
      label = paste0(
        "obs = ", round(observed, 3),
        "\nnull = ", round(null_mean, 3),
        "\np = ", signif(p_value, 3)
      )
    )
  
  ann_df <- plot_df %>%
    group_by(metric) %>%
    summarise(
      xmax = max(value, na.rm = TRUE),
      ymax = Inf,
      .groups = "drop"
    ) %>%
    left_join(obs_plot_df %>% select(metric, label), by = "metric")
  
  p_perm <- ggplot(plot_df, aes(x = value)) +
    geom_histogram(bins = 20, fill = "#BFD7EA", color = "white") +
    geom_vline(
      data = obs_plot_df,
      aes(xintercept = observed),
      color = "#D55E00",
      linewidth = 1
    ) +
    geom_text(
      data = ann_df,
      aes(x = xmax, y = ymax, label = label),
      hjust = 1.05,
      vjust = 1.2,
      size = 3.5,
      color = "#D55E00",
      inherit.aes = FALSE
    ) +
    facet_wrap(~ metric, scales = "free", ncol = 1) +
    labs(
      title = "Permutation validation against comparative light background",
      subtitle = "Null pool excludes strict conserved core genes",
      x = "Permutation value",
      y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0)
    )
  
  ggsave(
    file.path(PATH_CONFIG$results_dir, "VALIDATION_permutation_panel_light.png"),
    p_perm, width = 7.5, height = 9, dpi = 300
  )
  ggsave(
    file.path(PATH_CONFIG$results_dir, "VALIDATION_permutation_panel_light.pdf"),
    p_perm, width = 7.5, height = 9
  )
  
  # ---------------------------------------------------------------------------
  # 10. SUMMARY
  # ---------------------------------------------------------------------------
  cat("  [PV-L9] Writing summary...\n")
  
  summary_lines <- c(
    "PERMUTATION VALIDATION AGAINST COMPARATIVE LIGHT BACKGROUND (FIXED)",
    "",
    paste("Background size:", nrow(bg2)),
    paste("Strict core size:", length(strict_core)),
    paste("Strict core represented in background:", k),
    paste("Null pool size:", N_null),
    paste("Permutations:", B)
  )
  
  writeLines(
    summary_lines,
    file.path(PATH_CONFIG$results_dir, "VALIDATION_permutation_summary_light.txt")
  )
  
  assign("validation_background_light", bg2, envir = .GlobalEnv)
  assign("validation_obs_light", obs_vs_null, envir = .GlobalEnv)
  assign("validation_perm_light", perm_df, envir = .GlobalEnv)
  
  cat("\n")
  cat("  ✓ PERMUTATION VALIDATION COMPLETE\n\n")
  
  invisible(list(
    background = bg2,
    observed_vs_null = obs_vs_null,
    permutations = perm_df,
    diagnostics = diag_df
  ))
}