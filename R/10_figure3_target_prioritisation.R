################################################################################
# FIGURE 3: TARGET PRIORITISATION AND TRANSLATIONAL LANDSCAPE
################################################################################

module_figure3_target_prioritisation <- function() {
  
  cat("\n")
  cat("======================================================================\n")
  cat("FIGURE 3: TARGET PRIORITISATION AND TRANSLATIONAL LANDSCAPE\n")
  cat("======================================================================\n\n")
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(tidyr)
    library(stringr)
  })
  
  # ---------------------------------------------------------------------------
  # 1. DIRECTORIES
  # ---------------------------------------------------------------------------
  results_dir <- if (exists("PATH_CONFIG", envir = .GlobalEnv)) {
    get("PATH_CONFIG", envir = .GlobalEnv)$results_dir
  } else {
    "results"
  }
  
  fig_dir <- file.path(results_dir, "FIG3_panels")
  tab_dir <- file.path(results_dir, "FIG3_tables")
  
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("  Results dir: ", normalizePath(results_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("  FIG3 panel dir: ", normalizePath(fig_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("  FIG3 table dir: ", normalizePath(tab_dir, winslash = "/", mustWork = FALSE), "\n\n", sep = "")
  
  save_panel <- function(plot_obj, filename_base, width, height) {
    png_path <- file.path(fig_dir, paste0(filename_base, ".png"))
    pdf_path <- file.path(fig_dir, paste0(filename_base, ".pdf"))
    
    ggsave(png_path, plot_obj, width = width, height = height, dpi = 300)
    ggsave(pdf_path, plot_obj, width = width, height = height)
    
    cat("    Saved: ", normalizePath(png_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
    cat("    Saved: ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  }
  
  # ---------------------------------------------------------------------------
  # 2. INPUT FILES
  # ---------------------------------------------------------------------------
  shortlist_file <- file.path(results_dir, "shortlist_overall_inclusive_top15_cutoff.csv")
  scored_file <- file.path(results_dir, "target_prioritisation_scored.csv")
  master_file <- file.path(results_dir, "target_master.csv")
  
  req_files <- c(scored_file, master_file)
  missing_files <- req_files[!file.exists(req_files)]
  
  if (length(missing_files) > 0) {
    stop(
      "Missing required file(s) for Figure 3:\n",
      paste(" - ", missing_files, collapse = "\n")
    )
  }
  
  cat("  [FIG3-A] Loading inputs...\n")
  
  target_scored <- read_csv(scored_file, show_col_types = FALSE)
  target_master <- read_csv(master_file, show_col_types = FALSE)
  
  if (file.exists(shortlist_file)) {
    shortlist_overall <- read_csv(shortlist_file, show_col_types = FALSE)
  } else {
    cat("    NOTE: shortlist_overall_inclusive_top15_cutoff.csv not found.\n")
    cat("          Reconstructing inclusive shortlist from target_prioritisation_scored.csv\n")
    
    top_n <- if (exists("PARAM_CONFIG", envir = .GlobalEnv)) {
      get("PARAM_CONFIG", envir = .GlobalEnv)$top_n_targets
    } else {
      15
    }
    
    target_scored <- target_scored %>%
      arrange(desc(consensus_priority_score), desc(conservation_score), target_symbol)
    
    cutoff_score <- target_scored$consensus_priority_score[min(top_n, nrow(target_scored))]
    
    shortlist_overall <- target_scored %>%
      filter(consensus_priority_score >= cutoff_score) %>%
      arrange(desc(consensus_priority_score), desc(conservation_score), target_symbol)
  }
  
  # ---------------------------------------------------------------------------
  # 3. PANEL A DATA
  # ---------------------------------------------------------------------------
  cat("  [FIG3-B] Preparing top target ranking data...\n")
  
  fig3a_df <- shortlist_overall %>%
    arrange(desc(consensus_priority_score), desc(conservation_score), target_symbol) %>%
    mutate(
      target_symbol = factor(target_symbol, levels = rev(target_symbol))
    ) %>%
    select(
      rank, target_symbol, consensus_priority_score,
      conservation_score, n_drugs, max_drug_phase, any_approved
    )
  
  write_csv(fig3a_df, file.path(tab_dir, "FIG3A_top_ranked_targets_data.csv"))
  
  # ---------------------------------------------------------------------------
  # 4. PANEL B DATA
  # ---------------------------------------------------------------------------
  cat("  [FIG3-C] Preparing binary evidence tile data...\n")
  
  fig3b_table <- shortlist_overall %>%
    mutate(
      max_drug_phase_num = suppressWarnings(as.numeric(max_drug_phase)),
      
      approved_or_phase2 = as.numeric(
        (any_approved %in% TRUE) |
          (!is.na(max_drug_phase_num) & max_drug_phase_num >= 2)
      ),
      
      oncology_linked = as.numeric(any_onco_drug %in% TRUE),
      breast_linked   = as.numeric(any_breast_drug %in% TRUE)
    ) %>%
    transmute(
      target_symbol,
      approved_or_phase2 = approved_or_phase2,
      oncology_linked = oncology_linked,
      breast_linked = breast_linked,
      any_drug = as.numeric(n_drugs > 0),
      sm_tractable = as.numeric(any_sm %in% TRUE)
    )
  
  write_csv(fig3b_table, file.path(tab_dir, "FIG3B_target_component_table.csv"))
  
  metric_labels_b <- c(
    approved_or_phase2 = "Approved/\nphase >=2",
    oncology_linked = "Oncology-\nlinked",
    breast_linked = "Breast-\nlinked",
    any_drug = "Any\ndrug",
    sm_tractable = "SM\ntractable"
  )
  
  fig3b_df <- fig3b_table %>%
    pivot_longer(
      cols = -target_symbol,
      names_to = "metric",
      values_to = "present"
    ) %>%
    mutate(
      metric = factor(metric, levels = names(metric_labels_b), labels = metric_labels_b),
      target_symbol = factor(target_symbol, levels = rev(fig3b_table$target_symbol)),
      present = factor(present, levels = c(0, 1), labels = c("No", "Yes"))
    )
  
  write_csv(fig3b_df, file.path(tab_dir, "FIG3B_binary_evidence_tile_data.csv"))
  
  # ---------------------------------------------------------------------------
  # 5. PANEL C DATA
  # ---------------------------------------------------------------------------
  cat("  [FIG3-D] Preparing translational evidence tier data...\n")
  
  fig3c_target_level <- target_master %>%
    mutate(
      max_drug_phase_num = dplyr::coalesce(
        suppressWarnings(as.numeric(max_drug_phase)), 0
      ),
      n_drugs_num = dplyr::coalesce(
        suppressWarnings(as.numeric(n_drugs)), 0
      ),
      
      clinical_leverage =
        (any_approved %in% TRUE) |
        (max_drug_phase_num >= 2),
      
      retained_drug_evidence =
        n_drugs_num > 0,
      
      oncology_drug_evidence =
        any_onco_drug %in% TRUE,
      
      pharmacologically_tractable =
        (has_tractability_data %in% TRUE) |
        (any_druggable %in% TRUE) |
        (any_sm %in% TRUE) |
        (any_ab %in% TRUE) |
        (any_pr %in% TRUE) |
        (any_oc %in% TRUE),
      
      translational_tier = case_when(
        clinical_leverage ~
          "Clinical leverage (approved/phase >=2)",
        
        !clinical_leverage & oncology_drug_evidence ~
          "Early oncology drug evidence only",
        
        !clinical_leverage & !oncology_drug_evidence & retained_drug_evidence ~
          "Other drug evidence only",
        
        !retained_drug_evidence & pharmacologically_tractable ~
          "Tractable, no drug evidence",
        
        TRUE ~
          "Biology only"
      )
    )
  
  write_csv(fig3c_target_level, file.path(tab_dir, "FIG3C_target_level_tiers.csv"))
  
  fig3c_df <- fig3c_target_level %>%
    count(translational_tier, name = "count") %>%
    mutate(
      percent = 100 * count / sum(count),
      label_text = paste0(count, " (", sprintf("%.1f", percent), "%)")
    )
  
  preferred_levels <- c(
    "Biology only",
    "Tractable, no drug evidence",
    "Other drug evidence only",
    "Early oncology drug evidence only",
    "Clinical leverage (approved/phase >=2)"
  )
  
  fig3c_df <- fig3c_df %>%
    mutate(
      translational_tier = factor(
        translational_tier,
        levels = preferred_levels[preferred_levels %in% translational_tier]
      )
    ) %>%
    arrange(translational_tier)
  
  write_csv(fig3c_df, file.path(tab_dir, "FIG3C_translational_gap_tiers_data.csv"))
  
  # ---------------------------------------------------------------------------
  # 6. COMMON THEME
  # ---------------------------------------------------------------------------
  fig_theme <- theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      axis.text.y = element_text(face = "bold")
    )
  
  # ---------------------------------------------------------------------------
  # 7. PANEL A PLOT
  # ---------------------------------------------------------------------------
  cat("  [FIG3-E] Plotting Panel A...\n")
  
  pA <- ggplot(fig3a_df, aes(x = consensus_priority_score, y = target_symbol)) +
    geom_col(width = 0.7, fill = "black") +
    geom_text(
      aes(label = round(consensus_priority_score, 2)),
      hjust = -0.15,
      size = 3.5
    ) +
    labs(
      title = "A. Highest-priority conserved targets",
      x = "Consensus prioritisation score",
      y = NULL
    ) +
    fig_theme +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(15, 35, 15, 15)) +
    expand_limits(x = max(fig3a_df$consensus_priority_score, na.rm = TRUE) * 1.15)
  
  save_panel(pA, "FIG3A_top_ranked_targets", width = 8.8, height = 5.8)
  
  # ---------------------------------------------------------------------------
  # 8. PANEL B PLOT
  # ---------------------------------------------------------------------------
  cat("  [FIG3-F] Plotting Panel B...\n")
  
  pB <- ggplot(fig3b_df, aes(x = metric, y = target_symbol, fill = present)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_manual(
      values = c("No" = "grey90", "Yes" = "black"),
      name = "Evidence"
    ) +
    labs(
      title = "B. Translational evidence landscape",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid = element_blank(),
      legend.position = "top"
    )
  
  save_panel(pB, "FIG3B_translational_evidence_landscape", width = 8.8, height = 5.8)
  
  # ---------------------------------------------------------------------------
  # 9. PANEL C PLOT
  # ---------------------------------------------------------------------------
  cat("  [FIG3-G] Plotting Panel C...\n")
  
  pC <- ggplot(fig3c_df, aes(x = count, y = translational_tier)) +
    geom_col(width = 0.65, fill = "black") +
    geom_text(
      aes(label = label_text),
      hjust = -0.05,
      size = 4
    ) +
    labs(
      title = "C. Translational evidence tier distribution",
      x = "Number of conserved targets",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold")
    ) +
    expand_limits(x = max(fig3c_df$count, na.rm = TRUE) * 1.12)
  
  save_panel(pC, "FIG3C_translational_gap_tiers", width = 8.2, height = 4.8)
  
  # ---------------------------------------------------------------------------
  # 10. DONE
  # ---------------------------------------------------------------------------
  cat("\n")
  cat("  ✓ FIGURE 3 COMPLETE\n\n")
  cat("  Panel files saved in:\n")
  cat("   ", normalizePath(fig_dir, winslash = "/", mustWork = FALSE), "\n")
  cat("  Plot-ready tables saved in:\n")
  cat("   ", normalizePath(tab_dir, winslash = "/", mustWork = FALSE), "\n\n")
  
  invisible(list(
    fig3a_df = fig3a_df,
    fig3b_df = fig3b_df,
    fig3c_target_level = fig3c_target_level,
    fig3c_df = fig3c_df,
    pA = pA,
    pB = pB,
    pC = pC
  ))
}