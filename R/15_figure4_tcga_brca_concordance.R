#15_figure4_tcga_brca_concordance.R
################################################################################
# FIGURE 4: TCGA-BRCA EXTERNAL CONCORDANCE FIGURE
################################################################################

figure4_tcga_brca_concordance <- function(save_pdf = FALSE) {
  
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("FIGURE 4: TCGA-BRCA EXTERNAL CONCORDANCE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(stringr)
    library(tibble)
  })
  
  # ---------------------------------------------------------------------------
  # 1. PATHS / HELPERS
  # ---------------------------------------------------------------------------
  if (exists("PATH_CONFIG", envir = .GlobalEnv)) {
    results_dir <- get("PATH_CONFIG", envir = .GlobalEnv)$results_dir
  } else {
    results_dir <- "results"
  }
  
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  
  save_plot_multi <- function(plot_obj, filename_base, width = 8, height = 6) {
    
    filename_base <- stringr::str_remove(filename_base, "\\.png$")
    filename_base <- stringr::str_remove(filename_base, "\\.pdf$")
    
    png_path <- file.path(results_dir, paste0(filename_base, ".png"))
    pdf_path <- file.path(results_dir, paste0(filename_base, ".pdf"))
    
    ggplot2::ggsave(
      filename = png_path,
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      device = "png",
      bg = "white"
    )
    
    cat("  ✓ Saved: ", normalizePath(png_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
    
    if (isTRUE(save_pdf)) {
      pdf_status <- tryCatch(
        {
          ggplot2::ggsave(
            filename = pdf_path,
            plot = plot_obj,
            width = width,
            height = height,
            units = "in",
            device = grDevices::cairo_pdf,
            bg = "white"
          )
          TRUE
        },
        error = function(e) {
          warning(
            "Could not save PDF file: ",
            normalizePath(pdf_path, winslash = "/", mustWork = FALSE),
            "\nReason: ", conditionMessage(e),
            "\nPNG output was saved successfully. Close any open PDF viewer or shorten the project path, then rerun with save_pdf = TRUE if PDF output is required.",
            call. = FALSE
          )
          FALSE
        }
      )
      
      if (isTRUE(pdf_status)) {
        cat("  ✓ Saved: ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
      }
    }
    
    invisible(png_path)
  }
  
  save_path <- function(name) file.path(results_dir, name)
  
  # ---------------------------------------------------------------------------
  # 2. INPUT FILES
  # ---------------------------------------------------------------------------
  req_files <- c(
    file.path(results_dir, "TCGA_BRCA_core_gene_results.csv"),
    file.path(results_dir, "TCGA_BRCA_concordance_summary.csv")
  )
  
  missing_files <- req_files[!file.exists(req_files)]
  if (length(missing_files) > 0) {
    stop(
      "Missing required file(s) for TCGA figure:\n",
      paste(" - ", missing_files, collapse = "\n"),
      call. = FALSE
    )
  }
  
  cat("  [FIG4-A] Loading TCGA concordance tables...\n")
  
  core_res <- readr::read_csv(
    file.path(results_dir, "TCGA_BRCA_core_gene_results.csv"),
    show_col_types = FALSE
  )
  
  concordance_summary <- readr::read_csv(
    file.path(results_dir, "TCGA_BRCA_concordance_summary.csv"),
    show_col_types = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 3. CLEAN / DERIVED COLUMNS
  # ---------------------------------------------------------------------------
  cat("  [FIG4-B] Preparing plotting tables...\n")
  
  required_core_cols <- c(
    "gene_symbol",
    "core_direction",
    "tested_in_tcga",
    "logFC",
    "adj.P.Val",
    "significant",
    "directionally_concordant",
    "significant_and_concordant"
  )
  
  missing_core_cols <- setdiff(required_core_cols, colnames(core_res))
  if (length(missing_core_cols) > 0) {
    stop(
      "TCGA_BRCA_core_gene_results.csv is missing required column(s): ",
      paste(missing_core_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  core_plot <- core_res %>%
    mutate(
      core_direction = factor(core_direction, levels = c("UP", "DOWN")),
      tested_in_tcga = coalesce(as.logical(tested_in_tcga), FALSE),
      significant = coalesce(as.logical(significant), FALSE),
      directionally_concordant = coalesce(as.logical(directionally_concordant), FALSE),
      significant_and_concordant = coalesce(as.logical(significant_and_concordant), FALSE),
      logFC = suppressWarnings(as.numeric(logFC)),
      adj.P.Val = suppressWarnings(as.numeric(adj.P.Val)),
      neglog10_fdr = ifelse(!is.na(adj.P.Val) & adj.P.Val > 0, -log10(adj.P.Val), NA_real_),
      effect_zone = case_when(
        is.na(logFC) ~ "Inside",
        logFC <= -1 ~ "Outside_DOWN",
        logFC >= 1 ~ "Outside_UP",
        TRUE ~ "Inside"
      )
    )
  
  # Labels are selected reproducibly rather than hard-coded:
  # - all tested DOWN genes, because the DOWN core is small
  # - the five most significant UP genes with |log2FC| >= 1,
  #   prioritising significant + concordant genes
  up_labels <- core_plot %>%
    filter(
      tested_in_tcga %in% TRUE,
      core_direction == "UP",
      !is.na(logFC),
      abs(logFC) >= 1,
      significant_and_concordant %in% TRUE
    ) %>%
    arrange(adj.P.Val, desc(abs(logFC)), gene_symbol) %>%
    slice_head(n = 5)
  
  down_labels <- core_plot %>%
    filter(
      tested_in_tcga %in% TRUE,
      core_direction == "DOWN",
      !is.na(logFC)
    )
  
  label_tbl <- bind_rows(up_labels, down_labels) %>%
    distinct(gene_symbol, .keep_all = TRUE)
  
  readr::write_csv(core_plot, save_path("FIG4_TCGA_A_core_concordance_dotplot_data.csv"))
  readr::write_csv(label_tbl, save_path("FIG4_TCGA_A_core_concordance_dotplot_labels.csv"))
  
  # Summary barplot data: only significant + concordant
  required_summary_cols <- c("core_direction", "n_significant_and_concordant")
  missing_summary_cols <- setdiff(required_summary_cols, colnames(concordance_summary))
  if (length(missing_summary_cols) > 0) {
    stop(
      "TCGA_BRCA_concordance_summary.csv is missing required column(s): ",
      paste(missing_summary_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  summary_sigconc <- concordance_summary %>%
    filter(core_direction %in% c("UP", "DOWN")) %>%
    transmute(
      core_direction = factor(core_direction, levels = c("UP", "DOWN")),
      metric = "Significant and concordant",
      count = suppressWarnings(as.numeric(n_significant_and_concordant))
    )
  
  readr::write_csv(summary_sigconc, save_path("FIG4_TCGA_B_sigconc_barplot_data.csv"))
  
  # ---------------------------------------------------------------------------
  # 4. COMMON THEME / COLORS
  # ---------------------------------------------------------------------------
  col_up <- "#D55E00"
  col_down <- "blue"
  
  fig_theme <- theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold")
    )
  
  # ---------------------------------------------------------------------------
  # 5. PANEL A: DOT PLOT
  # ---------------------------------------------------------------------------
  cat("  [FIG4-C] Building Panel A dot plot...\n")
  
  p_a <- ggplot(
    core_plot %>% filter(tested_in_tcga %in% TRUE),
    aes(x = logFC, y = neglog10_fdr)
  ) +
    geom_vline(xintercept = -1, linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_vline(xintercept = 0,  linetype = "dashed", linewidth = 0.4, color = "black") +
    geom_vline(xintercept = 1,  linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "grey45") +
    
    geom_point(
      data = core_plot %>% filter(tested_in_tcga %in% TRUE, effect_zone == "Inside"),
      aes(x = logFC, y = neglog10_fdr),
      color = "grey65",
      size = 3,
      alpha = 0.8,
      show.legend = FALSE
    ) +
    
    geom_point(
      data = core_plot %>% filter(tested_in_tcga %in% TRUE, effect_zone != "Inside"),
      aes(x = logFC, y = neglog10_fdr, color = core_direction),
      size = 3,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    
    geom_text(
      data = label_tbl,
      aes(
        x = logFC,
        y = neglog10_fdr,
        label = gene_symbol
      ),
      size = 3,
      vjust = -0.7,
      color = "grey15",
      check_overlap = TRUE,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    
    scale_color_manual(
      values = c("UP" = col_up, "DOWN" = col_down),
      guide = "none"
    ) +
    
    labs(
      title = "TCGA-BRCA external concordance of the strict conserved core",
      x = "TCGA-BRCA log2FC (Primary Tumor vs Normal Tissue)",
      y = expression(-log[10]("adjusted p-value"))
    ) +
    
    fig_theme +
    theme(
      legend.position = "none"
    )
  
  save_plot_multi(
    p_a,
    "FIG4_TCGA_A_core_concordance_dotplot",
    width = 8.5,
    height = 6.5
  )
  
  # ---------------------------------------------------------------------------
  # 6. PANEL B: BARPLOT (ONLY SIGNIFICANT + CONCORDANT)
  # ---------------------------------------------------------------------------
  cat("  [FIG4-D] Building Panel B barplot...\n")
  
  p_b <- ggplot(summary_sigconc, aes(x = metric, y = count, fill = core_direction)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "grey25", linewidth = 0.2) +
    geom_text(
      aes(label = count),
      position = position_dodge(width = 0.7),
      vjust = -0.25,
      size = 4
    ) +
    scale_fill_manual(
      values = c("UP" = col_up, "DOWN" = col_down),
      name = "Conserved core\ndirection"
    ) +
    labs(
      title = "TCGA-BRCA significant and concordant genes",
      x = NULL,
      y = "Number of genes"
    ) +
    fig_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "right"
    ) +
    expand_limits(y = max(summary_sigconc$count, na.rm = TRUE) + 5)
  
  save_plot_multi(
    p_b,
    "FIG4_TCGA_B_sigconc_barplot",
    width = 6.5,
    height = 5
  )
  
  # ---------------------------------------------------------------------------
  # 7. TEXT SUMMARY
  # ---------------------------------------------------------------------------
  cat("  [FIG4-E] Writing figure summary...\n")
  
  n_tested <- sum(core_plot$tested_in_tcga %in% TRUE, na.rm = TRUE)
  n_sig <- sum(core_plot$significant %in% TRUE, na.rm = TRUE)
  n_conc <- sum(core_plot$directionally_concordant %in% TRUE, na.rm = TRUE)
  n_sig_conc <- sum(core_plot$significant_and_concordant %in% TRUE, na.rm = TRUE)
  
  txt <- c(
    "TCGA-BRCA FIGURE SUMMARY",
    paste("Strict conserved core genes tested:", n_tested),
    paste("Significant in TCGA-BRCA:", n_sig),
    paste("Directionally concordant:", n_conc),
    paste("Significant and concordant:", n_sig_conc),
    paste("PDF export:", save_pdf)
  )
  
  writeLines(txt, save_path("FIG4_TCGA_summary.txt"))
  
  assign("fig_tcga_core_plot_table", core_plot, envir = .GlobalEnv)
  assign("fig_tcga_summary_sigconc", summary_sigconc, envir = .GlobalEnv)
  assign("fig_tcga_panel_a", p_a, envir = .GlobalEnv)
  assign("fig_tcga_panel_b", p_b, envir = .GlobalEnv)
  
  cat("  ✓ FIGURE 4 COMPLETE\n\n")
  
  invisible(list(
    core_plot = core_plot,
    label_tbl = label_tbl,
    summary_sigconc = summary_sigconc,
    p_a = p_a,
    p_b = p_b
  ))
}