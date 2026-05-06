################################################################################
# MODULE 17: FIGURE 5 - TCGA-BRCA SIGNATURE VALIDATION — NO COX VERSION
#
# Reads outputs from:
#   16_module_tcga_brca_survival_signature_validation.R
#
# Main Figure 5 panels:
#   A. Conserved-core ssGSEA score by PAM50 subtype
#   B. Kaplan-Meier overall survival, High vs Low conserved-core signature
#
# This version intentionally removes:
#   - Cox forest panel.
#   - Cox source data.
#   - Secondary endpoint panels.
#   - PAM50-adjusted sensitivity figures.
#   - PH diagnostic figures.
#
# Expected Module 16 outputs:
#   TCGA_signature_dataset.csv
#   TCGA_KM_fit.rds
#   TCGA_KM_summary.csv
#   TCGA_subtype_score_association.csv
################################################################################

figure5_tcga_survival_signature_validation <- function(
    save_pdf = TRUE,
    dpi = 300
) {
  
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("FIGURE 5: TCGA-BRCA SIGNATURE VALIDATION — NO COX\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tibble)
    library(stringr)
    library(ggplot2)
    library(survival)
  })
  
  if (!requireNamespace("survminer", quietly = TRUE)) {
    stop("survminer is required for Figure 5B. Install it before running Module 17.", call. = FALSE)
  }
  
  # ---------------------------------------------------------------------------
  # 1. Paths and helpers
  # ---------------------------------------------------------------------------
  results_dir <- if (exists("PATH_CONFIG", envir = .GlobalEnv)) {
    get("PATH_CONFIG", envir = .GlobalEnv)$results_dir
  } else {
    "results"
  }
  
  surv_dir <- file.path(results_dir, "SURVIVAL")
  fig_dir <- file.path(results_dir, "FIG5_panels")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  
  read_surv_csv <- function(filename, required = FALSE) {
    path <- file.path(surv_dir, filename)
    if (!file.exists(path)) {
      if (isTRUE(required)) stop("Missing required Figure 5 input: ", path, call. = FALSE)
      cat("  NOTE: missing optional input: ", filename, "\n", sep = "")
      return(NULL)
    }
    readr::read_csv(path, show_col_types = FALSE)
  }
  
  read_surv_rds <- function(filename, required = FALSE) {
    path <- file.path(surv_dir, filename)
    if (!file.exists(path)) {
      if (isTRUE(required)) stop("Missing required Figure 5 input: ", path, call. = FALSE)
      cat("  NOTE: missing optional input: ", filename, "\n", sep = "")
      return(NULL)
    }
    readRDS(path)
  }
  
  save_fig <- function(plot_obj, filename_base, width = 7, height = 5) {
    png_path <- file.path(fig_dir, paste0(filename_base, ".png"))
    pdf_path <- file.path(fig_dir, paste0(filename_base, ".pdf"))
    
    ggplot2::ggsave(
      filename = png_path,
      plot = plot_obj,
      width = width,
      height = height,
      dpi = dpi,
      device = "png",
      bg = "white"
    )
    cat("  ✓ Saved: ", normalizePath(png_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
    
    if (isTRUE(save_pdf)) {
      tryCatch(
        {
          ggplot2::ggsave(
            filename = pdf_path,
            plot = plot_obj,
            width = width,
            height = height,
            device = grDevices::cairo_pdf,
            bg = "white"
          )
          cat("  ✓ Saved: ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
        },
        error = function(e) {
          warning("Could not save PDF for ", filename_base, ": ", conditionMessage(e), call. = FALSE)
        }
      )
    }
    
    invisible(png_path)
  }
  
  save_survminer <- function(ggsurv_obj, filename_base, width = 7, height = 7) {
    png_path <- file.path(fig_dir, paste0(filename_base, ".png"))
    pdf_path <- file.path(fig_dir, paste0(filename_base, ".pdf"))
    
    grDevices::png(png_path, width = width * dpi, height = height * dpi, res = dpi)
    print(ggsurv_obj)
    grDevices::dev.off()
    cat("  ✓ Saved: ", normalizePath(png_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
    
    if (isTRUE(save_pdf)) {
      tryCatch(
        {
          grDevices::cairo_pdf(pdf_path, width = width, height = height)
          print(ggsurv_obj)
          grDevices::dev.off()
          cat("  ✓ Saved: ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
        },
        error = function(e) {
          if (grDevices::dev.cur() != 1) grDevices::dev.off()
          warning("Could not save PDF for ", filename_base, ": ", conditionMessage(e), call. = FALSE)
        }
      )
    }
    
    invisible(png_path)
  }
  
  theme_fig <- function() {
    ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        axis.title = ggplot2::element_text(size = 10),
        axis.text = ggplot2::element_text(size = 9),
        legend.title = ggplot2::element_text(size = 9, face = "bold"),
        legend.text = ggplot2::element_text(size = 8),
        plot.margin = ggplot2::margin(8, 8, 8, 8)
      )
  }
  
  # ---------------------------------------------------------------------------
  # 2. Load Module 16 outputs
  # ---------------------------------------------------------------------------
  cat("  [17A] Loading Module 16 no-Cox outputs...\n")
  
  df <- read_surv_csv("TCGA_signature_dataset.csv", required = TRUE)
  km_fit <- read_surv_rds("TCGA_KM_fit.rds", required = FALSE)
  km_summary <- read_surv_csv("TCGA_KM_summary.csv", required = FALSE)
  subtype_tbl <- read_surv_csv("TCGA_subtype_score_association.csv", required = FALSE)
  
  if (!"group" %in% names(df)) stop("TCGA_signature_dataset.csv must contain group.", call. = FALSE)
  if (!all(c("OS.time", "OS") %in% names(df))) {
    stop("TCGA_signature_dataset.csv must contain OS.time and OS.", call. = FALSE)
  }
  
  score_col <- dplyr::case_when(
    "conserved_core_ssgsea" %in% names(df) ~ "conserved_core_ssgsea",
    "ssGSEA" %in% names(df) ~ "ssGSEA",
    TRUE ~ NA_character_
  )
  if (is.na(score_col)) {
    stop("Could not find conserved-core score column in TCGA_signature_dataset.csv.", call. = FALSE)
  }
  
  id_col <- dplyr::case_when(
    "patient_id" %in% names(df) ~ "patient_id",
    "barcode" %in% names(df) ~ "barcode",
    TRUE ~ NA_character_
  )
  if (is.na(id_col)) id_col <- NULL
  
  df <- df %>%
    dplyr::mutate(
      group = factor(group, levels = c("Low", "High")),
      OS = suppressWarnings(as.numeric(OS)),
      OS.time = suppressWarnings(as.numeric(OS.time))
    )
  
  # ---------------------------------------------------------------------------
  # 3. Panel A: ssGSEA by PAM50 subtype
  # ---------------------------------------------------------------------------
  cat("  [17B] Building Panel A PAM50 subtype plot...\n")
  
  pam50_colours <- c(
    "LumA" = "#4575B4",
    "LumB" = "#74ADD1",
    "Her2" = "#FEE090",
    "Basal" = "#D73027",
    "Normal" = "#A6D96A"
  )
  
  if ("PAM50" %in% names(df) && sum(!is.na(df$PAM50) & df$PAM50 != "") >= 20) {
    pam_df <- df %>%
      dplyr::filter(!is.na(PAM50), PAM50 != "", !is.na(.data[[score_col]])) %>%
      dplyr::mutate(
        PAM50 = factor(PAM50, levels = c("LumA", "LumB", "Normal", "Her2", "Basal")),
        signature_score = .data[[score_col]]
      )
    
    if (!is.null(subtype_tbl) && "p_value" %in% names(subtype_tbl) && nrow(subtype_tbl) > 0) {
      kw_p <- subtype_tbl$p_value[1]
    } else {
      kw_p <- stats::kruskal.test(signature_score ~ PAM50, data = pam_df)$p.value
    }
    
    fig5a <- ggplot2::ggplot(pam_df, ggplot2::aes(x = PAM50, y = signature_score, fill = PAM50)) +
      ggplot2::geom_violin(alpha = 0.65, trim = TRUE, colour = "grey40", linewidth = 0.4) +
      ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", colour = "grey20", linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = pam50_colours, guide = "none") +
      ggplot2::labs(
        x = "PAM50 subtype",
        y = "Conserved-core ssGSEA score",
        title = "Conserved-core signature score by PAM50 subtype"
      ) +
      theme_fig()
    
    save_fig(fig5a, "Figure5A_subtype_violin", width = 7, height = 5)
  } else {
    pam_df <- NULL
    kw_p <- NA_real_
    cat("      PAM50 data unavailable or insufficient; skipping Panel A.\n")
  }
  
  # ---------------------------------------------------------------------------
  # 4. Panel B: Kaplan-Meier overall survival
  # ---------------------------------------------------------------------------
  cat("  [17C] Building Panel B Kaplan-Meier plot...\n")
  
  if (is.null(km_fit)) {
    km_fit <- survival::survfit(survival::Surv(OS.time, OS) ~ group, data = df)
  }
  
  if (!is.null(km_summary) && "logrank_p" %in% names(km_summary) && nrow(km_summary) > 0) {
    lr_p <- km_summary$logrank_p[1]
  } else {
    lr <- survival::survdiff(survival::Surv(OS.time, OS) ~ group, data = df)
    lr_p <- 1 - stats::pchisq(lr$chisq, df = length(lr$n) - 1)
  }
  
  table_theme <- survminer::theme_cleantable() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 9, hjust = 0),
      axis.title.x = ggplot2::element_text(size = 8),
      axis.text.x  = ggplot2::element_text(size = 8, colour = "black"),
      axis.text.y  = ggplot2::element_text(size = 8),
      axis.ticks.x = ggplot2::element_line(colour = "black"),
      axis.line.x  = ggplot2::element_line(colour = "black")
    )
  
  km_plot <- survminer::ggsurvplot(
    km_fit,
    data = df,
    pval = FALSE,
    conf.int = TRUE,
    risk.table = "absolute",
    risk.table.title = "Number at risk",
    risk.table.height = 0.18,
    risk.table.col = "strata",
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    break.time.by = 2000,
    palette = c("#4575B4", "#D73027"),
    legend.title = "Conserved-core signature",
    legend.labs = c("Low", "High"),
    xlab = "Time (days)",
    ylab = "Overall survival probability",
    title = "TCGA-BRCA conserved-core signature and overall survival",
    ggtheme = theme_fig(),
    tables.theme = table_theme
  )
  
  save_survminer(km_plot, "Figure5B_KM", width = 7, height = 7)
  
  # ---------------------------------------------------------------------------
  # 5. Source data and finish
  # ---------------------------------------------------------------------------
  cat("  [17D] Saving Figure 5 source data...\n")
  
  source_cols <- c(
    id_col,
    "conserved_core_ssgsea", "conserved_core_ssgsea_z",
    "group",
    "OS.time", "OS",
    "PAM50", "ER_status",
    "stage", "stage_clean",
    "clinical_grade", "subtype_grade",
    "age_at_diagnosis", "age_years", "age_10y"
  )
  source_cols <- unique(source_cols[!is.na(source_cols) & source_cols %in% names(df)])
  
  source_data <- df %>% dplyr::select(dplyr::all_of(source_cols))
  readr::write_csv(source_data, file.path(fig_dir, "Figure5_source_data.csv"))
  cat("  ✓ Saved: ", normalizePath(file.path(fig_dir, "Figure5_source_data.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  
  if (exists("pam_df", inherits = FALSE) && !is.null(pam_df) && nrow(pam_df) > 0) {
    pam_source <- pam_df %>%
      dplyr::mutate(kruskal_wallis_p = kw_p)
    readr::write_csv(pam_source, file.path(fig_dir, "Figure5A_subtype_violin_source_data.csv"))
    cat("  ✓ Saved: ", normalizePath(file.path(fig_dir, "Figure5A_subtype_violin_source_data.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  }
  
  if (!is.null(km_summary) && nrow(km_summary) > 0) {
    readr::write_csv(km_summary, file.path(fig_dir, "Figure5B_KM_source_stats.csv"))
    cat("  ✓ Saved: ", normalizePath(file.path(fig_dir, "Figure5B_KM_source_stats.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  }
  
  cat("\n")
  cat("  ✓ FIGURE 5 NO COX COMPLETE\n\n")
  cat("  Output directory:\n")
  cat("   - ", normalizePath(fig_dir, winslash = "/", mustWork = FALSE), "\n\n", sep = "")
  
  invisible(list(
    dataset = df,
    km_fit = km_fit,
    subtype_plot_data = if (exists("pam_df", inherits = FALSE)) pam_df else NULL,
    subtype_summary = subtype_tbl
  ))
}
