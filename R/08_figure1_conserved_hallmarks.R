################################################################################
# FIGURE 1: CROSS-SPECIES CONSERVED HALLMARK PROGRAMS
#
# Panels:
#   A. Dog conserved hallmark enrichment (GSEA; x = NES)
#   B. Cat conserved hallmark enrichment (ORA; x = GeneRatio)
#
# Uses:
#   - results/DOG_Hallmark_fgsea_simple_plus.csv
#   - results/CAT_Hallmark_ORA_UP_universeFixed.csv
#   - results/CAT_Hallmark_ORA_DOWN_universeFixed.csv
#   - results/CONSERVED_HALLMARKS_summary.csv
#
# Outputs:
#   - results/FIG1_panels/FIG1A_dog_conserved_hallmarks_dotplot.{png,pdf}
#   - results/FIG1_panels/FIG1B_cat_conserved_hallmarks_dotplot.{png,pdf}
#   - results/FIG1_tables/FIG1A_dog_conserved_hallmarks_data.csv
#   - results/FIG1_tables/FIG1B_cat_conserved_hallmarks_data.csv
################################################################################

module_figure1_hallmarks <- function() {
  
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("FIGURE 1: CROSS-SPECIES CONSERVED HALLMARK PROGRAMS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(ggplot2)
    library(tools)
  })
  
  # ---------------------------------------------------------------------------
  # 1. DIRECTORIES
  # ---------------------------------------------------------------------------
  fig_dir <- file.path(PATH_CONFIG$results_dir, "FIG1_panels")
  tab_dir <- file.path(PATH_CONFIG$results_dir, "FIG1_tables")
  
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---------------------------------------------------------------------------
  # 2. INPUT FILES
  # ---------------------------------------------------------------------------
  required_files <- c(
    file.path(PATH_CONFIG$results_dir, "DOG_Hallmark_fgsea_simple_plus.csv"),
    file.path(PATH_CONFIG$results_dir, "CAT_Hallmark_ORA_UP_universeFixed.csv"),
    file.path(PATH_CONFIG$results_dir, "CAT_Hallmark_ORA_DOWN_universeFixed.csv"),
    file.path(PATH_CONFIG$results_dir, "CONSERVED_HALLMARKS_summary.csv")
  )
  
  missing_files <- required_files[!file.exists(required_files)]
  
  if (length(missing_files) > 0) {
    stop(
      "Missing required file(s) for Figure 1:\n",
      paste(" - ", missing_files, collapse = "\n")
    )
  }
  
  cat("  [FIG1-A] Loading input files...\n")
  
  dog_fg <- read_csv(
    file.path(PATH_CONFIG$results_dir, "DOG_Hallmark_fgsea_simple_plus.csv"),
    show_col_types = FALSE
  )
  
  cat_up <- read_csv(
    file.path(PATH_CONFIG$results_dir, "CAT_Hallmark_ORA_UP_universeFixed.csv"),
    show_col_types = FALSE
  )
  
  cat_down <- read_csv(
    file.path(PATH_CONFIG$results_dir, "CAT_Hallmark_ORA_DOWN_universeFixed.csv"),
    show_col_types = FALSE
  )
  
  cons_broad <- read_csv(
    file.path(PATH_CONFIG$results_dir, "CONSERVED_HALLMARKS_summary.csv"),
    show_col_types = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 3. HELPERS
  # ---------------------------------------------------------------------------
  clean_hallmark_label <- function(x) {
    x <- gsub("^HALLMARK_", "", x)
    x <- gsub("_", " ", x)
    x <- toTitleCase(tolower(x))
    x <- gsub("E2f", "E2F", x)
    x <- gsub("G2m", "G2M", x)
    x <- gsub("Mtorc1", "MTORC1", x)
    x
  }
  
  parse_gene_ratio <- function(x) {
    sapply(x, function(z) {
      if (is.na(z) || z == "") return(NA_real_)
      parts <- strsplit(z, "/")[[1]]
      if (length(parts) != 2) return(NA_real_)
      as.numeric(parts[1]) / as.numeric(parts[2])
    })
  }
  
  save_plot_multi <- function(plot_obj, filename_base, width = 8, height = 6) {
    ggsave(
      filename = file.path(fig_dir, paste0(filename_base, ".pdf")),
      plot = plot_obj, width = width, height = height, units = "in"
    )
    ggsave(
      filename = file.path(fig_dir, paste0(filename_base, ".png")),
      plot = plot_obj, width = width, height = height, units = "in", dpi = 300
    )
  }
  
  # ---------------------------------------------------------------------------
  # 4. CONSERVED HALLMARK ORDER
  # ---------------------------------------------------------------------------
  cat("  [FIG1-B] Defining hallmark order...\n")
  
  hallmark_order <- cons_broad %>%
    distinct(hallmark, direction) %>%
    mutate(direction_order = ifelse(direction == "UP", 1, 2)) %>%
    arrange(direction_order, hallmark) %>%
    pull(hallmark)
  
  if (length(hallmark_order) == 0) {
    stop("No conserved hallmarks found in CONSERVED_HALLMARKS_summary.csv")
  }
  
  cat("      Conserved hallmarks:", length(hallmark_order), "\n")
  
  # ---------------------------------------------------------------------------
  # 5. PANEL A DATA (DOG)
  # ---------------------------------------------------------------------------
  cat("  [FIG1-C] Preparing dog panel data...\n")
  
  fig1a_df <- dog_fg %>%
    filter(pathway %in% hallmark_order) %>%
    mutate(
      hallmark = pathway,
      hallmark_clean = clean_hallmark_label(hallmark),
      direction = ifelse(NES > 0, "UP", "DOWN"),
      direction = factor(direction, levels = c("UP", "DOWN")),
      neglog10_padj = -log10(pmax(padj, 1e-300)),
      fdr_label = paste0("FDR ", formatC(padj, format = "e", digits = 1))
    )
  
  fig1a_df <- fig1a_df %>%
    mutate(hallmark_clean = factor(hallmark_clean, levels = rev(clean_hallmark_label(hallmark_order))))
  
  write_csv(fig1a_df, file.path(tab_dir, "FIG1A_dog_conserved_hallmarks_data.csv"))
  
  # ---------------------------------------------------------------------------
  # 6. PANEL B DATA (CAT)
  # ---------------------------------------------------------------------------
  cat("  [FIG1-D] Preparing cat panel data...\n")
  
  cat_up2 <- cat_up %>%
    mutate(direction = "UP") %>%
    select(ID, Description, GeneRatio, Count, p.adjust, direction)
  
  cat_down2 <- cat_down %>%
    mutate(direction = "DOWN") %>%
    select(ID, Description, GeneRatio, Count, p.adjust, direction)
  
  conserved_pairs <- cons_broad %>%
    distinct(hallmark, direction)
  
  fig1b_df <- bind_rows(cat_up2, cat_down2) %>%
    mutate(
      hallmark = ID
    ) %>%
    inner_join(conserved_pairs, by = c("hallmark", "direction")) %>%
    mutate(
      hallmark_clean = clean_hallmark_label(hallmark),
      gene_ratio_num = parse_gene_ratio(GeneRatio),
      neglog10_padj = -log10(pmax(p.adjust, 1e-300)),
      fdr_label = paste0("FDR ", formatC(p.adjust, format = "e", digits = 1))
    )
  
  fig1b_df <- fig1b_df %>%
    mutate(
      hallmark_clean = factor(
        hallmark_clean,
        levels = rev(clean_hallmark_label(hallmark_order))
      )
    )
  
  write_csv(fig1b_df, file.path(tab_dir, "FIG1B_cat_conserved_hallmarks_data.csv"))
  # ---------------------------------------------------------------------------
  # 7. PANEL A PLOT (DOG)
  # ---------------------------------------------------------------------------
  cat("  [FIG1-E] Plotting dog panel...\n")
  
  p_a <- ggplot(
    fig1a_df,
    aes(
      x = NES,
      y = hallmark_clean,
      color = direction,
      size = neglog10_padj
    )
  ) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey35") +
    geom_point() +
    geom_text(
      aes(label = fdr_label),
      nudge_x = 0.18,
      hjust = 0,
      size = 4,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c("UP" = "red", "DOWN" = "blue"),
      drop = FALSE,
      name = "Direction"
    ) +
    scale_size_continuous(name = expression(-log[10](FDR))) +
    guides(
      color = guide_legend(
        override.aes = list(size = 5)
      )
    ) +
    coord_cartesian(
      xlim = c(min(fig1a_df$NES, na.rm = TRUE) - 0.35,
               max(fig1a_df$NES, na.rm = TRUE) + 1.25),
      clip = "off"
    ) +
    labs(
      title = "A. Dog conserved hallmark enrichment",
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      plot.margin = margin(10, 95, 10, 10),
      axis.text.y = element_text(face = "bold")
    )
  
  save_plot_multi(p_a, "FIG1A_dog_conserved_hallmarks_dotplot", width = 9, height = 6)
  
  # ---------------------------------------------------------------------------
  # 8. PANEL B PLOT (CAT)
  # ---------------------------------------------------------------------------
  cat("  [FIG1-F] Plotting cat panel...\n")
  
  p_b <- ggplot(
    fig1b_df,
    aes(
      x = gene_ratio_num,
      y = hallmark_clean,
      size = Count,
      color = neglog10_padj
    )
  ) +
    geom_point() +
    geom_text(
      aes(label = fdr_label),
      nudge_x = 0.012,
      hjust = 0,
      size = 4,
      show.legend = FALSE
    ) +
    scale_color_gradient(
      low = "blue",
      high = "red",
      name = expression(-log[10](FDR))
    ) +
    scale_size_continuous(name = "Gene count") +
    coord_cartesian(
      xlim = c(min(fig1b_df$gene_ratio_num, na.rm = TRUE) - 0.02,
               max(fig1b_df$gene_ratio_num, na.rm = TRUE) + 0.08),
      clip = "off"
    ) +
    labs(
      title = "B. Cat conserved hallmark enrichment",
      x = "GeneRatio",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      plot.margin = margin(10, 65, 10, 10),
      axis.text.y = element_text(face = "bold")
    )
  
  save_plot_multi(p_b, "FIG1B_cat_conserved_hallmarks_dotplot", width = 9, height = 6)
  
  # ---------------------------------------------------------------------------
  # 9. DONE
  # ---------------------------------------------------------------------------
  cat("\n")
  cat("  ✓ FIGURE 1 COMPLETE\n\n")
  cat("  Panel files saved in:\n")
  cat("   ", normalizePath(fig_dir, winslash = "/", mustWork = FALSE), "\n")
  cat("  Plot-ready tables saved in:\n")
  cat("   ", normalizePath(tab_dir, winslash = "/", mustWork = FALSE), "\n\n")
  
  invisible(list(
    fig1a_df = fig1a_df,
    fig1b_df = fig1b_df,
    pA = p_a,
    pB = p_b
  ))
}