#R/09_figure2_overlap_support.R
################################################################################
# FIGURE 2 + TABLE 1: GENE-LEVEL OVERLAP SUPPORT ACROSS CONSERVED HALLMARKS
################################################################################

module_figure2_overlap_support <- function() {
  
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("FIGURE 2 + TABLE 1: GENE-LEVEL OVERLAP SUPPORT\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(ggplot2)
    library(tidyr)
    library(grid)
  })
  
  # ---------------------------------------------------------------------------
  # 1. CHECK PATH_CONFIG
  # ---------------------------------------------------------------------------
  if (!exists("PATH_CONFIG")) {
    stop("PATH_CONFIG not found in environment.")
  }
  if (is.null(PATH_CONFIG$results_dir) || !nzchar(PATH_CONFIG$results_dir)) {
    stop("PATH_CONFIG$results_dir is missing or empty.")
  }
  
  # ---------------------------------------------------------------------------
  # 2. DIRECTORIES
  # ---------------------------------------------------------------------------
  fig_dir <- file.path(PATH_CONFIG$results_dir, "FIG2_panels")
  tab_dir <- file.path(PATH_CONFIG$results_dir, "FIG2_tables")
  
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---------------------------------------------------------------------------
  # 3. FILE CHECKS
  # ---------------------------------------------------------------------------
  req_files <- c(
    file.path(PATH_CONFIG$results_dir, "DOG_Hallmark_fgsea_simple_plus.csv"),
    file.path(PATH_CONFIG$results_dir, "CAT_Hallmark_ORA_UP_universeFixed.csv"),
    file.path(PATH_CONFIG$results_dir, "CAT_Hallmark_ORA_DOWN_universeFixed.csv"),
    file.path(PATH_CONFIG$results_dir, "CONSERVED_HALLMARKS_summary.csv"),
    file.path(PATH_CONFIG$results_dir, "CONSERVED_HALLMARKS_strict_summary.csv"),
    file.path(PATH_CONFIG$results_dir, "CONSERVED_CORE_TARGETS_by_hallmark.csv"),
    file.path(PATH_CONFIG$results_dir, "DOG_LE_vs_CAT_UP_Fisher_byHallmark.csv"),
    file.path(PATH_CONFIG$results_dir, "DOG_LE_vs_CAT_DOWN_Fisher_byHallmark.csv")
  )
  
  missing_files <- req_files[!file.exists(req_files)]
  if (length(missing_files) > 0) {
    stop(
      "Missing required file(s) for Figure 2/Table 1:\n",
      paste(" - ", missing_files, collapse = "\n")
    )
  }
  
  # ---------------------------------------------------------------------------
  # 4. HELPERS
  # ---------------------------------------------------------------------------
  parse_gene_ratio <- function(x) {
    sapply(x, function(z) {
      if (is.na(z) || z == "") return(NA_real_)
      parts <- strsplit(z, "/")[[1]]
      if (length(parts) != 2) return(NA_real_)
      as.numeric(parts[1]) / as.numeric(parts[2])
    })
  }
  
  clean_hallmark <- function(x) {
    x %>%
      as.character() %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ")
  }
  
  # ---------------------------------------------------------------------------
  # 5. LOAD FILES
  # ---------------------------------------------------------------------------
  cat("  [FIG2-A] Loading inputs...\n")
  
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
  
  cons_strict <- read_csv(
    file.path(PATH_CONFIG$results_dir, "CONSERVED_HALLMARKS_strict_summary.csv"),
    show_col_types = FALSE
  )
  
  core_by_h <- read_csv(
    file.path(PATH_CONFIG$results_dir, "CONSERVED_CORE_TARGETS_by_hallmark.csv"),
    show_col_types = FALSE
  )
  
  fish_up <- read_csv(
    file.path(PATH_CONFIG$results_dir, "DOG_LE_vs_CAT_UP_Fisher_byHallmark.csv"),
    show_col_types = FALSE
  )
  
  fish_down <- read_csv(
    file.path(PATH_CONFIG$results_dir, "DOG_LE_vs_CAT_DOWN_Fisher_byHallmark.csv"),
    show_col_types = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 6. PREP COMMON ORDER
  # ---------------------------------------------------------------------------
  cat("  [FIG2-B] Defining hallmark order...\n")
  
  hallmark_order <- cons_broad %>%
    distinct(hallmark, direction) %>%
    mutate(direction_order = ifelse(direction == "UP", 1, 2)) %>%
    arrange(direction_order, hallmark) %>%
    pull(hallmark)
  
  if (length(hallmark_order) == 0) {
    stop("No conserved hallmarks found.")
  }
  
  hallmark_order_clean <- clean_hallmark(hallmark_order)
  
  # ---------------------------------------------------------------------------
  # 7. DOG SUMMARY
  # ---------------------------------------------------------------------------
  cat("  [FIG2-C] Building dog hallmark summary...\n")
  
  dog_tbl <- dog_fg %>%
    filter(pathway %in% hallmark_order) %>%
    transmute(
      hallmark = pathway,
      dog_NES = NES,
      dog_pval = pval,
      dog_padj = padj,
      dog_size = size
    )
  
  # ---------------------------------------------------------------------------
  # 8. CAT SUMMARY
  # ---------------------------------------------------------------------------
  cat("  [FIG2-D] Building cat hallmark summary...\n")
  
  cat_up_tbl <- cat_up %>%
    transmute(
      hallmark = ID,
      direction = "UP",
      cat_GeneRatio = GeneRatio,
      cat_GeneRatio_excel = paste0("'", GeneRatio),
      cat_GeneRatio_num = parse_gene_ratio(GeneRatio),
      cat_Count = Count,
      cat_pval = pvalue,
      cat_padj = p.adjust
    )
  
  cat_down_tbl <- cat_down %>%
    transmute(
      hallmark = ID,
      direction = "DOWN",
      cat_GeneRatio = GeneRatio,
      cat_GeneRatio_excel = paste0("'", GeneRatio),
      cat_GeneRatio_num = parse_gene_ratio(GeneRatio),
      cat_Count = Count,
      cat_pval = pvalue,
      cat_padj = p.adjust
    )
  
  cat_tbl <- bind_rows(cat_up_tbl, cat_down_tbl) %>%
    filter(hallmark %in% hallmark_order)
  
  # ---------------------------------------------------------------------------
  # 9. FISHER SUMMARY
  # ---------------------------------------------------------------------------
  cat("  [FIG2-E] Building overlap statistics table...\n")
  
  fish_up_tbl <- fish_up %>%
    transmute(
      hallmark,
      direction = "UP",
      overlap,
      A_size,
      B_size,
      U_size,
      jaccard,
      odds_ratio,
      fisher_p = p_value,
      fisher_padj = padj
    )
  
  fish_down_tbl <- fish_down %>%
    transmute(
      hallmark,
      direction = "DOWN",
      overlap,
      A_size,
      B_size,
      U_size,
      jaccard,
      odds_ratio,
      fisher_p = p_value,
      fisher_padj = padj
    )
  
  fish_tbl <- bind_rows(fish_up_tbl, fish_down_tbl) %>%
    filter(hallmark %in% hallmark_order)
  
  # ---------------------------------------------------------------------------
  # 10. LAYER SUMMARY
  # ---------------------------------------------------------------------------
  cat("  [FIG2-F] Building hallmark layer summary...\n")
  
  broad_layer <- cons_broad %>%
    transmute(hallmark, direction, broad = 1L)
  
  strict_layer <- cons_strict %>%
    filter(class == "strict_overlap") %>%
    transmute(hallmark, direction, strict_overlap = 1L)
  
  fisher_layer <- cons_strict %>%
    filter(class == "fisher_significant") %>%
    transmute(hallmark, direction, fisher_significant = 1L)
  
  layer_tbl <- cons_broad %>%
    distinct(hallmark, direction) %>%
    left_join(broad_layer, by = c("hallmark", "direction")) %>%
    left_join(strict_layer, by = c("hallmark", "direction")) %>%
    left_join(fisher_layer, by = c("hallmark", "direction")) %>%
    mutate(
      broad = coalesce(broad, 0L),
      strict_overlap = coalesce(strict_overlap, 0L),
      fisher_significant = coalesce(fisher_significant, 0L)
    )
  
  # ---------------------------------------------------------------------------
  # 11. HALLMARK-WISE CORE GENE SUMMARY
  # ---------------------------------------------------------------------------
  cat("  [FIG2-G] Building hallmark-wise conserved gene summary...\n")
  
  required_core_cols <- c("hallmark", "direction", "n_dog_le_universe", "n_conserved_core", "genes")
  missing_core_cols <- setdiff(required_core_cols, names(core_by_h))
  
  if (length(missing_core_cols) > 0) {
    stop(
      "CONSERVED_CORE_TARGETS_by_hallmark.csv is missing required columns:\n",
      paste(" - ", missing_core_cols, collapse = "\n")
    )
  }
  
  core_strict_tbl <- core_by_h %>%
    filter(hallmark %in% hallmark_order) %>%
    transmute(
      hallmark,
      direction,
      n_dog_le_universe,
      n_conserved_core,
      genes
    )
  
  # ---------------------------------------------------------------------------
  # 12. TABLE 1 BACKBONE
  # ---------------------------------------------------------------------------
  cat("  [FIG2-H] Building Table 1 backbone...\n")
  
  table1_df <- cons_broad %>%
    distinct(hallmark, direction) %>%
    left_join(dog_tbl, by = "hallmark") %>%
    left_join(cat_tbl, by = c("hallmark", "direction")) %>%
    left_join(fish_tbl, by = c("hallmark", "direction")) %>%
    left_join(layer_tbl, by = c("hallmark", "direction")) %>%
    left_join(core_strict_tbl, by = c("hallmark", "direction")) %>%
    mutate(
      broad = ifelse(is.na(broad), 0L, broad),
      strict_overlap = ifelse(is.na(strict_overlap), 0L, strict_overlap),
      fisher_significant = ifelse(is.na(fisher_significant), 0L, fisher_significant),
      hallmark_clean = clean_hallmark(hallmark)
    ) %>%
    select(
      hallmark, hallmark_clean, direction,
      dog_NES, dog_padj, dog_size,
      cat_GeneRatio, cat_GeneRatio_excel, cat_GeneRatio_num, cat_Count, cat_padj,
      overlap, jaccard, odds_ratio, fisher_p, fisher_padj,
      n_dog_le_universe, n_conserved_core,
      broad, strict_overlap, fisher_significant,
      genes
    ) %>%
    arrange(match(hallmark, hallmark_order))
  
  write_csv(table1_df, file.path(PATH_CONFIG$results_dir, "TABLE1_conserved_hallmark_programs.csv"))
  write_csv(table1_df, file.path(tab_dir, "TABLE1_conserved_hallmark_programs.csv"))
  
  # ---------------------------------------------------------------------------
  # 13. FIGURE 2 PANEL DATA
  # ---------------------------------------------------------------------------
  cat("  [FIG2-I] Writing panel-ready data tables...\n")
  
  fig2a_df <- table1_df %>%
    transmute(
      hallmark,
      hallmark_clean,
      direction,
      fisher_padj,
      neglog10_fisher_padj = -log10(pmax(fisher_padj, 1e-300)),
      overlap
    ) %>%
    mutate(
      hallmark_clean = factor(hallmark_clean, levels = rev(hallmark_order_clean))
    )
  
  fig2b_df <- table1_df %>%
    transmute(
      hallmark,
      hallmark_clean,
      direction,
      overlap,
      jaccard,
      fisher_padj
    ) %>%
    mutate(
      hallmark_clean = factor(hallmark_clean, levels = rev(hallmark_order_clean))
    )
  
  fig2c_df <- table1_df %>%
    transmute(
      hallmark,
      hallmark_clean,
      direction,
      n_dog_le_universe,
      n_conserved_core
    ) %>%
    pivot_longer(
      cols = c(n_dog_le_universe, n_conserved_core),
      names_to = "metric",
      values_to = "n_genes"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c("n_dog_le_universe", "n_conserved_core"),
        labels = c("Dog leading-edge", "Conserved core")
      ),
      hallmark_clean = factor(hallmark_clean, levels = rev(hallmark_order_clean))
    )
  

  write_csv(fig2a_df, file.path(tab_dir, "FIG2A_fisher_support_data.csv"))
  write_csv(fig2b_df, file.path(tab_dir, "FIG2B_overlap_jaccard_data.csv"))
  write_csv(fig2c_df, file.path(tab_dir, "FIG2C_gene_count_support_data.csv"))
  
  # ---------------------------------------------------------------------------
  # 14. THEME
  # ---------------------------------------------------------------------------
  fig_theme <- theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "top",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold"),
      legend.spacing.x = unit(6, "pt"),
      legend.key.width = unit(0.8, "cm"),
      legend.key.height = unit(0.5, "cm"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      axis.text.y = element_text(face = "bold")
    )
  
  # ---------------------------------------------------------------------------
  # 15. PANEL A: Fisher support
  # ---------------------------------------------------------------------------
  cat("  [FIG2-J] Plotting Panel A...\n")
  
  pA <- ggplot(fig2a_df, aes(x = neglog10_fisher_padj, y = hallmark_clean, color = direction, size = overlap)) +
    geom_point(alpha = 0.9) +
    scale_color_manual(values = c("UP" = "#D55E00", "DOWN" = "#0072B2"), name = "Direction") +
    scale_size_continuous(name = "Overlap size", range = c(4, 10)) +
    guides(
      size = guide_legend(order = 1),
      color = guide_legend(
        order = 2,
        override.aes = list(size = 5)
      )
    ) +
    labs(
      title = "A. Fisher-supported overlap across conserved hallmarks",
      x = expression(-log[10]("adjusted Fisher p-value")),
      y = NULL
    ) +
    fig_theme +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(15, 50, 15, 15)) +
    expand_limits(x = max(fig2a_df$neglog10_fisher_padj, na.rm = TRUE) * 1.12)
  
  ggsave(file.path(fig_dir, "FIG2A_fisher_support.png"), pA, width = 9.2, height = 5.0, dpi = 300)
  ggsave(file.path(fig_dir, "FIG2A_fisher_support.pdf"), pA, width = 9.2, height = 5.0)
  
  # ---------------------------------------------------------------------------
  # 16. PANEL B: Jaccard / overlap
  # ---------------------------------------------------------------------------
  cat("  [FIG2-K] Plotting Panel B...\n")
  
  pB <- ggplot(fig2b_df, aes(x = jaccard, y = hallmark_clean, color = direction, size = overlap)) +
    geom_point(alpha = 0.9) +
    scale_color_manual(values = c("UP" = "#D55E00", "DOWN" = "#0072B2"), name = "Direction") +
    scale_size_continuous(name = "Overlap size", range = c(4, 10)) +
    guides(
      size = guide_legend(order = 1),
      color = guide_legend(
        order = 2,
        override.aes = list(size = 5)
      )
    ) +
    labs(
      title = "B. Overlap magnitude across conserved hallmarks",
      x = "Jaccard index",
      y = NULL
    ) +
    fig_theme +
    theme(plot.margin = margin(15, 50, 15, 15))
  
  ggsave(file.path(fig_dir, "FIG2B_jaccard_overlap.png"), pB, width = 9.4, height = 5.0, dpi = 300)
  ggsave(file.path(fig_dir, "FIG2B_jaccard_overlap.pdf"), pB, width = 9.4, height = 5.0)
  
  # ---------------------------------------------------------------------------
  # 17. PANEL C: Dog LE vs conserved core counts
  # ---------------------------------------------------------------------------
  cat("  [FIG2-L] Plotting Panel C...\n")
  
  pC <- ggplot(fig2c_df, aes(x = metric, y = n_genes, fill = direction)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ hallmark_clean, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = c("UP" = "#D55E00", "DOWN" = "#0072B2"), name = "Direction") +
    labs(
      title = "C. Hallmark-wise conserved gene support",
      x = NULL,
      y = "Number of genes"
    ) +
    fig_theme +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(size = 9)
    )
  
  ggsave(file.path(fig_dir, "FIG2C_gene_support_counts.png"), pC, width = 10.0, height = 6.5, dpi = 300)
  ggsave(file.path(fig_dir, "FIG2C_gene_support_counts.pdf"), pC, width = 10.0, height = 6.5)
  
  # ---------------------------------------------------------------------------
  # 18. DONE
  # ---------------------------------------------------------------------------
  cat("\n")
  cat("  ✓ FIGURE 2 + TABLE 1 COMPLETE\n\n")
  cat("  Main table:\n")
  cat("   ", file.path(PATH_CONFIG$results_dir, "TABLE1_conserved_hallmark_programs.csv"), "\n")
  cat("  Panels:\n")
  cat("   ", normalizePath(fig_dir, winslash = "/", mustWork = FALSE), "\n")
  cat("  Plot-ready tables:\n")
  cat("   ", normalizePath(tab_dir, winslash = "/", mustWork = FALSE), "\n\n")
  
  invisible(list(
    table1_df = table1_df,
    fig2a_df = fig2a_df,
    fig2b_df = fig2b_df,
    fig2c_df = fig2c_df
  ))
}