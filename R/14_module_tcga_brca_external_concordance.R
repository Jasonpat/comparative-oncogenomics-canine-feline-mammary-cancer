#14_module_tcga_brca_external_concordance_local.R
################################################################################
# STEP 14: TCGA-BRCA EXTERNAL HUMAN CONCORDANCE (LOCAL FILE VERSION)
################################################################################

module_tcga_brca_external_concordance_local <- function(
    tcga_dir = NULL,
    use_cached_query = FALSE
) {
  
  # Allows the function to run either as:
  # module_tcga_brca_external_concordance_local(tcga_dir = "...")
  # or directly from the pipeline if PATH_CONFIG$tcga_dir exists.
  if (is.null(tcga_dir) && exists("PATH_CONFIG", envir = .GlobalEnv)) {
    tcga_dir <- get("PATH_CONFIG", envir = .GlobalEnv)$tcga_dir
  }
  
  if (is.null(tcga_dir) || length(tcga_dir) == 0 || is.na(tcga_dir) || !nzchar(tcga_dir)) {
    stop("tcga_dir must be provided or PATH_CONFIG$tcga_dir must exist.", call. = FALSE)
  }
  
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("STEP 14: TCGA-BRCA EXTERNAL HUMAN CONCORDANCE (LOCAL FILE VERSION)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  suppressPackageStartupMessages({
    library(TCGAbiolinks)
    library(readr)
    library(dplyr)
    library(stringr)
    library(tibble)
    library(edgeR)
    library(limma)
    library(tools)
  })
  
  # ---------------------------------------------------------------------------
  # 1. PATHS / HELPERS
  # ---------------------------------------------------------------------------
  if (exists("PATH_CONFIG", envir = .GlobalEnv)) {
    results_dir <- get("PATH_CONFIG", envir = .GlobalEnv)$results_dir
    cache_dir   <- get("PATH_CONFIG", envir = .GlobalEnv)$cache_dir
  } else {
    results_dir <- "results"
    cache_dir   <- "cache"
  }
  
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  save_result_local <- function(obj, filename, type = c("csv", "rds", "txt")) {
    type <- match.arg(type)
    out <- file.path(results_dir, filename)
    
    if (type == "csv") readr::write_csv(obj, out)
    if (type == "rds") saveRDS(obj, out)
    if (type == "txt") writeLines(as.character(obj), out)
    
    cat("  ✓ Saved: ", normalizePath(out, winslash = "/", mustWork = FALSE), "\n", sep = "")
    invisible(out)
  }
  
  save_any <- function(obj, filename, type = c("csv", "rds", "txt")) {
    type <- match.arg(type)
    if (exists("save_result", envir = .GlobalEnv)) {
      get("save_result", envir = .GlobalEnv)(obj, filename, type)
    } else {
      save_result_local(obj, filename, type)
    }
  }
  
  up_file   <- file.path(results_dir, "Conserved_Core_UP_genes_strict.txt")
  down_file <- file.path(results_dir, "Conserved_Core_DOWN_genes_strict.txt")
  
  if (!file.exists(up_file)) {
    stop("Missing required file: ", up_file, call. = FALSE)
  }
  
  if (!dir.exists(tcga_dir)) {
    stop("TCGA directory does not exist: ", tcga_dir, call. = FALSE)
  }
  
  # ---------------------------------------------------------------------------
  # 2. LOAD STRICT CORE GENES
  # ---------------------------------------------------------------------------
  cat("  [14A] Loading strict conserved core gene sets...\n")
  
  up_genes <- readLines(up_file, warn = FALSE)
  up_genes <- unique(trimws(up_genes))
  up_genes <- up_genes[up_genes != ""]
  
  if (file.exists(down_file)) {
    down_genes <- readLines(down_file, warn = FALSE)
    down_genes <- unique(trimws(down_genes))
    down_genes <- down_genes[down_genes != ""]
  } else {
    down_genes <- character(0)
    cat("      NOTE: No DOWN strict core file found. Continuing with UP genes only.\n")
  }
  
  core_tbl <- bind_rows(
    tibble(gene_symbol = up_genes, core_direction = "UP"),
    tibble(gene_symbol = down_genes, core_direction = "DOWN")
  ) %>%
    distinct(gene_symbol, .keep_all = TRUE)
  
  save_any(core_tbl, "TCGA_BRCA_core_gene_input_table.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # 3. QUERY METADATA FROM GDC
  # ---------------------------------------------------------------------------
  cat("  [14B] Querying GDC metadata for TCGA-BRCA...\n")
  
  query_cache <- file.path(cache_dir, "TCGA_BRCA_query_results.rds")
  
  if (use_cached_query && file.exists(query_cache)) {
    meta <- readRDS(query_cache)
  } else {
    query <- GDCquery(
      project = "TCGA-BRCA",
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c("Primary Tumor", "Solid Tissue Normal")
    )
    meta <- getResults(query)
    saveRDS(meta, query_cache)
  }
  
  required_meta_cols <- c("file_id", "file_name", "sample_type", "cases.submitter_id")
  missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))
  if (length(missing_meta_cols) > 0) {
    stop(
      "Missing required metadata column(s): ",
      paste(missing_meta_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # 4. FIND LOCAL FILES
  # ---------------------------------------------------------------------------
  cat("  [14C] Finding local STAR-count files...\n")
  
  files <- list.files(
    tcga_dir,
    pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    stop("No STAR-count TSV files found under: ", tcga_dir, call. = FALSE)
  }
  
  local_tbl <- tibble(
    full_path = files,
    file_id = basename(dirname(files)),
    local_file_name = basename(files)
  )
  
  meta2 <- meta %>%
    filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal")) %>%
    select(file_id, file_name, sample_type, cases.submitter_id) %>%
    distinct()
  
  sample_tbl <- meta2 %>%
    inner_join(local_tbl, by = "file_id")
  
  if (nrow(sample_tbl) == 0) {
    stop("No overlap between local files and queried TCGA metadata.", call. = FALSE)
  }
  
  cat("      Local files matched to metadata:", nrow(sample_tbl), "\n")
  print(table(sample_tbl$sample_type))
  
  save_any(sample_tbl, "TCGA_BRCA_local_sample_table.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # 5. READ ONE FILE TO DEFINE GENE REFERENCE
  # ---------------------------------------------------------------------------
  cat("  [14D] Reading gene reference from first file...\n")
  
  first_df <- readr::read_tsv(
    sample_tbl$full_path[1],
    skip = 1,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  needed_cols <- c("gene_id", "gene_name", "gene_type", "unstranded")
  missing_needed <- setdiff(needed_cols, colnames(first_df))
  if (length(missing_needed) > 0) {
    stop(
      "Missing required column(s) in STAR-count file: ",
      paste(missing_needed, collapse = ", "),
      call. = FALSE
    )
  }
  
  gene_ref <- first_df %>%
    filter(!str_starts(gene_id, "^N_")) %>%
    transmute(
      gene_id = as.character(gene_id),
      gene_id_clean = sub("\\..*$", "", gene_id),
      gene_symbol = as.character(gene_name),
      gene_type = as.character(gene_type)
    )
  
  # ---------------------------------------------------------------------------
  # 6. BUILD COUNT MATRIX
  # ---------------------------------------------------------------------------
  cat("  [14E] Building count matrix from local files...\n")
  
  read_one_counts <- function(path) {
    df <- readr::read_tsv(
      path,
      skip = 1,
      show_col_types = FALSE,
      progress = FALSE
    )
    
    df %>%
      filter(!str_starts(gene_id, "^N_")) %>%
      transmute(
        gene_id = as.character(gene_id),
        unstranded = as.numeric(unstranded)
      )
  }
  
  count_list <- lapply(sample_tbl$full_path, read_one_counts)
  
  ref_gene_id <- count_list[[1]]$gene_id
  same_order <- vapply(count_list, function(x) identical(x$gene_id, ref_gene_id), logical(1))
  if (!all(same_order)) {
    stop("Not all files have identical gene_id order. Need more defensive merging.", call. = FALSE)
  }
  
  count_mat <- do.call(
    cbind,
    lapply(count_list, function(x) x$unstranded)
  )
  
  rownames(count_mat) <- ref_gene_id
  colnames(count_mat) <- paste0(
    sample_tbl$cases.submitter_id,
    "__",
    ifelse(sample_tbl$sample_type == "Primary Tumor", "Tumor", "Normal"),
    "__",
    sample_tbl$file_id
  )
  
  # ---------------------------------------------------------------------------
  # 7. DIFFERENTIAL EXPRESSION: TUMOR VS NORMAL
  # ---------------------------------------------------------------------------
  cat("  [14F] Running limma-voom differential expression...\n")
  
  group <- factor(sample_tbl$sample_type, levels = c("Solid Tissue Normal", "Primary Tumor"))
  
  dge <- edgeR::DGEList(counts = count_mat)
  keep <- edgeR::filterByExpr(dge, group = group)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group)
  v <- limma::voom(dge, design = design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  
  de_tbl <- limma::topTable(
    fit,
    coef = "groupPrimary Tumor",
    number = Inf,
    sort.by = "P"
  ) %>%
    tibble::rownames_to_column("gene_id") %>%
    as_tibble() %>%
    left_join(gene_ref, by = "gene_id") %>%
    mutate(
      tcga_direction = case_when(
        logFC > 0 ~ "UP",
        logFC < 0 ~ "DOWN",
        TRUE ~ "UNCHANGED"
      ),
      significant = adj.P.Val < 0.05
    ) %>%
    arrange(adj.P.Val, desc(abs(logFC)))
  
  save_any(de_tbl, "TCGA_BRCA_DE_full_table.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # 8. COLLAPSE TO BEST SYMBOL ENTRY
  # ---------------------------------------------------------------------------
  cat("  [14G] Collapsing TCGA results to best symbol-level entries...\n")
  
  de_symbol_best <- de_tbl %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    arrange(adj.P.Val, desc(abs(logFC)), gene_symbol) %>%
    distinct(gene_symbol, .keep_all = TRUE)
  
  core_results <- core_tbl %>%
    left_join(de_symbol_best, by = "gene_symbol") %>%
    mutate(
      tested_in_tcga = !is.na(logFC),
      directionally_concordant = case_when(
        core_direction == "UP" & logFC > 0 ~ TRUE,
        core_direction == "DOWN" & logFC < 0 ~ TRUE,
        core_direction %in% c("UP", "DOWN") & !is.na(logFC) ~ FALSE,
        TRUE ~ NA
      ),
      significant_and_concordant = significant %in% TRUE & directionally_concordant %in% TRUE
    ) %>%
    arrange(
      core_direction,
      desc(significant_and_concordant),
      adj.P.Val,
      desc(abs(logFC))
    )
  
  save_any(core_results, "TCGA_BRCA_core_gene_results.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # 9. SUMMARY
  # ---------------------------------------------------------------------------
  cat("  [14H] Summarising concordance...\n")
  
  summary_tbl <- core_results %>%
    group_by(core_direction) %>%
    summarise(
      n_core_genes = n(),
      n_tested_in_tcga = sum(tested_in_tcga %in% TRUE, na.rm = TRUE),
      n_significant = sum(significant %in% TRUE, na.rm = TRUE),
      n_directionally_concordant = sum(directionally_concordant %in% TRUE, na.rm = TRUE),
      n_significant_and_concordant = sum(significant_and_concordant %in% TRUE, na.rm = TRUE),
      .groups = "drop"
    )
  
  overall_summary <- tibble(
    core_direction = "ALL",
    n_core_genes = nrow(core_results),
    n_tested_in_tcga = sum(core_results$tested_in_tcga %in% TRUE, na.rm = TRUE),
    n_significant = sum(core_results$significant %in% TRUE, na.rm = TRUE),
    n_directionally_concordant = sum(core_results$directionally_concordant %in% TRUE, na.rm = TRUE),
    n_significant_and_concordant = sum(core_results$significant_and_concordant %in% TRUE, na.rm = TRUE)
  )
  
  concordance_summary <- bind_rows(summary_tbl, overall_summary)
  save_any(concordance_summary, "TCGA_BRCA_concordance_summary.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # 10. HEATMAP-FRIENDLY DATA
  # ---------------------------------------------------------------------------
  cat("  [14I] Exporting core-gene expression matrix...\n")
  
  voom_mat <- as.data.frame(v$E) %>%
    tibble::rownames_to_column("gene_id") %>%
    left_join(gene_ref, by = "gene_id") %>%
    filter(gene_symbol %in% core_tbl$gene_symbol) %>%
    select(gene_symbol, everything(), -gene_id_clean, -gene_type, -gene_id) %>%
    group_by(gene_symbol) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop") %>%
    left_join(core_tbl, by = "gene_symbol") %>%
    relocate(core_direction, .after = gene_symbol)
  
  save_any(voom_mat, "TCGA_BRCA_core_heatmap_data.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # 11. TEXT SUMMARY
  # ---------------------------------------------------------------------------
  txt <- c(
    "TCGA-BRCA EXTERNAL HUMAN CONCORDANCE SUMMARY",
    paste("Tumor samples:", sum(sample_tbl$sample_type == "Primary Tumor")),
    paste("Normal samples:", sum(sample_tbl$sample_type == "Solid Tissue Normal")),
    paste("Matched local STAR-count files:", nrow(sample_tbl)),
    paste("Strict conserved core genes tested:", nrow(core_results)),
    paste("Genes significant in TCGA-BRCA:", sum(core_results$significant %in% TRUE, na.rm = TRUE)),
    paste("Genes directionally concordant:", sum(core_results$directionally_concordant %in% TRUE, na.rm = TRUE)),
    paste("Genes significant and concordant:", sum(core_results$significant_and_concordant %in% TRUE, na.rm = TRUE))
  )
  
  save_any(txt, "TCGA_BRCA_summary.txt", "txt")
  
  # ---------------------------------------------------------------------------
  # 12. EXPORT OBJECTS
  # ---------------------------------------------------------------------------
  assign("tcga_brca_sample_table", sample_tbl, envir = .GlobalEnv)
  assign("tcga_brca_de", de_tbl, envir = .GlobalEnv)
  assign("tcga_brca_core_results", core_results, envir = .GlobalEnv)
  assign("tcga_brca_concordance_summary", concordance_summary, envir = .GlobalEnv)
  
  cat("  ✓ STEP 14 COMPLETE\n\n")
  
  invisible(list(
    sample_table = sample_tbl,
    de_table = de_tbl,
    core_results = core_results,
    concordance_summary = concordance_summary,
    heatmap_data = voom_mat
  ))
}
