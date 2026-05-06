################################################################################
# MODULE 16: TCGA-BRCA PATIENT-LEVEL SIGNATURE VALIDATION — NO COX VERSION
#
# Purpose:
#   Validates the strict conserved UP-core as a patient-level TCGA-BRCA tumour
#   signature using local STAR-count files, ssGSEA scoring, Kaplan-Meier OS
#   visualization/statistics, and PAM50 subtype annotation.
#
# This version intentionally removes all Cox analyses from the main workflow.
#
# Main analyses retained:
#   - Local TCGA STAR-count matching to GDC metadata.
#   - Patient-level primary tumour expression matrix construction.
#   - TMM + limma-voom normalization.
#   - ssGSEA score of the strict conserved UP-core.
#   - Gene-overlap check.
#   - Kaplan-Meier overall survival analysis using median signature split.
#   - PAM50 subtype association using Kruskal-Wallis test.
#   - Stage-cleaning QC for metadata transparency.
#
# Analyses intentionally removed:
#   - Cox proportional hazards models.
#   - PH assumption checks.
#   - Individual target-level Cox models.
#   - Secondary endpoint Cox analyses.
#   - PAM50-adjusted Cox sensitivity models.
#   - DOWN-core or UP-minus-DOWN exploratory survival models.
#
# Outputs written to results/SURVIVAL/:
#   TCGA_survival_local_sample_table.csv
#   TCGA_signature_dataset.csv
#   TCGA_signature_gene_overlap.csv
#   TCGA_stage_clean_qc.csv
#   TCGA_KM_fit.rds
#   TCGA_KM_summary.csv
#   TCGA_subtype_score_association.csv
#   TCGA_survival_method_parameters.csv
#   TCGA_survival_summary.txt
################################################################################

module_tcga_brca_survival_signature_validation <- function(
    tcga_dir = NULL,
    use_cached_query = TRUE,
    score_method = "ssgsea",
    group_cutoff = "median",
    seed = 1
) {
  
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("MODULE 16: TCGA-BRCA PATIENT-LEVEL SIGNATURE VALIDATION — NO COX\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tibble)
    library(tidyr)
    library(stringr)
    library(edgeR)
    library(limma)
    library(survival)
  })
  
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("TCGAbiolinks is required for TCGA metadata/survival/subtype annotation.", call. = FALSE)
  }
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA is required for ssGSEA scoring. Install it before running Module 16.", call. = FALSE)
  }
  
  set.seed(seed)
  
  # ---------------------------------------------------------------------------
  # 1. Paths and helpers
  # ---------------------------------------------------------------------------
  results_dir <- if (exists("PATH_CONFIG", envir = .GlobalEnv)) {
    get("PATH_CONFIG", envir = .GlobalEnv)$results_dir
  } else {
    "results"
  }
  
  cache_dir <- if (exists("PATH_CONFIG", envir = .GlobalEnv)) {
    get("PATH_CONFIG", envir = .GlobalEnv)$cache_dir
  } else {
    "cache"
  }
  
  if (is.null(tcga_dir)) {
    tcga_dir <- if (exists("PATH_CONFIG", envir = .GlobalEnv)) {
      get("PATH_CONFIG", envir = .GlobalEnv)$tcga_dir
    } else {
      NULL
    }
  }
  
  if (is.null(tcga_dir) || is.na(tcga_dir) || !nzchar(tcga_dir)) {
    stop("tcga_dir must be supplied or PATH_CONFIG$tcga_dir must be defined.", call. = FALSE)
  }
  if (!dir.exists(tcga_dir)) {
    stop("TCGA directory does not exist: ", tcga_dir, call. = FALSE)
  }
  
  surv_dir <- file.path(results_dir, "SURVIVAL")
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(surv_dir, showWarnings = FALSE, recursive = TRUE)
  
  save_surv <- function(obj, filename, type = c("csv", "rds", "txt")) {
    type <- match.arg(type)
    out <- file.path(surv_dir, filename)
    if (type == "csv") readr::write_csv(obj, out)
    if (type == "rds") saveRDS(obj, out)
    if (type == "txt") writeLines(as.character(obj), out)
    cat("  ✓ Saved: ", normalizePath(out, winslash = "/", mustWork = FALSE), "\n", sep = "")
    invisible(out)
  }
  
  clean_patient_barcode <- function(x) {
    x <- as.character(x)
    substr(x, 1, 12)
  }
  
  clean_sample_barcode <- function(x) {
    x <- as.character(x)
    substr(x, 1, 16)
  }
  
  find_first_col <- function(df, candidates) {
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) return(NULL)
    hit[1]
  }
  
  clean_tcga_stage <- function(x) {
    y <- as.character(x)
    y <- stringr::str_to_upper(stringr::str_squish(y))
    
    y[y %in% c(
      "", "NA", "NAN",
      "NOT REPORTED", "NOT AVAILABLE",
      "[NOT REPORTED]", "[NOT AVAILABLE]"
    )] <- NA_character_
    
    major <- stringr::str_match(
      y,
      "^STAGE\\s+(IV|III|II|I)(?:[A-C])?(?:\\b|\\s|$)"
    )[, 2]
    
    factor(major, levels = c("I", "II", "III", "IV"))
  }
  
  if (!identical(score_method, "ssgsea")) {
    stop("Only score_method='ssgsea' is currently implemented.", call. = FALSE)
  }
  
  # ---------------------------------------------------------------------------
  # 2. Load strict conserved UP-core genes
  # ---------------------------------------------------------------------------
  cat("  [16A] Loading strict conserved UP-core genes...\n")
  
  core_file <- file.path(results_dir, "Conserved_Core_UP_genes_strict.txt")
  if (!file.exists(core_file)) {
    stop("Missing required file: ", core_file, "\nRun the conserved core module first.", call. = FALSE)
  }
  
  core_genes <- readLines(core_file, warn = FALSE)
  core_genes <- unique(trimws(core_genes))
  core_genes <- core_genes[core_genes != ""]
  
  if (length(core_genes) == 0) {
    stop("Conserved_Core_UP_genes_strict.txt is empty.", call. = FALSE)
  }
  
  cat("      Strict UP-core genes:", length(core_genes), "\n")
  
  # ---------------------------------------------------------------------------
  # 3. Query/load TCGA metadata and match local STAR-count files
  # ---------------------------------------------------------------------------
  cat("  [16B] Matching local TCGA STAR-count files to GDC metadata...\n")
  
  query_cache <- file.path(cache_dir, "TCGA_BRCA_query_results.rds")
  
  if (isTRUE(use_cached_query) && file.exists(query_cache)) {
    meta <- readRDS(query_cache)
  } else {
    query <- TCGAbiolinks::GDCquery(
      project = "TCGA-BRCA",
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c("Primary Tumor", "Solid Tissue Normal")
    )
    meta <- TCGAbiolinks::getResults(query)
    saveRDS(meta, query_cache)
  }
  
  required_meta_cols <- c("file_id", "file_name", "sample_type", "cases.submitter_id")
  missing_meta_cols <- setdiff(required_meta_cols, names(meta))
  if (length(missing_meta_cols) > 0) {
    stop("TCGA metadata missing required column(s): ", paste(missing_meta_cols, collapse = ", "), call. = FALSE)
  }
  
  local_files <- list.files(
    tcga_dir,
    pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(local_files) == 0) {
    stop("No local STAR-count TSV files found under: ", tcga_dir, call. = FALSE)
  }
  
  local_tbl <- tibble::tibble(
    full_path = local_files,
    file_id = basename(dirname(local_files)),
    local_file_name = basename(local_files)
  )
  
  sample_tbl <- meta %>%
    dplyr::filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal")) %>%
    dplyr::select(file_id, file_name, sample_type, cases.submitter_id) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(local_tbl, by = "file_id") %>%
    dplyr::mutate(
      patient_id = clean_patient_barcode(cases.submitter_id),
      sample_barcode = clean_sample_barcode(cases.submitter_id),
      is_primary_tumor = sample_type == "Primary Tumor"
    )
  
  if (nrow(sample_tbl) == 0) {
    stop("No overlap between local files and queried TCGA metadata.", call. = FALSE)
  }
  
  cat("      Local files matched:", nrow(sample_tbl), "\n")
  cat("      Primary tumour files:", sum(sample_tbl$is_primary_tumor), "\n")
  cat("      Normal files:", sum(sample_tbl$sample_type == "Solid Tissue Normal"), "\n")
  
  save_surv(sample_tbl, "TCGA_survival_local_sample_table.csv", "csv")
  
  tumor_tbl <- sample_tbl %>%
    dplyr::filter(is_primary_tumor) %>%
    dplyr::arrange(patient_id, sample_barcode, file_id) %>%
    dplyr::distinct(patient_id, .keep_all = TRUE)
  
  if (nrow(tumor_tbl) < 50) {
    stop("Too few primary tumour patients after patient-level collapse: ", nrow(tumor_tbl), call. = FALSE)
  }
  
  cat("      Patient-level primary tumours retained:", nrow(tumor_tbl), "\n")
  
  # ---------------------------------------------------------------------------
  # 4. Build count matrix from local STAR-count files
  # ---------------------------------------------------------------------------
  cat("  [16C] Reading local STAR-count files and building tumour count matrix...\n")
  
  expr_mat_cache <- file.path(cache_dir, "TCGA_BRCA_surv_expr_mat.rds")
  
  if (isTRUE(use_cached_query) && file.exists(expr_mat_cache)) {
    cat("      Loading cached expression matrix (skipping [16C] and [16D])...\n")
    expr_mat <- readRDS(expr_mat_cache)
    cat("      Symbol-level expression genes:", nrow(expr_mat), "\n")
    cat("      Patients:", ncol(expr_mat), "\n")
  } else {
    read_star_counts <- function(path) {
      df <- readr::read_tsv(
        path,
        skip = 1,
        show_col_types = FALSE,
        progress = FALSE
      )
      
      required_cols <- c("gene_id", "gene_name", "gene_type", "unstranded")
      missing_cols <- setdiff(required_cols, names(df))
      if (length(missing_cols) > 0) {
        stop(
          "STAR-count file missing column(s): ", paste(missing_cols, collapse = ", "),
          "\nFile: ", path,
          call. = FALSE
        )
      }
      
      df %>%
        dplyr::filter(!stringr::str_starts(gene_id, "^N_")) %>%
        dplyr::transmute(
          gene_id = as.character(gene_id),
          gene_id_clean = sub("\\..*$", "", gene_id),
          gene_symbol = as.character(gene_name),
          gene_type = as.character(gene_type),
          unstranded = suppressWarnings(as.numeric(unstranded))
        )
    }
    
    first_df <- read_star_counts(tumor_tbl$full_path[1])
    gene_ref <- first_df %>%
      dplyr::select(gene_id, gene_id_clean, gene_symbol, gene_type)
    
    count_list <- lapply(tumor_tbl$full_path, read_star_counts)
    ref_ids <- count_list[[1]]$gene_id
    same_order <- vapply(count_list, function(x) identical(x$gene_id, ref_ids), logical(1))
    
    if (!all(same_order)) {
      cat("      Gene order differs across files. Using defensive merge by gene_id...\n")
      
      merged_counts <- count_list[[1]] %>%
        dplyr::select(gene_id, sample_1 = unstranded)
      
      if (length(count_list) > 1) {
        for (i in 2:length(count_list)) {
          tmp <- count_list[[i]] %>%
            dplyr::select(gene_id, !!paste0("sample_", i) := unstranded)
          merged_counts <- merged_counts %>% dplyr::inner_join(tmp, by = "gene_id")
        }
      }
      
      count_mat <- merged_counts %>%
        dplyr::select(-gene_id) %>%
        as.data.frame() %>%
        as.matrix()
      rownames(count_mat) <- merged_counts$gene_id
      gene_ref <- gene_ref %>% dplyr::filter(gene_id %in% rownames(count_mat))
    } else {
      count_mat <- do.call(cbind, lapply(count_list, function(x) x$unstranded))
      rownames(count_mat) <- ref_ids
    }
    
    colnames(count_mat) <- tumor_tbl$patient_id
    
    gene_ref <- gene_ref %>%
      dplyr::distinct(gene_id, .keep_all = TRUE) %>%
      dplyr::right_join(tibble::tibble(gene_id = rownames(count_mat)), by = "gene_id")
    
    if (!identical(gene_ref$gene_id, rownames(count_mat))) {
      idx <- match(rownames(count_mat), gene_ref$gene_id)
      gene_ref <- gene_ref[idx, ]
    }
    
    cat("      Count matrix genes:", nrow(count_mat), "\n")
    cat("      Count matrix patients:", ncol(count_mat), "\n")
    
    # -------------------------------------------------------------------------
    # 5. TMM + voom normalization
    # -------------------------------------------------------------------------
    cat("  [16D] Normalising tumour expression with TMM + voom...\n")
    
    dge <- edgeR::DGEList(counts = count_mat)
    keep <- edgeR::filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge)
    
    gene_ref_keep <- gene_ref[keep, ]
    
    design <- matrix(
      1,
      nrow = ncol(dge),
      ncol = 1,
      dimnames = list(colnames(dge), "intercept")
    )
    
    v <- limma::voom(dge, design = design, plot = FALSE)
    
    expr_tbl <- as.data.frame(v$E) %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::left_join(gene_ref_keep, by = "gene_id") %>%
      dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
      dplyr::group_by(gene_symbol) %>%
      dplyr::summarise(dplyr::across(where(is.numeric), mean), .groups = "drop")
    
    expr_mat <- expr_tbl %>%
      tibble::column_to_rownames("gene_symbol") %>%
      as.matrix()
    
    cat("      Symbol-level expression genes:", nrow(expr_mat), "\n")
    cat("      Patients:", ncol(expr_mat), "\n")
    
    saveRDS(expr_mat, expr_mat_cache)
    cat("      Expression matrix cached for future runs.\n")
  }
  
  # ---------------------------------------------------------------------------
  # 6. Core-gene overlap and signature scoring
  # ---------------------------------------------------------------------------
  cat("  [16E] Computing conserved-core signature scores...\n")
  
  common_genes <- intersect(core_genes, rownames(expr_mat))
  missing_genes <- setdiff(core_genes, rownames(expr_mat))
  
  overlap_tbl <- tibble::tibble(
    gene_symbol = core_genes,
    detected_in_tcga = gene_symbol %in% common_genes
  )
  save_surv(overlap_tbl, "TCGA_signature_gene_overlap.csv", "csv")
  
  cat("      Matched core genes:", length(common_genes), "/", length(core_genes), "\n")
  
  if (length(common_genes) < 10) {
    stop("Too few core genes matched in TCGA expression matrix: ", length(common_genes), call. = FALSE)
  }
  
  gene_set_list <- list(CONSERVED_UP_CORE = common_genes)
  
  ssgsea_scores <- tryCatch({
    params <- GSVA::ssgseaParam(exprData = expr_mat, geneSets = gene_set_list, normalize = TRUE)
    as.numeric(GSVA::gsva(params)[1, ])
  }, error = function(e1) {
    cat("      New GSVA API failed; trying legacy gsva(..., method='ssgsea')...\n")
    tryCatch({
      as.numeric(GSVA::gsva(
        expr_mat,
        gene_set_list,
        method = "ssgsea",
        kcdf = "Gaussian",
        verbose = FALSE
      )[1, ])
    }, error = function(e2) {
      stop(
        "ssGSEA scoring failed. New API error: ", conditionMessage(e1),
        " | Legacy API error: ", conditionMessage(e2),
        call. = FALSE
      )
    })
  })
  
  names(ssgsea_scores) <- colnames(expr_mat)
  
  score_tbl <- tibble::tibble(
    patient_id = colnames(expr_mat),
    conserved_core_ssgsea = as.numeric(ssgsea_scores[colnames(expr_mat)])
  )
  
  # ---------------------------------------------------------------------------
  # 7. TCGA OS endpoint via GDCquery_clinic
  # ---------------------------------------------------------------------------
  cat("  [16F] Loading TCGA OS endpoint...\n")
  
  surv_cache <- file.path(cache_dir, "TCGA_BRCA_survival.rds")
  
  if (isTRUE(use_cached_query) && file.exists(surv_cache)) {
    surv_raw <- readRDS(surv_cache)
  } else {
    surv_raw <- TCGAbiolinks::GDCquery_clinic("TCGA-BRCA", type = "clinical")
    saveRDS(surv_raw, surv_cache)
  }
  
  surv2 <- surv_raw %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      patient_id = clean_patient_barcode(bcr_patient_barcode),
      OS.time = suppressWarnings(as.numeric(
        dplyr::if_else(
          !is.na(days_to_death) & days_to_death > 0,
          days_to_death,
          days_to_last_follow_up
        )
      )),
      OS = dplyr::case_when(
        tolower(vital_status) %in% c("dead", "deceased") ~ 1,
        tolower(vital_status) %in% c("alive", "living")  ~ 0,
        TRUE ~ NA_real_
      ),
      age_at_diagnosis = suppressWarnings(as.numeric(age_at_diagnosis)),
      age_years = age_at_diagnosis / 365.25,
      age_10y = age_years / 10
    ) %>%
    dplyr::select(patient_id, OS.time, OS, age_at_diagnosis, age_years, age_10y) %>%
    dplyr::distinct(patient_id, .keep_all = TRUE)
  
  # ---------------------------------------------------------------------------
  # 8. PAM50 subtype metadata and clinical covariates
  # ---------------------------------------------------------------------------
  cat("  [16G] Loading PAM50/subtype metadata and clinical stage...\n")
  
  stage_grade <- surv_raw %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      patient_id = clean_patient_barcode(bcr_patient_barcode),
      stage = as.character(ajcc_pathologic_stage),
      stage_clean = clean_tcga_stage(ajcc_pathologic_stage),
      clinical_grade = as.character(tumor_grade)
    ) %>%
    dplyr::select(patient_id, stage, stage_clean, clinical_grade) %>%
    dplyr::distinct(patient_id, .keep_all = TRUE)
  
  stage_clean_qc <- stage_grade %>%
    dplyr::count(stage, stage_clean, name = "n") %>%
    dplyr::arrange(stage_clean, stage)
  save_surv(stage_clean_qc, "TCGA_stage_clean_qc.csv", "csv")
  
  cat("      Raw clinical stage distribution:\n")
  print(table(stage_grade$stage, useNA = "ifany"))
  cat("      Cleaned clinical stage distribution:\n")
  print(table(stage_grade$stage_clean, useNA = "ifany"))
  
  subtype_cache <- file.path(cache_dir, "TCGA_BRCA_subtypes.rds")
  subtype_tbl <- NULL
  
  subtype_raw <- tryCatch({
    if (isTRUE(use_cached_query) && file.exists(subtype_cache)) {
      readRDS(subtype_cache)
    } else {
      tmp <- TCGAbiolinks::TCGAquery_subtype("BRCA")
      saveRDS(tmp, subtype_cache)
      tmp
    }
  }, error = function(e) {
    warning("Could not retrieve TCGA-BRCA subtype metadata: ", conditionMessage(e), call. = FALSE)
    NULL
  })
  
  if (!is.null(subtype_raw)) {
    subtype_raw <- tibble::as_tibble(subtype_raw)
    subtype_patient_col <- find_first_col(subtype_raw, c("patient", "Patient", "bcr_patient_barcode", "submitter_id"))
    pam50_col <- find_first_col(subtype_raw, c("BRCA_Subtype_PAM50", "PAM50", "Subtype_mRNA", "subtype", "Subtype"))
    er_col <- find_first_col(subtype_raw, c("ER.Status", "ER_status", "ER", "er_status_by_ihc"))
    grade_col <- find_first_col(subtype_raw, c("grade", "Grade", "histological_grade"))
    
    if (!is.null(subtype_patient_col)) {
      subtype_tbl <- subtype_raw %>%
        dplyr::mutate(patient_id = clean_patient_barcode(.data[[subtype_patient_col]])) %>%
        dplyr::transmute(
          patient_id = patient_id,
          PAM50 = if (!is.null(pam50_col)) as.character(.data[[pam50_col]]) else NA_character_,
          ER_status = if (!is.null(er_col)) as.character(.data[[er_col]]) else NA_character_,
          subtype_grade = if (!is.null(grade_col)) as.character(.data[[grade_col]]) else NA_character_
        ) %>%
        dplyr::distinct(patient_id, .keep_all = TRUE)
    }
  }
  
  if (is.null(subtype_tbl)) {
    subtype_tbl <- tibble::tibble(
      patient_id = score_tbl$patient_id,
      PAM50 = NA_character_,
      ER_status = NA_character_,
      subtype_grade = NA_character_
    ) %>% dplyr::distinct(patient_id, .keep_all = TRUE)
  }
  
  # ---------------------------------------------------------------------------
  # 9. Final analysis dataset
  # ---------------------------------------------------------------------------
  cat("  [16H] Building patient-level survival/signature dataset...\n")
  
  final_df <- score_tbl %>%
    dplyr::left_join(surv2, by = "patient_id") %>%
    dplyr::left_join(subtype_tbl, by = "patient_id") %>%
    dplyr::left_join(stage_grade, by = "patient_id") %>%
    dplyr::filter(!is.na(OS.time), !is.na(OS), OS.time > 0) %>%
    dplyr::mutate(
      conserved_core_ssgsea_z = as.numeric(scale(conserved_core_ssgsea)),
      group = if (group_cutoff == "median") {
        ifelse(conserved_core_ssgsea >= stats::median(conserved_core_ssgsea, na.rm = TRUE), "High", "Low")
      } else {
        stop("Only group_cutoff='median' is currently implemented.", call. = FALSE)
      },
      group = factor(group, levels = c("Low", "High")),
      PAM50 = dplyr::case_when(
        is.na(PAM50) ~ NA_character_,
        stringr::str_detect(tolower(PAM50), "luminal a|luma") ~ "LumA",
        stringr::str_detect(tolower(PAM50), "luminal b|lumb") ~ "LumB",
        stringr::str_detect(tolower(PAM50), "her2") ~ "Her2",
        stringr::str_detect(tolower(PAM50), "basal") ~ "Basal",
        stringr::str_detect(tolower(PAM50), "normal") ~ "Normal",
        TRUE ~ as.character(PAM50)
      )
    )
  
  if (nrow(final_df) < 50) {
    stop("Too few patients with usable OS data after merging: ", nrow(final_df), call. = FALSE)
  }
  
  cat("      Patients with OS data:", nrow(final_df), "\n")
  cat("      OS events:", sum(final_df$OS == 1, na.rm = TRUE), "\n")
  cat("      High/Low groups:", paste(names(table(final_df$group)), as.integer(table(final_df$group)), collapse = "; "), "\n")
  cat("      Stage_clean distribution in final_df:\n")
  print(table(final_df$stage_clean, useNA = "ifany"))
  
  save_surv(final_df, "TCGA_signature_dataset.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # 10. Kaplan-Meier OS summary only
  # ---------------------------------------------------------------------------
  cat("  [16I] Running Kaplan-Meier OS summary...\n")
  
  km_fit <- survival::survfit(survival::Surv(OS.time, OS) ~ group, data = final_df)
  save_surv(km_fit, "TCGA_KM_fit.rds", "rds")
  
  lr <- survival::survdiff(survival::Surv(OS.time, OS) ~ group, data = final_df)
  lr_p <- 1 - stats::pchisq(lr$chisq, df = length(lr$n) - 1)
  
  km_summary <- tibble::tibble(
    comparison = "High vs Low conserved-core ssGSEA median split",
    n_low = as.integer(table(final_df$group)["Low"]),
    n_high = as.integer(table(final_df$group)["High"]),
    os_events = sum(final_df$OS == 1, na.rm = TRUE),
    logrank_chisq = as.numeric(lr$chisq),
    logrank_p = lr_p
  )
  save_surv(km_summary, "TCGA_KM_summary.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # 11. PAM50 subtype score association
  # ---------------------------------------------------------------------------
  cat("  [16J] Testing PAM50 subtype-score association...\n")
  
  subtype_dat <- final_df %>%
    dplyr::filter(!is.na(PAM50), PAM50 != "")
  
  if (nrow(subtype_dat) >= 20 && dplyr::n_distinct(subtype_dat$PAM50) >= 2) {
    kw <- stats::kruskal.test(conserved_core_ssgsea ~ PAM50, data = subtype_dat)
    
    subtype_summary <- subtype_dat %>%
      dplyr::group_by(PAM50) %>%
      dplyr::summarise(
        n = dplyr::n(),
        mean_ssgsea = mean(conserved_core_ssgsea, na.rm = TRUE),
        median_ssgsea = median(conserved_core_ssgsea, na.rm = TRUE),
        sd_ssgsea = stats::sd(conserved_core_ssgsea, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        test = "Kruskal-Wallis",
        statistic = as.numeric(kw$statistic),
        p_value = as.numeric(kw$p.value)
      )
  } else {
    subtype_summary <- tibble::tibble(
      PAM50 = character(), n = integer(), mean_ssgsea = numeric(),
      median_ssgsea = numeric(), sd_ssgsea = numeric(),
      test = character(), statistic = numeric(), p_value = numeric()
    )
  }
  save_surv(subtype_summary, "TCGA_subtype_score_association.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # 12. Method metadata and summary
  # ---------------------------------------------------------------------------
  method_parameters <- tibble::tibble(
    parameter = c(
      "tcga_dir",
      "local_STAR_count_files_matched",
      "patient_level_primary_tumours_retained",
      "patients_retained_for_OS_analysis",
      "OS_events",
      "OS_filtering_rule",
      "normalisation",
      "primary_signature_score",
      "group_cutoff",
      "primary_endpoint",
      "PAM50_source",
      "stage_cleaning",
      "core_genes_input",
      "core_genes_matched_in_TCGA",
      "seed"
    ),
    value = c(
      normalizePath(tcga_dir, winslash = "/", mustWork = FALSE),
      as.character(nrow(sample_tbl)),
      as.character(nrow(tumor_tbl)),
      as.character(nrow(final_df)),
      as.character(sum(final_df$OS == 1, na.rm = TRUE)),
      "Patients retained if matched primary tumour expression data, non-missing OS status, non-missing OS time, and OS.time > 0",
      "TMM normalisation followed by limma-voom on patient-level primary tumour counts",
      "ssGSEA score of strict conserved UP-core genes using GSVA",
      group_cutoff,
      "Overall survival (OS)",
      "TCGAbiolinks::TCGAquery_subtype('BRCA')",
      "Major AJCC stage extracted using clean_tcga_stage(); Stage IV is explicitly protected from Stage I misclassification",
      basename(core_file),
      paste0(length(common_genes), "/", length(core_genes)),
      as.character(seed)
    )
  )
  save_surv(method_parameters, "TCGA_survival_method_parameters.csv", "csv")
  
  summary_lines <- c(
    "TCGA-BRCA PATIENT-LEVEL SIGNATURE VALIDATION SUMMARY — NO COX",
    "",
    paste("Local STAR-count files matched:", nrow(sample_tbl)),
    paste("Patient-level primary tumours retained:", nrow(tumor_tbl)),
    paste("Patients with OS data:", nrow(final_df)),
    paste("Patients excluded from OS analysis after primary-tumour collapse:", nrow(tumor_tbl) - nrow(final_df)),
    paste("OS events:", sum(final_df$OS == 1, na.rm = TRUE)),
    paste("Strict conserved UP-core genes:", length(core_genes)),
    paste("Matched core genes in TCGA:", length(common_genes)),
    paste("Missing core genes in TCGA:", length(missing_genes)),
    paste("Median split groups:", paste(names(table(final_df$group)), as.integer(table(final_df$group)), collapse = "; ")),
    paste("Stage_clean distribution:", paste(names(table(final_df$stage_clean, useNA = "ifany")), as.integer(table(final_df$stage_clean, useNA = "ifany")), collapse = "; ")),
    paste("Log-rank p:", signif(lr_p, 4))
  )
  save_surv(summary_lines, "TCGA_survival_summary.txt", "txt")
  
  assign("tcga_survival_df", final_df, envir = .GlobalEnv)
  assign("tcga_km_os", km_fit, envir = .GlobalEnv)
  assign("tcga_survival_signature_overlap", overlap_tbl, envir = .GlobalEnv)
  assign("tcga_stage_clean_qc", stage_clean_qc, envir = .GlobalEnv)
  
  cat("\n")
  cat("  ✓ MODULE 16 NO COX COMPLETE\n\n")
  cat("  Main outputs:\n")
  cat("   - ", normalizePath(file.path(surv_dir, "TCGA_signature_dataset.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("   - ", normalizePath(file.path(surv_dir, "TCGA_KM_summary.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("   - ", normalizePath(file.path(surv_dir, "TCGA_subtype_score_association.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("   - ", normalizePath(file.path(surv_dir, "TCGA_survival_summary.txt"), winslash = "/", mustWork = FALSE), "\n\n", sep = "")
  
  invisible(list(
    dataset = final_df,
    km_fit = km_fit,
    km_summary = km_summary,
    subtype_association = subtype_summary,
    gene_overlap = overlap_tbl,
    stage_clean_qc = stage_clean_qc,
    method_parameters = method_parameters
  ))
}

