################################################################################
# STEP 1: DOG ANALYSIS, BLOCK-AWARE LIMMA & HALLMARK GSEA
# Primary Q1-oriented implementation:
#   - default: keep all classifiable dog samples
#   - model repeated tumor/normal IDs with limma::duplicateCorrelation when present
#   - optional: complete-pairs-only mode via PARAM_CONFIG$dog_analysis_mode
#
# Main outputs:
#   - DOG_mapping_1to1_highconf.rds
#   - DOG_expr_humanENSG_1to1.rds
#   - DOG_sample_structure_QC.csv
#   - DOG_DE_limma_blockaware_all_full.csv / .rds
#   - DOG_DE_limma_paired_full.csv / DOG_DE_limma_paired.rds (backward-compatible names)
#   - DOG_rank_tstat.rds
#   - DOG_Hallmark_fgsea.rds
#   - DOG_Hallmark_fgsea_simple_plus.csv
################################################################################

module_dog_gsea <- function() {
  suppressPackageStartupMessages({
    library(limma)
    library(biomaRt)
    library(fgsea)
    library(msigdbr)
    library(BiocParallel)
    library(dplyr)
    library(tibble)
  })

  cat("  [1A] Loading dog expression data...\n")

  if (!file.exists(PATH_CONFIG$dog_expr)) {
    stop("Dog expression file not found: ", PATH_CONFIG$dog_expr)
  }

  expr <- read.csv(PATH_CONFIG$dog_expr, row.names = 1, check.names = FALSE)

  if (nrow(expr) == 0 || ncol(expr) == 0) {
    stop("Dog expression matrix is empty.")
  }

  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"

  cat("      Genes:", nrow(expr), "| Samples:", ncol(expr), "\n")

  # ---------------------------------------------------------------------------
  # 1B. DOG-TO-HUMAN MAPPING (biomaRt, cached)
  # ---------------------------------------------------------------------------
  cat("  [1B] Querying biomaRt for dog-to-human orthologs...\n")

  genes <- rownames(expr)

  mapping <- load_cached("DOG_biomart_mapping_raw", function() {
    dog <- useEnsembl(
      biomart = "genes",
      dataset = "clfamiliaris_gene_ensembl",
      version = PARAM_CONFIG$ensembl_version
    )

    getBM(
      attributes = c(
        "ensembl_gene_id",
        "hsapiens_homolog_ensembl_gene",
        "hsapiens_homolog_associated_gene_name",
        "hsapiens_homolog_orthology_type",
        "hsapiens_homolog_orthology_confidence"
      ),
      filters = "ensembl_gene_id",
      values = genes,
      mart = dog
    )
  })

  if (is.null(mapping) || nrow(mapping) == 0) {
    stop("biomaRt mapping returned no rows.")
  }

  colnames(mapping) <- c(
    "dog_ensembl", "human_ensembl", "human_symbol",
    "orthology_type", "orthology_conf"
  )

  mapping_1to1 <- mapping %>%
    dplyr::filter(
      !is.na(human_symbol), human_symbol != "",
      !is.na(human_ensembl), human_ensembl != "",
      orthology_type == PARAM_CONFIG$ortholog_type_filter,
      orthology_conf == PARAM_CONFIG$ortholog_conf_filter
    ) %>%
    dplyr::distinct(dog_ensembl, .keep_all = TRUE)

  cat("      Input dog genes:", length(genes), "\n")
  cat("      Total mapped rows:", nrow(mapping), "\n")
  cat("      1:1 high-confidence mapped genes:", nrow(mapping_1to1), "\n")

  if (nrow(mapping_1to1) == 0) {
    stop("No 1:1 high-confidence dog-to-human mappings remained after filtering.")
  }

  save_result(mapping_1to1, "DOG_mapping_1to1_highconf.rds", "rds")

  # ---------------------------------------------------------------------------
  # 1C. MAP EXPRESSION TO HUMAN ENSG SPACE
  # ---------------------------------------------------------------------------
  cat("  [1C] Converting expression matrix to human ENSG space...\n")

  keep <- intersect(rownames(expr), mapping_1to1$dog_ensembl)
  expr_1to1 <- expr[keep, , drop = FALSE]

  dog_to_human_ensg <- setNames(mapping_1to1$human_ensembl, mapping_1to1$dog_ensembl)
  new_names <- dog_to_human_ensg[rownames(expr_1to1)]

  expr_1to1 <- expr_1to1[!is.na(new_names), , drop = FALSE]
  rownames(expr_1to1) <- new_names[!is.na(new_names)]

  dup_ensg <- sum(duplicated(rownames(expr_1to1)))
  cat("      Final matrix:", nrow(expr_1to1), "genes ×", ncol(expr_1to1), "samples\n")
  cat("      Duplicate human ENSG rownames:", dup_ensg, "\n")

  if (nrow(expr_1to1) == 0) {
    stop("Expression matrix became empty after mapping to human ENSG space.")
  }

  if (dup_ensg > 0) {
    warning("Duplicate human ENSG rownames detected after 1:1 filtering. Collapsing duplicates by mean.")
    expr_df <- as.data.frame(expr_1to1)
    expr_df$human_ensembl <- rownames(expr_1to1)
    expr_df <- stats::aggregate(. ~ human_ensembl, data = expr_df, FUN = mean)
    rownames(expr_df) <- expr_df$human_ensembl
    expr_df$human_ensembl <- NULL
    expr_1to1 <- as.matrix(expr_df)
    storage.mode(expr_1to1) <- "numeric"
  }

  save_result(expr_1to1, "DOG_expr_humanENSG_1to1.rds", "rds")

  # ---------------------------------------------------------------------------
  # 1D. SAMPLE PARSING + PARTIALLY PAIRED STRUCTURE QC
  # ---------------------------------------------------------------------------
  cat("  [1D] Parsing sample groups and checking paired/unpaired structure...\n")

  cn <- colnames(expr_1to1)
  group <- ifelse(
    grepl("-tumor$", cn), "tumor",
    ifelse(grepl("-normal$", cn), "normal", NA)
  )

  if (any(is.na(group))) {
    bad_samples <- cn[is.na(group)]
    stop(
      "Some samples could not be classified as tumor/normal from column names: ",
      paste(bad_samples, collapse = ", ")
    )
  }

  grp <- factor(group, levels = c("normal", "tumor"))
  block <- factor(sub("-(tumor|normal)$", "", cn))

  tab_id <- as.data.frame.matrix(table(block, grp))
  if (!"normal" %in% names(tab_id)) tab_id$normal <- 0L
  if (!"tumor" %in% names(tab_id)) tab_id$tumor <- 0L
  tab_id$sample_id <- rownames(tab_id)
  tab_id <- tab_id %>%
    dplyr::select(sample_id, normal, tumor) %>%
    dplyr::mutate(
      total_samples_for_id = normal + tumor,
      complete_tumor_normal_pair = normal > 0 & tumor > 0
    )

  n_complete_pairs <- sum(tab_id$complete_tumor_normal_pair)
  n_singleton_ids <- sum(tab_id$total_samples_for_id == 1)

  cat("      Samples:", length(grp), "\n")
  cat("      Tumor samples:", sum(grp == "tumor"), "\n")
  cat("      Normal samples:", sum(grp == "normal"), "\n")
  cat("      Unique IDs:", length(unique(block)), "\n")
  cat("      IDs with both normal & tumor:", n_complete_pairs, "\n")
  cat("      Singleton IDs:", n_singleton_ids, "\n")

  save_result(tab_id, "DOG_sample_structure_QC.csv", "csv")

  if (length(unique(grp)) < 2) {
    stop("Both normal and tumor samples are required for differential expression.")
  }

  dog_analysis_mode <- if (!is.null(PARAM_CONFIG$dog_analysis_mode)) {
    PARAM_CONFIG$dog_analysis_mode
  } else {
    "all_blockaware"
  }

  if (!dog_analysis_mode %in% c("all_blockaware", "complete_pairs_only")) {
    stop("Unsupported PARAM_CONFIG$dog_analysis_mode: ", dog_analysis_mode)
  }

  if (dog_analysis_mode == "complete_pairs_only") {
    complete_ids <- tab_id$sample_id[tab_id$complete_tumor_normal_pair]
    if (length(complete_ids) == 0) {
      stop("No complete tumor-normal pairs detected; complete_pairs_only cannot proceed.")
    }
    keep_samples <- block %in% complete_ids
    expr_1to1 <- expr_1to1[, keep_samples, drop = FALSE]
    grp <- droplevels(grp[keep_samples])
    block <- droplevels(block[keep_samples])
    cat("      Analysis mode: complete_pairs_only\n")
    cat("      Retained samples:", ncol(expr_1to1), "\n")
  } else {
    cat("      Analysis mode: all_blockaware\n")
    cat("      Retained samples:", ncol(expr_1to1), "\n")
  }

  # Recompute paired count after any optional filtering
  tab_used <- table(block, grp)
  n_repeated_ids_used <- sum(rowSums(tab_used) > 1)
  n_complete_pairs_used <- sum(rowSums(tab_used > 0) == 2)

  # ---------------------------------------------------------------------------
  # 1E. DIFFERENTIAL EXPRESSION (ALL-SAMPLE BLOCK-AWARE LIMMA BY DEFAULT)
  # ---------------------------------------------------------------------------
  cat("  [1E] Running limma differential expression...\n")
 
  # log2(FPKM + 1) transformation applied prior to limma.
  # Input data are FPKM values (continuous, already length-normalised).
  # log2(x+1) with trend=FALSE (default eBayes) is the standard approach for
  # pre-normalised microarray-style data; voom is not appropriate here because
  # voom expects raw counts and models mean-variance from count data.
  # Methods text: "FPKM values were log2(x+1) transformed prior to limma
  # differential expression analysis."
  Y <- log2(expr_1to1 + 1)
  design <- model.matrix(~ grp)

  use_block <- n_repeated_ids_used > 0 && n_complete_pairs_used > 0

  if (use_block) {
    corfit <- duplicateCorrelation(Y, design, block = block)
    cat("      duplicateCorrelation consensus:", round(corfit$consensus, 4), "\n")
    fit <- lmFit(Y, design, block = block, correlation = corfit$consensus)
  } else {
    corfit <- NULL
    cat("      No repeated paired IDs detected/used; running ordinary limma model.\n")
    fit <- lmFit(Y, design)
  }

  fit <- eBayes(fit)
  deg <- topTable(fit, coef = "grptumor", number = Inf, sort.by = "P")

  sig_genes_fdr <- sum(deg$adj.P.Val < 0.05, na.rm = TRUE)
  sig_genes_fdr_fc <- sum(abs(deg$logFC) >= 1 & deg$adj.P.Val < 0.05, na.rm = TRUE)

  cat("      Significant genes (FDR < 0.05):", sig_genes_fdr, "\n")
  cat("      Significant genes (|logFC| >= 1 & FDR < 0.05):", sig_genes_fdr_fc, "\n")

  save_result(deg, paste0("DOG_DE_limma_", dog_analysis_mode, "_full.csv"), "csv_rownames")
  save_result(deg, paste0("DOG_DE_limma_", dog_analysis_mode, ".rds"), "rds")

  # Backward-compatible names for downstream scripts / older manuscript drafts.
  # These are no longer semantically "paired-only"; methods should call this block-aware limma.
  save_result(deg, "DOG_DE_limma_paired_full.csv", "csv_rownames")
  save_result(deg, "DOG_DE_limma_paired.rds", "rds")

  deg_sig <- subset(deg, abs(logFC) >= 1 & adj.P.Val < 0.05)
  save_result(deg_sig, paste0("DOG_DE_limma_", dog_analysis_mode, "_sig_FDR0.05_logFC1.csv"), "csv_rownames")
  save_result(deg_sig, "DOG_DE_limma_sig_FDR0.05_logFC1.csv", "csv_rownames")

  # ---------------------------------------------------------------------------
  # 1F. HALLMARK GSEA
  # ---------------------------------------------------------------------------
  cat("  [1F] Running hallmark GSEA...\n")

  dog_rank <- deg$t
  names(dog_rank) <- rownames(deg)
  dog_rank <- dog_rank[!is.na(dog_rank)]
  dog_rank <- sort(dog_rank, decreasing = TRUE)

  if (length(dog_rank) == 0) {
    stop("dog_rank is empty after removing NA values.")
  }

  save_result(dog_rank, "DOG_rank_tstat.rds", "rds")

  gsea_min_size <- if (!is.null(PARAM_CONFIG$dog_gsea_min_size)) PARAM_CONFIG$dog_gsea_min_size else 15
  gsea_max_size <- if (!is.null(PARAM_CONFIG$dog_gsea_max_size)) PARAM_CONFIG$dog_gsea_max_size else 500
  gsea_nperm_simple <- if (!is.null(PARAM_CONFIG$dog_gsea_nperm_simple)) PARAM_CONFIG$dog_gsea_nperm_simple else 10000
  gsea_seed <- if (exists("EXECUTION_CONFIG", envir = .GlobalEnv) && !is.null(EXECUTION_CONFIG$random_seed)) {
    EXECUTION_CONFIG$random_seed
  } else {
    1
  }

  hall <- msigdbr(species = "Homo sapiens", collection = "H")
  hall_sets <- split(hall$ensembl_gene, hall$gs_name)
  hall_sets <- lapply(hall_sets, function(s) intersect(s, names(dog_rank)))
  hall_sets <- hall_sets[lengths(hall_sets) >= gsea_min_size]

  if (length(hall_sets) == 0) {
    stop("No Hallmark gene sets remained after intersecting with dog ranked universe.")
  }

  # Reproducibility:
  # fgsea() defaults to the multilevel implementation. We explicitly set
  # nPermSimple and use serial BiocParallel execution to make the preliminary
  # permutation stage and borderline p-values more reproducible across systems.
  set.seed(gsea_seed)

  fg <- fgsea::fgsea(
    pathways = hall_sets,
    stats = dog_rank,
    minSize = gsea_min_size,
    maxSize = gsea_max_size,
    scoreType = "std",
    nPermSimple = gsea_nperm_simple,
    BPPARAM = BiocParallel::SerialParam()
  )

  fg <- fg[order(fg$padj, -abs(fg$NES)), ]

  sig_hallmarks <- sum(fg$padj < PARAM_CONFIG$dog_hallmark_padj, na.rm = TRUE)
  cat("      Significant hallmarks (FDR <", PARAM_CONFIG$dog_hallmark_padj, "): ", sig_hallmarks, "\n", sep = "")

  save_result(fg, "DOG_Hallmark_fgsea.rds", "rds")

  fg_simple_plus <- fg[, c("pathway", "NES", "padj", "size", "leadingEdge", "pval")]
  fg_simple_plus$leadingEdge <- vapply(
    fg_simple_plus$leadingEdge,
    function(x) paste(x, collapse = ";"),
    character(1)
  )
  save_result(fg_simple_plus, "DOG_Hallmark_fgsea_simple_plus.csv", "csv")

  # ---------------------------------------------------------------------------
  # 1G. EXTRACT SIGNIFICANT HALLMARK PATHWAYS / LEADING-EDGE GENE LISTS
  # ---------------------------------------------------------------------------
  cat("  [1G] Exporting significant hallmark pathway and leading-edge summaries...\n")

  dog_sig <- subset(fg, padj < PARAM_CONFIG$dog_hallmark_padj)
  dog_up_paths <- dog_sig$pathway[dog_sig$NES > 0]
  dog_down_paths <- dog_sig$pathway[dog_sig$NES < 0]

  dog_up_genes <- unique(unlist(dog_sig$leadingEdge[dog_sig$NES > 0]))
  dog_down_genes <- unique(unlist(dog_sig$leadingEdge[dog_sig$NES < 0]))

  save_result(dog_up_paths, "DOG_up_pathways_padj0.05.txt", "txt")
  save_result(dog_down_paths, "DOG_down_pathways_padj0.05.txt", "txt")
  save_result(dog_up_genes, "DOG_up_leadingEdge_ENSG.txt", "txt")
  save_result(dog_down_genes, "DOG_down_leadingEdge_ENSG.txt", "txt")

  # ---------------------------------------------------------------------------
  # 1H. OPTIONAL PROGENy ANALYSIS
  # ---------------------------------------------------------------------------
  run_progeny <- isTRUE(PARAM_CONFIG$run_progeny)
  prog <- NULL
  prog_res <- NULL

  if (run_progeny) {
    cat("  [1H] Running optional PROGENy pathway activity analysis...\n")

    if (!requireNamespace("progeny", quietly = TRUE)) {
      stop("Package 'progeny' is required because PARAM_CONFIG$run_progeny = TRUE.")
    }

    ensg_to_sym <- setNames(mapping_1to1$human_symbol, mapping_1to1$human_ensembl)
    sym <- ensg_to_sym[rownames(Y)]
    keep_sym <- !is.na(sym) & sym != ""

    Ysym <- Y[keep_sym, , drop = FALSE]
    rownames(Ysym) <- sym[keep_sym]

    if (any(duplicated(rownames(Ysym)))) {
      df_tmp <- as.data.frame(Ysym)
      df_tmp$Symbol <- rownames(Ysym)
      Ysym2 <- stats::aggregate(. ~ Symbol, data = df_tmp, FUN = mean)
      rownames(Ysym2) <- Ysym2$Symbol
      Ysym2$Symbol <- NULL
      Ysym <- as.matrix(Ysym2)
    }

    Ysym <- as.matrix(Ysym)
    storage.mode(Ysym) <- "numeric"

    prog <- progeny::progeny(Ysym, organism = "Human", scale = TRUE, top = 1000)
    design_p <- model.matrix(~ grp)

    if (use_block) {
      corfit_p <- duplicateCorrelation(t(prog), design_p, block = block)
      fit_p <- lmFit(t(prog), design_p, block = block, correlation = corfit_p$consensus)
    } else {
      fit_p <- lmFit(t(prog), design_p)
    }

    fit_p <- eBayes(fit_p)
    prog_res <- topTable(fit_p, coef = "grptumor", number = Inf, sort.by = "P")

    save_result(prog, "DOG_PROGENy_scores.rds", "rds")
    save_result(prog_res, "DOG_PROGENy_tumor_vs_normal.rds", "rds")
    save_result(prog_res, "DOG_PROGENy_tumor_vs_normal.csv", "csv_rownames")
  } else {
    cat("  [1H] Skipping PROGENy (PARAM_CONFIG$run_progeny = FALSE).\n")
  }

  # ---------------------------------------------------------------------------
  # 1I. EXPORT TO GLOBAL ENV + RETURN
  # ---------------------------------------------------------------------------
  assign("mapping_1to1", mapping_1to1, envir = .GlobalEnv)
  assign("expr_1to1", expr_1to1, envir = .GlobalEnv)
  assign("dog_sample_qc", tab_id, envir = .GlobalEnv)
  assign("deg", deg, envir = .GlobalEnv)
  assign("dog_rank", dog_rank, envir = .GlobalEnv)
  assign("fg", fg, envir = .GlobalEnv)
  assign("prog", prog, envir = .GlobalEnv)
  assign("prog_res", prog_res, envir = .GlobalEnv)

  cat("  ✓ STEP 1 COMPLETE\n\n")

  invisible(list(
    mapping_1to1 = mapping_1to1,
    expr_1to1 = expr_1to1,
    dog_sample_qc = tab_id,
    deg = deg,
    dog_rank = dog_rank,
    fg = fg,
    prog = prog,
    prog_res = prog_res
  ))
}
