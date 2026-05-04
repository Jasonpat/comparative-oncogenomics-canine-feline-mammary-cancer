################################################################################
# STEP 4: EXTRACT CONSERVED CORE GENES
# Based on:
#   - conserved hallmarks from Step 3
#   - dog leading-edge genes from fgsea
#   - fixed universe from Step 2
#
# Logic:
#   - For each conserved hallmark:
#       dog leading-edge symbols ∩ universe_hs ∩ cat_up_u
#   - Core UP genes = union across conserved UP hallmarks
#   - Core DOWN genes = union across conserved DOWN hallmarks
#
# Outputs:
#   - DOG_leadingEdge_symbols_conserved_UP.rds
#   - DOG_leadingEdge_symbols_conserved_DOWN.rds
#   - Conserved_Core_UP_genes_strict.txt
#   - Conserved_Core_DOWN_genes_strict.txt
#   - Conserved_UP_core_gene_frequency.csv
#   - Conserved_DOWN_core_gene_frequency.csv
################################################################################

module_conserved_targets <- function() {
  
  cat("  [4A] Checking required upstream objects...\n")
  
  required_objs <- c("fg", "cons_up", "cons_down", "mapping_1to1", "cat_up_u", "cat_down_u", "universe_hs")
  missing_objs <- required_objs[!vapply(required_objs, exists, logical(1), envir = .GlobalEnv)]
  
  if (length(missing_objs) > 0) {
    stop(
      "Missing required object(s) for conserved target extraction: ",
      paste(missing_objs, collapse = ", "),
      ". Run previous modules first."
    )
  }
  
  fg <- get("fg", envir = .GlobalEnv)
  cons_up <- get("cons_up", envir = .GlobalEnv)
  cons_down <- get("cons_down", envir = .GlobalEnv)
  mapping_1to1 <- get("mapping_1to1", envir = .GlobalEnv)
  cat_up_u <- get("cat_up_u", envir = .GlobalEnv)
  cat_down_u <- get("cat_down_u", envir = .GlobalEnv)
  universe_hs <- get("universe_hs", envir = .GlobalEnv)
  
  # ---------------------------------------------------------------------------
  # 4B. BUILD ENSG -> SYMBOL MAPPING
  # ---------------------------------------------------------------------------
  cat("  [4B] Building ENSG-to-symbol mapping...\n")
  
  ensg_to_sym <- setNames(mapping_1to1$human_symbol, mapping_1to1$human_ensembl)
  
  # ---------------------------------------------------------------------------
  # 4C. EXTRACT DOG LEADING-EDGE SYMBOLS PER CONSERVED HALLMARK
  # ---------------------------------------------------------------------------
  cat("  [4C] Extracting dog leading-edge symbols per conserved hallmark...\n")
  
  get_dog_le_symbols <- function(pathway_name) {
    idx <- which(fg$pathway == pathway_name)
    
    if (length(idx) == 0) return(character(0))
    
    le_ensg <- fg$leadingEdge[[idx[1]]]
    
    if (is.null(le_ensg) || length(le_ensg) == 0) return(character(0))
    
    le_sym <- unname(ensg_to_sym[le_ensg])
    le_sym <- le_sym[!is.na(le_sym) & le_sym != ""]
    unique(le_sym)
  }
  
  dog_le_up <- setNames(lapply(cons_up, get_dog_le_symbols), cons_up)
  dog_le_down <- setNames(lapply(cons_down, get_dog_le_symbols), cons_down)
  
  # Restrict to fixed universe
  dog_le_up_u <- lapply(dog_le_up, function(x) intersect(x, universe_hs))
  dog_le_down_u <- lapply(dog_le_down, function(x) intersect(x, universe_hs))
  
  cat("      Conserved UP hallmarks:", length(dog_le_up_u), "\n")
  cat("      Conserved DOWN hallmarks:", length(dog_le_down_u), "\n")
  
  if (length(dog_le_up_u) > 0) {
    cat("      UP leading-edge sizes:", paste(sapply(dog_le_up_u, length), collapse = ", "), "\n")
  }
  if (length(dog_le_down_u) > 0) {
    cat("      DOWN leading-edge sizes:", paste(sapply(dog_le_down_u, length), collapse = ", "), "\n")
  }
  
  save_result(dog_le_up, "DOG_leadingEdge_symbols_conserved_UP.rds", "rds")
  save_result(dog_le_down, "DOG_leadingEdge_symbols_conserved_DOWN.rds", "rds")
  
  # ---------------------------------------------------------------------------
  # 4D. CORE GENE EXTRACTION (STRICT, HALLMARK-WISE)
  # ---------------------------------------------------------------------------
  cat("  [4D] Extracting strict conserved core genes hallmark-by-hallmark...\n")
  
  conserved_up_u <- lapply(cons_up, function(h) {
    intersect(dog_le_up_u[[h]], cat_up_u)
  })
  names(conserved_up_u) <- cons_up
  
  conserved_down_u <- lapply(cons_down, function(h) {
    intersect(dog_le_down_u[[h]], cat_down_u)
  })
  names(conserved_down_u) <- cons_down
  
  core_up_genes_u <- sort(unique(unlist(conserved_up_u)))
  core_down_genes_u <- sort(unique(unlist(conserved_down_u)))
  
  cat("      Core UP genes:", length(core_up_genes_u), "\n")
  cat("      Core DOWN genes:", length(core_down_genes_u), "\n")
  
  if (length(core_up_genes_u) > 0) {
    cat("      Example UP core genes:", paste(head(core_up_genes_u, 10), collapse = ", "), "\n")
  }
  if (length(core_down_genes_u) > 0) {
    cat("      Example DOWN core genes:", paste(head(core_down_genes_u, 10), collapse = ", "), "\n")
  }
  
  save_result(core_up_genes_u, "Conserved_Core_UP_genes_strict.txt", "txt")
  save_result(core_down_genes_u, "Conserved_Core_DOWN_genes_strict.txt", "txt")
  
  # Backward-compatible export name
  save_result(core_up_genes_u, "core_conserved_targets.txt", "txt")
  
  # ---------------------------------------------------------------------------
  # 4E. CORE GENE FREQUENCY ACROSS CONSERVED HALLMARKS
  # ---------------------------------------------------------------------------
  cat("  [4E] Calculating core gene frequency across conserved hallmarks...\n")
  
  if (length(conserved_up_u) > 0 && sum(lengths(conserved_up_u)) > 0) {
    gene_freq_up <- sort(table(unlist(conserved_up_u)), decreasing = TRUE)
    gene_freq_up_df <- data.frame(
      gene = names(gene_freq_up),
      freq = as.integer(gene_freq_up),
      stringsAsFactors = FALSE
    )
  } else {
    gene_freq_up_df <- data.frame(
      gene = character(),
      freq = integer(),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(conserved_down_u) > 0 && sum(lengths(conserved_down_u)) > 0) {
    gene_freq_down <- sort(table(unlist(conserved_down_u)), decreasing = TRUE)
    gene_freq_down_df <- data.frame(
      gene = names(gene_freq_down),
      freq = as.integer(gene_freq_down),
      stringsAsFactors = FALSE
    )
  } else {
    gene_freq_down_df <- data.frame(
      gene = character(),
      freq = integer(),
      stringsAsFactors = FALSE
    )
  }
  
  save_result(gene_freq_up_df, "Conserved_UP_core_gene_frequency.csv", "csv")
  save_result(gene_freq_down_df, "Conserved_DOWN_core_gene_frequency.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # 4F. SIMPLE SUMMARY TABLE
  # ---------------------------------------------------------------------------
  cat("  [4F] Writing hallmark-wise conserved target summary...\n")
  
  summary_up <- if (length(conserved_up_u) > 0) {
    data.frame(
      hallmark = names(conserved_up_u),
      direction = "UP",
      n_dog_le_universe = sapply(dog_le_up_u, length),
      n_conserved_core = sapply(conserved_up_u, length),
      genes = vapply(conserved_up_u, function(x) paste(x, collapse = ";"), character(1)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      hallmark = character(),
      direction = character(),
      n_dog_le_universe = integer(),
      n_conserved_core = integer(),
      genes = character(),
      stringsAsFactors = FALSE
    )
  }
  
  summary_down <- if (length(conserved_down_u) > 0) {
    data.frame(
      hallmark = names(conserved_down_u),
      direction = "DOWN",
      n_dog_le_universe = sapply(dog_le_down_u, length),
      n_conserved_core = sapply(conserved_down_u, length),
      genes = vapply(conserved_down_u, function(x) paste(x, collapse = ";"), character(1)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      hallmark = character(),
      direction = character(),
      n_dog_le_universe = integer(),
      n_conserved_core = integer(),
      genes = character(),
      stringsAsFactors = FALSE
    )
  }
  
  core_summary <- dplyr::bind_rows(summary_up, summary_down)
  save_result(core_summary, "CONSERVED_CORE_TARGETS_by_hallmark.csv", "csv")
  # ---------------------------------------------------------------------------
  # 4G. HALLMARK-WISE DOG LEADING-EDGE vs CAT DIRECTIONAL FISHER TESTS
  # ---------------------------------------------------------------------------
  cat("  [4G] Calculating Fisher overlap support for conserved hallmarks...\n")
  
  make_fisher_table <- function(dog_le_list, cat_set, direction_label, universe_hs) {
    
    if (length(dog_le_list) == 0) {
      return(data.frame(
        hallmark = character(),
        direction = character(),
        overlap = integer(),
        A_size = integer(),
        B_size = integer(),
        U_size = integer(),
        jaccard = numeric(),
        odds_ratio = numeric(),
        p_value = numeric(),
        stringsAsFactors = FALSE
      ))
    }
    
    out <- lapply(names(dog_le_list), function(h) {
      
      A <- intersect(unique(dog_le_list[[h]]), universe_hs)
      B <- intersect(unique(cat_set), universe_hs)
      U <- unique(universe_hs)
      
      a <- length(intersect(A, B))
      b <- length(setdiff(A, B))
      c <- length(setdiff(B, A))
      d <- length(setdiff(U, union(A, B)))
      
      mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      
      ft <- fisher.test(mat, alternative = "greater")
      
      union_size <- length(union(A, B))
      
      data.frame(
        hallmark = h,
        direction = direction_label,
        overlap = a,
        A_size = length(A),
        B_size = length(B),
        U_size = length(U),
        jaccard = ifelse(union_size > 0, a / union_size, NA_real_),
        odds_ratio = unname(ft$estimate),
        p_value = ft$p.value,
        stringsAsFactors = FALSE
      )
    })
    
    dplyr::bind_rows(out)
  }
  
  fish_up <- make_fisher_table(
    dog_le_list = dog_le_up_u,
    cat_set = cat_up_u,
    direction_label = "UP",
    universe_hs = universe_hs
  )
  
  fish_down <- make_fisher_table(
    dog_le_list = dog_le_down_u,
    cat_set = cat_down_u,
    direction_label = "DOWN",
    universe_hs = universe_hs
  )
  
  fish_all <- dplyr::bind_rows(fish_up, fish_down) %>%
    dplyr::mutate(
      padj = p.adjust(p_value, method = "BH")
    )
  
  fish_up <- fish_all %>% dplyr::filter(direction == "UP")
  fish_down <- fish_all %>% dplyr::filter(direction == "DOWN")
  
  save_result(fish_up, "DOG_LE_vs_CAT_UP_Fisher_byHallmark.csv", "csv")
  save_result(fish_down, "DOG_LE_vs_CAT_DOWN_Fisher_byHallmark.csv", "csv")
  
  strict_summary <- dplyr::bind_rows(
    fish_all %>%
      dplyr::filter(overlap > 0) %>%
      dplyr::transmute(hallmark, direction, class = "strict_overlap"),
    
    fish_all %>%
      dplyr::filter(!is.na(padj), padj < 0.05) %>%
      dplyr::transmute(hallmark, direction, class = "fisher_significant")
  ) %>%
    dplyr::distinct()
  
  save_result(strict_summary, "CONSERVED_HALLMARKS_strict_summary.csv", "csv")
  # ---------------------------------------------------------------------------
  # 4Η. EXPORT TO GLOBAL ENV + RETURN
  # ---------------------------------------------------------------------------
  assign("dog_le_up", dog_le_up, envir = .GlobalEnv)
  assign("dog_le_down", dog_le_down, envir = .GlobalEnv)
  assign("dog_le_up_u", dog_le_up_u, envir = .GlobalEnv)
  assign("dog_le_down_u", dog_le_down_u, envir = .GlobalEnv)
  
  assign("conserved_up_u", conserved_up_u, envir = .GlobalEnv)
  assign("conserved_down_u", conserved_down_u, envir = .GlobalEnv)
  
  assign("core_up_genes_u", core_up_genes_u, envir = .GlobalEnv)
  assign("core_down_genes_u", core_down_genes_u, envir = .GlobalEnv)
  
  assign("gene_freq_up_df", gene_freq_up_df, envir = .GlobalEnv)
  assign("gene_freq_down_df", gene_freq_down_df, envir = .GlobalEnv)
  
  cat("  ✓ STEP 4 COMPLETE\n\n")
  
  invisible(list(
    dog_le_up = dog_le_up,
    dog_le_down = dog_le_down,
    dog_le_up_u = dog_le_up_u,
    dog_le_down_u = dog_le_down_u,
    conserved_up_u = conserved_up_u,
    conserved_down_u = conserved_down_u,
    core_up_genes_u = core_up_genes_u,
    core_down_genes_u = core_down_genes_u,
    gene_freq_up_df = gene_freq_up_df,
    gene_freq_down_df = gene_freq_down_df,
    core_summary = core_summary
  ))
}