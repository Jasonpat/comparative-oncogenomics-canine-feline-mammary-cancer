################################################################################
# STEP 2: CAT GENE LISTS, HUMAN MAPPING & HALLMARK ORA (fixed universe)
# Based on original cat scripts, cleaned and merged
#
# Outputs:
#   - CAT_mapping_ens_to_humanSymbol_1to1_highconf.rds
#   - CAT_up_mapped_humanSymbols.txt
#   - CAT_down_mapped_humanSymbols.txt
#   - CAT_up_universeFixed.txt
#   - CAT_down_universeFixed.txt
#   - CAT_Hallmark_ORA_UP_universeFixed.csv
#   - CAT_Hallmark_ORA_DOWN_universeFixed.csv
################################################################################

module_cat_hallmarks <- function() {
  
  # ---------------------------------------------------------------------------
  # 2A. INPUT CHECKS
  # ---------------------------------------------------------------------------
  cat("  [2A] Checking cat input files...\n")
  
  if (!file.exists(PATH_CONFIG$cat_up)) {
    stop("Cat UP file not found: ", PATH_CONFIG$cat_up)
  }
  if (!file.exists(PATH_CONFIG$cat_down)) {
    stop("Cat DOWN file not found: ", PATH_CONFIG$cat_down)
  }
  
  if (!exists("mapping_1to1", envir = .GlobalEnv)) {
    stop("Dog mapping_1to1 not found. Run module_dog_gsea() first.")
  }
  
  mapping_1to1 <- get("mapping_1to1", envir = .GlobalEnv)
  
  # ---------------------------------------------------------------------------
  # 2B. READ & PARSE CAT GENE LISTS
  # ---------------------------------------------------------------------------
  cat("  [2B] Reading cat gene lists...\n")
  
  read_gene_list <- function(file) {
    x <- readLines(file, warn = FALSE)
    x <- paste(x, collapse = " ")
    x <- gsub("[\\r\\n\\t]", " ", x)
    # perl=TRUE is required so that \s is recognised as the whitespace class.
    # Without it, R's default ERE engine treats \s as a literal 's' character,
    # causing space-delimited gene lists to not be split correctly.
    x <- strsplit(x, "[,;\\s]+", perl = TRUE)[[1]]
    x <- trimws(x)
    x <- x[x != ""]
    unique(x)
  }
  
  cat_up_raw <- read_gene_list(PATH_CONFIG$cat_up)
  cat_down_raw <- read_gene_list(PATH_CONFIG$cat_down)
  
  if (length(cat_up_raw) == 0) stop("CAT_up file produced an empty gene list.")
  if (length(cat_down_raw) == 0) stop("CAT_down file produced an empty gene list.")
  
  cat("      UP genes:", length(cat_up_raw), "\n")
  cat("      DOWN genes:", length(cat_down_raw), "\n")
  
  # Split symbols vs ENSFCAG IDs
  is_ensfcag <- function(x) grepl("^ENSFCAG", x, ignore.case = FALSE)
  
  cat_up_sym <- cat_up_raw[!is_ensfcag(cat_up_raw)]
  cat_up_ens <- cat_up_raw[ is_ensfcag(cat_up_raw)]
  
  cat_down_sym <- cat_down_raw[!is_ensfcag(cat_down_raw)]
  cat_down_ens <- cat_down_raw[ is_ensfcag(cat_down_raw)]
  
  cat("      UP (symbols:", length(cat_up_sym), "| ENSFCAG:", length(cat_up_ens), ")\n")
  cat("      DOWN (symbols:", length(cat_down_sym), "| ENSFCAG:", length(cat_down_ens), ")\n")
  
  # ---------------------------------------------------------------------------
  # 2C. MAP CAT ENSFCAG IDs TO HUMAN SYMBOLS (biomaRt, cached)
  # ---------------------------------------------------------------------------
  cat("  [2C] Mapping cat ENSFCAG IDs to human symbols...\n")
  
  ens_to_map <- unique(c(cat_up_ens, cat_down_ens))
  
  if (length(ens_to_map) > 0) {
    map_cat <- load_cached("CAT_biomart_mapping_raw", function() {
      
      cat_mart <- useEnsembl(
        biomart = "genes",
        dataset = "fcatus_gene_ensembl",
        version = PARAM_CONFIG$ensembl_version
      )
      
      getBM(
        attributes = c(
          "ensembl_gene_id",
          "hsapiens_homolog_associated_gene_name",
          "hsapiens_homolog_orthology_type",
          "hsapiens_homolog_orthology_confidence"
        ),
        filters = "ensembl_gene_id",
        values = ens_to_map,
        mart = cat_mart
      )
    })
    
    if (is.null(map_cat) || nrow(map_cat) == 0) {
      warning("Cat biomaRt mapping returned no rows. Only symbol inputs will be retained.")
      map_cat <- data.frame(
        cat_ens = character(),
        human_symbol = character(),
        orthology_type = character(),
        orthology_conf = numeric(),
        stringsAsFactors = FALSE
      )
    } else {
      colnames(map_cat) <- c("cat_ens", "human_symbol", "orthology_type", "orthology_conf")
    }
  } else {
    map_cat <- data.frame(
      cat_ens = character(),
      human_symbol = character(),
      orthology_type = character(),
      orthology_conf = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  map_cat_1to1 <- map_cat %>%
    dplyr::filter(
      orthology_type == PARAM_CONFIG$ortholog_type_filter,
      orthology_conf == PARAM_CONFIG$ortholog_conf_filter,
      !is.na(human_symbol), human_symbol != ""
    ) %>%
    dplyr::distinct(cat_ens, .keep_all = TRUE)
  
  cat("      Raw mapping rows:", nrow(map_cat), "\n")
  cat("      1:1 high-confidence mapped cat IDs:", nrow(map_cat_1to1), "\n")
  
  cat_to_hs <- setNames(map_cat_1to1$human_symbol, map_cat_1to1$cat_ens)
  
  save_result(map_cat_1to1, "CAT_mapping_ens_to_humanSymbol_1to1_highconf.rds", "rds")
  
  # ---------------------------------------------------------------------------
  # 2D. COMBINE CAT SYMBOLS + MAPPED ENSFCAG IDs
  # ---------------------------------------------------------------------------
  cat("  [2D] Combining symbols with mapped orthologs...\n")
  
  cat_up_mapped <- unique(c(cat_up_sym, unname(cat_to_hs[cat_up_ens])))
  cat_down_mapped <- unique(c(cat_down_sym, unname(cat_to_hs[cat_down_ens])))
  
  cat_up_mapped <- cat_up_mapped[!is.na(cat_up_mapped) & cat_up_mapped != ""]
  cat_down_mapped <- cat_down_mapped[!is.na(cat_down_mapped) & cat_down_mapped != ""]
  
  cat("      Final mapped UP genes:", length(cat_up_mapped), "\n")
  cat("      Final mapped DOWN genes:", length(cat_down_mapped), "\n")
  
  save_result(cat_up_mapped, "CAT_up_mapped_humanSymbols.txt", "txt")
  save_result(cat_down_mapped, "CAT_down_mapped_humanSymbols.txt", "txt")
  
  # ---------------------------------------------------------------------------
  # 2E. DEFINE FIXED HUMAN UNIVERSE
  # universe_hs = dog human symbol universe ∩ hallmark genes
  # ---------------------------------------------------------------------------
  cat("  [2E] Defining fixed human universe for ORA...\n")
  
  hall_sym <- msigdbr(species = "Homo sapiens", collection = "H")
  term2gene <- hall_sym[, c("gs_name", "gene_symbol")]
  
  hall_genes <- unique(term2gene$gene_symbol)
  
  dog_universe_hs <- sort(unique(mapping_1to1$human_symbol))
  dog_universe_hs <- dog_universe_hs[!is.na(dog_universe_hs) & dog_universe_hs != ""]
  
  universe_hs <- sort(intersect(dog_universe_hs, hall_genes))
  
  if (length(universe_hs) == 0) {
    stop("Fixed universe_hs is empty.")
  }
  
  cat("      Fixed universe size:", length(universe_hs), "\n")
  
  # Restrict cat lists to the same universe
  cat_up_u <- intersect(cat_up_mapped, universe_hs)
  cat_down_u <- intersect(cat_down_mapped, universe_hs)
  
  cat("      UP genes before universe filter:", length(cat_up_mapped), "\n")
  cat("      UP genes after universe filter:", length(cat_up_u), "\n")
  cat("      DOWN genes before universe filter:", length(cat_down_mapped), "\n")
  cat("      DOWN genes after universe filter:", length(cat_down_u), "\n")
  
  if (length(cat_up_u) == 0) warning("cat_up_u is empty after universe restriction.")
  if (length(cat_down_u) == 0) warning("cat_down_u is empty after universe restriction.")
  
  save_result(cat_up_u, "CAT_up_universeFixed.txt", "txt")
  save_result(cat_down_u, "CAT_down_universeFixed.txt", "txt")
  
  # ---------------------------------------------------------------------------
  # 2F. HALLMARK ORA WITH FIXED UNIVERSE
  # ---------------------------------------------------------------------------
  cat("  [2F] Running Hallmark ORA with fixed universe...\n")
  
  # Reviewer-safe ORA settings:
  # pvalueCutoff/qvalueCutoff are set to 1 so clusterProfiler returns the
  # complete tested Hallmark result table. Significance is applied downstream
  # using PARAM_CONFIG$cat_ora_padj, preventing hidden pre-filtering by package
  # defaults and preserving reproducibility if thresholds are changed later.
  enr_up_u <- enricher(
    gene          = cat_up_u,
    TERM2GENE     = term2gene,
    universe      = universe_hs,
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    pAdjustMethod = "BH"
  )
  
  enr_down_u <- enricher(
    gene          = cat_down_u,
    TERM2GENE     = term2gene,
    universe      = universe_hs,
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    pAdjustMethod = "BH"
  )
  
  up_tab_u <- as.data.frame(enr_up_u)
  down_tab_u <- as.data.frame(enr_down_u)
  
  if (nrow(up_tab_u) == 0) {
    up_tab_u <- data.frame()
    cat("      No enriched UP hallmarks returned.\n")
  }
  if (nrow(down_tab_u) == 0) {
    down_tab_u <- data.frame()
    cat("      No enriched DOWN hallmarks returned.\n")
  }
  
  sig_up <- if (nrow(up_tab_u) > 0) sum(up_tab_u$p.adjust < PARAM_CONFIG$cat_ora_padj, na.rm = TRUE) else 0
  sig_down <- if (nrow(down_tab_u) > 0) sum(down_tab_u$p.adjust < PARAM_CONFIG$cat_ora_padj, na.rm = TRUE) else 0
  
  cat("      Significant UP hallmarks (FDR <", PARAM_CONFIG$cat_ora_padj, "): ", sig_up, "\n", sep = "")
  cat("      Significant DOWN hallmarks (FDR <", PARAM_CONFIG$cat_ora_padj, "): ", sig_down, "\n", sep = "")
  
  save_result(up_tab_u, "CAT_Hallmark_ORA_UP_universeFixed.csv", "csv")
  save_result(down_tab_u, "CAT_Hallmark_ORA_DOWN_universeFixed.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # 2G. EXPORT TO GLOBAL ENV + RETURN
  # ---------------------------------------------------------------------------
  assign("map_cat_1to1", map_cat_1to1, envir = .GlobalEnv)
  assign("cat_up_mapped", cat_up_mapped, envir = .GlobalEnv)
  assign("cat_down_mapped", cat_down_mapped, envir = .GlobalEnv)
  assign("cat_up_u", cat_up_u, envir = .GlobalEnv)
  assign("cat_down_u", cat_down_u, envir = .GlobalEnv)
  assign("up_tab_u", up_tab_u, envir = .GlobalEnv)
  assign("down_tab_u", down_tab_u, envir = .GlobalEnv)
  assign("universe_hs", universe_hs, envir = .GlobalEnv)
  assign("term2gene", term2gene, envir = .GlobalEnv)
  
  cat("  ✓ STEP 2 COMPLETE\n\n")
  
  invisible(list(
    map_cat_1to1 = map_cat_1to1,
    cat_up_mapped = cat_up_mapped,
    cat_down_mapped = cat_down_mapped,
    cat_up_u = cat_up_u,
    cat_down_u = cat_down_u,
    up_tab_u = up_tab_u,
    down_tab_u = down_tab_u,
    universe_hs = universe_hs,
    term2gene = term2gene
  ))
}