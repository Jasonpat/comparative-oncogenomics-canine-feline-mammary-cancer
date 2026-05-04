################################################################################
# OPTIONAL VALIDATION: NETWORK CONTEXT OF PRIORITISED CONSERVED TARGETS
# STRING PPI + igraph centrality analysis with one-sided empirical permutation
# testing and Benjamini-Hochberg correction.
#
# Fixes:
#   - robust STRING edge filtering before igraph::graph_from_data_frame()
#   - removes edges with NA from/to/weight
#   - writes NETWORK_STRING_edges_for_igraph.csv for debugging
#   - uses exact top-N prioritised targets by priority_rank
#   - no ggraph dependency
################################################################################

module_network_validation_prioritised_targets <- function(
    score_threshold = 700,
    string_version = "12.0",
    species = 9606,
    top_n = NULL,
    n_permutations = 10000,
    seed = 1,
    save_pdf = FALSE
) {
  
  cat("\n")
  cat("======================================================================\n")
  cat("NETWORK VALIDATION: PRIORITISED TARGETS WITHIN THE CONSERVED CORE\n")
  cat("======================================================================\n\n")
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tibble)
    library(stringr)
    library(ggplot2)
    library(igraph)
    library(STRINGdb)
  })
  
  set.seed(seed)
  options(timeout = max(600, getOption("timeout")))
  
  # ---------------------------------------------------------------------------
  # 1. Paths and parameters
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
  
  top_n <- if (is.null(top_n)) {
    if (exists("PARAM_CONFIG", envir = .GlobalEnv)) {
      get("PARAM_CONFIG", envir = .GlobalEnv)$top_n_targets
    } else {
      15
    }
  } else {
    top_n
  }
  
  network_dir <- file.path(results_dir, "NETWORK_validation")
  plot_dir <- file.path(network_dir, "plots")
  table_dir <- file.path(network_dir, "tables")
  string_cache_dir <- file.path(cache_dir, "STRING_cache")
  
  dir.create(network_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(string_cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  save_table <- function(x, filename) {
    out_path <- file.path(table_dir, filename)
    readr::write_csv(x, out_path)
    cat("  ✓ Saved: ", normalizePath(out_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
    invisible(out_path)
  }
  
  save_plot <- function(p, filename_base, width = 7.5, height = 5.5) {
    filename_base <- stringr::str_remove(filename_base, "\\.png$")
    filename_base <- stringr::str_remove(filename_base, "\\.pdf$")
    
    png_path <- file.path(plot_dir, paste0(filename_base, ".png"))
    pdf_path <- file.path(plot_dir, paste0(filename_base, ".pdf"))
    
    ggplot2::ggsave(
      filename = png_path,
      plot = p,
      width = width,
      height = height,
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
            plot = p,
            width = width,
            height = height,
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
            "\nPNG output was saved successfully.",
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
  
  as_bool_safe <- function(x) {
    if (is.logical(x)) return(dplyr::coalesce(x, FALSE))
    if (is.numeric(x)) return(dplyr::coalesce(x != 0, FALSE))
    x_chr <- tolower(trimws(as.character(x)))
    out <- x_chr %in% c("true", "t", "1", "yes", "y")
    out[is.na(x_chr)] <- FALSE
    out
  }
  
  # ---------------------------------------------------------------------------
  # 2. Load prioritisation and target master
  # ---------------------------------------------------------------------------
  
  cat("  [N1] Loading prioritisation and conserved-core data...\n")
  
  scored_file <- file.path(results_dir, "target_prioritisation_scored.csv")
  master_file <- file.path(results_dir, "target_master.csv")
  
  if (!file.exists(scored_file)) {
    stop("Missing target_prioritisation_scored.csv. Run Step 6 first.", call. = FALSE)
  }
  if (!file.exists(master_file)) {
    stop("Missing target_master.csv. Run Step 5 first.", call. = FALSE)
  }
  
  target_scored <- readr::read_csv(scored_file, show_col_types = FALSE)
  target_master <- readr::read_csv(master_file, show_col_types = FALSE)
  
  if (!"target_symbol" %in% names(target_scored)) {
    stop("target_prioritisation_scored.csv must contain target_symbol.", call. = FALSE)
  }
  if (!"target_symbol" %in% names(target_master)) {
    stop("target_master.csv must contain target_symbol.", call. = FALSE)
  }
  
  if (!"priority_rank" %in% names(target_scored)) {
    cat("      NOTE: priority_rank not found; reconstructing deterministic ordering.\n")
    
    if (all(c("consensus_priority_score", "conservation_score", "max_drug_phase", "n_drugs") %in% names(target_scored))) {
      target_scored <- target_scored %>%
        arrange(
          desc(suppressWarnings(as.numeric(consensus_priority_score))),
          desc(suppressWarnings(as.numeric(conservation_score))),
          desc(suppressWarnings(as.numeric(max_drug_phase))),
          desc(suppressWarnings(as.numeric(n_drugs))),
          target_symbol
        )
    } else if ("consensus_priority_score" %in% names(target_scored)) {
      target_scored <- target_scored %>%
        arrange(desc(suppressWarnings(as.numeric(consensus_priority_score))), target_symbol)
    } else {
      target_scored <- target_scored %>%
        arrange(target_symbol)
    }
    
    target_scored <- target_scored %>%
      mutate(
        priority_rank = row_number(),
        score_rank = if ("consensus_priority_score" %in% names(.)) {
          dense_rank(desc(suppressWarnings(as.numeric(consensus_priority_score))))
        } else {
          row_number()
        }
      )
  }
  
  if (!"score_rank" %in% names(target_scored)) {
    target_scored <- target_scored %>%
      mutate(
        score_rank = if ("consensus_priority_score" %in% names(.)) {
          dense_rank(desc(suppressWarnings(as.numeric(consensus_priority_score))))
        } else {
          priority_rank
        }
      )
  }
  
  for (col in c("consensus_priority_score", "conservation_score", "max_drug_phase", "n_drugs", "n_indications")) {
    if (!col %in% names(target_scored)) target_scored[[col]] <- 0
  }
  for (col in c("any_sm", "any_approved", "any_onco_drug", "any_breast_drug")) {
    if (!col %in% names(target_scored)) target_scored[[col]] <- FALSE
  }
  
  target_scored <- target_scored %>%
    mutate(
      target_symbol = as.character(target_symbol),
      priority_rank = suppressWarnings(as.numeric(priority_rank)),
      score_rank = suppressWarnings(as.numeric(score_rank)),
      consensus_priority_score = suppressWarnings(as.numeric(consensus_priority_score)),
      conservation_score = suppressWarnings(as.numeric(conservation_score)),
      max_drug_phase = suppressWarnings(as.numeric(max_drug_phase)),
      n_drugs = suppressWarnings(as.numeric(n_drugs)),
      n_indications = suppressWarnings(as.numeric(n_indications)),
      any_sm = as_bool_safe(any_sm),
      any_approved = as_bool_safe(any_approved),
      any_onco_drug = as_bool_safe(any_onco_drug),
      any_breast_drug = as_bool_safe(any_breast_drug)
    ) %>%
    arrange(priority_rank, target_symbol) %>%
    distinct(target_symbol, .keep_all = TRUE)
  
  core_symbols <- sort(unique(as.character(target_master$target_symbol)))
  core_symbols <- core_symbols[!is.na(core_symbols) & core_symbols != ""]
  
  prioritised_symbols <- target_scored %>%
    arrange(priority_rank, target_symbol) %>%
    slice_head(n = top_n) %>%
    pull(target_symbol) %>%
    unique()
  
  prioritised_label <- paste0("Prioritised top ", top_n)
  other_label <- "Other conserved core"
  
  cat("      Conserved-core targets:", length(core_symbols), "\n")
  cat("      Prioritised targets used:", length(prioritised_symbols), "\n")
  cat("      STRING score threshold:", score_threshold, "\n")
  
  priority_table <- target_scored %>%
    mutate(
      network_group = ifelse(
        target_symbol %in% prioritised_symbols,
        prioritised_label,
        other_label
      )
    )
  
  save_table(priority_table, "NETWORK_priority_table_used.csv")
  
  # ---------------------------------------------------------------------------
  # 3. Method parameters
  # ---------------------------------------------------------------------------
  
  method_parameters <- tibble(
    parameter = c(
      "string_version",
      "species",
      "score_threshold",
      "top_n",
      "n_permutations",
      "seed",
      "statistical_test",
      "p_adjustment_method",
      "selection_rule",
      "null_model",
      "primary_metrics",
      "pdf_export"
    ),
    value = c(
      as.character(string_version),
      as.character(species),
      as.character(score_threshold),
      as.character(top_n),
      as.character(n_permutations),
      as.character(seed),
      "One-sided empirical permutation enrichment test; p = proportion of random sets with mean >= observed mean",
      "Benjamini-Hochberg FDR correction across tested network metrics",
      "Exact top-N targets by priority_rank from target_prioritisation_scored.csv",
      "Random target sets of equal size sampled without replacement from the mapped conserved-core universe",
      "degree; weighted_degree; betweenness; harmonic_centrality; eigenvector_centrality",
      as.character(save_pdf)
    )
  )
  
  save_table(method_parameters, "NETWORK_method_parameters.csv")
  
  # ---------------------------------------------------------------------------
  # 4. Map symbols to STRING IDs
  # ---------------------------------------------------------------------------
  
  cat("  [N2] Mapping conserved-core targets to STRING IDs...\n")
  
  string_db <- STRINGdb::STRINGdb$new(
    version = string_version,
    species = species,
    score_threshold = score_threshold,
    input_directory = string_cache_dir
  )
  
  input_genes <- data.frame(
    target_symbol = core_symbols,
    stringsAsFactors = FALSE
  )
  
  string_map <- string_db$map(
    input_genes,
    "target_symbol",
    removeUnmappedRows = FALSE
  ) %>%
    as_tibble() %>%
    rename(string_id = STRING_id) %>%
    mutate(
      mapped_to_string = !is.na(string_id) & string_id != "",
      network_group = ifelse(
        target_symbol %in% prioritised_symbols,
        prioritised_label,
        other_label
      )
    )
  
  save_table(string_map, "NETWORK_STRING_mapping.csv")
  
  mapped <- string_map %>%
    filter(mapped_to_string) %>%
    arrange(target_symbol, string_id) %>%
    distinct(target_symbol, .keep_all = TRUE) %>%
    select(target_symbol, string_id, network_group)
  
  if (nrow(mapped) == 0) {
    stop("No conserved-core targets mapped to STRING IDs.", call. = FALSE)
  }
  
  cat("      Mapped to STRING:", nrow(mapped), "/", length(core_symbols), "\n")
  
  # ---------------------------------------------------------------------------
  # 5. Retrieve STRING interactions and build graph
  # ---------------------------------------------------------------------------
  
  cat("  [N3] Retrieving STRING interactions and building core network...\n")

  # STRINGdb::get_interactions() has a bug in v2.22.0: when the aliases/info
  # files are already cached, it internally calls graph_from_data_frame() and
  # crashes before returning any data. We bypass it by reading the
  # protein.links file directly — which is methodologically identical (same
  # STRING v12.0 data, same species, same score threshold).
  #
  # The protein.links file is NOT downloaded by STRINGdb$map() — only aliases
  # and info are. So we download it explicitly here if it is missing.
  #
  # Methods text: "STRING v12.0 PPI data (Homo sapiens, taxon 9606) were
  # accessed via the STRINGdb Bioconductor package (v2.22.0) and filtered at
  # combined score >= <score_threshold>."

  links_file <- file.path(
    string_cache_dir,
    paste0(species, ".protein.links.v", string_version, ".txt.gz")
  )

  if (!file.exists(links_file)) {
    cat("      protein.links file not in cache — downloading from STRING...\n")
    links_url <- paste0(
      "https://stringdb-downloads.org/download/protein.links.v",
      string_version, "/", species,
      ".protein.links.v", string_version, ".txt.gz"
    )
    tryCatch(
      {
        utils::download.file(
          url      = links_url,
          destfile = links_file,
          mode     = "wb",
          quiet    = FALSE
        )
        cat("      Download complete.\n")
      },
      error = function(e) {
        stop(
          "Could not download STRING links file.\n",
          "URL tried: ", links_url, "\n",
          "Error: ", conditionMessage(e), "\n",
          "Alternative: manually place the file at:\n  ", links_file,
          call. = FALSE
        )
      }
    )
  } else {
    cat("      Using cached STRING links file.\n")
  }

  cat("      Reading STRING links (this may take ~30s)...\n")
  links_all <- readr::read_delim(
    links_file,
    delim      = " ",
    col_types  = readr::cols(
      protein1       = readr::col_character(),
      protein2       = readr::col_character(),
      combined_score = readr::col_double()
    ),
    progress = FALSE
  )

  cat("      Total STRING links loaded:", nrow(links_all), "\n")

  # Filter to interactions between our mapped core targets at score threshold
  mapped_ids <- unique(mapped$string_id)   # format: "9606.ENSP..."

  interactions_raw <- links_all %>%
    filter(
      protein1 %in% mapped_ids,
      protein2 %in% mapped_ids,
      combined_score >= score_threshold
    ) %>%
    rename(from = protein1, to = protein2)

  cat("      Interactions between core targets (score >=",
      score_threshold, "):", nrow(interactions_raw), "\n")

  save_table(interactions_raw, "NETWORK_STRING_interactions_raw.csv")

  # Build symbol lookup
  id_to_symbol <- setNames(mapped$target_symbol, mapped$string_id)

  if (nrow(interactions_raw) == 0) {
    warning(
      "No STRING interactions found between core targets at score >= ",
      score_threshold, ". Try lowering score_threshold (e.g. 400).",
      call. = FALSE
    )
    edge_tbl <- tibble(
      from_symbol    = character(),
      to_symbol      = character(),
      combined_score = numeric()
    )
  } else {
    edge_tbl <- interactions_raw %>%
      transmute(
        from_string    = as.character(from),
        to_string      = as.character(to),
        from_symbol    = unname(id_to_symbol[from_string]),
        to_symbol      = unname(id_to_symbol[to_string]),
        combined_score = suppressWarnings(as.numeric(combined_score))
      ) %>%
      filter(
        !is.na(from_symbol), from_symbol != "",
        !is.na(to_symbol),   to_symbol   != "",
        from_symbol != to_symbol,
        !is.na(combined_score), is.finite(combined_score)
      ) %>%
      mutate(
        a = pmin(from_symbol, to_symbol),
        b = pmax(from_symbol, to_symbol)
      ) %>%
      group_by(a, b) %>%
      summarise(combined_score = max(combined_score, na.rm = TRUE), .groups = "drop") %>%
      transmute(from_symbol = a, to_symbol = b, combined_score = combined_score)

    cat("      Resolved unique undirected edges:", nrow(edge_tbl), "\n")
  }

  save_table(edge_tbl, "NETWORK_STRING_edges_core.csv")

  edge_for_graph <- edge_tbl %>%
    transmute(
      from   = as.character(from_symbol),
      to     = as.character(to_symbol),
      weight = suppressWarnings(as.numeric(combined_score)) / 1000
    ) %>%
    filter(
      !is.na(from), from != "",
      !is.na(to),   to   != "",
      from != to,
      !is.na(weight), is.finite(weight)
    ) %>%
    distinct(from, to, .keep_all = TRUE)

  save_table(edge_for_graph, "NETWORK_STRING_edges_for_igraph.csv")

  # ── Build vertex table (slim join to avoid fan-out duplicates) ──────────────
  priority_slim <- priority_table %>%
    select(
      target_symbol, priority_rank, score_rank, consensus_priority_score,
      conservation_score, max_drug_phase, n_drugs, n_indications,
      any_sm, any_approved, any_onco_drug, any_breast_drug
    ) %>%
    distinct(target_symbol, .keep_all = TRUE)

  node_tbl <- mapped %>%
    select(target_symbol, string_id) %>%
    left_join(priority_slim, by = "target_symbol") %>%
    mutate(
      name = as.character(target_symbol),
      network_group = ifelse(
        target_symbol %in% prioritised_symbols,
        prioritised_label,
        other_label
      )
    ) %>%
    filter(!is.na(name), name != "") %>%
    distinct(name, .keep_all = TRUE)

  if (nrow(node_tbl) == 0) {
    stop("Vertex table is empty after cleaning. Check STRING mapping.", call. = FALSE)
  }

  # Drop edges whose endpoints are absent from vertex table
  edge_symbols       <- unique(c(edge_for_graph$from, edge_for_graph$to))
  missing_from_nodes <- setdiff(edge_symbols, node_tbl$name)
  if (length(missing_from_nodes) > 0) {
    cat("      NOTE: removing", length(missing_from_nodes),
        "edge endpoint(s) absent from vertex table.\n")
    edge_for_graph <- edge_for_graph %>%
      filter(from %in% node_tbl$name, to %in% node_tbl$name)
    save_table(edge_for_graph, "NETWORK_STRING_edges_for_igraph.csv")
  }

  # Final safety filter before igraph
  edge_for_graph <- edge_for_graph %>%
    filter(
      !is.na(from), from != "",
      !is.na(to),   to   != "",
      from != to,
      !is.na(weight), is.finite(weight)
    ) %>%
    distinct(from, to, .keep_all = TRUE)

  vertex_df <- node_tbl %>%
    select(name, everything()) %>%
    distinct(name, .keep_all = TRUE)

  cat("      Final nodes:", nrow(vertex_df),
      "| Final edges:", nrow(edge_for_graph), "\n")

  g <- igraph::graph_from_data_frame(
    d        = edge_for_graph,
    directed = FALSE,
    vertices = vertex_df
  )

  cat("      Network nodes:", igraph::vcount(g), "\n")
  cat("      Network edges:", igraph::ecount(g), "\n")
  
  # ---------------------------------------------------------------------------
  # 6. Node-level centrality metrics
  # ---------------------------------------------------------------------------
  
  cat("  [N4] Calculating node-level network metrics...\n")
  
  if (igraph::ecount(g) > 0) {
    E(g)$weight <- ifelse(is.na(E(g)$weight), 1, E(g)$weight)
    
    degree_raw <- igraph::degree(g, mode = "all", normalized = FALSE)
    strength_raw <- igraph::strength(g, mode = "all", weights = E(g)$weight)
    betweenness_raw <- igraph::betweenness(g, directed = FALSE, normalized = TRUE, weights = NA)
    
    harmonic_raw <- igraph::harmonic_centrality(
      g,
      mode = "all",
      weights = NA,
      normalized = TRUE
    )
    
    eigen_raw <- tryCatch(
      igraph::eigen_centrality(g, directed = FALSE, weights = E(g)$weight)$vector,
      error = function(e) {
        out <- rep(NA_real_, igraph::vcount(g))
        names(out) <- V(g)$name
        out
      }
    )
  } else {
    degree_raw <- rep(0, igraph::vcount(g)); names(degree_raw) <- V(g)$name
    strength_raw <- degree_raw
    betweenness_raw <- degree_raw
    harmonic_raw <- degree_raw
    eigen_raw <- degree_raw
  }
  
  comp <- igraph::components(g)
  comp_size <- comp$csize[comp$membership]
  
  node_metrics <- tibble(
    target_symbol = V(g)$name,
    network_group = V(g)$network_group,
    priority_rank = suppressWarnings(as.numeric(V(g)$priority_rank)),
    score_rank = suppressWarnings(as.numeric(V(g)$score_rank)),
    consensus_priority_score = suppressWarnings(as.numeric(V(g)$consensus_priority_score)),
    conservation_score = suppressWarnings(as.numeric(V(g)$conservation_score)),
    max_drug_phase = suppressWarnings(as.numeric(V(g)$max_drug_phase)),
    n_drugs = suppressWarnings(as.numeric(V(g)$n_drugs)),
    n_indications = suppressWarnings(as.numeric(V(g)$n_indications)),
    any_sm = V(g)$any_sm %in% TRUE,
    any_approved = V(g)$any_approved %in% TRUE,
    any_onco_drug = V(g)$any_onco_drug %in% TRUE,
    any_breast_drug = V(g)$any_breast_drug %in% TRUE,
    degree = as.numeric(degree_raw[V(g)$name]),
    weighted_degree = as.numeric(strength_raw[V(g)$name]),
    betweenness = as.numeric(betweenness_raw[V(g)$name]),
    harmonic_centrality = as.numeric(harmonic_raw[V(g)$name]),
    eigenvector_centrality = as.numeric(eigen_raw[V(g)$name]),
    component_id = as.integer(comp$membership),
    component_size = as.integer(comp_size)
  ) %>%
    arrange(priority_rank, target_symbol)
  
  save_table(node_metrics, "NETWORK_node_metrics.csv")
  
  # ---------------------------------------------------------------------------
  # 7. Group summaries and one-sided permutation tests
  # ---------------------------------------------------------------------------
  
  cat("  [N5] Comparing prioritised targets with random conserved-core target sets...\n")
  
  group_summary <- node_metrics %>%
    group_by(network_group) %>%
    summarise(
      n_targets = n(),
      mean_degree = mean(degree, na.rm = TRUE),
      median_degree = median(degree, na.rm = TRUE),
      mean_weighted_degree = mean(weighted_degree, na.rm = TRUE),
      median_weighted_degree = median(weighted_degree, na.rm = TRUE),
      mean_betweenness = mean(betweenness, na.rm = TRUE),
      median_betweenness = median(betweenness, na.rm = TRUE),
      mean_harmonic_centrality = mean(harmonic_centrality, na.rm = TRUE),
      median_harmonic_centrality = median(harmonic_centrality, na.rm = TRUE),
      mean_eigenvector_centrality = mean(eigenvector_centrality, na.rm = TRUE),
      median_eigenvector_centrality = median(eigenvector_centrality, na.rm = TRUE),
      .groups = "drop"
    )
  
  save_table(group_summary, "NETWORK_group_summary.csv")
  
  prioritised_set <- node_metrics %>%
    filter(network_group == prioritised_label) %>%
    pull(target_symbol)
  
  universe_set <- node_metrics$target_symbol
  
  if (length(prioritised_set) == 0) {
    stop("No prioritised targets were represented in the mapped STRING network.", call. = FALSE)
  }
  
  metric_vector <- function(metric_name) {
    x <- node_metrics[[metric_name]]
    names(x) <- node_metrics$target_symbol
    x
  }
  
  tested_metrics <- c(
    "degree",
    "weighted_degree",
    "betweenness",
    "harmonic_centrality",
    "eigenvector_centrality"
  )
  
  observed_metrics <- tibble(
    metric = tested_metrics,
    observed_mean = vapply(
      tested_metrics,
      function(m) mean(metric_vector(m)[prioritised_set], na.rm = TRUE),
      numeric(1)
    )
  )
  
  permute_metric <- function(metric_name, B = 10000) {
    x <- metric_vector(metric_name)
    replicate(B, {
      sampled <- sample(universe_set, length(prioritised_set), replace = FALSE)
      mean(x[sampled], na.rm = TRUE)
    })
  }
  
  # One-sided upper-tail empirical p-value.
  # Directional hypothesis: prioritised conserved targets are more central than
  # random same-sized target sets from the mapped conserved-core universe.
  perm_results_raw <- lapply(tested_metrics, function(m) {
    null_vals <- permute_metric(m, B = n_permutations)
    obs <- observed_metrics$observed_mean[observed_metrics$metric == m]
    n_valid <- sum(!is.na(null_vals))
    null_mean <- mean(null_vals, na.rm = TRUE)
    null_sd <- stats::sd(null_vals, na.rm = TRUE)
    
    tibble(
      metric = m,
      test_statistic = "mean centrality of prioritised target set",
      test_direction = "greater_or_equal",
      observed_mean = obs,
      null_mean = null_mean,
      null_sd = null_sd,
      effect_delta = obs - null_mean,
      effect_ratio = ifelse(null_mean > 0, obs / null_mean, NA_real_),
      empirical_p_greater_equal = (sum(null_vals >= obs, na.rm = TRUE) + 1) / (n_valid + 1),
      empirical_z = ifelse(null_sd > 0, (obs - null_mean) / null_sd, NA_real_),
      n_permutations = n_permutations,
      null_n_valid = n_valid
    )
  }) %>%
    bind_rows()
  
  perm_results <- perm_results_raw %>%
    mutate(
      p_adj_BH = stats::p.adjust(empirical_p_greater_equal, method = "BH"),
      passes_raw_0.05 = empirical_p_greater_equal < 0.05,
      passes_BH_0.05 = p_adj_BH < 0.05,
      interpretation = case_when(
        passes_BH_0.05 ~ "Supported after BH correction",
        passes_raw_0.05 ~ "Supported before BH correction only",
        empirical_p_greater_equal < 0.10 ~ "Weak trend",
        TRUE ~ "Not enriched"
      )
    )
  
  save_table(perm_results, "NETWORK_permutation_summary.csv")
  save_table(perm_results, "NETWORK_permutation_summary_with_BH.csv")
  
  manuscript_perm_table <- perm_results %>%
    transmute(
      metric,
      observed_mean = round(observed_mean, 4),
      null_mean = round(null_mean, 4),
      empirical_one_sided_p = signif(empirical_p_greater_equal, 3),
      BH_adjusted_p = signif(p_adj_BH, 3),
      empirical_z = round(empirical_z, 3),
      interpretation
    )
  
  save_table(manuscript_perm_table, "NETWORK_permutation_summary_manuscript_ready.csv")
  
  network_summary <- tibble(
    string_version = string_version,
    string_species = species,
    string_score_threshold = score_threshold,
    n_core_targets = length(core_symbols),
    n_string_mapped_targets = nrow(mapped),
    n_prioritised_targets = length(prioritised_symbols),
    n_prioritised_targets_mapped = length(prioritised_set),
    n_network_nodes = igraph::vcount(g),
    n_network_edges = igraph::ecount(g),
    network_density = ifelse(igraph::vcount(g) > 1, igraph::edge_density(g, loops = FALSE), NA_real_),
    n_components = igraph::components(g)$no,
    largest_component_size = ifelse(length(igraph::components(g)$csize) > 0, max(igraph::components(g)$csize), NA_integer_),
    seed = seed,
    n_permutations = n_permutations,
    empirical_test = "one-sided greater/equal",
    p_adjustment_method = "Benjamini-Hochberg",
    n_metrics_tested = length(tested_metrics),
    n_metrics_raw_p_lt_0.05 = sum(perm_results$passes_raw_0.05, na.rm = TRUE),
    n_metrics_BH_p_lt_0.05 = sum(perm_results$passes_BH_0.05, na.rm = TRUE)
  )
  
  save_table(network_summary, "NETWORK_summary.csv")
  
  # ---------------------------------------------------------------------------
  # 8. Plots
  # ---------------------------------------------------------------------------
  
  cat("  [N6] Plotting network validation outputs...\n")
  
  # Color scheme:
  # - red: prioritised top 15 / BH-significant metrics
  # - grey: remaining conserved core / non-significant metrics
  col_prioritised     <- "#B2182B"
  col_other           <- "#BDBDBD"
  col_supported       <- "#B2182B"
  col_not_supported   <- "#D9D9D9"
  col_edge            <- "#C9C9C9"
  
  group_cols <- c("Other conserved core" = col_other)
  group_cols[prioritised_label] <- col_prioritised
  
  common_plot_theme <- theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "plain"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  # ---------------------------------------------------------------------------
  # 8A. Permutation support barplot
  # ---------------------------------------------------------------------------
  
  metric_order <- c(
    "Degree",
    "Weighted degree",
    "Betweenness",
    "Harmonic centrality",
    "Eigenvector centrality"
  )
  
  p_perm <- perm_results %>%
    mutate(
      metric_label = case_when(
        metric == "degree" ~ "Degree",
        metric == "weighted_degree" ~ "Weighted degree",
        metric == "betweenness" ~ "Betweenness",
        metric == "harmonic_centrality" ~ "Harmonic centrality",
        metric == "eigenvector_centrality" ~ "Eigenvector centrality",
        TRUE ~ metric
      ),
      metric_label = factor(metric_label, levels = rev(metric_order)),
      neg_log10_BH = -log10(p_adj_BH),
      BH_support = ifelse(
        p_adj_BH < 0.05,
        "BH < 0.05",
        "BH \u2265 0.05"
      ),
      BH_support = factor(
        BH_support,
        levels = c("BH < 0.05", "BH \u2265 0.05")
      )
    ) %>%
    ggplot(aes(x = neg_log10_BH, y = metric_label, fill = BH_support)) +
    geom_col(
      width = 0.65,
      color = "grey25",
      linewidth = 0.2
    ) +
    geom_vline(
      xintercept = -log10(0.05),
      linetype = "dashed",
      color = "grey30"
    ) +
    scale_fill_manual(
      values = c(
        "BH < 0.05" = col_supported,
        "BH \u2265 0.05" = col_not_supported
      ),
      breaks = c("BH < 0.05", "BH \u2265 0.05"),
      labels = c("BH < 0.05", "BH \u2265 0.05"),
      name = "Support"
    ) +
    labs(
      title = "Network enrichment of prioritised targets",
      x = expression(-log[10](BH-adjusted~p)),
      y = NULL
    ) +
    common_plot_theme
  
  save_plot(
    p_perm,
    "NETWORK_BH_adjusted_permutation_pvalues",
    width = 7.2,
    height = 4.8
  )
  
  # ---------------------------------------------------------------------------
  # 8B. STRING network layout
  # ---------------------------------------------------------------------------
  
  if (igraph::vcount(g) > 0) {
    set.seed(seed)
    
    layout_mat <- igraph::layout_with_fr(
      g,
      weights = if (igraph::ecount(g) > 0) E(g)$weight else NULL,
      niter = 2000
    )
    
    layout_df <- tibble(
      target_symbol = V(g)$name,
      x = as.numeric(scale(layout_mat[, 1])),
      y = as.numeric(scale(layout_mat[, 2])),
      network_group = V(g)$network_group,
      degree = as.numeric(degree_raw[V(g)$name])
    )
    
    if (igraph::ecount(g) > 0) {
      edge_ends <- igraph::ends(g, E(g), names = TRUE)
      
      edge_plot_df <- tibble(
        from = edge_ends[, 1],
        to = edge_ends[, 2],
        weight = E(g)$weight
      ) %>%
        left_join(
          layout_df %>%
            select(from = target_symbol, x_from = x, y_from = y),
          by = "from"
        ) %>%
        left_join(
          layout_df %>%
            select(to = target_symbol, x_to = x, y_to = y),
          by = "to"
        )
    } else {
      edge_plot_df <- tibble(
        from = character(),
        to = character(),
        x_from = numeric(),
        y_from = numeric(),
        x_to = numeric(),
        y_to = numeric(),
        weight = numeric()
      )
    }
    
    label_df <- layout_df %>%
      filter(network_group == prioritised_label)
    
    p_network <- ggplot() +
      geom_segment(
        data = edge_plot_df,
        aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
        color = col_edge,
        alpha = 0.22,
        linewidth = 0.25
      ) +
      geom_point(
        data = layout_df,
        aes(x = x, y = y, fill = network_group, size = degree),
        shape = 21,
        color = "white",
        stroke = 0.35,
        alpha = 0.95
      ) +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(x = x, y = y, label = target_symbol),
        size = 3,
        fontface = "bold",
        color = "grey15",
        box.padding = 0.35,
        point.padding = 0.2,
        segment.color = "grey50",
        segment.size = 0.25,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = seed,
        show.legend = FALSE
      ) +
      scale_fill_manual(
        values = group_cols,
        name = NULL
      ) +
      scale_size_continuous(
        range = c(2.2, 7),
        name = "Degree"
      ) +
      labs(
        title = "Conserved-core interaction network",
        subtitle = paste0("Prioritised top ", top_n, " targets are highlighted and labelled"),
        x = NULL,
        y = NULL
      ) +
      coord_cartesian(clip = "off") +
      theme_void(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", color = "grey10"),
        plot.subtitle = element_text(color = "grey25"),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(10, 20, 10, 10)
      )
    
    save_plot(
      p_network,
      "NETWORK_STRING_core_prioritised_layout",
      width = 8.8,
      height = 6.6
    )
    
    save_table(layout_df, "NETWORK_plot_layout_nodes.csv")
    save_table(edge_plot_df, "NETWORK_plot_layout_edges.csv")
  }
  # ---------------------------------------------------------------------------
  # 9. Export graph object and finish
  # ---------------------------------------------------------------------------
  
  saveRDS(g, file.path(network_dir, "NETWORK_STRING_core_graph.rds"))
  
  assign("network_graph", g, envir = .GlobalEnv)
  assign("network_node_metrics", node_metrics, envir = .GlobalEnv)
  assign("network_group_summary", group_summary, envir = .GlobalEnv)
  assign("network_permutation_summary", perm_results, envir = .GlobalEnv)
  assign("network_summary", network_summary, envir = .GlobalEnv)
  assign("network_method_parameters", method_parameters, envir = .GlobalEnv)
  
  cat("\n")
  cat("  ✓ NETWORK VALIDATION COMPLETE\n\n")
  cat("  Main outputs:\n")
  cat("   - ", normalizePath(file.path(table_dir, "NETWORK_node_metrics.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("   - ", normalizePath(file.path(table_dir, "NETWORK_group_summary.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("   - ", normalizePath(file.path(table_dir, "NETWORK_permutation_summary_with_BH.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("   - ", normalizePath(file.path(table_dir, "NETWORK_STRING_edges_for_igraph.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("   - ", normalizePath(file.path(table_dir, "NETWORK_method_parameters.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("   - ", normalizePath(plot_dir, winslash = "/", mustWork = FALSE), "\n\n", sep = "")
  
  invisible(list(
    graph = g,
    string_mapping = string_map,
    edges = edge_tbl,
    edges_for_graph = edge_for_graph,
    node_metrics = node_metrics,
    group_summary = group_summary,
    permutation_summary = perm_results,
    network_summary = network_summary,
    method_parameters = method_parameters
  ))
}
