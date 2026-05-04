################################################################################
# STEP 5B: PROBE-LINKED AND DRUG-LEVEL INDICATION EXPANSION
#
# Purpose:
#   Adds a secondary Open Targets evidence layer:
#     target -> chemicalProbes -> drug
#     drug   -> indications
#
#   It also expands indications for primary target-level drugs by querying the
#   drug-level indications page. This helps recover drug-level disease evidence
#   that may not be returned in target-level drugAndClinicalCandidates output.
#
# Probe inclusion criteria:
#   1. isHighQuality = TRUE
#   2. associated drug maximumClinicalStage >= min_probe_phase
#      UNKNOWN / NULL probe-drug stages are excluded from the probe-linked layer.
#
# Important:
#   - Primary target-level evidence is retained for audit.
#   - Drug-level/probe-linked evidence is retained separately for audit.
#   - Integrated evidence is built from the union of:
#       primary target-level rows
#       primary drug-level indication expansion
#       high-quality probe-linked drug-level indication expansion
#   - overwrite_target_master = TRUE replaces the main Step-5 target_master.csv
#     with clean integrated columns so that Step 6 and Figure 3 use the updated
#     integrated evidence.
#
# Outputs:
#   - STEP5B_method_parameters.csv
#   - OT_chemical_probe_links_all_unfiltered.csv
#   - OT_chemical_probe_links_high_quality_only.csv
#   - OT_probe_drug_phase_check.csv
#   - OT_probe_drugs_excluded_by_phase_filter.csv
#   - OT_chemical_probe_drug_links.csv
#   - OT_target_drug_links_for_deep_indication_query.csv
#   - OT_drug_level_indications_raw.csv
#   - OT_target_drug_deep_indication_rows.csv
#   - OT_primary_target_level_rows_for_integration.csv
#   - OT_integrated_drug_indication_rows.csv
#   - target_deep_indication_summary.csv
#   - target_master_integrated_deep_indications.csv
#   - QC_targets_with_phase_gain_from_deep_indications.csv
#   - QC_targets_with_breast_gain_from_deep_indications.csv
#   - QC_TTK_deep_indication_rows.csv, if TTK is present
#
# Usage:
#   source("05b_module_probe_linked_indications.R")
#   res_5b <- module_probe_linked_indications(overwrite_target_master = FALSE)
#
#   # Review QC files, then if satisfied:
#   res_5b <- module_probe_linked_indications(overwrite_target_master = TRUE)
################################################################################

module_probe_linked_indications <- function(
    overwrite_target_master = FALSE,
    min_probe_phase = 1,
    sleep_sec = 0.35,
    force_refresh = FALSE
) {
  
  cat("\n")
  cat("======================================================================\n")
  cat("STEP 5B: PROBE-LINKED AND DRUG-LEVEL INDICATION EXPANSION\n")
  cat("======================================================================\n\n")
  cat(
    "  Probe filter: isHighQuality = TRUE AND maximumClinicalStage >=",
    min_probe_phase, "\n\n"
  )
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)
    library(tibble)
    library(httr)
    library(jsonlite)
    library(purrr)
  })
  
  # ---------------------------------------------------------------------------
  # 1. Paths
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
  
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  ot_probe_cache_dir <- file.path(cache_dir, "OT_probe_links")
  ot_drug_cache_dir  <- file.path(cache_dir, "OT_drug_indications")
  
  dir.create(ot_probe_cache_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(ot_drug_cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  save_result_local <- function(x, filename) {
    out <- file.path(results_dir, filename)
    readr::write_csv(x, out)
    cat(
      "  \u2713 Saved: ",
      normalizePath(out, winslash = "/", mustWork = FALSE),
      "\n",
      sep = ""
    )
    invisible(out)
  }
  
  # ---------------------------------------------------------------------------
  # 2. Helper functions
  # ---------------------------------------------------------------------------
  
  safe_chr <- function(x, default = NA_character_) {
    if (is.null(x) || length(x) == 0) return(default)
    if (is.list(x)) x <- x[[1]]
    if (is.null(x) || length(x) == 0 || is.na(x)) return(default)
    as.character(x)
  }
  
  safe_num <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0) return(default)
    if (is.list(x)) x <- x[[1]]
    out <- suppressWarnings(as.numeric(x))
    ifelse(is.na(out), default, out)
  }
  
  safe_bool <- function(x, default = FALSE) {
    if (is.null(x) || length(x) == 0) return(default)
    if (is.list(x)) x <- x[[1]]
    if (is.logical(x)) return(isTRUE(x))
    x_chr <- tolower(trimws(as.character(x)))
    x_chr %in% c("true", "t", "1", "yes", "y")
  }
  
  col_or <- function(df, candidates, default = NA_character_) {
    for (c in candidates) {
      if (c %in% names(df)) return(df[[c]])
    }
    rep(default, nrow(df))
  }
  
  first_non_missing <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0) NA_character_ else x[1]
  }
  
  ensure_columns <- function(df, defaults) {
    for (nm in names(defaults)) {
      if (!nm %in% names(df)) df[[nm]] <- defaults[[nm]]
    }
    df
  }
  
  stage_to_numeric <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA_real_)
    
    x_chr <- toupper(trimws(as.character(x)))
    x_chr[x_chr %in% c("", "NA", "NAN", "NULL", "NONE", "UNKNOWN")] <- NA_character_
    
    out <- rep(NA_real_, length(x_chr))
    
    # Approved/marketed/regulatory terms.
    out[grepl("APPROVAL|APPROVED|REGULATORY|MARKETED", x_chr)] <- 4
    
    # Phase III and mixed II/III.
    out[is.na(out) & grepl(
      "PHASE[_ -]*(3|III)|PHASE[_ -]*(2|II)[_/ -]*(3|III)|PHASE[_ -]*II[_/ -]*III",
      x_chr
    )] <- 3
    
    # Phase II and mixed I/II, including phase 1b/2.
    out[is.na(out) & grepl(
      "PHASE[_ -]*(2|II)|PHASE[_ -]*(1|I|1B|IB)[_/ -]*(2|II)|PHASE[_ -]*I[_/ -]*II",
      x_chr
    )] <- 2
    
    # Phase I.
    out[is.na(out) & grepl(
      "PHASE[_ -]*(1|I)|EARLY[_ -]*PHASE[_ -]*(1|I)",
      x_chr
    )] <- 1
    
    # Preclinical / phase 0.
    out[is.na(out) & grepl(
      "PHASE[_ -]*0|PRECLINICAL|PRE[_ -]*CLINICAL|IND",
      x_chr
    )] <- 0
    
    numeric_direct <- suppressWarnings(as.numeric(x_chr))
    idx <- is.na(out) & !is.na(numeric_direct)
    out[idx] <- numeric_direct[idx]
    
    out
  }
  
  disease_key <- function(id, name) {
    id <- as.character(id)
    name <- as.character(name)
    
    out <- ifelse(
      !is.na(id) & id != "",
      paste0("ID:", id),
      paste0("NAME:", str_to_lower(str_squish(name)))
    )
    
    out[is.na(out) | out %in% c("NAME:", "ID:")] <- NA_character_
    out
  }
  
  # Keep these terms synchronized with the rules in 05_module_drug_targets.R.
  oncology_regex <- paste(
    c(
      "cancer", "carcinoma", "sarcoma", "tumou?r", "neoplasm",
      "leukemia", "leukaemia", "lymphoma", "melanoma", "glioma",
      "blastoma", "myeloma", "malignan"
    ),
    collapse = "|"
  )
  
  breast_regex <- "breast|mammary"
  
  classify_indications <- function(df) {
    df <- ensure_columns(df, list(
      disease_name = NA_character_,
      disease_from_source = NA_character_,
      indication_source = NA_character_,
      disease_id = NA_character_
    ))
    
    df %>%
      mutate(
        disease_text = str_to_lower(str_squish(paste(
          coalesce(as.character(disease_name), ""),
          coalesce(as.character(disease_from_source), ""),
          coalesce(as.character(indication_source), ""),
          sep = " "
        ))),
        is_oncology_drug_indication = str_detect(disease_text, oncology_regex),
        is_breast_drug_indication = str_detect(disease_text, breast_regex),
        disease_key = disease_key(
          disease_id,
          coalesce(as.character(disease_name), as.character(disease_from_source))
        )
      )
  }
  
  OT_URL <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  ot_post <- function(query, variables = list()) {
    res <- httr::POST(
      OT_URL,
      body = list(query = query, variables = variables),
      encode = "json",
      httr::add_headers(
        `Content-Type` = "application/json",
        `Accept` = "application/json",
        `User-Agent` = "R httr OpenTargets probe-linked indication client"
      ),
      httr::timeout(90)
    )
    
    txt <- httr::content(res, "text", encoding = "UTF-8")
    
    if (httr::status_code(res) != 200) {
      return(list(ok = FALSE, status = httr::status_code(res), text = txt))
    }
    
    j <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
    
    if (!is.null(j$errors)) {
      return(list(
        ok = FALSE,
        status = 200,
        graphql_errors = j$errors,
        text = txt,
        data = j
      ))
    }
    
    list(ok = TRUE, data = j)
  }
  
  method_parameters <- tibble(
    parameter = c(
      "primary_layer",
      "secondary_layer",
      "probe_filter_isHighQuality",
      "probe_filter_minimum_drug_phase",
      "phase_encoding",
      "oncology_regex",
      "breast_regex",
      "overwrite_target_master",
      "force_refresh"
    ),
    value = c(
      "Open Targets target-level drug/clinical-candidate output from Step 5",
      "Open Targets chemicalProbes -> drug-level indications expansion",
      "TRUE",
      as.character(min_probe_phase),
      "Phase I/II encoded as 2; approved/marketed encoded as 4",
      oncology_regex,
      breast_regex,
      as.character(overwrite_target_master),
      as.character(force_refresh)
    )
  )
  save_result_local(method_parameters, "STEP5B_method_parameters.csv")
  
  # ---------------------------------------------------------------------------
  # 3. Load target master and primary drug-evidence table
  # ---------------------------------------------------------------------------
  
  cat("  [5B-1] Loading target and primary drug evidence tables...\n")
  
  target_master_file <- file.path(results_dir, "target_master.csv")
  
  if (exists("target_master", envir = .GlobalEnv)) {
    target_master <- get("target_master", envir = .GlobalEnv)
  } else if (file.exists(target_master_file)) {
    target_master <- readr::read_csv(target_master_file, show_col_types = FALSE)
  } else {
    stop("target_master.csv not found. Run Step 5 first.", call. = FALSE)
  }
  
  if (!all(c("target_symbol", "target_ensembl") %in% names(target_master))) {
    stop("target_master must contain target_symbol and target_ensembl.", call. = FALSE)
  }
  
  target_master <- ensure_columns(target_master, list(
    n_drugs = 0,
    n_indications = 0,
    n_oncology_drug_indications = 0,
    n_breast_drug_indications = 0,
    max_drug_phase = 0,
    any_approved = FALSE,
    any_onco_drug = FALSE,
    any_breast_drug = FALSE,
    n_rows_with_phase = 0,
    n_rows_without_phase = 0
  ))
  
  primary_files <- c(
    file.path(results_dir, "drug_target_table_with_ATC_and_fallback.csv"),
    file.path(results_dir, "drug_target_table_filtered.csv"),
    file.path(results_dir, "drug_target_table_raw_OpenTargets.csv")
  )
  primary_file <- primary_files[file.exists(primary_files)][1]
  
  if (is.na(primary_file) || length(primary_file) == 0) {
    warning(
      "No primary drug evidence table found. Step 5B will use probe-linked drugs only.",
      call. = FALSE
    )
    primary_drug <- tibble()
  } else {
    primary_drug <- readr::read_csv(primary_file, show_col_types = FALSE)
    cat("      Primary drug table: ", basename(primary_file), "\n", sep = "")
  }
  
  target_table <- target_master %>%
    select(target_symbol, target_ensembl) %>%
    filter(
      !is.na(target_symbol), target_symbol != "",
      !is.na(target_ensembl), target_ensembl != ""
    ) %>%
    distinct()
  
  cat("      Targets for probe query:", nrow(target_table), "\n")
  
  # ---------------------------------------------------------------------------
  # 4. Query Open Targets chemical probes
  # ---------------------------------------------------------------------------
  
  cat("  [5B-2] Querying Open Targets chemical probes...\n")
  
  query_target_probes <- function(ensembl_id) {
    query <- "
      query getTargetChemicalProbes($ensemblId: String!) {
        target(ensemblId: $ensemblId) {
          id
          approvedSymbol
          approvedName
          chemicalProbes {
            id
            drugId
            mechanismOfAction
            isHighQuality
            probesDrugsScore
            scoreInCells
            scoreInOrganisms
          }
        }
      }
    "
    
    ot_post(query, variables = list(ensemblId = ensembl_id))
  }
  
  get_probes_cached <- function(ensembl_id) {
    f <- file.path(ot_probe_cache_dir, paste0(ensembl_id, "_probes.json"))
    
    if (file.exists(f) && !isTRUE(force_refresh)) {
      j <- jsonlite::fromJSON(readr::read_file(f), simplifyVector = FALSE)
      return(list(ok = TRUE, data = j, from_cache = TRUE))
    }
    
    out <- query_target_probes(ensembl_id)
    
    if (isTRUE(out$ok)) {
      writeLines(
        jsonlite::toJSON(out$data, auto_unbox = TRUE, pretty = TRUE),
        con = f
      )
    }
    
    Sys.sleep(sleep_sec)
    out
  }
  
  extract_probe_links <- function(x, ensembl_id) {
    if (!isTRUE(x$ok)) return(NULL)
    
    t <- x$data$data$target
    if (is.null(t)) return(NULL)
    
    probes <- t$chemicalProbes
    if (is.null(probes) || length(probes) == 0) return(NULL)
    
    out <- vector("list", length(probes))
    
    for (i in seq_along(probes)) {
      p <- probes[[i]]
      
      out[[i]] <- tibble(
        target_ensembl = ensembl_id,
        target_symbol = safe_chr(t$approvedSymbol),
        target_name = safe_chr(t$approvedName),
        probe_id = safe_chr(p$id),
        probe_drug_chembl = safe_chr(p$drugId),
        probe_mechanism_of_action = safe_chr(p$mechanismOfAction),
        probe_is_high_quality = safe_bool(p$isHighQuality, default = FALSE),
        probe_probesdrugs_score = safe_num(p$probesDrugsScore),
        probe_score_in_cells = safe_num(p$scoreInCells),
        probe_score_in_organisms = safe_num(p$scoreInOrganisms)
      )
    }
    
    bind_rows(out)
  }
  
  probe_raw <- lapply(target_table$target_ensembl, get_probes_cached)
  names(probe_raw) <- target_table$target_ensembl
  
  probe_links_all <- bind_rows(Map(extract_probe_links, probe_raw, names(probe_raw)))
  
  if (is.null(probe_links_all) || nrow(probe_links_all) == 0) {
    probe_links_all <- tibble(
      target_ensembl = character(),
      target_symbol = character(),
      target_name = character(),
      probe_id = character(),
      probe_drug_chembl = character(),
      probe_mechanism_of_action = character(),
      probe_is_high_quality = logical(),
      probe_probesdrugs_score = numeric(),
      probe_score_in_cells = numeric(),
      probe_score_in_organisms = numeric()
    )
  }
  
  save_result_local(probe_links_all, "OT_chemical_probe_links_all_unfiltered.csv")
  
  probe_links_hq <- probe_links_all %>%
    filter(
      probe_is_high_quality %in% TRUE,
      !is.na(probe_drug_chembl), probe_drug_chembl != ""
    ) %>%
    distinct()
  
  save_result_local(probe_links_hq, "OT_chemical_probe_links_high_quality_only.csv")
  
  cat("      Total probe records found:", nrow(probe_links_all), "\n")
  cat("      After high-quality filter:", nrow(probe_links_hq), "\n")
  
  # ---------------------------------------------------------------------------
  # 5. Query clinical stage for high-quality probe-linked drugs
  # ---------------------------------------------------------------------------
  
  cat("  [5B-3] Querying clinical stage for high-quality probe-linked drugs...\n")
  
  unique_probe_drugs <- sort(unique(probe_links_hq$probe_drug_chembl))
  cat("      Unique high-quality probe drugs to check:", length(unique_probe_drugs), "\n")
  
  query_drug_phase <- function(chembl_id) {
    query_current <- "
      query getDrugPhase($chemblId: String!) {
        drug(chemblId: $chemblId) {
          id
          name
          maximumClinicalStage
        }
      }
    "
    
    out <- ot_post(query_current, variables = list(chemblId = chembl_id))
    
    if (!isTRUE(out$ok)) {
      query_legacy <- "
        query getDrugPhaseLegacy($chemblId: String!) {
          drug(chemblId: $chemblId) {
            id
            name
            maximumClinicalTrialPhase
          }
        }
      "
      out <- ot_post(query_legacy, variables = list(chemblId = chembl_id))
    }
    
    out
  }
  
  get_drug_phase_cached <- function(chembl_id) {
    # Reuse drug-indication cache when available, because it contains drug-level
    # maximum clinical stage.
    f_ind <- file.path(ot_drug_cache_dir, paste0(chembl_id, "_drug_indications.json"))
    f_phase <- file.path(ot_drug_cache_dir, paste0(chembl_id, "_drug_phase.json"))
    
    if (file.exists(f_ind) && !isTRUE(force_refresh)) {
      j <- jsonlite::fromJSON(readr::read_file(f_ind), simplifyVector = FALSE)
      d <- j$data$drug
      
      phase_raw <- safe_chr(d$maximumClinicalStage)
      if (is.na(phase_raw)) phase_raw <- safe_chr(d$maximumClinicalTrialPhase)
      
      return(tibble(
        drug_chembl = chembl_id,
        drug_name = safe_chr(d$name),
        max_phase_raw = phase_raw,
        max_phase_num = stage_to_numeric(phase_raw)[1],
        phase_query_source = "cached_drug_indications"
      ))
    }
    
    if (file.exists(f_phase) && !isTRUE(force_refresh)) {
      j <- jsonlite::fromJSON(readr::read_file(f_phase), simplifyVector = FALSE)
      d <- j$data$drug
      
      phase_raw <- safe_chr(d$maximumClinicalStage)
      if (is.na(phase_raw)) phase_raw <- safe_chr(d$maximumClinicalTrialPhase)
      
      return(tibble(
        drug_chembl = chembl_id,
        drug_name = safe_chr(d$name),
        max_phase_raw = phase_raw,
        max_phase_num = stage_to_numeric(phase_raw)[1],
        phase_query_source = "cached_drug_phase"
      ))
    }
    
    out <- query_drug_phase(chembl_id)
    Sys.sleep(sleep_sec)
    
    if (!isTRUE(out$ok)) {
      return(tibble(
        drug_chembl = chembl_id,
        drug_name = NA_character_,
        max_phase_raw = NA_character_,
        max_phase_num = NA_real_,
        phase_query_source = "query_failed"
      ))
    }
    
    writeLines(
      jsonlite::toJSON(out$data, auto_unbox = TRUE, pretty = TRUE),
      con = f_phase
    )
    
    d <- out$data$data$drug
    phase_raw <- safe_chr(d$maximumClinicalStage)
    if (is.na(phase_raw)) phase_raw <- safe_chr(d$maximumClinicalTrialPhase)
    
    tibble(
      drug_chembl = chembl_id,
      drug_name = safe_chr(d$name),
      max_phase_raw = phase_raw,
      max_phase_num = stage_to_numeric(phase_raw)[1],
      phase_query_source = "live_drug_phase_query"
    )
  }
  
  probe_drug_phases <- bind_rows(lapply(unique_probe_drugs, get_drug_phase_cached))
  save_result_local(probe_drug_phases, "OT_probe_drug_phase_check.csv")
  
  probe_drugs_clinical <- probe_drug_phases %>%
    filter(!is.na(max_phase_num), max_phase_num >= min_probe_phase)
  
  probe_drugs_excluded <- probe_drug_phases %>%
    filter(is.na(max_phase_num) | max_phase_num < min_probe_phase)
  
  save_result_local(probe_drugs_excluded, "OT_probe_drugs_excluded_by_phase_filter.csv")
  
  cat("      Probe drugs retained after phase filter:", nrow(probe_drugs_clinical), "\n")
  if (nrow(probe_drugs_clinical) > 0) {
    cat("      Retained probe-linked drugs:", paste(probe_drugs_clinical$drug_name, collapse = ", "), "\n")
  }
  
  probe_links <- probe_links_hq %>%
    filter(probe_drug_chembl %in% probe_drugs_clinical$drug_chembl) %>%
    left_join(
      probe_drug_phases %>%
        select(drug_chembl, drug_name, max_phase_raw, max_phase_num),
      by = c("probe_drug_chembl" = "drug_chembl")
    ) %>%
    rename(drug_chembl = probe_drug_chembl) %>%
    distinct()
  
  save_result_local(probe_links, "OT_chemical_probe_drug_links.csv")
  
  cat("      Probe-linked target-drug records retained:", nrow(probe_links), "\n")
  cat("      Unique retained probe-linked drugs:", n_distinct(probe_links$drug_chembl), "\n")
  
  # ---------------------------------------------------------------------------
  # 6. Build target-drug link table for drug-level indication expansion
  # ---------------------------------------------------------------------------
  
  cat("  [5B-4] Building target-drug link table for deep indication expansion...\n")
  
  primary_links <- if (nrow(primary_drug) > 0 &&
                       all(c("target_symbol", "drug_chembl") %in% names(primary_drug))) {
    tibble(
      target_symbol = as.character(col_or(primary_drug, "target_symbol")),
      target_ensembl = as.character(col_or(primary_drug, "target_ensembl")),
      drug_chembl = as.character(col_or(primary_drug, "drug_chembl")),
      drug_name_primary = as.character(col_or(primary_drug, "drug_name")),
      link_source = "primary_target_level"
    ) %>%
      filter(!is.na(drug_chembl), drug_chembl != "") %>%
      distinct()
  } else {
    tibble(
      target_symbol = character(),
      target_ensembl = character(),
      drug_chembl = character(),
      drug_name_primary = character(),
      link_source = character()
    )
  }
  
  probe_target_links <- probe_links %>%
    transmute(
      target_symbol = as.character(target_symbol),
      target_ensembl = as.character(target_ensembl),
      drug_chembl = as.character(drug_chembl),
      drug_name_primary = as.character(drug_name),
      link_source = "chemical_probe_linked"
    ) %>%
    distinct()
  
  target_drug_links <- bind_rows(primary_links, probe_target_links) %>%
    filter(
      !is.na(target_symbol), target_symbol != "",
      !is.na(drug_chembl), drug_chembl != ""
    ) %>%
    group_by(target_symbol, target_ensembl, drug_chembl) %>%
    summarise(
      link_source = paste(sort(unique(link_source)), collapse = "|"),
      drug_name_primary = first_non_missing(drug_name_primary),
      .groups = "drop"
    )
  
  save_result_local(target_drug_links, "OT_target_drug_links_for_deep_indication_query.csv")
  
  cat("      Target-drug links:", nrow(target_drug_links), "\n")
  cat("      Unique drugs to query:", n_distinct(target_drug_links$drug_chembl), "\n")
  
  # ---------------------------------------------------------------------------
  # 7. Query drug-level indications
  # ---------------------------------------------------------------------------
  
  cat("  [5B-5] Querying Open Targets drug-level indications...\n")
  
  query_drug_indications <- function(chembl_id) {
    query_current <- "
      query getDrugIndications($chemblId: String!) {
        drug(chemblId: $chemblId) {
          id
          name
          drugType
          maximumClinicalStage
          indications {
            count
            rows {
              maxClinicalStage
              disease {
                id
                name
              }
            }
          }
        }
      }
    "
    
    out <- ot_post(query_current, variables = list(chemblId = chembl_id))
    
    if (!isTRUE(out$ok)) {
      query_legacy <- "
        query getDrugIndicationsLegacy($chemblId: String!) {
          drug(chemblId: $chemblId) {
            id
            name
            drugType
            maximumClinicalTrialPhase
            indications {
              count
              rows {
                maxPhaseForIndication
                disease {
                  id
                  name
                }
              }
            }
          }
        }
      "
      out <- ot_post(query_legacy, variables = list(chemblId = chembl_id))
    }
    
    out
  }
  
  get_drug_indications_cached <- function(chembl_id) {
    f <- file.path(ot_drug_cache_dir, paste0(chembl_id, "_drug_indications.json"))
    
    if (file.exists(f) && !isTRUE(force_refresh)) {
      j <- jsonlite::fromJSON(readr::read_file(f), simplifyVector = FALSE)
      return(list(ok = TRUE, data = j, from_cache = TRUE))
    }
    
    out <- query_drug_indications(chembl_id)
    
    if (isTRUE(out$ok)) {
      writeLines(
        jsonlite::toJSON(out$data, auto_unbox = TRUE, pretty = TRUE),
        con = f
      )
    }
    
    Sys.sleep(sleep_sec)
    out
  }
  
  extract_drug_indications <- function(x, chembl_id) {
    empty_row <- tibble(
      drug_chembl = chembl_id,
      drug_query_ok = FALSE,
      drug_query_error = NA_character_,
      drug_name = NA_character_,
      drug_type = NA_character_,
      drug_max_phase = NA_real_,
      drug_max_phase_raw = NA_character_,
      phase_evidence = NA_real_,
      phase_evidence_raw = NA_character_,
      disease_id = NA_character_,
      disease_name = NA_character_,
      disease_from_source = NA_character_,
      indication_source = "drug_level_indications"
    )
    
    if (!isTRUE(x$ok)) {
      empty_row$drug_query_error <- if (!is.null(x$text)) substr(x$text, 1, 500) else "GraphQL error"
      return(empty_row)
    }
    
    d <- x$data$data$drug
    if (is.null(d)) {
      empty_row$drug_query_ok <- TRUE
      return(empty_row)
    }
    
    drug_stage_raw <- safe_chr(d$maximumClinicalStage)
    if (is.na(drug_stage_raw)) drug_stage_raw <- safe_chr(d$maximumClinicalTrialPhase)
    drug_stage_num <- stage_to_numeric(drug_stage_raw)[1]
    
    rows <- d$indications$rows
    
    if (is.null(rows) || length(rows) == 0) {
      return(tibble(
        drug_chembl = chembl_id,
        drug_query_ok = TRUE,
        drug_query_error = NA_character_,
        drug_name = safe_chr(d$name),
        drug_type = safe_chr(d$drugType),
        drug_max_phase = drug_stage_num,
        drug_max_phase_raw = drug_stage_raw,
        phase_evidence = drug_stage_num,
        phase_evidence_raw = drug_stage_raw,
        disease_id = NA_character_,
        disease_name = NA_character_,
        disease_from_source = NA_character_,
        indication_source = "drug_level_indications"
      ))
    }
    
    out <- vector("list", length(rows))
    
    for (i in seq_along(rows)) {
      r <- rows[[i]]
      disease <- r$disease
      
      row_stage_raw <- safe_chr(r$maxClinicalStage)
      if (is.na(row_stage_raw)) row_stage_raw <- safe_chr(r$maxPhaseForIndication)
      row_stage_num <- stage_to_numeric(row_stage_raw)[1]
      
      phase_for_record <- suppressWarnings(max(c(row_stage_num, drug_stage_num), na.rm = TRUE))
      if (is.infinite(phase_for_record)) phase_for_record <- NA_real_
      
      out[[i]] <- tibble(
        drug_chembl = chembl_id,
        drug_query_ok = TRUE,
        drug_query_error = NA_character_,
        drug_name = safe_chr(d$name),
        drug_type = safe_chr(d$drugType),
        drug_max_phase = drug_stage_num,
        drug_max_phase_raw = drug_stage_raw,
        phase_evidence = phase_for_record,
        phase_evidence_raw = row_stage_raw,
        disease_id = safe_chr(disease$id),
        disease_name = safe_chr(disease$name),
        disease_from_source = safe_chr(disease$name),
        indication_source = "drug_level_indications"
      )
    }
    
    bind_rows(out)
  }
  
  unique_drugs <- sort(unique(target_drug_links$drug_chembl))
  drug_raw <- lapply(unique_drugs, get_drug_indications_cached)
  names(drug_raw) <- unique_drugs
  
  drug_indications <- bind_rows(Map(extract_drug_indications, drug_raw, names(drug_raw))) %>%
    classify_indications() %>%
    distinct()
  
  save_result_local(drug_indications, "OT_drug_level_indications_raw.csv")
  cat("      Drug-level indication rows:", nrow(drug_indications), "\n")
  
  # ---------------------------------------------------------------------------
  # 8. Join target-drug links to drug-level indications
  # ---------------------------------------------------------------------------
  
  cat("  [5B-6] Expanding target-drug links with drug-level indications...\n")
  
  deep_rows <- target_drug_links %>%
    left_join(drug_indications, by = "drug_chembl", relationship = "many-to-many") %>%
    mutate(
      evidence_layer = case_when(
        str_detect(link_source, "primary_target_level") &
          str_detect(link_source, "chemical_probe_linked") ~
          "primary_and_probe_linked_drug_indications",
        str_detect(link_source, "chemical_probe_linked") ~
          "probe_linked_drug_indications",
        str_detect(link_source, "primary_target_level") ~
          "primary_drug_level_indications",
        TRUE ~ "drug_level_indications"
      ),
      drug_name = coalesce(as.character(drug_name), as.character(drug_name_primary)),
      drug_is_approved = !is.na(phase_evidence) & phase_evidence >= 4
    ) %>%
    classify_indications() %>%
    distinct()
  
  save_result_local(deep_rows, "OT_target_drug_deep_indication_rows.csv")
  
  # ---------------------------------------------------------------------------
  # 9. Primary target-level rows for true integrated union
  # ---------------------------------------------------------------------------
  
  cat("  [5B-7] Preparing primary target-level rows for integrated summary...\n")
  
  if (nrow(primary_drug) > 0 && all(c("target_symbol", "drug_chembl") %in% names(primary_drug))) {
    primary_phase_a <- suppressWarnings(as.numeric(col_or(
      primary_drug,
      c("phase_for_filter", "phase_evidence", "max_drug_phase"),
      NA_real_
    )))
    primary_phase_b <- suppressWarnings(as.numeric(col_or(
      primary_drug,
      c("drug_max_phase", "max_drug_phase"),
      NA_real_
    )))
    
    primary_phase <- pmax(primary_phase_a, primary_phase_b, na.rm = TRUE)
    primary_phase[is.infinite(primary_phase)] <- NA_real_
    
    primary_target_rows <- tibble(
      target_symbol = as.character(col_or(primary_drug, "target_symbol")),
      target_ensembl = as.character(col_or(primary_drug, "target_ensembl")),
      drug_chembl = as.character(col_or(primary_drug, "drug_chembl")),
      drug_name = as.character(col_or(primary_drug, "drug_name")),
      drug_type = as.character(col_or(primary_drug, "drug_type")),
      drug_max_phase = suppressWarnings(as.numeric(col_or(primary_drug, "drug_max_phase", NA_real_))),
      drug_max_phase_raw = as.character(col_or(primary_drug, "drug_max_phase_raw")),
      phase_evidence = primary_phase,
      phase_evidence_raw = as.character(col_or(primary_drug, c("phase_evidence_raw", "drug_max_phase_raw"))),
      disease_id = as.character(col_or(primary_drug, "disease_id")),
      disease_name = as.character(col_or(primary_drug, "disease_name")),
      disease_from_source = as.character(col_or(primary_drug, c("disease_from_source", "disease_name"))),
      drug_is_approved = col_or(primary_drug, "drug_is_approved", FALSE) %in% TRUE,
      link_source = "primary_target_level",
      evidence_layer = "primary_target_level_rows",
      indication_source = "target_level_drugAndClinicalCandidates"
    ) %>%
      filter(!is.na(target_symbol), target_symbol != "") %>%
      classify_indications() %>%
      distinct()
  } else {
    primary_target_rows <- tibble(
      target_symbol = character(),
      target_ensembl = character(),
      drug_chembl = character(),
      drug_name = character(),
      drug_type = character(),
      drug_max_phase = numeric(),
      drug_max_phase_raw = character(),
      phase_evidence = numeric(),
      phase_evidence_raw = character(),
      disease_id = character(),
      disease_name = character(),
      disease_from_source = character(),
      drug_is_approved = logical(),
      link_source = character(),
      evidence_layer = character(),
      indication_source = character(),
      disease_text = character(),
      is_oncology_drug_indication = logical(),
      is_breast_drug_indication = logical(),
      disease_key = character()
    )
  }
  
  save_result_local(primary_target_rows, "OT_primary_target_level_rows_for_integration.csv")
  
  integrated_rows <- bind_rows(primary_target_rows, deep_rows) %>%
    distinct()
  
  save_result_local(integrated_rows, "OT_integrated_drug_indication_rows.csv")
  
  # ---------------------------------------------------------------------------
  # 10. Summarise at target level
  # ---------------------------------------------------------------------------
  
  cat("  [5B-8] Summarising evidence at target level...\n")
  
  summarise_evidence <- function(df, suffix) {
    if (nrow(df) == 0) {
      return(target_table %>% mutate(
        !!paste0("n_drug_records_", suffix) := 0L,
        !!paste0("n_indication_records_", suffix) := 0L,
        !!paste0("n_oncology_indications_", suffix) := 0L,
        !!paste0("n_breast_indications_", suffix) := 0L,
        !!paste0("max_phase_", suffix) := 0,
        !!paste0("any_approved_", suffix) := FALSE,
        !!paste0("any_oncology_", suffix) := FALSE,
        !!paste0("any_breast_", suffix) := FALSE,
        !!paste0("n_rows_with_phase_", suffix) := 0L,
        !!paste0("n_rows_without_phase_", suffix) := 0L
      ))
    }
    
    df %>%
      group_by(target_symbol, target_ensembl) %>%
      summarise(
        !!paste0("n_drug_records_", suffix) :=
          n_distinct(drug_chembl[!is.na(drug_chembl) & drug_chembl != ""]),
        
        !!paste0("n_indication_records_", suffix) :=
          n_distinct(disease_key[!is.na(disease_key) & disease_key != ""]),
        
        !!paste0("n_oncology_indications_", suffix) :=
          n_distinct(disease_key[
            is_oncology_drug_indication %in% TRUE &
              !is.na(disease_key) & disease_key != ""
          ]),
        
        !!paste0("n_breast_indications_", suffix) :=
          n_distinct(disease_key[
            is_breast_drug_indication %in% TRUE &
              !is.na(disease_key) & disease_key != ""
          ]),
        
        !!paste0("max_phase_", suffix) :=
          ifelse(all(is.na(phase_evidence)), 0, max(phase_evidence, na.rm = TRUE)),
        
        !!paste0("any_approved_", suffix) :=
          any(drug_is_approved %in% TRUE | (!is.na(phase_evidence) & phase_evidence >= 4)),
        
        !!paste0("any_oncology_", suffix) :=
          any(is_oncology_drug_indication %in% TRUE),
        
        !!paste0("any_breast_", suffix) :=
          any(is_breast_drug_indication %in% TRUE),
        
        !!paste0("n_rows_with_phase_", suffix) :=
          sum(!is.na(phase_evidence)),
        
        !!paste0("n_rows_without_phase_", suffix) :=
          sum(is.na(phase_evidence)),
        
        .groups = "drop"
      )
  }
  
  summary_primary_targetlevel <- summarise_evidence(primary_target_rows, "primary_targetlevel")
  summary_primary_druglevel <- summarise_evidence(
    deep_rows %>% filter(str_detect(link_source, "primary_target_level")),
    "primary_druglevel"
  )
  summary_probe_linked <- summarise_evidence(
    deep_rows %>% filter(str_detect(link_source, "chemical_probe_linked")),
    "probe_linked"
  )
  summary_integrated <- summarise_evidence(integrated_rows, "integrated")
  
  deep_summary <- target_table %>%
    left_join(summary_primary_targetlevel, by = c("target_symbol", "target_ensembl")) %>%
    left_join(summary_primary_druglevel, by = c("target_symbol", "target_ensembl")) %>%
    left_join(summary_probe_linked, by = c("target_symbol", "target_ensembl")) %>%
    left_join(summary_integrated, by = c("target_symbol", "target_ensembl")) %>%
    mutate(
      across(starts_with("n_"), ~ coalesce(as.integer(.x), 0L)),
      across(starts_with("max_phase_"), ~ coalesce(as.numeric(.x), 0)),
      across(starts_with("any_"), ~ coalesce(.x %in% TRUE, FALSE))
    )
  
  save_result_local(deep_summary, "target_deep_indication_summary.csv")
  
  # ---------------------------------------------------------------------------
  # 11. Integrated target_master table
  # ---------------------------------------------------------------------------
  
  cat("  [5B-9] Creating integrated target_master table...\n")
  
  target_master_integrated <- target_master %>%
    mutate(
      n_drugs_primary_original = coalesce(suppressWarnings(as.numeric(n_drugs)), 0),
      n_indications_primary_original = coalesce(suppressWarnings(as.numeric(n_indications)), 0),
      n_oncology_drug_indications_primary_original =
        coalesce(suppressWarnings(as.numeric(n_oncology_drug_indications)), 0),
      n_breast_drug_indications_primary_original =
        coalesce(suppressWarnings(as.numeric(n_breast_drug_indications)), 0),
      max_drug_phase_primary_original = coalesce(suppressWarnings(as.numeric(max_drug_phase)), 0),
      any_approved_primary_original = any_approved %in% TRUE,
      any_onco_drug_primary_original = any_onco_drug %in% TRUE,
      any_breast_drug_primary_original = any_breast_drug %in% TRUE,
      n_rows_with_phase_primary_original = coalesce(suppressWarnings(as.numeric(n_rows_with_phase)), 0),
      n_rows_without_phase_primary_original = coalesce(suppressWarnings(as.numeric(n_rows_without_phase)), 0)
    ) %>%
    left_join(deep_summary, by = c("target_symbol", "target_ensembl"))
  
  get_col <- function(col, default) {
    if (col %in% names(target_master_integrated)) {
      target_master_integrated[[col]]
    } else {
      rep(default, nrow(target_master_integrated))
    }
  }
  
  target_master_integrated <- target_master_integrated %>%
    mutate(
      n_drugs_integrated = pmax(
        n_drugs_primary_original,
        coalesce(as.numeric(get_col("n_drug_records_integrated", 0)), 0),
        na.rm = TRUE
      ),
      
      n_indications_integrated = pmax(
        n_indications_primary_original,
        coalesce(as.numeric(get_col("n_indication_records_integrated", 0)), 0),
        na.rm = TRUE
      ),
      
      n_oncology_drug_indications_integrated = pmax(
        n_oncology_drug_indications_primary_original,
        coalesce(as.numeric(get_col("n_oncology_indications_integrated", 0)), 0),
        na.rm = TRUE
      ),
      
      n_breast_drug_indications_integrated = pmax(
        n_breast_drug_indications_primary_original,
        coalesce(as.numeric(get_col("n_breast_indications_integrated", 0)), 0),
        na.rm = TRUE
      ),
      
      max_drug_phase_integrated = pmax(
        max_drug_phase_primary_original,
        coalesce(as.numeric(get_col("max_phase_integrated", 0)), 0),
        na.rm = TRUE
      ),
      
      any_approved_integrated =
        any_approved_primary_original | (get_col("any_approved_integrated", FALSE) %in% TRUE),
      
      any_onco_drug_integrated =
        any_onco_drug_primary_original | (get_col("any_oncology_integrated", FALSE) %in% TRUE),
      
      any_breast_drug_integrated =
        any_breast_drug_primary_original | (get_col("any_breast_integrated", FALSE) %in% TRUE),
      
      n_rows_with_phase_integrated = pmax(
        n_rows_with_phase_primary_original,
        coalesce(as.numeric(get_col("n_rows_with_phase_integrated", 0)), 0),
        na.rm = TRUE
      ),
      
      n_rows_without_phase_integrated = pmax(
        n_rows_without_phase_primary_original,
        coalesce(as.numeric(get_col("n_rows_without_phase_integrated", 0)), 0),
        na.rm = TRUE
      ),
      
      has_probe_linked_drug_evidence =
        coalesce(as.numeric(get_col("n_drug_records_probe_linked", 0)), 0) > 0,
      
      has_probe_linked_oncology_evidence =
        get_col("any_oncology_probe_linked", FALSE) %in% TRUE,
      
      has_probe_linked_breast_evidence =
        get_col("any_breast_probe_linked", FALSE) %in% TRUE,
      
      has_primary_vs_integrated_phase_change =
        max_drug_phase_integrated > max_drug_phase_primary_original,
      
      has_primary_vs_integrated_oncology_change =
        any_onco_drug_integrated & !any_onco_drug_primary_original,
      
      has_primary_vs_integrated_breast_change =
        any_breast_drug_integrated & !any_breast_drug_primary_original
    )
  
  save_result_local(
    target_master_integrated,
    "target_master_integrated_deep_indications.csv"
  )
  
  # ---------------------------------------------------------------------------
  # 12. Optional overwrite of target_master.csv with clean integrated columns
  # ---------------------------------------------------------------------------
  
  if (isTRUE(overwrite_target_master)) {
    original_cols <- c(
      "target_symbol", "target_ensembl", "conservation_score",
      "has_tractability_data", "any_sm", "any_ab", "any_pr", "any_oc",
      "any_druggable", "protac_like",
      "n_drugs", "n_indications", "n_oncology_drug_indications",
      "n_breast_drug_indications", "max_drug_phase",
      "any_approved", "n_rows_with_phase", "n_rows_without_phase",
      "any_onco_drug", "any_breast_drug"
    )
    
    target_master_overwrite <- target_master_integrated %>%
      mutate(
        n_drugs = n_drugs_integrated,
        n_indications = n_indications_integrated,
        n_oncology_drug_indications = n_oncology_drug_indications_integrated,
        n_breast_drug_indications = n_breast_drug_indications_integrated,
        max_drug_phase = max_drug_phase_integrated,
        any_approved = any_approved_integrated,
        n_rows_with_phase = n_rows_with_phase_integrated,
        n_rows_without_phase = n_rows_without_phase_integrated,
        any_onco_drug = any_onco_drug_integrated,
        any_breast_drug = any_breast_drug_integrated
      ) %>%
      select(any_of(original_cols))
    
    save_result_local(target_master_overwrite, "target_master.csv")
    assign("target_master", target_master_overwrite, envir = .GlobalEnv)
    
    cat("      NOTE: target_master.csv overwritten with clean integrated columns.\n")
    cat("            Full audit columns retained in target_master_integrated_deep_indications.csv\n")
  } else {
    assign("target_master_integrated_deep_indications", target_master_integrated, envir = .GlobalEnv)
    cat("      NOTE: target_master.csv NOT overwritten.\n")
    cat("            Set overwrite_target_master = TRUE to use integrated evidence in Step 6/Figure 3.\n")
  }
  
  # ---------------------------------------------------------------------------
  # 13. QC outputs
  # ---------------------------------------------------------------------------
  
  qc_phase_changes <- target_master_integrated %>%
    filter(has_primary_vs_integrated_phase_change %in% TRUE) %>%
    select(any_of(c(
      "target_symbol",
      "max_drug_phase_primary_original",
      "max_drug_phase_integrated",
      "n_drug_records_probe_linked",
      "n_breast_indications_integrated",
      "any_breast_drug_integrated",
      "has_probe_linked_drug_evidence"
    ))) %>%
    arrange(desc(max_drug_phase_integrated), target_symbol)
  
  qc_breast_changes <- target_master_integrated %>%
    filter(has_primary_vs_integrated_breast_change %in% TRUE) %>%
    select(any_of(c(
      "target_symbol",
      "any_breast_drug_primary_original",
      "any_breast_drug_integrated",
      "n_breast_drug_indications_primary_original",
      "n_breast_drug_indications_integrated",
      "n_drug_records_probe_linked",
      "has_probe_linked_drug_evidence"
    ))) %>%
    arrange(desc(n_breast_drug_indications_integrated), target_symbol)
  
  qc_oncology_changes <- target_master_integrated %>%
    filter(has_primary_vs_integrated_oncology_change %in% TRUE) %>%
    select(any_of(c(
      "target_symbol",
      "any_onco_drug_primary_original",
      "any_onco_drug_integrated",
      "n_oncology_drug_indications_primary_original",
      "n_oncology_drug_indications_integrated",
      "n_drug_records_probe_linked",
      "has_probe_linked_drug_evidence"
    ))) %>%
    arrange(desc(n_oncology_drug_indications_integrated), target_symbol)
  
  save_result_local(qc_phase_changes, "QC_targets_with_phase_gain_from_deep_indications.csv")
  save_result_local(qc_breast_changes, "QC_targets_with_breast_gain_from_deep_indications.csv")
  save_result_local(qc_oncology_changes, "QC_targets_with_oncology_gain_from_deep_indications.csv")
  
  cat("\n")
  cat("  Integrated indication expansion summary:\n")
  cat("    Targets with phase gain:          ", nrow(qc_phase_changes), "\n")
  cat("    Targets with oncology gain:       ", nrow(qc_oncology_changes), "\n")
  cat("    Targets with breast/mammary gain: ", nrow(qc_breast_changes), "\n")
  
  if ("TTK" %in% target_master_integrated$target_symbol) {
    cat("\n")
    cat("  TTK integrated evidence check:\n")
    
    ttk_check <- target_master_integrated %>%
      filter(target_symbol == "TTK") %>%
      select(any_of(c(
        "target_symbol",
        "max_drug_phase_primary_original",
        "max_drug_phase_integrated",
        "any_onco_drug_primary_original",
        "any_onco_drug_integrated",
        "any_breast_drug_primary_original",
        "any_breast_drug_integrated",
        "n_drugs_primary_original",
        "n_drugs_integrated",
        "n_indications_primary_original",
        "n_indications_integrated",
        "n_breast_drug_indications_primary_original",
        "n_breast_drug_indications_integrated",
        "has_probe_linked_drug_evidence"
      )))
    
    print(ttk_check)
    
    save_result_local(
      integrated_rows %>% filter(target_symbol == "TTK"),
      "QC_TTK_deep_indication_rows.csv"
    )
  }
  
  assign("OT_probe_links", probe_links, envir = .GlobalEnv)
  assign("OT_target_drug_deep_indication_rows", deep_rows, envir = .GlobalEnv)
  assign("OT_integrated_drug_indication_rows", integrated_rows, envir = .GlobalEnv)
  assign("target_deep_indication_summary", deep_summary, envir = .GlobalEnv)
  
  cat("\n")
  cat("  \u2713 STEP 5B COMPLETE\n\n")
  
  invisible(list(
    method_parameters = method_parameters,
    probe_links_all = probe_links_all,
    probe_links_high_quality = probe_links_hq,
    probe_drug_phases = probe_drug_phases,
    probe_links = probe_links,
    target_drug_links = target_drug_links,
    drug_indications = drug_indications,
    deep_rows = deep_rows,
    primary_target_rows = primary_target_rows,
    integrated_rows = integrated_rows,
    deep_summary = deep_summary,
    target_master_integrated = target_master_integrated,
    qc_phase_changes = qc_phase_changes,
    qc_oncology_changes = qc_oncology_changes,
    qc_breast_changes = qc_breast_changes
  ))
}
