################################################################################
# STEP 5: DRUGGABILITY ASSESSMENT (Open Targets + ChEMBL)
################################################################################

module_drug_targets <- function() {
  
  cat("  [5A] Checking required upstream objects...\n")
  
  required_objs <- c("core_up_genes_u", "mapping_1to1")
  missing_objs <- required_objs[
    !vapply(required_objs, exists, logical(1), envir = .GlobalEnv)
  ]
  
  if (length(missing_objs) > 0) {
    stop(
      "Missing required object(s) for druggability step: ",
      paste(missing_objs, collapse = ", "),
      ". Run previous modules first.",
      call. = FALSE
    )
  }
  
  core_up_genes_u <- get("core_up_genes_u", envir = .GlobalEnv)
  mapping_1to1 <- get("mapping_1to1", envir = .GlobalEnv)
  
  has_conservation_objs <- all(vapply(
    c("cons_up", "dog_le_up_u", "cat_up_u"),
    exists,
    logical(1),
    envir = .GlobalEnv
  ))
  
  if (has_conservation_objs) {
    cons_up <- get("cons_up", envir = .GlobalEnv)
    dog_le_up_u <- get("dog_le_up_u", envir = .GlobalEnv)
    cat_up_u <- get("cat_up_u", envir = .GlobalEnv)
  }
  
  # ---------------------------------------------------------------------------
  # Helper functions
  # ---------------------------------------------------------------------------
  
  safe_chr <- function(x, default = NA_character_) {
    if (is.null(x) || length(x) == 0) return(default)
    x <- x[[1]]
    if (is.null(x) || length(x) == 0 || is.na(x)) return(default)
    as.character(x)
  }
  
  safe_bool <- function(x, default = FALSE) {
    if (is.null(x) || length(x) == 0) return(default)
    if (is.logical(x[[1]])) return(isTRUE(x[[1]]))
    
    x_chr <- tolower(trimws(as.character(x[[1]])))
    x_chr %in% c("true", "t", "1", "yes", "y")
  }
  
  stage_to_numeric <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA_real_)
    
    x_chr <- toupper(trimws(as.character(x)))
    x_chr[x_chr %in% c("", "NA", "NAN", "NULL", "NONE")] <- NA_character_
    
    out <- rep(NA_real_, length(x_chr))
    
    # Open Targets labels observed in drugAndClinicalCandidates:
    # APPROVAL, PHASE_3, PHASE_2_3, PHASE_2, PHASE_1_2, PHASE_1
    out[grepl("APPROVAL|APPROVED|REGULATORY|MARKETED", x_chr)] <- 4
    
    out[is.na(out) & grepl(
      "PHASE_3|PHASE_III|PHASE_2_3|PHASE_II_III|PHASE 3|PHASE III|PHASE 2/3|PHASE II/III",
      x_chr
    )] <- 3
    
    out[is.na(out) & grepl(
      "PHASE_2|PHASE_II|PHASE_1_2|PHASE_I_II|PHASE 2|PHASE II|PHASE 1/2|PHASE I/II",
      x_chr
    )] <- 2
    
    out[is.na(out) & grepl(
      "PHASE_1|PHASE_I|PHASE 1|PHASE I",
      x_chr
    )] <- 1
    
    out[is.na(out) & grepl(
      "PHASE_0|PHASE 0|PRECLINICAL|PRE_CLINICAL|PRE-CLINICAL|PRE CLINICAL",
      x_chr
    )] <- 0
    
    numeric_direct <- suppressWarnings(as.numeric(x_chr))
    idx <- is.na(out) & !is.na(numeric_direct)
    out[idx] <- numeric_direct[idx]
    
    out
  }
  
  extract_ot_error <- function(x) {
    if (isTRUE(x$ok)) return(NA_character_)
    
    if (!is.null(x$graphql_errors)) {
      msgs <- vapply(
        x$graphql_errors,
        function(e) {
          if (is.list(e) && !is.null(e$message)) {
            return(as.character(e$message))
          }
          paste(capture.output(str(e)), collapse = " | ")
        },
        character(1)
      )
      return(paste(msgs, collapse = " || "))
    }
    
    if (!is.null(x$text)) {
      return(substr(as.character(x$text), 1, 2000))
    }
    
    if (!is.null(x$status)) {
      return(paste0("HTTP status: ", x$status))
    }
    
    "Unknown Open Targets query failure"
  }
  
  as_logical_robust <- function(x) {
    if (is.logical(x)) return(ifelse(is.na(x), NA, x))
    if (is.numeric(x)) return(ifelse(is.na(x), NA, x != 0))
    
    x_chr <- tolower(trimws(as.character(x)))
    out <- rep(NA, length(x_chr))
    out[x_chr %in% c("true", "t", "1", "yes", "y")] <- TRUE
    out[x_chr %in% c("false", "f", "0", "no", "n")] <- FALSE
    out
  }
  
  # ---------------------------------------------------------------------------
  # Map strict conserved core symbols to human Ensembl IDs
  # ---------------------------------------------------------------------------
  
  cat("  [5B] Mapping core symbols to human ENSG...\n")
  
  map_clean <- mapping_1to1 %>%
    dplyr::select(human_symbol, human_ensembl) %>%
    dplyr::filter(
      !is.na(human_symbol), human_symbol != "",
      !is.na(human_ensembl), human_ensembl != ""
    ) %>%
    dplyr::distinct()
  
  collisions <- map_clean %>%
    dplyr::count(human_symbol) %>%
    dplyr::filter(n > 1)
  
  if (nrow(collisions) > 0) {
    cat(
      "      WARNING:",
      nrow(collisions),
      "symbol→ENSG collisions found; excluding them\n"
    )
    
    map_clean <- map_clean %>%
      dplyr::anti_join(collisions, by = "human_symbol")
  }
  
  sym_to_ensg <- setNames(map_clean$human_ensembl, map_clean$human_symbol)
  
  core_up_ensg <- unique(unname(sym_to_ensg[core_up_genes_u]))
  core_up_ensg <- core_up_ensg[!is.na(core_up_ensg) & core_up_ensg != ""]
  
  cat("      Core symbols:", length(core_up_genes_u), "\n")
  cat("      Mapped ENSG (strict 1:1):", length(core_up_ensg), "\n")
  
  if (length(core_up_ensg) == 0) {
    stop("No core symbols could be mapped to human ENSG.", call. = FALSE)
  }
  
  # ---------------------------------------------------------------------------
  # Open Targets GraphQL query
  # ---------------------------------------------------------------------------
  
  cat("  [5C] Preparing Open Targets queries...\n")
  
  OT_URL <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  # Fresh cache namespace for current Open Targets schema.
  # Use short cache folder names to avoid Windows MAX_PATH issues.
  ot_cache_dir <- file.path(PATH_CONFIG$cache_dir, "OTv3")
  chembl_cache_dir <- file.path(PATH_CONFIG$cache_dir, "ChEMBL")
  
  dir.create(ot_cache_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(chembl_cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (!dir.exists(ot_cache_dir)) {
    stop(
      "Could not create Open Targets cache directory: ",
      ot_cache_dir,
      "\nThis is usually a Windows path-length issue. Move the project to a shorter path, e.g. C:/omics_q1/",
      call. = FALSE
    )
  }
  
  if (!dir.exists(chembl_cache_dir)) {
    stop(
      "Could not create ChEMBL cache directory: ",
      chembl_cache_dir,
      call. = FALSE
    )
  }
  ot_post <- function(query, variables = list()) {
    res <- httr::POST(
      OT_URL,
      body = list(query = query, variables = variables),
      encode = "json",
      httr::add_headers(
        `Content-Type` = "application/json",
        `Accept` = "application/json",
        `User-Agent` = "R httr OpenTargets client"
      ),
      httr::timeout(60)
    )
    
    txt <- httr::content(res, "text", encoding = "UTF-8")
    
    if (httr::status_code(res) != 200) {
      return(list(
        ok = FALSE,
        status = httr::status_code(res),
        text = txt
      ))
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
  
  ot_query_target <- function(ensembl_id) {
    query <- "
      query targetDrugAndTractability($ensemblId: String!) {
        target(ensemblId: $ensemblId) {
          id
          approvedSymbol
          approvedName
          biotype
          tractability {
            label
            modality
            value
          }
          drugAndClinicalCandidates {
            count
            rows {
              id
              maxClinicalStage
              drug {
                id
                name
                drugType
                maximumClinicalStage
              }
              diseases {
                diseaseFromSource
                disease {
                  id
                  name
                }
              }
            }
          }
        }
      }
    "
    
    ot_post(query, variables = list(ensemblId = ensembl_id))
  }
  
  get_target_cached <- function(ensembl_id) {
    # Defensive cache-folder creation. This prevents:
    # cannot open file ... No such file or directory
    dir.create(ot_cache_dir, showWarnings = FALSE, recursive = TRUE)
    
    f <- file.path(ot_cache_dir, paste0(ensembl_id, ".json"))
    
    if (file.exists(f)) {
      j <- jsonlite::fromJSON(readr::read_file(f), simplifyVector = FALSE)
      return(list(ok = TRUE, data = j, from_cache = TRUE))
    }
    
    out <- ot_query_target(ensembl_id)
    
    if (isTRUE(out$ok)) {
      dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
      
      writeLines(
        jsonlite::toJSON(out$data, auto_unbox = TRUE, pretty = TRUE),
        con = f
      )
    }
    
    Sys.sleep(0.4)
    out
  }
  
  cat("  [5D] Querying Open Targets for clinical candidates and tractability...\n")
  cat("      Querying", length(core_up_ensg), "targets...\n")
  
  target_raw <- lapply(core_up_ensg, get_target_cached)
  names(target_raw) <- core_up_ensg
  
  target_ok <- target_raw[
    vapply(target_raw, function(x) isTRUE(x$ok), logical(1))
  ]
  
  target_fail <- target_raw[
    !vapply(target_raw, function(x) isTRUE(x$ok), logical(1))
  ]
  
  ot_query_qc <- data.frame(
    n_targets_queried = length(target_raw),
    n_targets_success = length(target_ok),
    n_targets_failed = length(target_fail),
    failure_pct = ifelse(
      length(target_raw) > 0,
      100 * length(target_fail) / length(target_raw),
      NA_real_
    ),
    stringsAsFactors = FALSE
  )
  
  ot_error_log <- data.frame(
    target_ensembl = names(target_raw),
    ok = vapply(target_raw, function(x) isTRUE(x$ok), logical(1)),
    error_message = vapply(target_raw, extract_ot_error, character(1)),
    stringsAsFactors = FALSE
  )
  
  save_result(ot_query_qc, "OT_query_QC.csv", "csv")
  save_result(ot_error_log, "OT_query_error_log.csv", "csv")
  
  cat("      Success:", length(target_ok), "| Failed:", length(target_fail), "\n")
  
  if (length(target_fail) > 0) {
    save_result(names(target_fail), "OT_failed_targets_ENSG.txt", "txt")
    saveRDS(
      target_fail,
      file.path(PATH_CONFIG$cache_dir, "OT_failed_targets_payloads.rds")
    )
  }
  
  if (length(target_ok) == 0) {
    stop(
      "Open Targets failed for all queried targets. ",
      "See results/OT_query_error_log.csv before using translational outputs.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Extract drug/clinical-candidate evidence
  # ---------------------------------------------------------------------------
  
  cat("  [5E] Extracting target-level clinical candidate evidence...\n")
  
  extract_drugs <- function(x, ensembl_id) {
    t <- x$data$data$target
    if (is.null(t)) return(NULL)
    
    clinical <- t$drugAndClinicalCandidates
    if (is.null(clinical)) return(NULL)
    
    rows <- clinical$rows
    if (is.null(rows) || length(rows) == 0) return(NULL)
    
    out <- list()
    
    for (i in seq_along(rows)) {
      r <- rows[[i]]
      drug <- r$drug
      
      drug_id <- safe_chr(drug$id)
      drug_name <- safe_chr(drug$name)
      drug_type <- safe_chr(drug$drugType)
      
      row_stage_chr <- safe_chr(r$maxClinicalStage)
      drug_stage_chr <- safe_chr(drug$maximumClinicalStage)
      
      row_stage_num <- stage_to_numeric(row_stage_chr)[1]
      drug_stage_num <- stage_to_numeric(drug_stage_chr)[1]
      
      phase_for_record <- suppressWarnings(max(
        c(row_stage_num, drug_stage_num),
        na.rm = TRUE
      ))
      
      if (is.infinite(phase_for_record)) {
        phase_for_record <- NA_real_
      }
      
      drug_is_approved <- isTRUE(!is.na(row_stage_num) && row_stage_num >= 4) ||
        isTRUE(!is.na(drug_stage_num) && drug_stage_num >= 4)
      
      diseases <- r$diseases
      
      if (is.null(diseases) || length(diseases) == 0) {
        out[[length(out) + 1]] <- data.frame(
          target_ensembl = ensembl_id,
          target_symbol = safe_chr(t$approvedSymbol),
          target_name = safe_chr(t$approvedName),
          phase_evidence = phase_for_record,
          phase_evidence_raw = row_stage_chr,
          moa = NA_character_,
          drug_chembl = drug_id,
          drug_name = drug_name,
          drug_type = drug_type,
          drug_is_approved = drug_is_approved,
          drug_max_phase = drug_stage_num,
          drug_max_phase_raw = drug_stage_chr,
          drug_year_first_approval = NA_real_,
          disease_id = NA_character_,
          disease_name = NA_character_,
          disease_from_source = NA_character_,
          stringsAsFactors = FALSE
        )
      } else {
        for (j in seq_along(diseases)) {
          d <- diseases[[j]]
          disease <- d$disease
          
          out[[length(out) + 1]] <- data.frame(
            target_ensembl = ensembl_id,
            target_symbol = safe_chr(t$approvedSymbol),
            target_name = safe_chr(t$approvedName),
            phase_evidence = phase_for_record,
            phase_evidence_raw = row_stage_chr,
            moa = NA_character_,
            drug_chembl = drug_id,
            drug_name = drug_name,
            drug_type = drug_type,
            drug_is_approved = drug_is_approved,
            drug_max_phase = drug_stage_num,
            drug_max_phase_raw = drug_stage_chr,
            drug_year_first_approval = NA_real_,
            disease_id = safe_chr(disease$id),
            disease_name = safe_chr(disease$name),
            disease_from_source = safe_chr(d$diseaseFromSource),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    dplyr::bind_rows(out)
  }
  
  drug_target <- dplyr::bind_rows(
    Map(extract_drugs, target_ok, names(target_ok))
  ) %>%
    dplyr::distinct()
  
  if (is.null(drug_target) || nrow(drug_target) == 0) {
    warning("No clinical drug/candidate evidence found for the queried targets.")
    
    drug_target <- data.frame(
      target_ensembl = character(),
      target_symbol = character(),
      target_name = character(),
      phase_evidence = numeric(),
      phase_evidence_raw = character(),
      moa = character(),
      drug_chembl = character(),
      drug_name = character(),
      drug_type = character(),
      drug_is_approved = logical(),
      drug_max_phase = numeric(),
      drug_max_phase_raw = character(),
      drug_year_first_approval = numeric(),
      disease_id = character(),
      disease_name = character(),
      disease_from_source = character(),
      stringsAsFactors = FALSE
    )
  }
  
  cat("      Raw clinical drug/candidate rows:", nrow(drug_target), "\n")
  save_result(drug_target, "drug_target_table_raw_OpenTargets.csv", "csv")
  
  phase_diagnostic <- drug_target %>%
    dplyr::mutate(
      phase_evidence_from_raw = stage_to_numeric(phase_evidence_raw),
      drug_max_phase_from_raw = stage_to_numeric(drug_max_phase_raw),
      phase_for_filter_preview = dplyr::case_when(
        !is.na(phase_evidence_from_raw) & !is.na(drug_max_phase_from_raw) ~
          pmax(phase_evidence_from_raw, drug_max_phase_from_raw),
        !is.na(phase_evidence_from_raw) ~ phase_evidence_from_raw,
        !is.na(drug_max_phase_from_raw) ~ drug_max_phase_from_raw,
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::summarise(
      n_raw_rows = dplyr::n(),
      n_with_drug_id = sum(!is.na(drug_chembl) & drug_chembl != ""),
      n_with_phase_evidence_raw = sum(!is.na(phase_evidence_raw) & phase_evidence_raw != ""),
      n_with_drug_max_phase_raw = sum(!is.na(drug_max_phase_raw) & drug_max_phase_raw != ""),
      n_with_parsed_phase = sum(!is.na(phase_for_filter_preview)),
      min_parsed_phase = suppressWarnings(min(phase_for_filter_preview, na.rm = TRUE)),
      max_parsed_phase = suppressWarnings(max(phase_for_filter_preview, na.rm = TRUE))
    ) %>%
    dplyr::mutate(
      min_parsed_phase = ifelse(is.infinite(min_parsed_phase), NA_real_, min_parsed_phase),
      max_parsed_phase = ifelse(is.infinite(max_parsed_phase), NA_real_, max_parsed_phase)
    )
  
  save_result(phase_diagnostic, "OT_phase_parsing_diagnostic.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # Filter drug evidence
  # ---------------------------------------------------------------------------
  
  cat("  [5F] Applying drug evidence filters...\n")
  
  drug_target_filt <- drug_target %>%
    dplyr::mutate(
      phase_evidence = dplyr::coalesce(
        suppressWarnings(as.numeric(phase_evidence)),
        stage_to_numeric(phase_evidence_raw)
      ),
      drug_max_phase = dplyr::coalesce(
        suppressWarnings(as.numeric(drug_max_phase)),
        stage_to_numeric(drug_max_phase_raw)
      ),
      phase_for_filter = dplyr::case_when(
        !is.na(phase_evidence) & !is.na(drug_max_phase) ~ pmax(phase_evidence, drug_max_phase),
        !is.na(phase_evidence) ~ phase_evidence,
        !is.na(drug_max_phase) ~ drug_max_phase,
        TRUE ~ NA_real_
      ),
      drug_is_approved = drug_is_approved | (!is.na(phase_for_filter) & phase_for_filter >= 4),
      has_phase_information = !is.na(phase_for_filter),
      has_valid_drug_id = !is.na(drug_chembl) & drug_chembl != ""
    ) %>%
    dplyr::filter(has_valid_drug_id)
  
  if (isTRUE(PARAM_CONFIG$drug_keep_approved_only)) {
    drug_target_filt <- drug_target_filt %>%
      dplyr::filter(drug_is_approved %in% TRUE)
  }
  
  # Important:
  # Missing-phase rows are retained as target-linked drug/candidate evidence.
  # Clinical maturity is assessed later using phase_for_filter and drug_is_approved.
  if (!isTRUE(PARAM_CONFIG$drug_keep_na_phase)) {
    cat("      NOTE: missing-phase drug rows are retained as drug-linked evidence.\n")
    cat("            Phase is used for clinical maturity, not for deleting evidence.\n")
  }
  
  cat("      Rows with valid drug IDs:", nrow(drug_target_filt), "\n")
  cat("      Rows with parsed phase:", sum(!is.na(drug_target_filt$phase_for_filter)), "\n")
  cat("      Distinct drugs:", dplyr::n_distinct(drug_target_filt$drug_chembl), "\n")
  
  # ---------------------------------------------------------------------------
  # Disease indication flags
  # ---------------------------------------------------------------------------
  
  cat("  [5G] Adding disease indication flags...\n")
  
  if (nrow(drug_target_filt) > 0) {
    drug_target_filt <- drug_target_filt %>%
      dplyr::mutate(
        disease_name_lc = tolower(dplyr::coalesce(disease_name, disease_from_source, "")),
        is_oncology_drug_indication = stringr::str_detect(
          disease_name_lc,
          "cancer|carcinoma|sarcoma|tumou?r|neoplasm|leukemia|lymphoma|melanoma|glioma|blastoma|myeloma|malignan"
        ),
        is_breast_drug_indication = stringr::str_detect(
          disease_name_lc,
          "breast|mammary"
        )
      ) %>%
      dplyr::select(-disease_name_lc)
  } else {
    drug_target_filt$is_oncology_drug_indication <- logical(0)
    drug_target_filt$is_breast_drug_indication <- logical(0)
  }
  
  # ---------------------------------------------------------------------------
  # ChEMBL ATC classification
  # ---------------------------------------------------------------------------
  
  cat("  [5H] Querying ChEMBL for ATC classifications...\n")
  
  CHEMBL_BASE <- "https://www.ebi.ac.uk/chembl/api/data/"
  
  chembl_get_atc <- function(chembl_id) {
    url <- paste0(
      CHEMBL_BASE,
      "atc_class?molecule_chembl_id=",
      chembl_id,
      "&limit=1000"
    )
    
    res <- httr::GET(
      url,
      httr::add_headers(`Accept` = "application/json"),
      httr::timeout(60)
    )
    
    if (httr::status_code(res) != 200) {
      return(NULL)
    }
    
    j <- jsonlite::fromJSON(
      httr::content(res, "text", encoding = "UTF-8")
    )
    
    if (is.null(j$atc_classes) || nrow(j$atc_classes) == 0) {
      return(NULL)
    }
    
    df <- as.data.frame(j$atc_classes, stringsAsFactors = FALSE)
    df$molecule_chembl_id <- chembl_id
    df
  }
  
  get_atc_cached <- function(chembl_id) {
    # Defensive cache-folder creation.
    dir.create(chembl_cache_dir, showWarnings = FALSE, recursive = TRUE)
    
    f <- file.path(chembl_cache_dir, paste0(chembl_id, "_atc.rds"))
    
    if (file.exists(f)) {
      return(readRDS(f))
    }
    
    out <- chembl_get_atc(chembl_id)
    
    dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
    saveRDS(out, f)
    
    Sys.sleep(0.2)
    out
  }
  
  unique_drugs <- if (nrow(drug_target_filt) > 0) {
    sort(unique(drug_target_filt$drug_chembl))
  } else {
    character(0)
  }
  
  atc_list <- lapply(unique_drugs, get_atc_cached)
  atc_list <- atc_list[!vapply(atc_list, is.null, logical(1))]
  
  if (length(atc_list) == 0) {
    cat("No ATC classifications found; using Open Targets drug_type as fallback\n")
    
    atc_df <- data.frame(
      molecule_chembl_id = character(),
      atc_code = character(),
      stringsAsFactors = FALSE
    )
  } else {
    atc_df <- dplyr::bind_rows(atc_list)
    
    if (!"molecule_chembl_id" %in% colnames(atc_df)) {
      atc_df$molecule_chembl_id <- character()
    }
    
    if (!"atc_code" %in% colnames(atc_df)) {
      atc_df$atc_code <- NA_character_
    }
  }
  
  drug_target_with_atc <- drug_target_filt %>%
    dplyr::left_join(
      atc_df,
      by = c("drug_chembl" = "molecule_chembl_id")
    ) %>%
    dplyr::mutate(
      atc_missing = is.na(atc_code) | atc_code == "",
      class_fallback = ifelse(atc_missing, drug_type, atc_code)
    )
  
  save_result(
    drug_target_with_atc,
    "drug_target_table_with_ATC_and_fallback.csv",
    "csv"
  )
  
  # ---------------------------------------------------------------------------
  # Extract tractability data
  # ---------------------------------------------------------------------------
  
  cat("  [5I] Extracting tractability data...\n")
  
  extract_tractability <- function(x, ensembl_id) {
    t <- x$data$data$target
    if (is.null(t)) return(NULL)
    
    tr <- t$tractability
    if (is.null(tr) || length(tr) == 0) return(NULL)
    
    out <- lapply(tr, function(z) {
      data.frame(
        target_ensembl = ensembl_id,
        target_symbol = safe_chr(t$approvedSymbol),
        label = safe_chr(z$label),
        modality = safe_chr(z$modality),
        tractability_value = safe_bool(z$value, default = NA),
        stringsAsFactors = FALSE
      )
    })
    
    dplyr::bind_rows(out)
  }
  
  tractability_table <- dplyr::bind_rows(
    Map(extract_tractability, target_ok, names(target_ok))
  )
  
  if (is.null(tractability_table) || nrow(tractability_table) == 0) {
    tractability_table <- data.frame(
      target_ensembl = character(),
      target_symbol = character(),
      label = character(),
      modality = character(),
      tractability_value = logical(),
      stringsAsFactors = FALSE
    )
  }
  
  save_result(tractability_table, "tractability_table_raw.csv", "csv")
  
  tract_summary_fixed <- tractability_table %>%
    dplyr::mutate(
      tractability_value = as_logical_robust(tractability_value),
      modality = as.character(modality)
    ) %>%
    dplyr::group_by(target_ensembl, target_symbol) %>%
    dplyr::summarise(
      has_tractability_data = any(!is.na(tractability_value)),
      any_sm = any(modality == "SM" & tractability_value %in% TRUE, na.rm = TRUE),
      any_ab = any(modality == "AB" & tractability_value %in% TRUE, na.rm = TRUE),
      any_pr = any(modality == "PR" & tractability_value %in% TRUE, na.rm = TRUE),
      any_oc = any(modality == "OC" & tractability_value %in% TRUE, na.rm = TRUE),
      any_druggable = any(tractability_value %in% TRUE, na.rm = TRUE),
      protac_like = any_pr,
      .groups = "drop"
    )
  
  save_result(
    tract_summary_fixed,
    "TRACTABILITY_summary_by_target_fixed.csv",
    "csv"
  )
  
  # ---------------------------------------------------------------------------
  # Conservation score
  # ---------------------------------------------------------------------------
  
  cat("  [5J] Building conservation score table...\n")
  
  scored_targets <- sort(unique(c(core_up_genes_u, drug_target_with_atc$target_symbol)))
  
  conservation_score <- data.frame(
    target_symbol = scored_targets,
    stringsAsFactors = FALSE
  )
  
  conservation_score$conservation_score <- 0
  
  if (has_conservation_objs && length(cons_up) > 0) {
    gene_freq_up <- table(unlist(lapply(cons_up, function(h) {
      intersect(dog_le_up_u[[h]], cat_up_u)
    })))
    
    if (length(gene_freq_up) > 0) {
      tmp <- as.numeric(gene_freq_up[conservation_score$target_symbol])
      tmp[is.na(tmp)] <- 0
      conservation_score$conservation_score <- tmp
    }
  }
  
  save_result(conservation_score, "conservation_score_targets.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # Target-level drug summary
  # ---------------------------------------------------------------------------
  
  cat("  [5K] Building target-level drug summary...\n")
  
  if (nrow(drug_target_with_atc) > 0) {
    target_drug_summary <- drug_target_with_atc %>%
      dplyr::group_by(target_symbol, target_ensembl) %>%
      dplyr::summarise(
        n_drugs = dplyr::n_distinct(drug_chembl[!is.na(drug_chembl) & drug_chembl != ""]),
        n_indications = dplyr::n_distinct(disease_id[!is.na(disease_id) & disease_id != ""]),
        n_oncology_drug_indications = dplyr::n_distinct(
          disease_id[
            is_oncology_drug_indication %in% TRUE &
              !is.na(disease_id) & disease_id != ""
          ]
        ),
        n_breast_drug_indications = dplyr::n_distinct(
          disease_id[
            is_breast_drug_indication %in% TRUE &
              !is.na(disease_id) & disease_id != ""
          ]
        ),
        max_drug_phase = ifelse(
          all(is.na(phase_for_filter)),
          NA_real_,
          max(phase_for_filter, na.rm = TRUE)
        ),
        any_approved = any(drug_is_approved %in% TRUE),
        n_rows_with_phase = sum(!is.na(phase_for_filter)),
        n_rows_without_phase = sum(is.na(phase_for_filter)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        max_drug_phase = ifelse(
          is.infinite(max_drug_phase),
          NA_real_,
          max_drug_phase
        ),
        any_onco_drug = n_oncology_drug_indications > 0,
        any_breast_drug = n_breast_drug_indications > 0
      )
  } else {
    target_drug_summary <- data.frame(
      target_symbol = character(),
      target_ensembl = character(),
      n_drugs = integer(),
      n_indications = integer(),
      n_oncology_drug_indications = integer(),
      n_breast_drug_indications = integer(),
      max_drug_phase = numeric(),
      any_approved = logical(),
      n_rows_with_phase = integer(),
      n_rows_without_phase = integer(),
      any_onco_drug = logical(),
      any_breast_drug = logical(),
      stringsAsFactors = FALSE
    )
  }
  
  save_result(
    target_drug_summary,
    "target_drug_summary_with_indication_counts.csv",
    "csv"
  )
  
  # ---------------------------------------------------------------------------
  # Integrated target master table
  # ---------------------------------------------------------------------------
  
  cat("  [5L] Building integrated target master table...\n")
  
  target_master <- data.frame(
    target_symbol = sort(unique(core_up_genes_u)),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(
      data.frame(
        target_symbol = names(sym_to_ensg),
        target_ensembl = unname(sym_to_ensg),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::distinct(),
      by = "target_symbol"
    ) %>%
    dplyr::left_join(conservation_score, by = "target_symbol") %>%
    dplyr::left_join(
      tract_summary_fixed %>%
        dplyr::select(
          target_symbol,
          target_ensembl,
          has_tractability_data,
          any_sm,
          any_ab,
          any_pr,
          any_oc,
          any_druggable,
          protac_like
        ),
      by = c("target_symbol", "target_ensembl")
    ) %>%
    dplyr::left_join(
      target_drug_summary,
      by = c("target_symbol", "target_ensembl")
    ) %>%
    dplyr::mutate(
      conservation_score = dplyr::coalesce(conservation_score, 0),
      has_tractability_data = dplyr::coalesce(has_tractability_data, FALSE),
      any_sm = dplyr::coalesce(any_sm, FALSE),
      any_ab = dplyr::coalesce(any_ab, FALSE),
      any_pr = dplyr::coalesce(any_pr, FALSE),
      any_oc = dplyr::coalesce(any_oc, FALSE),
      any_druggable = dplyr::coalesce(any_druggable, FALSE),
      protac_like = dplyr::coalesce(protac_like, FALSE),
      n_drugs = dplyr::coalesce(n_drugs, 0L),
      n_indications = dplyr::coalesce(n_indications, 0L),
      n_oncology_drug_indications = dplyr::coalesce(
        n_oncology_drug_indications,
        0L
      ),
      n_breast_drug_indications = dplyr::coalesce(
        n_breast_drug_indications,
        0L
      ),
      max_drug_phase = dplyr::coalesce(max_drug_phase, NA_real_),
      any_approved = dplyr::coalesce(any_approved, FALSE),
      n_rows_with_phase = dplyr::coalesce(n_rows_with_phase, 0L),
      n_rows_without_phase = dplyr::coalesce(n_rows_without_phase, 0L),
      any_onco_drug = dplyr::coalesce(any_onco_drug, FALSE),
      any_breast_drug = dplyr::coalesce(any_breast_drug, FALSE)
    )
  
  save_result(target_master, "target_master.csv", "csv")
  # Save primary target-level evidence snapshot before any Step 5B integration.
  # This file is used for the primary-vs-integrated prioritisation sensitivity analysis.
  save_result(target_master, "target_master_primary_targetlevel.csv", "csv")
  
  assign("target_master_primary_targetlevel", target_master, envir = .GlobalEnv)
  
  # ---------------------------------------------------------------------------
  # Export objects to GlobalEnv for downstream modules
  # ---------------------------------------------------------------------------
  
  assign("core_up_ensg", core_up_ensg, envir = .GlobalEnv)
  assign("target_raw", target_raw, envir = .GlobalEnv)
  assign("target_ok", target_ok, envir = .GlobalEnv)
  assign("target_fail", target_fail, envir = .GlobalEnv)
  assign("ot_query_qc", ot_query_qc, envir = .GlobalEnv)
  assign("ot_error_log", ot_error_log, envir = .GlobalEnv)
  
  assign("drug_target", drug_target, envir = .GlobalEnv)
  assign("drug_target_filt", drug_target_filt, envir = .GlobalEnv)
  assign("drug_target_with_atc", drug_target_with_atc, envir = .GlobalEnv)
  assign("tractability_table", tractability_table, envir = .GlobalEnv)
  assign("tract_summary_fixed", tract_summary_fixed, envir = .GlobalEnv)
  assign("conservation_score", conservation_score, envir = .GlobalEnv)
  assign("target_drug_summary", target_drug_summary, envir = .GlobalEnv)
  assign("target_master", target_master, envir = .GlobalEnv)
  
  cat("  ✓ STEP 5 COMPLETE\n\n")
  
  invisible(list(
    core_up_ensg = core_up_ensg,
    target_raw = target_raw,
    target_ok = target_ok,
    target_fail = target_fail,
    ot_query_qc = ot_query_qc,
    ot_error_log = ot_error_log,
    drug_target = drug_target,
    drug_target_filt = drug_target_filt,
    drug_target_with_atc = drug_target_with_atc,
    tractability_table = tractability_table,
    tract_summary_fixed = tract_summary_fixed,
    conservation_score = conservation_score,
    target_drug_summary = target_drug_summary,
    target_master = target_master
  ))
}