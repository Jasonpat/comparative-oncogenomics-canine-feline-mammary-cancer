################################################################################
# VALIDATION SUPPORT: BUILD LIGHT BACKGROUND TRACTABILITY UNIVERSE
#
# Hallmark-constrained background:
#   all dog leading-edge genes from all significant dog hallmarks mapped to human
#   targets, then annotated through Open Targets tractability and clinical/drug
#   candidate evidence.
#
# Reviewer/reproducibility fixes:
#   - uses current Open Targets target { tractability, drugAndClinicalCandidates }
#     schema, matching the corrected Step 5 logic
#   - uses a short shared cache folder, cache/OTv3, to reduce Windows path issues
#   - writes Open Targets QC and error logs for reproducibility/debugging
#   - keeps drug rows even when phase is missing; phase is used only to summarise
#     clinical maturity
#   - final background uniqueness is enforced by target_ensembl
################################################################################

module_build_background_tractability_universe_light <- function(
    padj_cutoff = 0.05,
    max_targets = NULL,
    ot_page_size = 200,      # retained for backward compatibility; not used by current OT schema
    ot_max_pages = 200,      # retained for backward compatibility; not used by current OT schema
    progress_every = 50,
    force_refresh_ot = FALSE
) {
  
  cat("\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat("BUILD LIGHT BACKGROUND TRACTABILITY UNIVERSE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(httr)
    library(jsonlite)
    library(tibble)
    library(stringr)
  })
  
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
  
  as_logical_robust <- function(x) {
    if (is.logical(x)) return(ifelse(is.na(x), NA, x))
    if (is.numeric(x)) return(ifelse(is.na(x), NA, x != 0))
    x_chr <- tolower(trimws(as.character(x)))
    out <- rep(NA, length(x_chr))
    out[x_chr %in% c("true", "t", "1", "yes", "y")] <- TRUE
    out[x_chr %in% c("false", "f", "0", "no", "n")] <- FALSE
    out
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
          if (is.list(e) && !is.null(e$message)) return(as.character(e$message))
          paste(capture.output(str(e)), collapse = " | ")
        },
        character(1)
      )
      return(paste(msgs, collapse = " || "))
    }
    
    if (!is.null(x$text)) return(substr(as.character(x$text), 1, 2000))
    if (!is.null(x$status)) return(paste0("HTTP status: ", x$status))
    
    "Unknown Open Targets query failure"
  }
  
  # ---------------------------------------------------------------------------
  # 1. File checks
  # ---------------------------------------------------------------------------
  req_files <- c(
    file.path(PATH_CONFIG$results_dir, "DOG_mapping_1to1_highconf.rds"),
    file.path(PATH_CONFIG$results_dir, "DOG_Hallmark_fgsea.rds")
  )
  
  missing_files <- req_files[!file.exists(req_files)]
  if (length(missing_files) > 0) {
    stop(
      "Missing required file(s) for light background universe build:\n",
      paste(" - ", missing_files, collapse = "\n"),
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # 2. Load inputs
  # ---------------------------------------------------------------------------
  cat("  [BG-L1] Loading dog mapping and hallmark results...\n")
  
  mapping_1to1 <- readRDS(file.path(PATH_CONFIG$results_dir, "DOG_mapping_1to1_highconf.rds"))
  fg <- readRDS(file.path(PATH_CONFIG$results_dir, "DOG_Hallmark_fgsea.rds"))
  
  map_tbl <- mapping_1to1 %>%
    dplyr::select(human_ensembl, human_symbol) %>%
    dplyr::filter(
      !is.na(human_ensembl), human_ensembl != "",
      !is.na(human_symbol), human_symbol != ""
    ) %>%
    dplyr::distinct()
  
  collisions <- map_tbl %>%
    dplyr::count(human_symbol) %>%
    dplyr::filter(n > 1)
  
  if (nrow(collisions) > 0) {
    cat("      WARNING:", nrow(collisions), "symbol→ENSG collisions found; excluding them\n")
    map_tbl <- map_tbl %>%
      dplyr::anti_join(collisions, by = "human_symbol")
  }
  
  ensg_to_sym_vec <- setNames(map_tbl$human_symbol, map_tbl$human_ensembl)
  
  # ---------------------------------------------------------------------------
  # 3. Build hallmark-supported background candidates
  # ---------------------------------------------------------------------------
  cat("  [BG-L2] Extracting genes from all significant dog hallmarks...\n")
  
  dog_sig <- fg %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff)
  
  if (nrow(dog_sig) == 0) {
    stop("No significant dog hallmarks found at padj < ", padj_cutoff, call. = FALSE)
  }
  
  all_le_ensg <- unique(unlist(dog_sig$leadingEdge))
  all_le_ensg <- all_le_ensg[!is.na(all_le_ensg) & all_le_ensg != ""]
  
  background_candidates <- tibble::tibble(
    target_ensembl = all_le_ensg,
    target_symbol_input = unname(ensg_to_sym_vec[all_le_ensg])
  ) %>%
    dplyr::filter(!is.na(target_symbol_input), target_symbol_input != "") %>%
    dplyr::distinct(target_ensembl, .keep_all = TRUE) %>%
    dplyr::arrange(target_symbol_input)
  
  if (!is.null(max_targets)) {
    background_candidates <- background_candidates %>%
      dplyr::slice_head(n = max_targets)
    cat("      Limiting targets to first", nrow(background_candidates), "for this run\n")
  }
  
  cat("      Significant dog hallmarks:", nrow(dog_sig), "\n")
  cat("      Leading-edge ENSG candidates:", length(all_le_ensg), "\n")
  cat("      Mapped background candidates:", nrow(background_candidates), "\n")
  
  if (nrow(background_candidates) == 0) {
    stop("No hallmark-supported background candidates found.", call. = FALSE)
  }
  
  # ---------------------------------------------------------------------------
  # 4. Open Targets helpers: current drugAndClinicalCandidates schema
  # ---------------------------------------------------------------------------
  cat("  [BG-L3] Preparing Open Targets queries...\n")
  
  OT_URL <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  # Short cache folder, matching corrected Step 5 and reducing Windows MAX_PATH risk.
  ot_cache_dir <- file.path(PATH_CONFIG$cache_dir, "OTv3")
  dir.create(ot_cache_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (!dir.exists(ot_cache_dir)) {
    stop(
      "Could not create Open Targets cache directory: ", ot_cache_dir,
      "\nThis is usually a Windows path-length issue. Move the project to a shorter path, e.g. C:/omics_q1/",
      call. = FALSE
    )
  }
  
  if (isTRUE(force_refresh_ot)) {
    cat("      force_refresh_ot=TRUE -> cached target files will be refreshed.\n")
  }
  
  # Parameters retained only to avoid breaking old calls.
  invisible(ot_page_size)
  invisible(ot_max_pages)
  
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
  
  get_target_cached <- function(ensembl_id, idx = NA_integer_, total_n = NA_integer_, force_refresh = FALSE) {
    dir.create(ot_cache_dir, showWarnings = FALSE, recursive = TRUE)
    f <- file.path(ot_cache_dir, paste0(ensembl_id, "_target.json"))
    meta_f <- file.path(ot_cache_dir, paste0(ensembl_id, "_meta.json"))
    
    if (file.exists(f) && !force_refresh) {
      if (!is.na(idx) && idx %% progress_every == 0) {
        cat("      Progress:", idx, "/", total_n, "(cached)\n")
      }
      j <- jsonlite::fromJSON(readr::read_file(f), simplifyVector = FALSE)
      return(list(ok = TRUE, data = j, from_cache = TRUE))
    }
    
    if (!is.na(idx) && idx %% progress_every == 0) {
      cat("      Progress:", idx, "/", total_n, "\n")
    }
    
    out <- ot_query_target(ensembl_id)
    
    if (isTRUE(out$ok)) {
      dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
      writeLines(
        jsonlite::toJSON(out$data, auto_unbox = TRUE, pretty = TRUE),
        con = f
      )
      
      meta_obj <- list(
        ensembl_id = ensembl_id,
        retrieved_at = as.character(Sys.time()),
        source = "Open Targets GraphQL",
        schema = "target { tractability, drugAndClinicalCandidates }"
      )
      writeLines(
        jsonlite::toJSON(meta_obj, auto_unbox = TRUE, pretty = TRUE),
        con = meta_f
      )
    }
    
    Sys.sleep(0.25)
    out
  }
  
  # ---------------------------------------------------------------------------
  # 5. Query Open Targets
  # ---------------------------------------------------------------------------
  cat("  [BG-L4] Querying Open Targets background universe...\n")
  
  ensg_vec <- background_candidates$target_ensembl
  total_n <- length(ensg_vec)
  
  target_raw <- vector("list", total_n)
  for (i in seq_along(ensg_vec)) {
    target_raw[[i]] <- get_target_cached(
      ensembl_id = ensg_vec[i],
      idx = i,
      total_n = total_n,
      force_refresh = force_refresh_ot
    )
  }
  names(target_raw) <- ensg_vec
  
  target_ok <- target_raw[vapply(target_raw, function(x) isTRUE(x$ok), logical(1))]
  target_fail <- target_raw[!vapply(target_raw, function(x) isTRUE(x$ok), logical(1))]
  
  ot_query_qc <- tibble::tibble(
    n_targets_queried = length(target_raw),
    n_targets_success = length(target_ok),
    n_targets_failed = length(target_fail),
    failure_pct = ifelse(
      length(target_raw) > 0,
      100 * length(target_fail) / length(target_raw),
      NA_real_
    )
  )
  
  ot_error_log <- tibble::tibble(
    target_ensembl = names(target_raw),
    ok = vapply(target_raw, function(x) isTRUE(x$ok), logical(1)),
    error_message = vapply(target_raw, extract_ot_error, character(1))
  )
  
  readr::write_csv(
    ot_query_qc,
    file.path(PATH_CONFIG$results_dir, "background_light_OT_query_QC.csv")
  )
  readr::write_csv(
    ot_error_log,
    file.path(PATH_CONFIG$results_dir, "background_light_OT_query_error_log.csv")
  )
  
  cat("      Success:", length(target_ok), "| Failed:", length(target_fail), "\n")
  
  if (length(target_fail) > 0) {
    writeLines(
      names(target_fail),
      file.path(PATH_CONFIG$results_dir, "background_light_OT_failed_targets.txt")
    )
  }
  
  if (length(target_ok) == 0) {
    stop(
      "Open Targets failed for all background targets. See results/background_light_OT_query_error_log.csv.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # 6. Extract tractability
  # ---------------------------------------------------------------------------
  cat("  [BG-L5] Extracting tractability annotations...\n")
  
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
    stop("No tractability annotations retrieved for light background universe.", call. = FALSE)
  }
  
  readr::write_csv(
    tractability_table,
    file.path(PATH_CONFIG$results_dir, "background_light_tractability_raw.csv")
  )
  
  background_tractability <- tractability_table %>%
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
  
  # ---------------------------------------------------------------------------
  # 7. Extract drug/clinical-candidate evidence summary
  # ---------------------------------------------------------------------------
  cat("  [BG-L6] Extracting background drug/candidate evidence summary...\n")
  
  extract_drug_rows <- function(x, ensembl_id) {
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
      
      phase_for_record <- suppressWarnings(max(c(row_stage_num, drug_stage_num), na.rm = TRUE))
      if (is.infinite(phase_for_record)) phase_for_record <- NA_real_
      
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
          drug_chembl = drug_id,
          drug_name = drug_name,
          drug_type = drug_type,
          drug_is_approved = drug_is_approved,
          drug_max_phase = drug_stage_num,
          drug_max_phase_raw = drug_stage_chr,
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
            drug_chembl = drug_id,
            drug_name = drug_name,
            drug_type = drug_type,
            drug_is_approved = drug_is_approved,
            drug_max_phase = drug_stage_num,
            drug_max_phase_raw = drug_stage_chr,
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
  
  background_drug_rows <- dplyr::bind_rows(
    Map(extract_drug_rows, target_ok, names(target_ok))
  ) %>%
    dplyr::distinct()
  
  if (is.null(background_drug_rows) || nrow(background_drug_rows) == 0) {
    background_drug_rows <- tibble::tibble(
      target_ensembl = character(),
      target_symbol = character(),
      target_name = character(),
      phase_evidence = numeric(),
      phase_evidence_raw = character(),
      drug_chembl = character(),
      drug_name = character(),
      drug_type = character(),
      drug_is_approved = logical(),
      drug_max_phase = numeric(),
      drug_max_phase_raw = character(),
      disease_id = character(),
      disease_name = character(),
      disease_from_source = character()
    )
  }
  
  readr::write_csv(
    background_drug_rows,
    file.path(PATH_CONFIG$results_dir, "comparative_background_drug_rows_light.csv")
  )
  
  if (nrow(background_drug_rows) == 0) {
    background_drug_summary <- tibble::tibble(
      target_ensembl = character(),
      target_symbol = character(),
      n_drugs = integer(),
      n_indications = integer(),
      n_oncology_drug_indications = integer(),
      n_breast_drug_indications = integer(),
      max_drug_phase = numeric(),
      any_approved = logical(),
      any_onco_drug = logical(),
      any_breast_drug = logical()
    )
  } else {
    background_drug_summary <- background_drug_rows %>%
      dplyr::mutate(
        phase_evidence = suppressWarnings(as.numeric(phase_evidence)),
        drug_max_phase = suppressWarnings(as.numeric(drug_max_phase)),
        phase_for_summary = dplyr::case_when(
          !is.na(phase_evidence) & !is.na(drug_max_phase) ~ pmax(phase_evidence, drug_max_phase),
          !is.na(phase_evidence) ~ phase_evidence,
          !is.na(drug_max_phase) ~ drug_max_phase,
          TRUE ~ NA_real_
        ),
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
      dplyr::group_by(target_ensembl, target_symbol) %>%
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
          all(is.na(phase_for_summary)),
          0,
          max(phase_for_summary, na.rm = TRUE)
        ),
        any_approved = any(drug_is_approved %in% TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        max_drug_phase = ifelse(is.infinite(max_drug_phase), 0, max_drug_phase),
        any_onco_drug = n_oncology_drug_indications > 0,
        any_breast_drug = n_breast_drug_indications > 0
      )
  }
  
  # ---------------------------------------------------------------------------
  # 8. Merge final background table
  # ---------------------------------------------------------------------------
  cat("  [BG-L7] Merging final background universe...\n")
  
  background_universe <- background_candidates %>%
    dplyr::left_join(background_tractability, by = "target_ensembl") %>%
    dplyr::left_join(
      background_drug_summary %>%
        dplyr::select(
          target_ensembl,
          target_symbol,
          n_drugs,
          n_indications,
          n_oncology_drug_indications,
          n_breast_drug_indications,
          max_drug_phase,
          any_approved,
          any_onco_drug,
          any_breast_drug
        ),
      by = "target_ensembl",
      suffix = c("_tract", "_drug")
    ) %>%
    dplyr::mutate(
      target_symbol = dplyr::coalesce(target_symbol_tract, target_symbol_drug, target_symbol_input),
      has_tractability_data = dplyr::coalesce(has_tractability_data, FALSE),
      any_sm = dplyr::coalesce(any_sm, FALSE),
      any_ab = dplyr::coalesce(any_ab, FALSE),
      any_pr = dplyr::coalesce(any_pr, FALSE),
      any_oc = dplyr::coalesce(any_oc, FALSE),
      any_druggable = dplyr::coalesce(any_druggable, FALSE),
      protac_like = dplyr::coalesce(protac_like, FALSE),
      n_drugs = dplyr::coalesce(n_drugs, 0L),
      n_indications = dplyr::coalesce(n_indications, 0L),
      n_oncology_drug_indications = dplyr::coalesce(n_oncology_drug_indications, 0L),
      n_breast_drug_indications = dplyr::coalesce(n_breast_drug_indications, 0L),
      max_drug_phase = dplyr::coalesce(max_drug_phase, 0),
      any_approved = dplyr::coalesce(any_approved, FALSE),
      any_onco_drug = dplyr::coalesce(any_onco_drug, FALSE),
      any_breast_drug = dplyr::coalesce(any_breast_drug, FALSE)
    ) %>%
    dplyr::select(
      target_symbol,
      target_ensembl,
      has_tractability_data,
      any_sm,
      any_ab,
      any_pr,
      any_oc,
      any_druggable,
      protac_like,
      n_drugs,
      n_indications,
      n_oncology_drug_indications,
      n_breast_drug_indications,
      max_drug_phase,
      any_approved,
      any_onco_drug,
      any_breast_drug
    ) %>%
    dplyr::filter(has_tractability_data) %>%
    dplyr::distinct(target_ensembl, .keep_all = TRUE) %>%
    dplyr::arrange(target_symbol)
  
  if (nrow(background_universe) == 0) {
    stop(
      "The final background universe is empty after filtering to has_tractability_data. ",
      "Check background_light_tractability_raw.csv and background_light_OT_query_QC.csv.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # 9. Save outputs
  # ---------------------------------------------------------------------------
  cat("  [BG-L8] Saving outputs...\n")
  
  readr::write_csv(
    background_universe,
    file.path(PATH_CONFIG$results_dir, "comparative_background_tractability_light.csv")
  )
  
  readr::write_csv(
    background_drug_summary,
    file.path(PATH_CONFIG$results_dir, "comparative_background_drug_summary_light.csv")
  )
  
  summary_lines <- c(
    "LIGHT BACKGROUND TRACTABILITY UNIVERSE SUMMARY",
    "",
    paste("Significant dog hallmarks:", nrow(dog_sig)),
    paste("Leading-edge ENSG candidates:", length(all_le_ensg)),
    paste("Mapped background candidates:", nrow(background_candidates)),
    paste("Successfully queried OT targets:", length(target_ok)),
    paste("Failed OT targets:", length(target_fail)),
    paste("Final tractability-assessed background size:", nrow(background_universe)),
    paste("SM-tractable in background:", sum(background_universe$any_sm %in% TRUE)),
    paste("Any druggable in background:", sum(background_universe$any_druggable %in% TRUE)),
    paste("Targets with n_drugs > 0:", sum(background_universe$n_drugs > 0, na.rm = TRUE)),
    paste("Any approved in background:", sum(background_universe$any_approved %in% TRUE)),
    paste("Any oncology-linked drug indication in background:", sum(background_universe$any_onco_drug %in% TRUE)),
    paste("Any breast/mammary-linked drug indication in background:", sum(background_universe$any_breast_drug %in% TRUE)),
    paste("Open Targets cache directory:", normalizePath(ot_cache_dir, winslash = "/", mustWork = FALSE)),
    paste("Open Targets schema:", "target { tractability, drugAndClinicalCandidates }")
  )
  
  writeLines(
    summary_lines,
    file.path(PATH_CONFIG$results_dir, "VALIDATION_background_light_summary.txt")
  )
  
  assign("background_tractability_light", background_tractability, envir = .GlobalEnv)
  assign("background_drug_rows_light", background_drug_rows, envir = .GlobalEnv)
  assign("background_drug_summary_light", background_drug_summary, envir = .GlobalEnv)
  assign("background_universe_light", background_universe, envir = .GlobalEnv)
  assign("background_light_ot_query_qc", ot_query_qc, envir = .GlobalEnv)
  assign("background_light_ot_error_log", ot_error_log, envir = .GlobalEnv)
  
  cat("\n")
  cat("  ✓ LIGHT BACKGROUND UNIVERSE COMPLETE\n\n")
  cat("  Output file:\n")
  cat("   ", file.path(PATH_CONFIG$results_dir, "comparative_background_tractability_light.csv"), "\n\n")
  
  invisible(list(
    background_universe = background_universe,
    background_tractability = background_tractability,
    background_drug_rows = background_drug_rows,
    background_drug_summary = background_drug_summary,
    ot_query_qc = ot_query_qc,
    ot_error_log = ot_error_log
  ))
}
