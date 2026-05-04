################################################################################
# STEP 3: IDENTIFY CONSERVED HALLMARKS
# Logic:
#   - dog significant UP hallmarks ∩ cat significant UP hallmarks = conserved UP
#   - dog significant DOWN hallmarks ∩ cat significant DOWN hallmarks = conserved DOWN
#
# Uses:
#   - fg        : Dog Hallmark fgsea results
#   - up_tab_u  : Cat Hallmark ORA UP results (fixed universe)
#   - down_tab_u: Cat Hallmark ORA DOWN results (fixed universe)
#
# Outputs:
#   - cons_up (vector)
#   - cons_down (vector)
#   - CONSERVED_HALLMARKS_summary.csv
################################################################################

module_conserved_hallmarks <- function() {
  
  cat("  [3A] Checking required upstream objects...\n")
  
  required_objs <- c("fg", "up_tab_u", "down_tab_u")
  missing_objs <- required_objs[!vapply(required_objs, exists, logical(1), envir = .GlobalEnv)]
  
  if (length(missing_objs) > 0) {
    stop(
      "Missing required object(s) for conserved hallmarks step: ",
      paste(missing_objs, collapse = ", "),
      ". Run previous modules first."
    )
  }
  
  fg <- get("fg", envir = .GlobalEnv)
  up_tab_u <- get("up_tab_u", envir = .GlobalEnv)
  down_tab_u <- get("down_tab_u", envir = .GlobalEnv)
  
  # ---------------------------------------------------------------------------
  # 3B. DOG SIGNIFICANT HALLMARKS BY DIRECTION
  # ---------------------------------------------------------------------------
  cat("  [3B] Filtering significant dog hallmarks by direction...\n")
  
  dog_sig <- fg %>%
    dplyr::filter(!is.na(padj), padj < PARAM_CONFIG$conservation_padj)
  
  dog_sig_up <- dog_sig$pathway[dog_sig$NES > 0]
  dog_sig_down <- dog_sig$pathway[dog_sig$NES < 0]
  
  cat("      Dog significant UP hallmarks:", length(dog_sig_up), "\n")
  cat("      Dog significant DOWN hallmarks:", length(dog_sig_down), "\n")
  
  # ---------------------------------------------------------------------------
  # 3C. CAT SIGNIFICANT HALLMARKS BY DIRECTION (fixed universe ORA)
  # ---------------------------------------------------------------------------
  cat("  [3C] Filtering significant cat hallmarks by direction...\n")
  
  if (nrow(up_tab_u) > 0) {
    cat_sig_up <- up_tab_u$ID[!is.na(up_tab_u$p.adjust) & up_tab_u$p.adjust < PARAM_CONFIG$cat_ora_padj]
  } else {
    cat_sig_up <- character(0)
  }
  
  if (nrow(down_tab_u) > 0) {
    cat_sig_down <- down_tab_u$ID[!is.na(down_tab_u$p.adjust) & down_tab_u$p.adjust < PARAM_CONFIG$cat_ora_padj]
  } else {
    cat_sig_down <- character(0)
  }
  
  cat("      Cat significant UP hallmarks:", length(cat_sig_up), "\n")
  cat("      Cat significant DOWN hallmarks:", length(cat_sig_down), "\n")
  
  # ---------------------------------------------------------------------------
  # 3D. CONSERVED HALLMARKS = significant in BOTH species, same direction
  # ---------------------------------------------------------------------------
  cat("  [3D] Finding conserved hallmarks (intersection in same direction)...\n")
  
  cons_up <- intersect(dog_sig_up, cat_sig_up)
  cons_down <- intersect(dog_sig_down, cat_sig_down)
  
  cat("      Conserved UP hallmarks:", length(cons_up), "\n")
  cat("      Conserved DOWN hallmarks:", length(cons_down), "\n")
  
  if (length(cons_up) > 0) {
    cat("      UP example(s):", paste(cons_up[seq_len(min(5, length(cons_up)))], collapse = ", "), "\n")
  }
  if (length(cons_down) > 0) {
    cat("      DOWN example(s):", paste(cons_down[seq_len(min(5, length(cons_down)))], collapse = ", "), "\n")
  }
  
  # ---------------------------------------------------------------------------
  # 3E. SUMMARY TABLE
  # ---------------------------------------------------------------------------
  cat("  [3E] Writing conserved hallmark summary table...\n")
  
  hall_summary <- dplyr::bind_rows(
    data.frame(
      hallmark = cons_up,
      direction = "UP",
      stringsAsFactors = FALSE
    ),
    data.frame(
      hallmark = cons_down,
      direction = "DOWN",
      stringsAsFactors = FALSE
    )
  )
  
  if (nrow(hall_summary) == 0) {
    hall_summary <- data.frame(
      hallmark = character(),
      direction = character(),
      stringsAsFactors = FALSE
    )
  }
  
  save_result(hall_summary, "CONSERVED_HALLMARKS_summary.csv", "csv")
  
  # ---------------------------------------------------------------------------
  # 3F. EXPORT TO GLOBAL ENV + RETURN
  # ---------------------------------------------------------------------------
  assign("dog_sig", dog_sig, envir = .GlobalEnv)
  assign("dog_sig_up", dog_sig_up, envir = .GlobalEnv)
  assign("dog_sig_down", dog_sig_down, envir = .GlobalEnv)
  
  assign("cat_sig_up", cat_sig_up, envir = .GlobalEnv)
  assign("cat_sig_down", cat_sig_down, envir = .GlobalEnv)
  
  assign("cons_up", cons_up, envir = .GlobalEnv)
  assign("cons_down", cons_down, envir = .GlobalEnv)
  
  cat("  ✓ STEP 3 COMPLETE\n\n")
  
  invisible(list(
    dog_sig = dog_sig,
    dog_sig_up = dog_sig_up,
    dog_sig_down = dog_sig_down,
    cat_sig_up = cat_sig_up,
    cat_sig_down = cat_sig_down,
    cons_up = cons_up,
    cons_down = cons_down,
    hall_summary = hall_summary
  ))
}