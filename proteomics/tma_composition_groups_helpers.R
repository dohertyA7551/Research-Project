# TMA-level RSF risk from metacluster abundances → patient composition groups.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

MC_PREFIX <- "abundance_metacluster_leiden_1.0_"

resolve_tma_comp_paths <- function() {
  list(
    proteomics = "/work_space/files/proteomics",
    features = "/work_space/files/Features/Run2_all/Manual_classes",
    spatial = "/work_space/files/spatial",
    transcriptomics = "/work_space/files/Transcriptomics",
    output = "/work_space/files/proteomics/results/tma_composition_groups"
  )
}

load_rsf_and_cutpoint <- function(p) {
  if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
    stop("Package 'randomForestSRC' is required.")
  }
  if (!requireNamespace("survminer", quietly = TRUE)) {
    stop("Package 'survminer' is required for the RSF cutpoint.")
  }
  rsf <- readRDS(file.path(p$proteomics, "rsf_fit_v2.rds"))
  surv <- readRDS(file.path(p$proteomics, "survival_df_with_risk.rds"))
  mc_vars <- names(sort(rsf$importance, decreasing = TRUE))
  cut_res <- survminer::surv_cutpoint(
    surv,
    time = "dfs_months",
    event = "dfs_event",
    variables = "rsf_risk",
    minprop = 0.15
  )
  cut_val <- cut_res$cutpoint$cutpoint
  list(rsf = rsf, surv = surv, mc_vars = mc_vars, cutpoint = cut_val)
}

load_tma_master_features <- function(p) {
  mf <- read.csv(
    file.path(p$features, "master_features.csv"),
    check.names = FALSE
  )
  clin_a <- read.csv(file.path(p$spatial, "ColonRO1_2020_RCSI_clinical.csv"))
  clin_b <- read.csv(file.path(p$spatial, "ColonRO1_2020_Taxonomy_clinical.csv"))
  clin_all <- bind_rows(clin_a, clin_b) %>%
    distinct(SS, Patient) %>%
    rename(tma_id = SS, patient_id = Patient)
  mf %>%
    left_join(clin_all, by = "tma_id")
}

score_tmas_with_rsf <- function(tma_df, rsf, mc_vars, cutpoint) {
  miss <- setdiff(mc_vars, names(tma_df))
  if (length(miss) > 0L) {
    stop("Missing RSF metacluster columns in master_features: ", paste(miss, collapse = ", "))
  }
  use <- tma_df %>%
    filter(if_all(all_of(mc_vars), ~ !is.na(.x)))
  if (nrow(use) == 0L) {
    stop("No TMA rows with complete RSF metacluster abundances.")
  }
  pred <- stats::predict(rsf, newdata = use[, mc_vars, drop = FALSE])$predicted
  use$rsf_tma_score <- pred
  use$tma_risk <- ifelse(use$rsf_tma_score >= cutpoint, "High", "Low")
  use
}

assign_patient_composition_group <- function(n_low, n_high, n_tma) {
  if (n_high >= 3L) {
    return("G4_high_3plus")
  }
  if (n_high == 0L) {
    return("G1_all_low")
  }
  if (n_high >= 2L) {
    return("G3_high_2plus")
  }
  if (n_low >= 2L) {
    return("G2_low_2plus")
  }
  "Other_mixed"
}

build_patient_composition_groups <- function(tma_scored) {
  tma_scored %>%
    filter(!is.na(.data$patient_id)) %>%
    group_by(.data$patient_id) %>%
    summarise(
      n_tma = n(),
      n_tma_low = sum(.data$tma_risk == "Low"),
      n_tma_high = sum(.data$tma_risk == "High"),
      mean_rsf_tma = mean(.data$rsf_tma_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      tma_comp_grp = assign_patient_composition_group(
        n_tma_low, n_tma_high, n_tma
      )
    ) %>%
    ungroup() %>%
    mutate(
      tma_comp_grp = factor(
        tma_comp_grp,
        levels = c(
          "G1_all_low",
          "G2_low_2plus",
          "G3_high_2plus",
          "G4_high_3plus",
          "Other_mixed"
        )
      )
    )
}

mc_rsf_direction_table <- function(surv, mc_vars) {
  bind_rows(lapply(mc_vars, function(col) {
    x <- surv[[col]]
    ct <- stats::cor.test(x, surv$rsf_risk, method = "spearman", exact = FALSE)
    tibble(
      metacluster = sub(paste0("^", MC_PREFIX), "", col),
      mc_column = col,
      spearman_rho_rsf = unname(ct$estimate),
      pval = ct$p.value,
      mean_in_low_risk = mean(x[surv$risk_grp == "Low"], na.rm = TRUE),
      mean_in_high_risk = mean(x[surv$risk_grp == "High"], na.rm = TRUE)
    )
  })) %>%
    arrange(desc(abs(.data$spearman_rho_rsf)))
}

run_tma_composition_group_assignment <- function(output_dir = NULL) {
  p <- resolve_tma_comp_paths()
  if (is.null(output_dir)) {
    output_dir <- p$output
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  res <- load_rsf_and_cutpoint(p)
  tma_raw <- load_tma_master_features(p)
  tma_scored <- score_tmas_with_rsf(tma_raw, res$rsf, res$mc_vars, res$cutpoint)
  patient_grp <- build_patient_composition_groups(tma_scored)

  patient_full <- patient_grp %>%
    left_join(
      res$surv %>%
        select(
          patient_id, dfs_months, dfs_event, rsf_risk, risk_grp,
          any_of(res$mc_vars)
        ),
      by = "patient_id"
    )

  mc_dir <- mc_rsf_direction_table(res$surv, res$mc_vars)

  readr::write_csv(
    tma_scored %>%
      select(tma_id, patient_id, rsf_tma_score, tma_risk, all_of(res$mc_vars)),
    file.path(output_dir, "tma_rsf_scores.csv")
  )
  readr::write_csv(patient_full, file.path(output_dir, "patient_tma_composition_groups.csv"))
  readr::write_csv(mc_dir, file.path(output_dir, "rsf_metacluster_direction.csv"))

  saveRDS(
    list(
      cutpoint = res$cutpoint,
      mc_vars = res$mc_vars,
      patient_groups = patient_full,
      tma_scores = tma_scored,
      mc_direction = mc_dir
    ),
    file.path(p$proteomics, "tma_composition_patient_groups.rds")
  )

  summary_lines <- c(
    "TMA composition groups (RSF on TMA metacluster abundances)",
    paste0("RSF model: rsf_fit_v2.rds (", length(res$mc_vars), " metaclusters)"),
    paste0("Patient-level RSF cutpoint (same as survival_df): ", signif(res$cutpoint, 4)),
    paste0("TMAs scored: ", nrow(tma_scored)),
    "",
    "TMA risk (High/Low per spot):",
    paste(capture.output(print(table(tma_scored$tma_risk))), collapse = "\n"),
    "",
    "Patient composition groups (mutually exclusive rules):",
    "  G1_all_low     = all TMAs Low",
    "  G2_low_2plus   = >=2 Low TMAs (mixed, not all low)",
    "  G3_high_2plus  = >=2 High TMAs (and <3 High)",
    "  G4_high_3plus  = >=3 High TMAs",
    "  Other_mixed    = does not meet above (often 2 TMAs: 1 Low + 1 High)",
    "",
    paste(capture.output(print(table(patient_grp$tma_comp_grp))), collapse = "\n"),
    "",
    "RSF metacluster direction (Spearman vs patient rsf_risk):",
    paste(
      mc_dir$metacluster,
      "rho=",
      signif(mc_dir$spearman_rho_rsf, 3),
      collapse = "; "
    )
  )
  writeLines(summary_lines, file.path(output_dir, "assignment_summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(
    patient_groups = patient_full,
    tma_scores = tma_scored,
    mc_direction = mc_dir,
    cutpoint = res$cutpoint
  ))
}
