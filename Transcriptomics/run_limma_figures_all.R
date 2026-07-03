# Generate logFC-focused volcano + top-gene heatmaps from limma DE results.
# Use: /work_space/envs/transcriptomics2/bin/Rscript run_limma_figures_all.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})
source("analysis_helpers.R")

resolve_paths <- function() {
  if (file.exists("stage2/normalized_genes.csv")) {
    list(
      base = ".",
      proteomics = "../proteomics",
      spatial = "../spatial",
      integrated = "../Integrated "
    )
  } else if (file.exists("Transcriptomics/stage2/normalized_genes.csv")) {
    list(
      base = "Transcriptomics",
      proteomics = "proteomics",
      spatial = "spatial",
      integrated = "Integrated "
    )
  } else {
    stop("Run from Transcriptomics/ or its parent directory.")
  }
}

run_cohort_figures <- function(
  cohort, expr_path, load_expr_fn, build_meta_fn, file_prefix, cohort_title,
  extra_copies = character(0),
  xmax_override = NULL,
  ymax_override = NULL
) {
  dirs <- setup_analysis_dirs(p$base, cohort)
  survival_orig <- file.path(p$proteomics, "survival_df_with_risk.rds")
  copies <- c(expression = expr_path, survival = survival_orig, extra_copies)
  w <- ensure_working_copies(copies, dirs$work_dir)
  survival_df <- readRDS(w["survival"])

  cat("\n=== ", cohort_title, " figures ===\n", sep = "")
  expr_mat <- filter_impute_expr(load_expr_fn(w["expression"]))
  sample_meta <- build_meta_fn(expr_mat, survival_df, w)

  results_path <- file.path(dirs$output_dir, paste0(file_prefix, "_limma_high_vs_low_all_genes.csv"))
  if (!file.exists(results_path)) {
    stop("Missing limma results: ", results_path, "\nRun run_limma_de_all.R first.")
  }
  results <- read_csv(results_path, show_col_types = FALSE)

  export_limma_figures(
    results = results,
    expr_mat = expr_mat,
    sample_meta = sample_meta,
    output_dir = dirs$output_dir,
    file_prefix = file_prefix,
    cohort_title = cohort_title,
    xmax_override = xmax_override,
    ymax_override = ymax_override
  )
}

p <- resolve_paths()

dirs_r <- setup_analysis_dirs(p$base, "retrospective")
dirs_c <- setup_analysis_dirs(p$base, "colossus")
results_r <- read_csv(
  file.path(dirs_r$output_dir, "Retrospective_limma_high_vs_low_all_genes.csv"),
  show_col_types = FALSE
)
results_c <- read_csv(
  file.path(dirs_c$output_dir, "Colossus_limma_high_vs_low_all_genes.csv"),
  show_col_types = FALSE
)
lim_r <- compute_limma_volcano_limits(results_r)
lim_c <- compute_limma_volcano_limits(results_c)
shared_xmax <- max(lim_r$xmax, lim_c$xmax)
shared_ymax <- max(lim_r$ymax, lim_c$ymax)
cat(
  "\nShared Retrospective/Colossus volcano axes: x = ±", round(shared_xmax, 2),
  ", y = 0–", round(shared_ymax, 2), "\n",
  sep = ""
)

run_cohort_figures(
  cohort = "stage2",
  expr_path = file.path(p$base, "stage2/normalized_genes.csv"),
  load_expr_fn = function(path) {
    prepare_expr_samples_x_genes(read_csv(path, show_col_types = FALSE), "patient_id")
  },
  build_meta_fn = function(expr_mat, survival_df, w) {
    clinical <- read_csv(
      ensure_working_copy(
        file.path(p$base, "stage2/clinical.csv"),
        file.path(p$base, "analysis_working/stage2")
      ),
      show_col_types = FALSE
    )
    mapping <- read_csv(
      ensure_working_copy(
        file.path(p$spatial, "master_patient_mapping.csv"),
        file.path(p$base, "analysis_working/stage2")
      ),
      show_col_types = FALSE
    )
    sample_ids <- colnames(expr_mat)
    data.frame(
      sample_id = sample_ids,
      patient_id = mapping$patient_id_g[
        match(clinical$r_code[match(sample_ids, clinical$patient_id)], mapping$r_code)
      ],
      stringsAsFactors = FALSE
    ) %>%
      left_join(survival_df %>% select(patient_id, risk_grp, rsf_risk), by = "patient_id")
  },
  file_prefix = "Stage2",
  cohort_title = "Stage 2 Transcriptomics"
)

run_cohort_figures(
  cohort = "retrospective",
  expr_path = file.path(p$base, "retrospective/Leuven_rlog_values.txt"),
  load_expr_fn = function(path) {
    prepare_expr_genes_x_samples(read_tsv(path, show_col_types = FALSE), "Geneid")
  },
  build_meta_fn = function(expr_mat, survival_df, w) {
    sample_ids <- colnames(expr_mat)
    data.frame(sample_id = sample_ids, patient_id = sample_ids, stringsAsFactors = FALSE) %>%
      left_join(survival_df %>% select(patient_id, risk_grp, rsf_risk), by = "patient_id")
  },
  file_prefix = "Retrospective",
  cohort_title = "Retrospective Transcriptomics",
  xmax_override = shared_xmax,
  ymax_override = shared_ymax
)

run_cohort_figures(
  cohort = "colossus",
  expr_path = file.path(p$base, "colossus/rna.csv"),
  load_expr_fn = function(path) {
    prepare_expr_samples_x_genes(read_csv(path, show_col_types = FALSE), "patient_id")
  },
  build_meta_fn = function(expr_mat, survival_df, w) {
    master_df <- read_csv(w["master"], show_col_types = FALSE)
    build_colossus_id_map(master_df, colnames(expr_mat)) %>%
      rename(sample_id = sample) %>%
      left_join(survival_df %>% select(patient_id, risk_grp, rsf_risk), by = "patient_id")
  },
  file_prefix = "Colossus",
  cohort_title = "Colossus Transcriptomics",
  extra_copies = c(
    master = file.path(p$integrated, "All_patinet_IDS_concatonated_masterdoc.csv")
  ),
  xmax_override = shared_xmax,
  ymax_override = shared_ymax
)

cat("\nLimma figures complete.\n")
