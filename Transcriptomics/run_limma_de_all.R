# Run limma High vs Low DE for all four transcriptomic cohorts.
# Use the transcriptomics2 R on /work_space (has limma, dplyr, readr, etc.):
#   /work_space/envs/transcriptomics2/bin/Rscript run_limma_de_all.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
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

p <- resolve_paths()
survival_orig <- file.path(p$proteomics, "survival_df_with_risk.rds")
stopifnot(file.exists(survival_orig))

# ---- Stage 2 ----
dirs <- setup_analysis_dirs(p$base, "stage2")
w <- ensure_working_copies(c(
  expression = file.path(p$base, "stage2/normalized_genes.csv"),
  clinical = file.path(p$base, "stage2/clinical.csv"),
  mapping = file.path(p$spatial, "master_patient_mapping.csv"),
  survival = survival_orig
), dirs$work_dir)
survival_df <- readRDS(w["survival"])
clinical <- read_csv(w["clinical"], show_col_types = FALSE)
mapping <- read_csv(w["mapping"], show_col_types = FALSE)
expr_mat <- filter_impute_expr(
  prepare_expr_samples_x_genes(read_csv(w["expression"], show_col_types = FALSE), "patient_id")
)
sample_ids <- colnames(expr_mat)
sample_meta <- data.frame(
  sample_id = sample_ids,
  patient_id = mapping$patient_id_g[
    match(clinical$r_code[match(sample_ids, clinical$patient_id)], mapping$r_code)
  ],
  stringsAsFactors = FALSE
) %>%
  left_join(survival_df %>% select(patient_id, risk_grp, rsf_risk), by = "patient_id")
cat("\n=== Stage 2 ===\n")
print(table(sample_meta$risk_grp, useNA = "always"))
s2 <- run_limma_high_vs_low(expr_mat, sample_meta, dirs$output_dir, "Stage2")
print(s2$top)
export_limma_figures(s2$results, expr_mat, sample_meta, dirs$output_dir, "Stage2", "Stage 2 Transcriptomics")

# ---- Retrospective ----
dirs <- setup_analysis_dirs(p$base, "retrospective")
w <- ensure_working_copies(c(
  expression = file.path(p$base, "retrospective/Leuven_rlog_values.txt"),
  survival = survival_orig
), dirs$work_dir)
survival_df <- readRDS(w["survival"])
expr_mat <- filter_impute_expr(
  prepare_expr_genes_x_samples(read_tsv(w["expression"], show_col_types = FALSE), "Geneid")
)
sample_ids <- colnames(expr_mat)
sample_meta <- data.frame(sample_id = sample_ids, patient_id = sample_ids, stringsAsFactors = FALSE) %>%
  left_join(survival_df %>% select(patient_id, risk_grp, rsf_risk), by = "patient_id")
cat("\n=== Retrospective ===\n")
print(table(sample_meta$risk_grp, useNA = "always"))
r <- run_limma_high_vs_low(expr_mat, sample_meta, dirs$output_dir, "Retrospective")
print(r$top)
export_limma_figures(r$results, expr_mat, sample_meta, dirs$output_dir, "Retrospective", "Retrospective Transcriptomics")

# ---- Colossus ----
dirs <- setup_analysis_dirs(p$base, "colossus")
w <- ensure_working_copies(c(
  expression = file.path(p$base, "colossus/rna.csv"),
  master = file.path(p$integrated, "All_patinet_IDS_concatonated_masterdoc.csv"),
  survival = survival_orig
), dirs$work_dir)
survival_df <- readRDS(w["survival"])
master_df <- read_csv(w["master"], show_col_types = FALSE)
expr_mat <- filter_impute_expr(
  prepare_expr_samples_x_genes(read_csv(w["expression"], show_col_types = FALSE), "patient_id")
)
sample_meta <- build_colossus_id_map(master_df, colnames(expr_mat)) %>%
  rename(sample_id = sample) %>%
  left_join(survival_df %>% select(patient_id, risk_grp, rsf_risk), by = "patient_id")
cat("\n=== Colossus ===\n")
print(table(sample_meta$risk_grp, useNA = "always"))
c <- run_limma_high_vs_low(expr_mat, sample_meta, dirs$output_dir, "Colossus")
print(c$top)
export_limma_figures(c$results, expr_mat, sample_meta, dirs$output_dir, "Colossus", "Colossus Transcriptomics")

# ---- Taxonomy ----
dirs <- setup_analysis_dirs(p$base, "taxonomy")
w <- ensure_working_copies(c(
  expression = file.path(p$base, "Taxonomy_calls/Taxonomy_manuela_with_rna_rlog.txt"),
  map = file.path(p$base, "Taxonomy_calls/Manuela_and_Belfast_RNA_classifications_CMS_CRIS.txt"),
  survival = survival_orig
), dirs$work_dir)
cat("\n=== Taxonomy (loading rlog) ===\n")
expr_mat <- load_taxonomy_expression_matrix(w["expression"])
sample_meta <- build_taxonomy_sample_meta(w["map"], w["survival"]) %>%
  filter(sample_id %in% colnames(expr_mat))
print(table(sample_meta$risk_grp[match(colnames(expr_mat), sample_meta$sample_id)], useNA = "always"))
t <- run_limma_high_vs_low(expr_mat, sample_meta, dirs$output_dir, "Taxonomy")
print(t$top)
export_limma_figures(t$results, expr_mat, sample_meta, dirs$output_dir, "Taxonomy", "Taxonomy RNA-Seq")

cat("\nLimma DE complete.\n")
