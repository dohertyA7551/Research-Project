


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(limma)
})
source("analysis_helpers.R")

resolve_paths <- function() {
  if (file.exists("stage2/normalized_genes.csv")) {
    list(base = ".", proteomics = "../proteomics", spatial = "../spatial", integrated = "../Integrated ")
  } else if (file.exists("Transcriptomics/stage2/normalized_genes.csv")) {
    list(base = "Transcriptomics", proteomics = "proteomics", spatial = "spatial", integrated = "Integrated ")
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
expr_mat_stage2 <- expr_mat; sample_meta_stage2 <- sample_meta  # <- ADDED

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
expr_mat_retro <- expr_mat; sample_meta_retro <- sample_meta  # <- ADDED

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
expr_mat_colossus <- expr_mat; sample_meta_colossus <- sample_meta  # <- ADDED

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

## ---- Legacy (pooled: Stage2 + Retrospective + Colossus, quantile-normalized) ----
common_genes <- Reduce(intersect, list(rownames(expr_mat_stage2), rownames(expr_mat_retro), rownames(expr_mat_colossus)))
cat("\nCommon genes across Legacy cohorts:", length(common_genes), "\n")

legacy_expr <- cbind(
  expr_mat_stage2[common_genes, ],
  expr_mat_retro[common_genes, ],
  expr_mat_colossus[common_genes, ]
)
legacy_expr_qn <- normalizeQuantiles(as.matrix(legacy_expr))

legacy_meta <- bind_rows(
  sample_meta_stage2   %>% mutate(cohort = "Stage2"),
  sample_meta_retro    %>% mutate(cohort = "Retrospective"),
  sample_meta_colossus %>% mutate(cohort = "Colossus")
) %>%
  filter(sample_id %in% colnames(legacy_expr_qn)) %>%
  distinct(sample_id, .keep_all = TRUE)

legacy_meta <- legacy_meta[match(colnames(legacy_expr_qn), legacy_meta$sample_id), ]
stopifnot(identical(legacy_meta$sample_id, colnames(legacy_expr_qn)))

cat("\n=== Legacy (pooled, all samples before NA-filtering) ===\n")
print(table(legacy_meta$risk_grp, legacy_meta$cohort, useNA = "always"))

# ---- NA-filter: keep only patients with a valid RSF risk_grp classification ----
# model.matrix() silently drops NA rows, which desyncs design from legacy_expr_qn
# (271 samples) and throws "row dimension of design doesn't match column dimension
# of data object" in lmFit(). Subset both metadata and expression matrix first.
legacy_meta_valid <- legacy_meta %>% filter(!is.na(risk_grp))
legacy_expr_qn_valid <- legacy_expr_qn[, legacy_meta_valid$sample_id]
stopifnot(identical(legacy_meta_valid$sample_id, colnames(legacy_expr_qn_valid)))

cat("\n=== Legacy (pooled, NA-filtered, n =", nrow(legacy_meta_valid), ") ===\n")
print(table(legacy_meta_valid$risk_grp, legacy_meta_valid$cohort, useNA = "always"))

legacy_meta_valid$risk_grp <- factor(legacy_meta_valid$risk_grp, levels = c("Low", "High"))
legacy_meta_valid$cohort <- factor(legacy_meta_valid$cohort)

design <- model.matrix(~ cohort + risk_grp, data = legacy_meta_valid)
fit <- lmFit(legacy_expr_qn_valid, design)
fit <- eBayes(fit)

legacy_limma_results <- topTable(fit, coef = "risk_grpHigh", number = Inf, sort.by = "P")
legacy_limma_results$gene <- rownames(legacy_limma_results)

dirs_legacy <- setup_analysis_dirs(p$base, "legacy_combined")
write_csv(legacy_limma_results, file.path(dirs_legacy$output_dir, "Legacy_QN_limma_high_vs_low_all_genes.csv"))

cat("\n=== Legacy limma top results ===\n")
print(head(legacy_limma_results))

# ---- Sanity checks ----
ncol(legacy_expr_qn_valid)          # should equal nrow(legacy_meta_valid), e.g. 89
nrow(legacy_meta_valid)             # should equal ncol(legacy_expr_qn_valid)
nrow(design)                        # should equal nrow(legacy_meta_valid)
sum(is.na(legacy_meta_valid$risk_grp))  # should be 0
sum(is.na(legacy_meta_valid$cohort))    # should be 0

# ---- Robust save: Legacy limma results (CSV + RDS + Rmd report) ----

# Build an explicit, guaranteed-to-exist output directory
legacy_output_dir <- file.path(p$base, "legacy_combined", "output")
dir.create(legacy_output_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(legacy_output_dir))  # fail loudly if this didn't work

# Save full results as CSV
csv_path <- file.path(legacy_output_dir, "Legacy_QN_limma_high_vs_low_all_genes.csv")
write_csv(legacy_limma_results, csv_path)
stopifnot(file.exists(csv_path))
cat("Saved CSV to:", normalizePath(csv_path), "\n")

# Also save as RDS (preserves exact R object, useful if you reload for figures later)
rds_path <- file.path(legacy_output_dir, "Legacy_QN_limma_high_vs_low_all_genes.rds")
saveRDS(legacy_limma_results, rds_path)
cat("Saved RDS to:", normalizePath(rds_path), "\n")

# Save top 5 up / top 5 down as a compact summary CSV too
legacy_limma_results$abs_logFC <- abs(legacy_limma_results$logFC)
top5_up <- legacy_limma_results %>% filter(logFC > 0) %>% arrange(adj.P.Val, desc(abs_logFC)) %>% head(5)
top5_down <- legacy_limma_results %>% filter(logFC < 0) %>% arrange(adj.P.Val, desc(abs_logFC)) %>% head(5)
write_csv(bind_rows(top5_up, top5_down), file.path(legacy_output_dir, "Legacy_limma_top5_up_down.csv"))