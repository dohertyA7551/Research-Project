# CMS-stratified transcriptomics: Low vs High risk within CMS2 and CMS4.
# Runs limma DE, GSEA (Hallmark + KEGG), PROGENy, and xCell for Taxonomy RNA-Seq.
#
#   cd /work_space/files/Transcriptomics
#   /work_space/envs/transcriptomics2/bin/Rscript run_cms_stratified_transcriptomics.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("cms_stratified_helpers.R")
source("progeny_helpers.R")
source("xcell_combine_helpers.R")

resolve_paths <- function() {
  if (file.exists("Taxonomy_calls/Taxonomy_manuela_with_rna_rlog.txt")) {
    list(
      base = ".",
      proteomics = "../proteomics",
      output_root = file.path("analysis_output", "cms_stratified")
    )
  } else if (file.exists("Transcriptomics/Taxonomy_calls/Taxonomy_manuela_with_rna_rlog.txt")) {
    list(
      base = "Transcriptomics",
      proteomics = "proteomics",
      output_root = file.path("Transcriptomics", "analysis_output", "cms_stratified")
    )
  } else {
    stop("Run from Transcriptomics/ or its parent directory.")
  }
}

p <- resolve_paths()
survival_path <- file.path(p$proteomics, "survival_df_with_risk.rds")
map_path <- file.path(
  p$base,
  "Taxonomy_calls",
  "Manuela_and_Belfast_RNA_classifications_CMS_CRIS.txt"
)
rlog_path <- file.path(p$base, "Taxonomy_calls/Taxonomy_manuela_with_rna_rlog.txt")
stopifnot(file.exists(survival_path), file.exists(map_path), file.exists(rlog_path))

cat("Loading Taxonomy expression and metadata...\n")
expr_mat <- load_taxonomy_expression_matrix(rlog_path)
sample_meta <- build_taxonomy_sample_meta_with_cms(map_path, survival_path) %>%
  dplyr::filter(.data$sample_id %in% colnames(expr_mat))

cat("\nTaxonomy CMS x risk (all samples with expression):\n")
print(table(sample_meta$cms_subtype, sample_meta$risk_grp, useNA = "ifany"))

survival_df <- readRDS(survival_path)
xcell_scores <- load_xcell_scores_taxonomy(
  list(base = p$base, proteomics = p$proteomics, spatial = NA, integrated = NA),
  survival_df
)

counts_list <- list()
results <- list()

for (cms_group in CMS_STRATIFIED_GROUPS) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("CMS subgroup:", cms_group, "\n")
  cat(strrep("=", 60), "\n")

  out_dir <- cms_output_dir(p$output_root, cms_group, "taxonomy")
  file_prefix <- paste0("Taxonomy_", cms_group)
  title <- paste0("Taxonomy RNA-Seq — ", cms_group)

  cms_counts <- count_cms_risk(sample_meta, cms_group)
  counts_list[[cms_group]] <- c(cms_counts, modality = "Taxonomy_RNA")
  cat(
    "Samples: Low=", cms_counts["n_low"],
    ", High=", cms_counts["n_high"], "\n",
    sep = ""
  )

  limma_gsea <- run_cms_taxonomy_limma_gsea(
    cms_group = cms_group,
    expr_mat = expr_mat,
    sample_meta = sample_meta,
    output_dir = out_dir,
    file_prefix = file_prefix,
    cohort_title = title
  )
  progeny_out <- run_cms_taxonomy_progeny(
    cms_group = cms_group,
    expr_mat = expr_mat,
    sample_meta = sample_meta,
    output_dir = out_dir,
    file_prefix = file_prefix,
    title = title
  )
  xcell_out <- run_cms_taxonomy_xcell(
    cms_group = cms_group,
    scores_long = xcell_scores,
    output_dir = out_dir,
    file_prefix = file_prefix,
    title = title
  )

  results[[cms_group]] <- list(
    limma = limma_gsea$limma,
    gsea = limma_gsea$gsea,
    progeny = progeny_out,
    xcell = xcell_out
  )
}

write_cms_counts_summary(
  counts_list,
  file.path(p$output_root, "taxonomy_sample_counts.csv")
)

cat("\nCMS-stratified Taxonomy transcriptomics complete.\n")
cat("Outputs: ", p$output_root, "/{CMS2,CMS4}/taxonomy/\n", sep = "")
