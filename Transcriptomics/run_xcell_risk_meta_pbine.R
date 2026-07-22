# Pbine meta-analysis of xCell2 High vs Low DFS risk (per cohort Wilcoxon -> combined).
#   cd /work_space/files/Transcriptomics
#   /work_space/envs/transcriptomics2/bin/Rscript run_xcell_risk_meta_pbine.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("analysis_helpers.R")
source("xcell_combine_helpers.R")
source("pbine_meta_helpers.R")

output_dir <- "analysis_output/meta"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- c("Stage2", "Retrospective", "Colossus", "Taxonomy")
diff_list <- list()
sample_counts <- list()

for (cn in cohorts) {
  scores <- load_xcell_scores_long(harmonize_cell_types = TRUE, cohorts = cn)
  n_low <- dplyr::n_distinct(scores$sample[scores$risk_grp == "Low"])
  n_high <- dplyr::n_distinct(scores$sample[scores$risk_grp == "High"])
  sample_counts[[cn]] <- c(Low = n_low, High = n_high)

  diff <- wilcox_enrichment_diff(scores)
  diff_list[[cn]] <- diff
  out_path <- file.path(output_dir, paste0(cn, "_xcell_high_vs_low_wilcoxon.csv"))
  readr::write_csv(
    diff %>% dplyr::select(-dplyr::any_of("padj")),
    out_path
  )
  message("Wrote ", out_path)
}

meta <- run_xcell_risk_meta_pbine(
  diff_list,
  output_dir = output_dir,
  file_prefix = "xcell_high_vs_low",
  p_cutoff = 0.05,
  method = "Fisher",
  min_cohorts = 2L
)

summary_lines <- c(
  "xCell2 High vs Low DFS risk — Pbine meta-analysis",
  "",
  paste0("Cohorts: ", paste(cohorts, collapse = ", ")),
  paste0(
    "Samples per cohort (Low / High): ",
    paste(
      names(sample_counts),
      vapply(sample_counts, function(x) paste0(x["Low"], " / ", x["High"]), character(1)),
      sep = "=",
      collapse = "; "
    )
  ),
  "",
  "Per-cohort Wilcoxon on harmonized xCell2 cell types (High vs Low).",
  "No multiple-testing correction (raw p-values at meta and per-cohort level).",
  "Cross-cohort synthesis: hierarchical Fisher (Stage2+Retro, Colossus+Taxonomy).",
  "",
  paste0("Cell types meta-tested: ", nrow(meta)),
  paste0("Significant at meta p < 0.05: ", sum(meta$meta_p < 0.05, na.rm = TRUE)),
  paste0("Sign-consistent across cohorts: ", sum(meta$sign_consistent %in% TRUE, na.rm = TRUE)),
  "",
  "Outputs:",
  "  {cohort}_xcell_high_vs_low_wilcoxon.csv",
  "  xcell_high_vs_low_meta_pbine.csv",
  "  xcell_high_vs_low_meta_pbine_dotplot.png",
  "  xcell_high_vs_low_meta_pbine_volcano.png"
)
writeLines(summary_lines, file.path(output_dir, "xcell_high_vs_low_meta_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

if (nrow(meta) > 0) {
  cat("\nTop meta hits (by meta p):\n")
  print(
    meta %>%
      dplyr::select(cell_type, mean_log2FC, meta_p, sign_consistent) %>%
      head(10)
  )
}

cat("\nDone.\n")
