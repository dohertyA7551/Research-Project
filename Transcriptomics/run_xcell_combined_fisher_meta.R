# Fisher meta-analysis: xCell2 High vs Low + G1 vs G4 across four cohorts.
#   cd /work_space/files/Transcriptomics
#   /work_space/envs/transcriptomics2/bin/Rscript run_xcell_combined_fisher_meta.R
#
# Step 1: hierarchical Fisher meta within each contrast (4 cohorts).
# Step 2: Fisher combine the two contrast-level meta p-values (8 cohort tests total).
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("analysis_helpers.R")
source("xcell_combine_helpers.R")
source("tma_g1_vs_g4_transcriptomics_helpers.R")
source("pbine_meta_helpers.R")

output_dir <- "analysis_output/meta"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- c("Stage2", "Retrospective", "Colossus", "Taxonomy")
comp_groups <- load_g1_g4_comp_groups()

hl_list <- list()
g1g4_list <- list()
hl_counts <- list()
g1_counts <- list()

for (cn in cohorts) {
  scores <- load_xcell_scores_long(harmonize_cell_types = TRUE, cohorts = cn)

  hl_diff <- wilcox_enrichment_diff(scores)
  hl_list[[cn]] <- hl_diff
  hl_counts[[cn]] <- c(
    Low = dplyr::n_distinct(scores$sample[scores$risk_grp == "Low"]),
    High = dplyr::n_distinct(scores$sample[scores$risk_grp == "High"])
  )
  write_csv(
    hl_diff %>% dplyr::select(-dplyr::any_of("padj")),
    file.path(output_dir, paste0(cn, "_xcell_high_vs_low_wilcoxon.csv"))
  )

  g1_scores <- scores %>%
    dplyr::left_join(comp_groups, by = "patient_id") %>%
    dplyr::filter(.data$tma_comp_grp %in% c(G1_GRP, G4_GRP))
  g1_diff <- wilcox_g1_vs_g4(
    g1_scores,
    feature_col = "cell_type",
    value_col = "enrichment",
    group_col = "tma_comp_grp"
  )
  g1g4_list[[cn]] <- g1_diff
  g1_counts[[cn]] <- c(
    G1 = dplyr::n_distinct(g1_scores$sample[g1_scores$tma_comp_grp == G1_GRP]),
    G4 = dplyr::n_distinct(g1_scores$sample[g1_scores$tma_comp_grp == G4_GRP])
  )
  write_csv(
    g1_diff %>% dplyr::select(-dplyr::any_of("padj")),
    file.path(output_dir, paste0(cn, "_xcell_g1_vs_g4_wilcoxon.csv"))
  )
}

out <- run_xcell_combined_hl_g1g4_fisher_meta(
  hl_diff_list = hl_list,
  g1g4_diff_list = g1g4_list,
  output_dir = output_dir,
  file_prefix = "xcell_highlow_plus_g1g4",
  p_cutoff = 0.05,
  min_cohorts = 2L
)

fmt_counts <- function(x) paste0(x[1], " / ", x[2])
summary_lines <- c(
  "xCell2 Fisher meta: High vs Low DFS risk + G1 vs G4 TMA composition",
  "",
  paste0("Cohorts: ", paste(cohorts, collapse = ", ")),
  paste0(
    "High/Low samples (Low / High): ",
    paste(names(hl_counts), vapply(hl_counts, fmt_counts, character(1)), sep = "=", collapse = "; ")
  ),
  paste0(
    "G1/G4 samples (G1 / G4): ",
    paste(names(g1_counts), vapply(g1_counts, fmt_counts, character(1)), sep = "=", collapse = "; ")
  ),
  "",
  "Per-cohort Wilcoxon on harmonized xCell2 scores (same score source for both contrasts).",
  "Contrast 1: High vs Low (proteomics DFS risk).",
  "Contrast 2: G1 (all TMAs low) vs G4 (>=3 TMAs high).",
  "",
  "Meta-analysis:",
  "  (a) Hierarchical Fisher across 4 cohorts per contrast",
  "      (Stage2+Retrospective, Colossus+Taxonomy, then combined).",
  "  (b) Fisher combine the two contrast-level meta p-values -> combined meta.",
  "",
  "Note: High/Low and G1/G4 are related but not identical groupings;",
  "within-cohort tests are correlated. Combined meta treats contrasts as",
  "complementary evidence, not fully independent replicates.",
  "",
  "No multiple-testing correction applied (raw meta p-values used).",
  "",
  paste0("Combined cell types tested: ", nrow(out$combined)),
  paste0("Combined meta p < 0.05: ", sum(out$combined$meta_p < 0.05, na.rm = TRUE)),
  paste0("High/Low-only meta p < 0.05: ", sum(out$high_vs_low$meta_p < 0.05, na.rm = TRUE)),
  paste0("G1/G4-only meta p < 0.05: ", sum(out$g1_vs_g4$meta_p < 0.05, na.rm = TRUE)),
  "",
  "Outputs:",
  "  {cohort}_xcell_high_vs_low_wilcoxon.csv",
  "  {cohort}_xcell_g1_vs_g4_wilcoxon.csv",
  "  xcell_high_vs_low_meta_fisher.csv + plots",
  "  xcell_g1_vs_g4_meta_fisher.csv + plots",
  "  xcell_highlow_plus_g1g4_meta_fisher.csv + plots"
)
writeLines(summary_lines, file.path(output_dir, "xcell_combined_fisher_meta_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

cat("\nTop combined meta hits:\n")
print(
  out$combined %>%
    dplyr::select(
      cell_type, mean_log2FC, meta_p,
      meta_p_high_vs_low, meta_p_g1_vs_g4
    ) %>%
    head(10)
)

cat("\nDone.\n")
