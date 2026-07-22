# Anderson-Darling mergeability checks for Hallmark GSEA (High vs Low) across cohorts.
#   /work_space/envs/transcriptomics2/bin/Rscript run_gsea_hallmark_ad_combined.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("gsea_helpers.R")
source("gsea_merge_helpers.R")

out <- run_gsea_hallmark_ad_and_merge_check(ad_alpha = 0.05, min_genes = 15L)

cat("\nHallmark merge check complete.\n")
cat("Pathways tested:", out$n_tested, "\n")
cat("Pass (raw logFC AD):", out$n_pass_raw, "\n")
cat("Pass (z-scored logFC AD):", out$n_pass_z, "\n")
cat("NES profile AD p:", signif(out$nes_profile_ad$ad_p, 4), "\n")
cat("\nSee:", file.path(out$output_dir, "hallmark_merge_summary.txt"), "\n")

if (out$n_pass_raw == 0L) {
  cat("\nNo pathways pass raw logFC AD — review distribution_checks/ before combining.\n")
} else {
  cat("\nSome pathways pass AD — see hallmark_ad_by_pathway.csv for which.\n")
}

cat("\nDone.\n")
