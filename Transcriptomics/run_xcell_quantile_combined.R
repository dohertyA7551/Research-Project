# Quantile-normalize harmonized xCell scores across cohorts, then combined High vs Low Wilcoxon.
#   /work_space/envs/transcriptomics2/bin/Rscript run_xcell_quantile_combined.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})
source("analysis_helpers.R")
source("xcell_combine_helpers.R")

out <- run_xcell_quantile_combined_analysis()

cat("\nTop combined hits (quantile-normalized):\n")
print(
  out$diff_combined %>%
    arrange(pval) %>%
    select(cell_type, log2FC, pval) %>%
    head(10)
)

cat("\nDone. Outputs in analysis_output/combined/Combined_QN_*\n")
