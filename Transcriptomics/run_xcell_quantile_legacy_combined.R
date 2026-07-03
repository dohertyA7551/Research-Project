# Quantile-normalize Stage 2 + Retrospective + Colossus only; Taxonomy stays separate.
#   /work_space/envs/transcriptomics2/bin/Rscript run_xcell_quantile_legacy_combined.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})
source("analysis_helpers.R")
source("xcell_combine_helpers.R")

out <- run_xcell_quantile_legacy_combined_analysis()

cat("\nTop legacy combined hits (QN, 3 cohorts):\n")
print(
  out$diff_combined %>%
    arrange(pval) %>%
    select(cell_type, log2FC, pval) %>%
    head(10)
)

cat("\nTaxonomy remains separate: analysis_output/taxonomy/\n")
cat("Done. Legacy combined outputs: analysis_output/legacy_combined/\n")
