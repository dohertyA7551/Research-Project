# Per-cell-type Anderson-Darling k-sample tests on xCell2 scores across cohorts.
# Compares each cell type's score distribution between datasets (not High vs Low).
#   /work_space/envs/transcriptomics2/bin/Rscript run_xcell_ad_combined.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})
source("analysis_helpers.R")
source("xcell_combine_helpers.R")

out <- run_xcell_ad_and_combined(ad_alpha = 0.05)

if (!out$can_combine) {
  cat("\nReview analysis_output/combined/xcell_ad_by_cell_type.csv\n")
} else {
  cat("\nPassing cell types written to xcell_ad_passing_cell_types.csv\n")
}

cat("\nDone.\n")
