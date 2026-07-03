# Dataset-level Anderson-Darling k-sample tests on xCell scores (4 cohorts);
# pool and rerun volcano + top-5 if the gate test passes.
#   /work_space/envs/transcriptomics2/bin/Rscript run_xcell_ad_combined.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})
source("analysis_helpers.R")
source("xcell_combine_helpers.R")

out <- run_xcell_ad_and_combined(
  ad_alpha = 0.05,
  combine_test = "zscore_by_cell_type_all_samples"
)

if (!out$can_combine) {
  cat("\nReview analysis_output/combined/xcell_ad_k_sample_dataset.csv\n")
  cat("To use raw scores as gate: combine_test = 'raw_all_samples'\n")
}

cat("\nDone.\n")
