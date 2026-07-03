# Density plots: xCell score distributions for selected cell types across 4 cohorts.
#   /work_space/envs/transcriptomics2/bin/Rscript run_xcell_density_plots.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})
source("analysis_helpers.R")
source("xcell_combine_helpers.R")

p <- resolve_transcriptomics_paths()
output_dir <- file.path(p$base, "analysis_output", "combined", "density_plots")
scores_long <- load_all_xcell_scores_long()

export_xcell_cohort_density_plots(
  scores_long,
  output_dir = output_dir,
  n_cell_types = 12L,
  z_score_within_cell_type = FALSE
)
export_xcell_cohort_density_plots(
  scores_long,
  output_dir = output_dir,
  n_cell_types = 12L,
  z_score_within_cell_type = TRUE
)

cat("\nDone. See ", output_dir, "\n", sep = "")
