# Generate figure panels for metacluster vs PROGENy/xCell correlation research.
# Uses existing CSV outputs in analysis_output/metacluster_correlation/.
#   /work_space/envs/transcriptomics2/bin/Rscript run_metacluster_correlation_figures.R
# Optional: Rscript run_metacluster_correlation_figures.R --scatter  (slow; re-runs PROGENy)
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})
source("metacluster_omics_correlation_helpers.R")
source("metacluster_omics_plot_helpers.R")

include_scatter <- "--scatter" %in% commandArgs(trailingOnly = TRUE)
fig_dir <- run_metacluster_correlation_figures(include_patient_scatter = include_scatter)
cat("\nFigures directory:", fig_dir, "\nDone.\n")
