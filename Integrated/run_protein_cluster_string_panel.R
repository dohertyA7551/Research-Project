#!/usr/bin/env Rscript
# Combined STRING network figure: top up/down clusters (High vs Low + G1 vs G4).
#
#   cd /work_space/files/Integrated
#   /work_space/envs/transcriptomics2/bin/Rscript run_protein_cluster_string_panel.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("epi_protein_cluster_helpers.R")
source("protein_cluster_string_helpers.R")

export_all_compartment_cluster_string_panels()

cat("\nDone. Figures under */string_panels/\n")
