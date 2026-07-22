#!/usr/bin/env Rscript
# Protein interactions that differ between Low and High RSF risk.
#
#   cd /work_space/files/Integrated
#   /work_space/envs/transcriptomics2/bin/Rscript run_protein_risk_interaction_diff.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("protein_risk_interaction_helpers.R")

export_all_compartment_risk_interactions()

cat("\nDone. Outputs under results/protein_risk_interaction_diff/\n")
