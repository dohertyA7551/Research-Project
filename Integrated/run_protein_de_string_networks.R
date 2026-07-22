#!/usr/bin/env Rscript
# DE protein STRING networks: High vs Low + G1 vs G4, all cell compartments.
#
#   cd /work_space/files/Integrated
#   /work_space/envs/transcriptomics2/bin/Rscript run_protein_de_string_networks.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("protein_de_string_helpers.R")

export_all_compartment_de_string_networks()

cat("\nDone. Outputs under results/protein_de_string_networks/\n")
