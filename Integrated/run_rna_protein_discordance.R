#!/usr/bin/env Rscript
# Patient-level RNA -> protein prediction and discordance (observed - predicted).
# Reads data from /work_space/files/* ; writes only to Integrated/results/.
#
#   cd /work_space/files/Integrated
#   /work_space/envs/transcriptomics2/bin/Rscript run_rna_protein_discordance.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("rna_protein_discordance_helpers.R")
run_rna_protein_discordance_all()
