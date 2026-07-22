#!/usr/bin/env Rscript
# Per-cohort xCell2 / PROGENy / KEGG GSEA (G1 vs G4), then distribution-checked merge + Pbine meta.
#   cd /work_space/files/Transcriptomics
#   /work_space/envs/transcriptomics2/bin/Rscript run_tma_g1_vs_g4_transcriptomics.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})
source("tma_g1_vs_g4_transcriptomics_helpers.R")
run_tma_g1_vs_g4_transcriptomics_all()
