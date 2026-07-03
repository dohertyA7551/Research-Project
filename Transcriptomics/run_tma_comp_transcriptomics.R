#!/usr/bin/env Rscript
# Transcriptomics xCell / PROGENy / GSEA by TMA composition group.
#   cd /work_space/files/Transcriptomics && /work_space/envs/transcriptomics2/bin/Rscript run_tma_comp_transcriptomics.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})
source("tma_comp_transcriptomics_helpers.R")
run_tma_comp_transcriptomics_all()
