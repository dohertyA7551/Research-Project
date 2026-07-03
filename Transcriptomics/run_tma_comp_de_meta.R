#!/usr/bin/env Rscript
# Compile internal per-cohort limma/GSEA with Pbine meta (no raw expression pooling).
#   cd /work_space/files/Transcriptomics && /work_space/envs/transcriptomics2/bin/Rscript run_tma_comp_de_meta.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("tma_comp_transcriptomics_helpers.R")
run_tma_comp_de_meta_all()
