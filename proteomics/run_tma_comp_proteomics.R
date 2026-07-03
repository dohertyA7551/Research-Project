#!/usr/bin/env Rscript
# Proteomics limma/Kruskal by TMA composition group.
#   cd /work_space/files/proteomics && /work_space/envs/transcriptomics2/bin/Rscript run_tma_comp_proteomics.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})
source("tma_comp_proteomics_helpers.R")
run_tma_comp_proteomics_limma()
run_tma_comp_proteomics_g1_vs_g4()
