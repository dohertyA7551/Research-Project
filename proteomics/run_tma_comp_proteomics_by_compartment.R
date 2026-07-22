#!/usr/bin/env Rscript
# Proteomics all TMAs Low vs all TMAs High — separate limma per cell compartment.
#   cd /work_space/files/proteomics
#   /work_space/envs/transcriptomics2/bin/Rscript run_tma_comp_proteomics_by_compartment.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})
source("tma_comp_proteomics_helpers.R")
run_tma_comp_proteomics_all_low_vs_all_high_by_compartment()
