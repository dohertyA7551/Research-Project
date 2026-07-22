#!/usr/bin/env Rscript
# Generate figure panels for RNA-protein discordance research.
# Uses existing CSV outputs in results/rna_protein_discordance/.
#
#   cd /work_space/files/Integrated
#   /work_space/envs/transcriptomics2/bin/Rscript run_discordance_figures.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("rna_protein_discordance_helpers.R")

fig_dir <- run_discordance_figures()
cat("\nFigures directory:", fig_dir, "\nDone.\n")
