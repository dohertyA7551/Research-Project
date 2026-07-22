# Four-panel figure: MC_13 vs MC_8 JAK-STAT regulons + STRING networks.
#   cd /work_space/files/Transcriptomics
#   /work_space/envs/transcriptomics2/bin/Rscript run_mc_jakstat_figure_panel.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(cowplot)
})
source("metacluster_jakstat_helpers.R")

out <- export_mc_jakstat_figure_panel()
cat("\nFigure panel written to:\n  ", out$png, "\n  ", out$pdf, "\n", sep = "")
