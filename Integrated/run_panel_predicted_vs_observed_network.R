#!/usr/bin/env Rscript
# Predicted (bulk RNA + STRING PPI) vs observed (epithelial proteomics co-abundance)
# for the TMA CODEX marker panel, High vs Low RSF risk.
#
#   cd /work_space/files/Integrated
#   /work_space/envs/transcriptomics2/bin/Rscript run_panel_predicted_vs_observed_network.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("panel_predicted_observed_network_helpers.R")

# ---- Config ----
CORR_THRESHOLD <- 0.3
STRING_SCORE <- 400L
HUB_TOP_N <- 10L
SEED <- 42L

out <- run_panel_predicted_vs_observed(
  corr_threshold = CORR_THRESHOLD,
  string_score = STRING_SCORE,
  hub_top_n = HUB_TOP_N,
  seed = SEED
)

cat("\nDone. Outputs under results/panel_predicted_vs_observed/\n")
