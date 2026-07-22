#!/usr/bin/env Rscript
# Epithelial proteomics: cluster markers and identify top up/down clusters (High vs Low).
# For both High vs Low and G1 vs G4, use run_protein_de_clusters_both_contrasts.R
#
#   cd /work_space/files/Integrated
#   /work_space/envs/transcriptomics2/bin/Rscript run_epi_protein_de_clusters.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("epi_protein_cluster_helpers.R")

N_CLUSTERS <- 6L
TOP_N <- 10L

run_epi_protein_de_clusters(
  n_clusters = N_CLUSTERS,
  top_n = TOP_N,
  seed = 42L
)

cat("\nDone. Outputs under results/epi_protein_high_vs_low_clusters/\n")
