#!/usr/bin/env Rscript
# Epithelial proteomics: top up/down marker clusters for High vs Low and G1 vs G4.
#
#   cd /work_space/files/Integrated
#   /work_space/envs/transcriptomics2/bin/Rscript run_protein_de_clusters_both_contrasts.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("epi_protein_cluster_helpers.R")

N_CLUSTERS <- 6L
TOP_N <- 10L

run_protein_de_clusters_both_contrasts(
  n_clusters = N_CLUSTERS,
  top_n = TOP_N,
  seed = 42L
)

cat("\nDone. Outputs under results/protein_de_clusters/\n")
