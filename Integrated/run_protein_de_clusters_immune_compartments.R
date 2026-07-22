#!/usr/bin/env Rscript
# Protein DE clusters across immune compartments (helperT, cytT, regT).
#
#   cd /work_space/files/Integrated
#   /work_space/envs/transcriptomics2/bin/Rscript run_protein_de_clusters_immune_compartments.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("epi_protein_cluster_helpers.R")

N_CLUSTERS <- 6L
TOP_N <- 10L

run_protein_de_clusters_immune_compartments(
  n_clusters = N_CLUSTERS,
  top_n = TOP_N,
  seed = 42L
)

cat("\nDone. Outputs under results/protein_de_clusters_immune/\n")
