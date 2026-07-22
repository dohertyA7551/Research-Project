# STRINGdb protein-protein interaction networks from limma High vs Low DE genes.
# Requires limma CSVs (run_limma_de_all.R) and STRINGdb (BiocManager::install("STRINGdb")).
#
#   cd /work_space/files/Transcriptomics
#   /work_space/envs/transcriptomics2/bin/Rscript run_string_ppi_high_vs_low.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("string_ppi_helpers.R")

out <- run_string_high_vs_low_all(
  fdr_cutoff = 0.05,
  fc_cutoff = 0.5,
  max_genes = 200L,
  score_threshold = 400L,
  string_version = "12.0",
  p_use = "auto"
)

cat("\nDone. Outputs under analysis_output/string_ppi/\n")
