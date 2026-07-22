# MC_13 (low risk) vs MC_8 (high risk): STRING marker networks + Hallmark regulon correlations.
#   cd /work_space/files/Transcriptomics
#   /work_space/envs/transcriptomics2/bin/Rscript run_mc_jakstat_string_regulon.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})
source("metacluster_jakstat_helpers.R")

out <- run_mc_jakstat_string_and_regulon()

cat("\nHallmark regulon Fisher meta (MC_13 & MC_8):\n")
if (nrow(out$regulon$meta) > 0L) {
  print(
    out$regulon$meta %>%
      arrange(fisher_padj) %>%
      select(metacluster, regulon_short, mean_rho, n_cohorts, fisher_p, fisher_padj) %>%
      head(20)
  )
}

cat("\nDone. Outputs under analysis_output/metacluster_correlation/jakstat_mc8_mc13/\n")
