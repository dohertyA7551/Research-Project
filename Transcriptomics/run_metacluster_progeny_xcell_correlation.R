# Metacluster abundance vs PROGENy / xCell (Spearman per cohort, Fisher meta).
#   /work_space/envs/transcriptomics2/bin/Rscript run_metacluster_progeny_xcell_correlation.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})
source("metacluster_omics_correlation_helpers.R")

out <- run_metacluster_progeny_xcell_correlation()

cat("\nPROGENy meta (top 10 by Fisher padj):\n")
if (nrow(out$progeny_meta) > 0L) {
  print(
    out$progeny_meta %>%
      arrange(fisher_padj) %>%
      select(metacluster, feature, mean_rho, n_cohorts, fisher_p, fisher_padj) %>%
      head(10)
  )
}

cat("\nxCell meta (top 10 by Fisher padj):\n")
if (nrow(out$xcell_meta) > 0L) {
  print(
    out$xcell_meta %>%
      arrange(fisher_padj) %>%
      select(metacluster, feature, mean_rho, n_cohorts, fisher_p, fisher_padj) %>%
      head(10)
  )
}

cat("\nDone.\n")
