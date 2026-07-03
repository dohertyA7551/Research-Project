# Class-specific (risk group) normalization + PROGENy for legacy cohorts and Taxonomy.
#   /work_space/envs/transcriptomics2/bin/Rscript run_progeny_legacy_taxonomy.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
})
source("progeny_helpers.R")

out <- run_progeny_legacy_and_taxonomy()

cat("\nLegacy top pathways:\n")
print(out$legacy$diff %>% arrange(pval) %>% select(pathway, log2FC, pval) %>% head(8))

cat("\nTaxonomy top pathways:\n")
print(out$taxonomy$diff %>% arrange(pval) %>% select(pathway, log2FC, pval) %>% head(8))

cat("\nDone.\n")
