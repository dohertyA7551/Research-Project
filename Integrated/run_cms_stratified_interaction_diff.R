#!/usr/bin/env Rscript
# Differential protein co-abundance interactions: Low vs High risk within CMS2 and CMS4.
#
#   cd /work_space/files/Integrated
#   /work_space/envs/transcriptomics2/bin/Rscript run_cms_stratified_interaction_diff.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("protein_risk_interaction_helpers.R")

export_all_cms_compartment_risk_interactions(
  cms_groups = CMS_STRATIFIED_GROUPS,
  min_low = 6L,
  min_high = 6L
)

cat("\nDone. Outputs under results/protein_risk_interaction_diff_cms/\n")
