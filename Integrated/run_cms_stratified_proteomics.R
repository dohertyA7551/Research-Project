# CMS-stratified proteomics: Low vs High risk within CMS2 and CMS4.
# Runs patient-level limma on all five CODEX compartments.
#
#   cd /work_space/files/Integrated
#   /work_space/envs/transcriptomics2/bin/Rscript run_cms_stratified_proteomics.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("cms_stratified_proteomics_helpers.R")

OUTPUT_ROOT <- file.path("results", "cms_stratified")
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

cat("CMS-stratified proteomics: Low vs High within CMS2 and CMS4\n")
cat("Compartments:", paste(ALL_PROTEIN_COMPARTMENTS$suffix, collapse = ", "), "\n\n")

cms_counts <- list()
all_results <- list()

for (cms_group in CMS_STRATIFIED_GROUPS) {
  cat(strrep("=", 60), "\n")
  cat("CMS subgroup:", cms_group, "\n")
  cat(strrep("=", 60), "\n")

  out_dir <- file.path(OUTPUT_ROOT, cms_group)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  prot_epi <- load_proteomics_with_cms(compartment_suffix = "NonImmuneEpithelium")
  aligned <- tryCatch(
    subset_protein_cms_risk(prot_epi, cms_group = cms_group),
    error = function(e) {
      cat("Skipping ", cms_group, ": ", conditionMessage(e), "\n", sep = "")
      NULL
    }
  )
  if (is.null(aligned)) next

  cms_counts[[cms_group]] <- data.frame(
    cms = cms_group,
    n_low = aligned$n_low,
    n_high = aligned$n_high,
    n_total = aligned$n_low + aligned$n_high,
    stringsAsFactors = FALSE
  )
  cat("Patients: Low=", aligned$n_low, ", High=", aligned$n_high, "\n", sep = "")

  cms_results <- list()
  for (i in seq_len(nrow(ALL_PROTEIN_COMPARTMENTS))) {
    suffix <- ALL_PROTEIN_COMPARTMENTS$suffix[i]
    label <- ALL_PROTEIN_COMPARTMENTS$label[i]
    cat("  ", label, " (", suffix, ")...\n", sep = "")
    cms_results[[suffix]] <- run_compartment_cms_limma(
      cms_group = cms_group,
      compartment_suffix = suffix,
      output_dir = out_dir
    )
  }
  all_results[[cms_group]] <- cms_results

  combined <- dplyr::bind_rows(lapply(cms_results, function(x) x$de))
  readr::write_csv(
    combined,
    file.path(out_dir, paste0(cms_group, "_all_compartments_limma.csv"))
  )
}

if (length(cms_counts) > 0) {
  counts_tbl <- dplyr::bind_rows(cms_counts)
  readr::write_csv(counts_tbl, file.path(OUTPUT_ROOT, "proteomics_patient_counts.csv"))
  cat("\nPatient counts:\n")
  print(counts_tbl)
}

cat("\nCMS-stratified proteomics complete.\n")
cat("Outputs: ", OUTPUT_ROOT, "/{CMS2,CMS4}/\n", sep = "")
