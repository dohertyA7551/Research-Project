#!/usr/bin/env Rscript
# Export patient-level RNA (genes x samples) and protein (markers x samples) for GeneExpressPred.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("rna_protein_discordance_helpers.R")

out_dir <- "/work_space/files/Integrated/results/geneexpresspred/input"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p <- resolve_discordance_paths()
prot <- load_patient_protein_matrix(p$merged_rds)
rna <- load_patient_rna_cohort_z()
pts <- intersect(colnames(rna$expr), colnames(prot$matrix))
message("Patients with RNA + protein: ", length(pts))

x_mat <- rna$expr[, pts, drop = FALSE]
y_mat <- prot$matrix[, pts, drop = FALSE]

write.table(
  x_mat,
  file.path(out_dir, "colon_codex_RNA.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE
)
write.table(
  y_mat,
  file.path(out_dir, "colon_codex_protein.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE
)
write_csv(
  data.frame(
    patient_id = pts,
    cohort = rna$meta$cohort[match(pts, rna$meta$patient_id)]
  ),
  file.path(out_dir, "patient_meta.csv")
)
writeLines(
  c(
    paste0("Patients: ", length(pts)),
    paste0("RNA genes: ", nrow(x_mat)),
    paste0("Protein features: ", nrow(y_mat))
  ),
  file.path(out_dir, "export_summary.txt")
)
cat("Exported to ", out_dir, "\n")
