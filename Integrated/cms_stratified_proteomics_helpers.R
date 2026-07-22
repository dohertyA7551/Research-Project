# CMS-stratified proteomics Low vs High risk (CMS2 and CMS4).
source("panel_predicted_observed_network_helpers.R")

CMS_STRATIFIED_GROUPS <- c("CMS2", "CMS4")

ALL_PROTEIN_COMPARTMENTS <- data.frame(
  suffix = c(
    "NonImmuneEpithelium",
    "NonImmuneStroma",
    "helperT",
    "cytT",
    "regT"
  ),
  label = c(
    "Non-immune epithelium (tumour)",
    "Stromal cells",
    "Helper T cells",
    "Cytotoxic T cells",
    "Regulatory T cells"
  ),
  stringsAsFactors = FALSE
)

load_proteomics_with_cms <- function(
  compartment_suffix = DEFAULT_COMPARTMENT_SUFFIX,
  merged_rds = NULL,
  survival_rds = NULL
) {
  paths <- resolve_panel_network_paths()
  if (is.null(merged_rds)) merged_rds <- paths$merged_rds
  if (is.null(survival_rds)) survival_rds <- paths$survival_rds

  merged <- readRDS(merged_rds)
  surv <- readRDS(survival_rds)
  pat <- intersect(merged$patient_id, surv$patient_id)
  merged <- merged[match(pat, merged$patient_id), , drop = FALSE]
  surv <- surv[match(pat, surv$patient_id), , drop = FALSE]

  pat_suffix <- paste0("_mean_", compartment_suffix, "$")
  epi_cols <- grep(pat_suffix, names(merged), value = TRUE, ignore.case = TRUE)
  if (length(epi_cols) == 0L) {
    stop("No proteomics columns for compartment suffix: ", compartment_suffix)
  }

  mat <- as.matrix(merged[, epi_cols, drop = FALSE])
  rownames(mat) <- merged$patient_id
  storage.mode(mat) <- "double"
  markers <- vapply(epi_cols, marker_from_compartment_col, character(1), compartment_suffix)
  mat <- t(mat)
  rownames(mat) <- markers

  meta <- data.frame(
    patient_id = merged$patient_id,
    risk_grp = surv$risk_grp,
    rsf_risk = surv$rsf_risk,
    cms_subtype = merged$cms_subtype,
    stringsAsFactors = FALSE
  )
  list(matrix = mat, meta = meta, markers = markers, feature_cols = epi_cols)
}

subset_protein_cms_risk <- function(
  prot,
  cms_group,
  min_patients = 20L,
  min_low = 3L,
  min_high = 3L
) {
  meta <- prot$meta %>%
    dplyr::filter(
      .data$cms_subtype == cms_group,
      .data$risk_grp %in% c("Low", "High")
    )
  n_low <- sum(meta$risk_grp == "Low")
  n_high <- sum(meta$risk_grp == "High")
  if (n_low < min_low || n_high < min_high) {
    stop(
      "Insufficient patients for ", cms_group,
      ": Low=", n_low, ", High=", n_high
    )
  }
  pts <- meta$patient_id
  if (length(pts) < min_patients) {
    warning(
      "Fewer than ", min_patients, " patients for ", cms_group,
      " (n=", length(pts), "); continuing."
    )
  }
  list(
    matrix = prot$matrix[, pts, drop = FALSE],
    meta = meta,
    markers = rownames(prot$matrix),
    patients = pts,
    cms_group = cms_group,
    n_low = n_low,
    n_high = n_high
  )
}

run_compartment_cms_limma <- function(
  cms_group,
  compartment_suffix,
  output_dir
) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("readr is required.")
  }
  prot <- load_proteomics_with_cms(compartment_suffix = compartment_suffix)
  aligned <- subset_protein_cms_risk(prot, cms_group = cms_group)
  de <- run_limma_high_vs_low_matrix(
    aligned$matrix,
    aligned$meta,
    feature_col = "marker"
  )
  label <- ALL_PROTEIN_COMPARTMENTS$label[
    match(compartment_suffix, ALL_PROTEIN_COMPARTMENTS$suffix)
  ]
  file_prefix <- paste0(cms_group, "_", compartment_suffix)

  de_df <- de$results %>%
    dplyr::mutate(
      marker = rownames(de$results),
      compartment = compartment_suffix,
      compartment_label = label,
      cms_subtype = cms_group
    )
  readr::write_csv(
    de_df,
    file.path(output_dir, paste0(file_prefix, "_limma_high_vs_low.csv"))
  )

  up <- de_df %>%
    dplyr::filter(.data$logFC > 0) %>%
    dplyr::arrange(dplyr::desc(.data$logFC)) %>%
    dplyr::slice_head(n = 10)
  down <- de_df %>%
    dplyr::filter(.data$logFC < 0) %>%
    dplyr::arrange(.data$logFC) %>%
    dplyr::slice_head(n = 10)
  readr::write_csv(
    dplyr::bind_rows(up, down),
    file.path(output_dir, paste0(file_prefix, "_limma_top10_logfc.csv"))
  )

  summary_lines <- c(
    paste0("Proteomics limma: ", cms_group, " — ", label),
    paste0("Patients: ", aligned$n_low + aligned$n_high,
           " (Low=", aligned$n_low, ", High=", aligned$n_high, ")"),
    paste0("Markers: ", nrow(de_df)),
    paste0("Nominal P <= 0.05: ", sum(de_df$P.Value <= 0.05, na.rm = TRUE)),
    paste0("BH FDR <= 0.05: ", sum(de_df$adj.P.Val <= 0.05, na.rm = TRUE)),
    "",
    "Top up in High (logFC):",
    paste0("  ", up$marker[1:min(5, nrow(up))],
           " (logFC=", signif(up$logFC[1:min(5, nrow(up))], 3), ")"),
    "",
    "Top down in High (logFC):",
    paste0("  ", down$marker[1:min(5, nrow(down))],
           " (logFC=", signif(down$logFC[1:min(5, nrow(down))], 3), ")")
  )
  writeLines(
    summary_lines,
    file.path(output_dir, paste0(file_prefix, "_summary.txt"))
  )

  invisible(list(de = de_df, aligned = aligned, summary = summary_lines))
}
