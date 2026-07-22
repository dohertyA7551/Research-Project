# CMS-stratified Low vs High risk helpers (CMS2 and CMS4 only).
# CMS3 is excluded from stratified comparisons (insufficient / no High risk for analysis).
source("analysis_helpers.R")

CMS_STRATIFIED_GROUPS <- c("CMS2", "CMS4")

resolve_cms_label <- function(belfast = NA_character_, manuela = NA_character_) {
  b <- trimws(as.character(belfast))
  m <- trimws(as.character(manuela))
  valid <- c("CMS1", "CMS2", "CMS3", "CMS4")
  out <- rep(NA_character_, length(b))
  b_ok <- !is.na(b) & b %in% valid
  m_ok <- !is.na(m) & m %in% valid
  out[b_ok] <- b[b_ok]
  out[!b_ok & m_ok] <- m[!b_ok & m_ok]
  if (length(out) == 1L) {
    return(out[[1L]])
  }
  out
}

load_proteomics_cms_table <- function(
  merged_rds = NULL,
  survival_rds = NULL
) {
  if (is.null(merged_rds)) {
    merged_rds <- file.path(
      dirname(dirname(normalizePath("."))),
      "proteomics",
      "RCSI_Taxonomy_clin_prot_merged_df.rds"
    )
    if (!file.exists(merged_rds)) {
      merged_rds <- file.path("proteomics", "RCSI_Taxonomy_clin_prot_merged_df.rds")
    }
    if (!file.exists(merged_rds)) {
      merged_rds <- file.path("../proteomics", "RCSI_Taxonomy_clin_prot_merged_df.rds")
    }
  }
  if (is.null(survival_rds)) {
    survival_rds <- sub("clin_prot_merged_df", "survival_df_with_risk", merged_rds)
    if (!file.exists(survival_rds)) {
      survival_rds <- file.path(dirname(merged_rds), "survival_df_with_risk.rds")
    }
  }
  merged <- readRDS(merged_rds)
  surv <- readRDS(survival_rds)
  merged %>%
    dplyr::transmute(
      patient_id = .data$patient_id,
      cms_subtype = .data$cms_subtype
    ) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      surv %>% dplyr::select(.data$patient_id, .data$risk_grp, .data$rsf_risk),
      by = "patient_id"
    )
}

load_taxonomy_cms_table <- function(
  map_path = NULL,
  survival_path = NULL
) {
  if (is.null(map_path)) {
    map_path <- file.path(
      "Taxonomy_calls",
      "Manuela_and_Belfast_RNA_classifications_CMS_CRIS.txt"
    )
    if (!file.exists(map_path)) {
      map_path <- file.path("Transcriptomics", map_path)
    }
  }
  if (is.null(survival_path)) {
    survival_path <- file.path("proteomics", "survival_df_with_risk.rds")
    if (!file.exists(survival_path)) {
      survival_path <- file.path("../proteomics", "survival_df_with_risk.rds")
    }
  }
  map_df <- readr::read_tsv(map_path, show_col_types = FALSE)
  surv <- readRDS(survival_path)
  map_df %>%
    dplyr::transmute(
      sample_id = .data$Patient,
      patient_id = .data$Code,
      cms_subtype = resolve_cms_label(.data$CMS_Belfast, .data$CMS_Manuela)
    ) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      surv %>% dplyr::select(.data$patient_id, .data$risk_grp, .data$rsf_risk),
      by = "patient_id"
    )
}

build_taxonomy_sample_meta_with_cms <- function(
  map_path = NULL,
  survival_path = NULL
) {
  load_taxonomy_cms_table(map_path = map_path, survival_path = survival_path) %>%
    dplyr::transmute(
      sample_id = .data$sample_id,
      patient_id = .data$patient_id,
      risk_grp = .data$risk_grp,
      rsf_risk = .data$rsf_risk,
      cms_subtype = .data$cms_subtype,
      cohort = "Taxonomy"
    )
}

count_cms_risk <- function(meta, cms_group) {
  sub <- meta %>%
    dplyr::filter(
      .data$cms_subtype == cms_group,
      .data$risk_grp %in% c("Low", "High")
    )
  c(
    cms = cms_group,
    n_total = nrow(sub),
    n_low = sum(sub$risk_grp == "Low"),
    n_high = sum(sub$risk_grp == "High")
  )
}

filter_sample_meta_cms_risk <- function(
  sample_meta,
  cms_group,
  min_low = 3L,
  min_high = 3L
) {
  meta <- sample_meta %>%
    dplyr::filter(
      .data$cms_subtype == cms_group,
      .data$risk_grp %in% c("Low", "High")
    )
  n_low <- sum(meta$risk_grp == "Low")
  n_high <- sum(meta$risk_grp == "High")
  if (n_low < min_low || n_high < min_high) {
    stop(
      "Insufficient samples for ", cms_group,
      ": Low=", n_low, ", High=", n_high,
      " (need >= ", min_low, " Low and >= ", min_high, " High)"
    )
  }
  meta
}

filter_patient_meta_cms_risk <- function(
  patient_meta,
  cms_group,
  min_low = 3L,
  min_high = 3L
) {
  meta <- patient_meta %>%
    dplyr::filter(
      .data$cms_subtype == cms_group,
      .data$risk_grp %in% c("Low", "High")
    )
  n_low <- sum(meta$risk_grp == "Low")
  n_high <- sum(meta$risk_grp == "High")
  if (n_low < min_low || n_high < min_high) {
    stop(
      "Insufficient patients for ", cms_group,
      ": Low=", n_low, ", High=", n_high,
      " (need >= ", min_low, " Low and >= ", min_high, " High)"
    )
  }
  meta
}

cms_output_dir <- function(base_dir, cms_group, subdir = NULL) {
  out <- file.path(base_dir, cms_group)
  if (!is.null(subdir)) {
    out <- file.path(out, subdir)
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

write_cms_counts_summary <- function(counts_list, out_path) {
  tbl <- dplyr::bind_rows(lapply(counts_list, function(x) {
    data.frame(
      cms = x[["cms"]],
      n_total = as.integer(x[["n_total"]]),
      n_low = as.integer(x[["n_low"]]),
      n_high = as.integer(x[["n_high"]]),
      modality = x[["modality"]],
      stringsAsFactors = FALSE
    )
  }))
  readr::write_csv(tbl, out_path)
  writeLines(
    c(
      "CMS-stratified Low vs High risk sample counts",
      paste0("Generated: ", Sys.time()),
      "",
      paste(capture.output(print(tbl, row.names = FALSE)), collapse = "\n")
    ),
    sub("\\.csv$", "_summary.txt", out_path)
  )
  invisible(tbl)
}

export_xcell_volcano <- function(
  diff,
  title,
  subtitle,
  out_path,
  fc_cutoff = 0.5,
  label_top = 12L
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for xCell volcano plots.")
  }
  plot_df <- diff %>%
    dplyr::filter(is.finite(.data$log2FC), !is.na(.data$pval)) %>%
    dplyr::mutate(
      neg_log10_p = -log10(pmax(.data$pval, .Machine$double.xmin)),
      sig = .data$pval <= 0.05 & abs(.data$log2FC) >= fc_cutoff
    )
  if (nrow(plot_df) == 0) {
    message("No xCell data for volcano: ", out_path)
    return(invisible(NULL))
  }
  axes <- volcano_axis_limits(plot_df, fc_cap = 5, xmax_ceiling = 2)
  label_df <- plot_df %>%
    dplyr::arrange(dplyr::desc(abs(.data$log2FC))) %>%
    dplyr::slice_head(n = label_top)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$log2FC, y = .data$neg_log10_p, colour = .data$sig)
  ) +
    ggplot2::geom_point(size = 2, alpha = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60")) +
    ggplot2::coord_cartesian(xlim = axes$xlim, ylim = axes$ylim) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Log2 fold change (High / Low)",
      y = expression(-log[10] * "(p-value)")
    ) +
    ggplot2::theme_bw(base_size = 11)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = .data$cell_type),
      size = 3,
      max.overlaps = Inf,
      colour = "black"
    )
  }

  ggplot2::ggsave(out_path, plot = p, width = 10, height = 7, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(out_path)
}

run_cms_taxonomy_limma_gsea <- function(
  cms_group,
  expr_mat,
  sample_meta,
  output_dir,
  file_prefix,
  cohort_title
) {
  source("gsea_helpers.R")
  meta <- filter_sample_meta_cms_risk(sample_meta, cms_group)
  expr <- expr_mat[, meta$sample_id, drop = FALSE]

  limma_out <- run_limma_high_vs_low(
    expr,
    meta,
    output_dir,
    file_prefix = file_prefix,
    min_group = 3L
  )
  export_limma_figures(
    limma_out$results,
    expr,
    meta,
    output_dir,
    file_prefix,
    cohort_title
  )

  gsea_out <- run_gsea_from_limma_csv(
    limma_csv = file.path(
      output_dir,
      paste0(file_prefix, "_limma_high_vs_low_all_genes.csv")
    ),
    output_dir = output_dir,
    file_prefix = file_prefix,
    cohort_title = cohort_title,
    collections = c("hallmark", "kegg")
  )

  invisible(list(limma = limma_out, gsea = gsea_out))
}

run_cms_taxonomy_progeny <- function(
  cms_group,
  expr_mat,
  sample_meta,
  output_dir,
  file_prefix,
  title
) {
  source("progeny_helpers.R")
  meta <- filter_sample_meta_cms_risk(sample_meta, cms_group)
  expr <- expr_mat[, meta$sample_id, drop = FALSE]
  run_progeny_risk_analysis(
    expr_mat = expr,
    sample_meta = meta,
    output_dir = output_dir,
    file_prefix = file_prefix,
    title = title,
    apply_class_norm = TRUE
  )
}

run_cms_taxonomy_xcell <- function(
  cms_group,
  scores_long,
  output_dir,
  file_prefix,
  title
) {
  cms_map <- scores_long %>%
    dplyr::distinct(.data$sample, .data$patient_id) %>%
    dplyr::left_join(
      load_taxonomy_cms_table() %>%
        dplyr::select(.data$sample_id, .data$cms_subtype),
      by = c("sample" = "sample_id")
    )
  scores <- scores_long %>%
    dplyr::left_join(
      cms_map %>% dplyr::select(.data$sample, .data$cms_subtype),
      by = "sample"
    ) %>%
    dplyr::filter(
      .data$cms_subtype == cms_group,
      .data$risk_grp %in% c("Low", "High")
    )
  n_low <- dplyr::n_distinct(scores$sample[scores$risk_grp == "Low"])
  n_high <- dplyr::n_distinct(scores$sample[scores$risk_grp == "High"])
  if (n_low < 3L || n_high < 3L) {
    stop(
      "Insufficient xCell samples for ", cms_group,
      ": Low=", n_low, ", High=", n_high
    )
  }

  diff <- wilcox_enrichment_diff(scores)
  readr::write_csv(
    diff,
    file.path(output_dir, paste0(file_prefix, "_xcell_wilcox.csv"))
  )
  export_top_log2fc_table(
    diff,
    file.path(output_dir, paste0(file_prefix, "_xcell_top10_logfc.csv")),
    n = 10
  )
  export_xcell_volcano(
    diff,
    title = paste0("xCell: High vs Low risk (", title, ")"),
    subtitle = paste0(cms_group, " | n = ", n_low + n_high, " (Low ", n_low, ", High ", n_high, ")"),
    out_path = file.path(output_dir, paste0(file_prefix, "_xcell_volcano.png"))
  )

  summary_lines <- c(
    paste0("xCell analysis: ", title),
    paste0("CMS subgroup: ", cms_group),
    paste0("Samples: ", n_low + n_high, " (Low=", n_low, ", High=", n_high, ")"),
    paste0(
      "Significant cell types (p <= 0.05, |log2FC| > 0.5): ",
      sum(diff$pval <= 0.05 & abs(diff$log2FC) > 0.5, na.rm = TRUE),
      " / ", nrow(diff)
    )
  )
  writeLines(summary_lines, file.path(output_dir, paste0(file_prefix, "_xcell_summary.txt")))
  cat(paste(summary_lines, collapse = "\n"), "\n")
  invisible(diff)
}
