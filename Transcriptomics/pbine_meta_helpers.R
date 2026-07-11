# Pbine p-value meta-analysis (Lin et al. 2022).
# Package source lives in workspace: vendor/Pbine/ (not installed system-wide).

PBINE_SOURCE <- "/work_space/files/Transcriptomics/vendor/Pbine/Pbine.r"

load_pbine <- function() {
  if (!exists("Pbine", mode = "function")) {
    if (!file.exists(PBINE_SOURCE)) {
      stop(
        "Pbine source not found at ", PBINE_SOURCE,
        ". Clone https://github.com/Yinchun-Lin/Pbine into vendor/Pbine."
      )
    }
    source(PBINE_SOURCE, local = FALSE)
  }
  invisible(TRUE)
}

sanitize_pvals_for_pbine <- function(pm, floor = 1e-6, ceiling = 1 - 1e-6) {
  pm <- as.matrix(pm)
  pm[!is.finite(pm)] <- NA_real_
  pm <- pmax(pmin(pm, ceiling), floor)
  pm
}

fisher_combine_pval_matrix <- function(pm) {
  pm <- as.matrix(pm)
  apply(pm, 1, function(row) {
    row <- row[is.finite(row)]
    if (length(row) < 2L) {
      return(NA_real_)
    }
    stats::pchisq(-2 * sum(log(row)), df = 2 * length(row), lower.tail = FALSE)
  })
}

combine_pval_matrix <- function(
  pm,
  method = "Int",
  p_floor = 1e-6,
  rowwise_fallback = TRUE,
  rowwise_max_rows = 500L
) {
  if (method == "Fisher") {
    return(fisher_combine_pval_matrix(sanitize_pvals_for_pbine(pm, floor = p_floor)))
  }
  pbine_combine_pval_matrix(
    pm,
    method = method,
    p_floor = p_floor,
    rowwise_fallback = rowwise_fallback,
    rowwise_max_rows = rowwise_max_rows
  )
}

estimate_pbine_sigma <- function(pm, shrink = 0.25) {
  pm <- as.matrix(pm)
  sigma <- stats::cor(pm, use = "pairwise.complete.obs")
  if (any(is.na(sigma))) {
    sigma[is.na(sigma)] <- 0
  }
  diag(sigma) <- 1
  m <- ncol(pm)
  if (m > 1L) {
    sigma <- (1 - shrink) * sigma + shrink * diag(m)
  }
  sigma
}

suppress_pbine_output <- function(expr) {
  out <- tempfile()
  sink(out)
  on.exit(sink(), add = TRUE)
  force(expr)
}

pbine_combine_row <- function(row, sigma, method = "Int") {
  ok <- is.finite(row)
  if (sum(ok) < 2L) {
    return(NA_real_)
  }
  sub <- matrix(row[ok], nrow = 1)
  sub_sigma <- sigma[ok, ok, drop = FALSE]
  suppress_pbine_output({
    if (sum(ok) == 2L) {
      Pbine(sub, sigma = sub_sigma[1, 2], method = method)[1]
    } else {
      Pbine(sub, sigma = sub_sigma, method = method)[1]
    }
  })
}

# Combine an n x m p-value matrix with Pbine (method = "Int" by default).
pbine_combine_pval_matrix <- function(
  pm,
  method = "Int",
  sigma = NULL,
  p_floor = 1e-6,
  rowwise_fallback = TRUE,
  rowwise_max_rows = 500L
) {
  load_pbine()
  pm <- sanitize_pvals_for_pbine(pm, floor = p_floor)
  if (ncol(pm) < 2L) {
    stop("Pbine requires at least two cohort p-value columns.")
  }

  if (ncol(pm) == 2L) {
    if (is.null(sigma)) {
      return(suppress_pbine_output(Pbine(pm, method = method)))
    }
    return(suppress_pbine_output(Pbine(pm, sigma = sigma, method = method)))
  }

  if (is.null(sigma)) {
    sigma <- estimate_pbine_sigma(pm)
  }

  bulk <- tryCatch(
    suppress_pbine_output(Pbine(pm, sigma = sigma, method = method)),
    error = function(e) e
  )
  if (!inherits(bulk, "error")) {
    return(bulk)
  }

  message(
    "Pbine bulk failed for ", ncol(pm), " cohorts (",
    conditionMessage(bulk), "); using hierarchical 2-cohort Pbine."
  )
  pbine_hierarchical_pvals(pm, method = method, p_floor = p_floor)
}

# Combine >2 cohort p-values via nested 2-cohort Pbine (microarray pair, RNA-seq pair, then final).
pbine_hierarchical_pvals <- function(
  pm,
  method = "Int",
  p_floor = 1e-6,
  platform_groups = NULL
) {
  pm <- sanitize_pvals_for_pbine(pm, floor = p_floor)
  cn <- colnames(pm)
  if (is.null(cn)) {
    cn <- paste0("C", seq_len(ncol(pm)))
    colnames(pm) <- cn
  }
  if (is.null(platform_groups)) {
    platform_groups <- list(
      microarray = intersect(c("Stage2", "Retrospective"), cn),
      rnaseq = intersect(c("Colossus", "Taxonomy"), cn)
    )
  }
  platform_groups <- platform_groups[
    vapply(platform_groups, function(cols) sum(cols %in% cn) >= 2L, logical(1))
  ]
  if (length(platform_groups) == 0L) {
    stop("Need at least one platform group with >=2 cohorts for hierarchical Pbine.")
  }

  stage1 <- lapply(platform_groups, function(cols) {
    sub <- pm[, cols[cols %in% cn], drop = FALSE]
    if (ncol(sub) == 2L) {
      combine_pval_matrix(sub, method = method, p_floor = p_floor, rowwise_fallback = FALSE)
    } else {
      pbine_hierarchical_pvals(sub, method = method, p_floor = p_floor, platform_groups = list(all = colnames(sub)))
    }
  })
  stage_mat <- do.call(cbind, stage1)
  colnames(stage_mat) <- names(stage1)
  if (ncol(stage_mat) == 1L) {
    return(stage_mat[, 1L])
  }
  combine_pval_matrix(
    stage_mat,
    method = method,
    p_floor = p_floor,
    rowwise_fallback = FALSE
  )
}

pbine_meta_limma_genes <- function(
  cohort_results,
  min_cohorts = 2L,
  method = "Int",
  rowwise_max_rows = 500L
) {
  cohort_results <- cohort_results[!vapply(cohort_results, is.null, logical(1))]
  if (length(cohort_results) < min_cohorts) {
    return(NULL)
  }

  cohort_names <- names(cohort_results)
  all_genes <- unique(unlist(lapply(cohort_results, function(df) df$gene)))

  pval_mat <- vapply(
    cohort_results,
    function(df) df$P.Value[match(all_genes, df$gene)],
    numeric(length(all_genes))
  )
  logfc_mat <- vapply(
    cohort_results,
    function(df) df$logFC[match(all_genes, df$gene)],
    numeric(length(all_genes))
  )
  if (is.vector(pval_mat)) {
    pval_mat <- matrix(pval_mat, ncol = 1)
    logfc_mat <- matrix(logfc_mat, ncol = 1)
  }
  colnames(pval_mat) <- cohort_names
  colnames(logfc_mat) <- cohort_names
  rownames(pval_mat) <- all_genes
  rownames(logfc_mat) <- all_genes

  n_cohorts <- rowSums(is.finite(pval_mat))
  keep <- n_cohorts >= min_cohorts
  if (!any(keep)) {
    return(NULL)
  }

  pm <- pval_mat[keep, , drop = FALSE]
  if (nrow(pm) > rowwise_max_rows) {
    message(
      "Gene-level Pbine skipped: ", nrow(pm),
      " genes exceeds rowwise_max_rows=", rowwise_max_rows, "."
    )
    return(NULL)
  }
  message(
    "Pbine gene meta: combining ", nrow(pm), " genes across ",
    ncol(pm), " cohort(s) (method=", method, ")"
  )
  meta_p <- tryCatch(
    pbine_combine_pval_matrix(
      pm,
      method = method,
      rowwise_max_rows = rowwise_max_rows
    ),
    error = function(e) {
      message("Skipping gene-level Pbine meta: ", conditionMessage(e))
      return(rep(NA_real_, nrow(pm)))
    }
  )

  mean_logFC <- rowMeans(logfc_mat[keep, , drop = FALSE], na.rm = TRUE)
  sign_consistent <- apply(logfc_mat[keep, , drop = FALSE], 1, function(row) {
    row <- row[is.finite(row)]
    if (length(row) < 2L) {
      return(NA)
    }
    sum(row > 0) == length(row) || sum(row < 0) == length(row)
  })

  meta <- data.frame(
    gene = rownames(pm),
    meta_p = meta_p,
    meta_padj = stats::p.adjust(meta_p, method = "BH"),
    n_cohorts = n_cohorts[keep],
    mean_logFC = mean_logFC,
    sign_consistent = sign_consistent,
    stringsAsFactors = FALSE
  )
  for (cn in cohort_names) {
    meta[[paste0("logFC_", cn)]] <- logfc_mat[keep, cn]
    meta[[paste0("P.Value_", cn)]] <- pval_mat[keep, cn]
  }
  meta[order(meta$meta_padj, meta$meta_p, na.last = TRUE), , drop = FALSE]
}

run_fgsea_meta_pbine <- function(
  fgsea_results_list,
  output_dir,
  file_prefix = "Transcriptomics",
  collection = "hallmark",
  fdr_cutoff = 0.05,
  method = "Int"
) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    library(readr)
  }

  cohort_names <- names(fgsea_results_list)
  if (length(cohort_names) < 2L) {
    stop("Pbine GSEA meta requires at least two cohort fgsea results.")
  }

  all_pathways <- unique(unlist(lapply(fgsea_results_list, function(df) df$pathway)))
  pval_mat <- vapply(
    fgsea_results_list,
    function(df) df$pval[match(all_pathways, df$pathway)],
    numeric(length(all_pathways))
  )
  nes_mat <- vapply(
    fgsea_results_list,
    function(df) df$NES[match(all_pathways, df$pathway)],
    numeric(length(all_pathways))
  )
  if (is.vector(pval_mat)) {
    pval_mat <- matrix(pval_mat, ncol = 1)
    nes_mat <- matrix(nes_mat, ncol = 1)
  }
  colnames(pval_mat) <- cohort_names
  colnames(nes_mat) <- cohort_names
  rownames(pval_mat) <- all_pathways

  meta_p <- pbine_combine_pval_matrix(
    pval_mat,
    method = method,
    rowwise_fallback = FALSE
  )

  meta <- data.frame(
    pathway = all_pathways,
    meta_p = meta_p,
    meta_padj = stats::p.adjust(meta_p, method = "BH"),
    n_cohorts = apply(pval_mat, 1, function(row) sum(is.finite(row))),
    mean_NES = rowMeans(nes_mat, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  for (cn in cohort_names) {
    meta[[paste0("pval_", cn)]] <- pval_mat[, cn]
    meta[[paste0("NES_", cn)]] <- nes_mat[, cn]
  }
  meta$pathway_label <- gsub("^HALLMARK_", "", meta$pathway)
  meta$pathway_label <- gsub("^KEGG_", "", meta$pathway_label)
  meta$pathway_label <- gsub("_", " ", meta$pathway_label)
  meta <- meta[order(meta$meta_padj, meta$meta_p), ]

  out_path <- file.path(output_dir, paste0(file_prefix, "_meta_", collection, "_pbine.csv"))
  readr::write_csv(meta, out_path)
  message("Wrote Pbine meta: ", out_path)

  sig <- meta %>% dplyr::filter(.data$meta_padj < fdr_cutoff)
  message("Meta pathways FDR < ", fdr_cutoff, ": ", nrow(sig))

  export_fgsea_meta_dotplot_pbine(
    meta, output_dir, file_prefix, collection, fdr_cutoff
  )
  export_fgsea_meta_heatmap_pbine(
    meta, output_dir, file_prefix, collection, fdr_cutoff
  )

  invisible(meta)
}

export_fgsea_meta_dotplot_pbine <- function(
  meta,
  output_dir,
  file_prefix,
  collection = "hallmark",
  fdr_cutoff = 0.05,
  top_n = 25,
  width = 10,
  height = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for meta dot plots.")
  }

  plot_df <- meta %>%
    dplyr::filter(.data$meta_padj < fdr_cutoff) %>%
    dplyr::arrange(.data$meta_padj, dplyr::desc(abs(.data$mean_NES))) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(plot_df) == 0) {
    plot_df <- meta %>%
      dplyr::arrange(.data$meta_padj, dplyr::desc(abs(.data$mean_NES))) %>%
      dplyr::slice_head(n = top_n)
  }

  plot_df$pathway_label <- factor(
    plot_df$pathway_label,
    levels = plot_df$pathway_label[order(plot_df$mean_NES, decreasing = TRUE)]
  )

  title <- paste0(file_prefix, " meta — Pbine ", collection)
  subtitle <- if (any(meta$meta_padj < fdr_cutoff, na.rm = TRUE)) {
    paste0("FDR < ", fdr_cutoff, "; top ", nrow(plot_df), " by |mean NES|")
  } else {
    paste0("No FDR < ", fdr_cutoff, "; top ", nrow(plot_df), " by rank")
  }

  if (is.null(height)) {
    height <- max(7, 0.32 * nrow(plot_df) + 2)
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$mean_NES,
      y = .data$pathway_label,
      colour = -log10(.data$meta_padj),
      size = .data$n_cohorts
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_gradient(low = "#4575b4", high = "#d73027") +
    ggplot2::scale_size(range = c(3, 7), breaks = sort(unique(plot_df$n_cohorts))) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Mean NES across cohorts",
      y = NULL,
      colour = expression(-log[10] * "(meta FDR)"),
      size = "Cohorts"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9)
    )

  png_path <- file.path(
    output_dir,
    paste0(file_prefix, "_meta_", collection, "_pbine_dotplot.png")
  )
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = 150)
  message("Wrote ", png_path)
  invisible(png_path)
}

export_fgsea_meta_heatmap_pbine <- function(
  meta,
  output_dir,
  file_prefix,
  collection = "hallmark",
  fdr_cutoff = 0.05,
  top_n = 30,
  width = 9,
  height = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for meta heatmaps.")
  }

  nes_cols <- grep("^NES_", names(meta), value = TRUE)
  if (length(nes_cols) == 0) {
    stop("No NES_<cohort> columns found in meta results.")
  }

  plot_df <- meta %>%
    dplyr::filter(.data$meta_padj < fdr_cutoff) %>%
    dplyr::arrange(.data$meta_padj) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(plot_df) == 0) {
    message("No meta pathways FDR < ", fdr_cutoff, "; skipping heatmap.")
    return(invisible(NULL))
  }

  long <- plot_df[, c("pathway_label", nes_cols), drop = FALSE]
  long <- stats::setNames(long, c("pathway_label", gsub("^NES_", "", nes_cols)))
  rownames(long) <- NULL
  long <- tidyr::pivot_longer(
    long,
    cols = -pathway_label,
    names_to = "cohort",
    values_to = "NES"
  )

  pathway_order <- plot_df$pathway_label[order(plot_df$mean_NES, decreasing = TRUE)]
  long$pathway_label <- factor(long$pathway_label, levels = rev(pathway_order))
  long$cohort <- factor(long$cohort, levels = gsub("^NES_", "", nes_cols))

  if (is.null(height)) {
    height <- max(6, 0.28 * length(pathway_order) + 2)
  }

  nes_range <- max(abs(long$NES[is.finite(long$NES)]), na.rm = TRUE)
  if (!is.finite(nes_range) || nes_range == 0) {
    nes_range <- 1
  }

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$cohort, y = .data$pathway_label, fill = .data$NES)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-nes_range, nes_range),
      name = "NES"
    ) +
    ggplot2::labs(
      title = paste0(file_prefix, " meta — per-cohort NES (", collection, ", Pbine)"),
      subtitle = paste0("Meta FDR < ", fdr_cutoff, "; top ", nrow(plot_df), " pathways"),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid = ggplot2::element_blank()
    )

  png_path <- file.path(
    output_dir,
    paste0(file_prefix, "_meta_", collection, "_pbine_heatmap.png")
  )
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = 150)
  message("Wrote ", png_path)
  invisible(png_path)
}

run_xcell_risk_meta_pbine <- function(
  diff_list,
  output_dir,
  file_prefix = "xcell_high_vs_low",
  p_cutoff = 0.05,
  method = "Int",
  min_cohorts = 2L,
  export_plots = TRUE
) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    library(readr)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    library(dplyr)
  }

  diff_list <- diff_list[!vapply(diff_list, is.null, logical(1))]
  cohort_names <- names(diff_list)
  if (length(cohort_names) < 2L) {
    stop("Pbine xCell meta requires at least two cohort Wilcoxon results.")
  }

  all_cell_types <- unique(unlist(lapply(diff_list, function(df) df$cell_type)))
  pval_mat <- vapply(
    diff_list,
    function(df) df$pval[match(all_cell_types, df$cell_type)],
    numeric(length(all_cell_types))
  )
  logfc_mat <- vapply(
    diff_list,
    function(df) df$log2FC[match(all_cell_types, df$cell_type)],
    numeric(length(all_cell_types))
  )
  if (is.vector(pval_mat)) {
    pval_mat <- matrix(pval_mat, ncol = 1)
    logfc_mat <- matrix(logfc_mat, ncol = 1)
  }
  colnames(pval_mat) <- cohort_names
  colnames(logfc_mat) <- cohort_names
  rownames(pval_mat) <- all_cell_types

  n_cohorts <- apply(pval_mat, 1, function(row) sum(is.finite(row)))
  keep <- n_cohorts >= min_cohorts
  if (!any(keep)) {
    stop("No cell types with p-values in >= ", min_cohorts, " cohorts.")
  }

  pm <- pval_mat[keep, , drop = FALSE]
  meta_p <- tryCatch(
    pbine_hierarchical_pvals(pm, method = method),
    error = function(e) {
      if (method != "Fisher") {
        message(
          "Pbine method=", method, " failed (", conditionMessage(e),
          "); retrying with Fisher."
        )
        return(pbine_hierarchical_pvals(pm, method = "Fisher"))
      }
      stop(e)
    }
  )

  mean_log2FC <- rowMeans(logfc_mat[keep, , drop = FALSE], na.rm = TRUE)
  sign_consistent <- apply(logfc_mat[keep, , drop = FALSE], 1, function(row) {
    row <- row[is.finite(row)]
    if (length(row) < 2L) {
      return(NA)
    }
    sum(row > 0) == length(row) || sum(row < 0) == length(row)
  })

  meta <- data.frame(
    cell_type = rownames(pm),
    meta_p = meta_p,
    n_cohorts = n_cohorts[keep],
    mean_log2FC = mean_log2FC,
    sign_consistent = sign_consistent,
    stringsAsFactors = FALSE
  )
  for (cn in cohort_names) {
    meta[[paste0("pval_", cn)]] <- pval_mat[keep, cn]
    meta[[paste0("log2FC_", cn)]] <- logfc_mat[keep, cn]
  }
  meta <- meta[order(meta$meta_p), ]

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(output_dir, paste0(file_prefix, "_meta_pbine.csv"))
  readr::write_csv(meta, out_path)
  message("Wrote Pbine meta: ", out_path)

  sig <- meta %>% dplyr::filter(.data$meta_p < p_cutoff)
  message("Meta cell types p < ", p_cutoff, ": ", nrow(sig))

  if (export_plots) {
    export_xcell_risk_meta_dotplot(meta, output_dir, file_prefix, p_cutoff)
    export_xcell_risk_meta_volcano(meta, output_dir, file_prefix, p_cutoff)
  }

  invisible(meta)
}

run_xcell_fisher_meta <- function(
  diff_list,
  output_dir,
  file_prefix = "xcell",
  p_cutoff = 0.05,
  min_cohorts = 2L,
  export_plots = TRUE
) {
  meta <- run_xcell_risk_meta_pbine(
    diff_list = diff_list,
    output_dir = output_dir,
    file_prefix = file_prefix,
    p_cutoff = p_cutoff,
    method = "Fisher",
    min_cohorts = min_cohorts,
    export_plots = export_plots
  )

  fisher_path <- file.path(output_dir, paste0(file_prefix, "_meta_fisher.csv"))
  pbine_path <- file.path(output_dir, paste0(file_prefix, "_meta_pbine.csv"))
  if (file.exists(pbine_path)) {
    file.copy(pbine_path, fisher_path, overwrite = TRUE)
    message("Wrote ", fisher_path)
  }
  invisible(meta)
}

run_xcell_combined_hl_g1g4_fisher_meta <- function(
  hl_diff_list,
  g1g4_diff_list,
  output_dir,
  file_prefix = "xcell_highlow_plus_g1g4",
  p_cutoff = 0.05,
  min_cohorts = 2L
) {
  if (!requireNamespace("readr", quietly = TRUE)) library(readr)
  if (!requireNamespace("dplyr", quietly = TRUE)) library(dplyr)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  meta_hl <- run_xcell_fisher_meta(
    hl_diff_list,
    output_dir = output_dir,
    file_prefix = "xcell_high_vs_low",
    p_cutoff = p_cutoff,
    min_cohorts = min_cohorts,
    export_plots = TRUE
  )
  meta_g1 <- run_xcell_fisher_meta(
    g1g4_diff_list,
    output_dir = output_dir,
    file_prefix = "xcell_g1_vs_g4",
    p_cutoff = p_cutoff,
    min_cohorts = min_cohorts,
    export_plots = TRUE
  )

  shared <- intersect(meta_hl$cell_type, meta_g1$cell_type)
  if (length(shared) == 0L) {
    stop("No shared cell types between High/Low and G1/G4 meta results.")
  }

  hl_p <- meta_hl$meta_p[match(shared, meta_hl$cell_type)]
  g1_p <- meta_g1$meta_p[match(shared, meta_g1$cell_type)]
  pm2 <- cbind(hl_p, g1_p)
  colnames(pm2) <- c("high_vs_low", "g1_vs_g4")
  rownames(pm2) <- shared

  combined_p <- fisher_combine_pval_matrix(pm2)
  mean_hl <- meta_hl$mean_log2FC[match(shared, meta_hl$cell_type)]
  mean_g1 <- meta_g1$mean_log2FC[match(shared, meta_g1$cell_type)]
  mean_combined <- rowMeans(cbind(mean_hl, mean_g1), na.rm = TRUE)

  combined <- data.frame(
    cell_type = shared,
    meta_p = combined_p,
    meta_p_high_vs_low = hl_p,
    meta_p_g1_vs_g4 = g1_p,
    mean_log2FC_high_vs_low = mean_hl,
    mean_log2FC_g1_vs_g4 = mean_g1,
    mean_log2FC = mean_combined,
    n_inputs = apply(pm2, 1, function(row) sum(is.finite(row))),
    stringsAsFactors = FALSE
  )
  combined <- combined[order(combined$meta_p), ]

  out_path <- file.path(output_dir, paste0(file_prefix, "_meta_fisher.csv"))
  readr::write_csv(combined, out_path)
  message("Wrote combined Fisher meta: ", out_path)
  message(
    "Combined meta p < ", p_cutoff, ": ",
    sum(combined$meta_p < p_cutoff, na.rm = TRUE)
  )

  export_xcell_risk_meta_dotplot(combined, output_dir, file_prefix, p_cutoff)
  export_xcell_risk_meta_volcano(combined, output_dir, file_prefix, p_cutoff)

  invisible(list(
    high_vs_low = meta_hl,
    g1_vs_g4 = meta_g1,
    combined = combined
  ))
}

export_xcell_risk_meta_dotplot <- function(
  meta,
  output_dir,
  file_prefix = "xcell_high_vs_low",
  p_cutoff = 0.05,
  top_n = 20,
  width = 10,
  height = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for meta dot plots.")
  }

  plot_df <- meta %>%
    dplyr::filter(.data$meta_p < p_cutoff) %>%
    dplyr::arrange(.data$meta_p, dplyr::desc(abs(.data$mean_log2FC))) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(plot_df) == 0) {
    plot_df <- meta %>%
      dplyr::arrange(.data$meta_p, dplyr::desc(abs(.data$mean_log2FC))) %>%
      dplyr::slice_head(n = top_n)
  }

  plot_df$cell_type <- factor(
    plot_df$cell_type,
    levels = plot_df$cell_type[order(plot_df$mean_log2FC, decreasing = TRUE)]
  )

  title <- paste0(file_prefix, " meta — Fisher xCell2")
  subtitle <- if (any(meta$meta_p < p_cutoff, na.rm = TRUE)) {
    paste0("p < ", p_cutoff, "; top ", nrow(plot_df), " by |mean log2FC|")
  } else {
    paste0("No p < ", p_cutoff, "; top ", nrow(plot_df), " by rank")
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$mean_log2FC,
      y = .data$cell_type,
      size = -log10(pmax(.data$meta_p, .Machine$double.xmin)),
      colour = .data$meta_p < p_cutoff
    )
  ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "#B2182B", `FALSE` = "grey55"),
      labels = c(`TRUE` = paste0("p < ", p_cutoff), `FALSE` = "Not sig."),
      name = NULL
    ) +
    ggplot2::scale_size_continuous(range = c(3, 10), name = expression(-log[10] * "(p)")) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Mean log2FC (High / Low) across cohorts",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (is.null(height)) {
    height <- max(5, 0.35 * nrow(plot_df) + 2)
  }
  png_path <- file.path(output_dir, paste0(file_prefix, "_meta_pbine_dotplot.png"))
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = 150)
  message("Wrote ", png_path)
  invisible(png_path)
}

export_xcell_risk_meta_volcano <- function(
  meta,
  output_dir,
  file_prefix = "xcell_high_vs_low",
  p_cutoff = 0.05,
  fc_cutoff = 0.5,
  width = 10,
  height = 7
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for meta volcano plots.")
  }

  plot_df <- meta %>%
    dplyr::filter(is.finite(.data$mean_log2FC), is.finite(.data$meta_p)) %>%
    dplyr::mutate(
      neg_log10_p = -log10(pmax(.data$meta_p, .Machine$double.xmin)),
      sig = .data$meta_p < p_cutoff & abs(.data$mean_log2FC) > fc_cutoff
    )
  if (nrow(plot_df) == 0) {
    message("No data for meta volcano.")
    return(invisible(NULL))
  }

  xmax <- min(3, max(1.25, stats::quantile(abs(plot_df$mean_log2FC), 0.98, na.rm = TRUE) * 1.15))
  ymax <- max(1.5, max(plot_df$neg_log10_p, na.rm = TRUE) * 1.08)
  label_df <- plot_df %>%
    dplyr::arrange(.data$meta_p, dplyr::desc(abs(.data$mean_log2FC))) %>%
    dplyr::slice_head(n = 12)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$mean_log2FC, y = .data$neg_log10_p, colour = .data$sig)
  ) +
    ggplot2::geom_point(size = 2.5, alpha = 0.75) +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", colour = "grey70") +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", colour = "grey70") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60")) +
    ggplot2::coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax)) +
    ggplot2::labs(
      title = paste0(file_prefix, " meta volcano — Fisher xCell2"),
      subtitle = paste0(
        "Stage2 + Retrospective + Colossus + Taxonomy | ",
        sum(plot_df$sig, na.rm = TRUE), " sig at p<", p_cutoff, " & |log2FC|>", fc_cutoff
      ),
      x = "Mean log2FC",
      y = expression(-log[10] * "(meta p-value)"),
      colour = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = .data$cell_type),
      size = 3,
      max.overlaps = 20,
      show.legend = FALSE
    )
  }

  png_path <- file.path(output_dir, paste0(file_prefix, "_meta_pbine_volcano.png"))
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = 150, bg = "white")
  message("Wrote ", png_path)
  invisible(png_path)
}
