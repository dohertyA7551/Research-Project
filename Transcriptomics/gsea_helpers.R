# Preranked GSEA (fgsea) helpers for limma High vs Low DE results.
# Source after base_dir is set (same pattern as analysis_helpers.R).

msigdb_pathways <- function(
  collection = "H",
  subcollection = NULL,
  species = "Homo sapiens"
) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package 'msigdbr' is required.")
  }
  df <- if (is.null(subcollection)) {
    msigdbr::msigdbr(species = species, collection = collection)
  } else {
    msigdbr::msigdbr(
      species = species,
      collection = collection,
      subcollection = subcollection
    )
  }
  df <- df %>% dplyr::select(gs_name, gene_symbol)
  split(df$gene_symbol, df$gs_name)
}

hallmark_pathways <- function(species = "Homo sapiens") {
  msigdb_pathways(collection = "H", species = species)
}

kegg_pathways <- function(species = "Homo sapiens") {
  msigdb_pathways(
    collection = "C2",
    subcollection = "CP:KEGG_LEGACY",
    species = species
  )
}

prepare_fgsea_stats <- function(limma_results, gene_col = "gene", stat_col = "logFC") {
  df <- limma_results
  df <- df[is.finite(df[[stat_col]]) & !is.na(df[[gene_col]]) & nzchar(df[[gene_col]]), ]
  df <- df[!duplicated(df[[gene_col]]), ]
  stats <- df[[stat_col]]
  names(stats) <- df[[gene_col]]
  stats <- sort(stats, decreasing = TRUE)
  stats
}

run_fgsea_preranked <- function(
  stats,
  pathways,
  minSize = 15,
  maxSize = 500,
  nPermSimple = 10000
) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required.")
  }
  fgsea::fgsea(
    pathways = pathways,
    stats    = stats,
    minSize  = minSize,
    maxSize  = maxSize,
    nPermSimple = nPermSimple
  ) %>%
    as.data.frame() %>%
    dplyr::arrange(dplyr::desc(NES))
}

export_fgsea_tables <- function(
  fgsea_out,
  output_dir,
  file_prefix,
  collection = "hallmark",
  fdr_cutoff = 0.05
) {
  if (!requireNamespace("readr", quietly = TRUE)) library(readr)

  out <- fgsea_out
  if ("leadingEdge" %in% names(out)) {
    out$leadingEdge <- vapply(out$leadingEdge, paste, collapse = ";", character(1))
  }

  all_path <- file.path(
    output_dir,
    paste0(file_prefix, "_fgsea_", collection, "_all.csv")
  )
  readr::write_csv(out, all_path)

  sig <- out %>% dplyr::filter(padj < fdr_cutoff)
  sig_path <- file.path(
    output_dir,
    paste0(file_prefix, "_fgsea_", collection, "_fdr05.csv")
  )
  readr::write_csv(sig, sig_path)

  message("Wrote ", all_path, " (", nrow(out), " pathways)")
  message("Wrote ", sig_path, " (", nrow(sig), " FDR < ", fdr_cutoff, ")")
  invisible(list(all = all_path, sig = sig_path, sig_df = sig))
}

export_fgsea_dotplot <- function(
  fgsea_out,
  output_dir,
  file_prefix,
  collection = "hallmark",
  cohort_title = NULL,
  fdr_cutoff = 0.05,
  top_n = 20,
  width = 10,
  height = 7
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for dot plots.")
  }

  plot_df <- fgsea_out %>%
    dplyr::filter(padj < fdr_cutoff) %>%
    dplyr::arrange(dplyr::desc(abs(NES))) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(plot_df) == 0) {
    plot_df <- fgsea_out %>%
      dplyr::arrange(padj, dplyr::desc(abs(NES))) %>%
      dplyr::slice_head(n = top_n)
  }

  plot_df$pathway <- gsub("^HALLMARK_", "", plot_df$pathway)
  plot_df$pathway <- gsub("^KEGG_", "", plot_df$pathway)
  plot_df$pathway <- gsub("_", " ", plot_df$pathway)
  plot_df$direction <- ifelse(plot_df$NES > 0, "High > Low", "High < Low")
  plot_df$pathway <- factor(
    plot_df$pathway,
    levels = plot_df$pathway[order(plot_df$NES, decreasing = TRUE)]
  )

  title <- paste0(
    if (!is.null(cohort_title)) cohort_title else file_prefix,
    " — fgsea ",
    collection,
    " (High vs Low risk)"
  )
  subtitle <- if (any(fgsea_out$padj < fdr_cutoff, na.rm = TRUE)) {
    paste0("FDR < ", fdr_cutoff, "; top ", nrow(plot_df), " by |NES|")
  } else {
    paste0("No FDR < ", fdr_cutoff, "; top ", nrow(plot_df), " by nominal rank")
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = NES, y = pathway, colour = -log10(padj), size = size)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_gradient(low = "#4575b4", high = "#d73027") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Normalized enrichment score (NES)",
      y = NULL,
      colour = expression(-log[10] * "(FDR)"),
      size = "Set size"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10)
    )

  png_path <- file.path(
    output_dir,
    paste0(file_prefix, "_fgsea_", collection, "_dotplot.png")
  )
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = 150)
  message("Wrote ", png_path)
  invisible(png_path)
}

run_gsea_from_limma_csv <- function(
  limma_csv,
  output_dir,
  file_prefix,
  cohort_title = NULL,
  collections = c("hallmark"),
  fdr_cutoff = 0.05
) {
  if (!requireNamespace("readr", quietly = TRUE)) library(readr)

  limma_results <- readr::read_csv(limma_csv, show_col_types = FALSE)
  stats <- prepare_fgsea_stats(limma_results)
  message(file_prefix, ": ranked ", length(stats), " genes for fgsea")

  pathway_map <- list(
    hallmark = hallmark_pathways(),
    kegg = kegg_pathways()
  )

  results <- list()
  for (coll in collections) {
    if (!coll %in% names(pathway_map)) {
      stop("Unknown collection: ", coll)
    }
    fgsea_out <- run_fgsea_preranked(stats, pathway_map[[coll]])
    export_fgsea_tables(fgsea_out, output_dir, file_prefix, coll, fdr_cutoff)
    export_fgsea_dotplot(
      fgsea_out, output_dir, file_prefix, coll, cohort_title, fdr_cutoff
    )
    sig_n <- sum(fgsea_out$padj < fdr_cutoff, na.rm = TRUE)
    message(file_prefix, " ", coll, ": ", sig_n, " pathways FDR < ", fdr_cutoff)
    results[[coll]] <- fgsea_out
  }
  results
}

run_fgsea_meta_fisher <- function(
  fgsea_results_list,
  output_dir,
  file_prefix = "Transcriptomics",
  collection = "hallmark",
  fdr_cutoff = 0.05
) {
  if (!requireNamespace("readr", quietly = TRUE)) library(readr)

  cohort_names <- names(fgsea_results_list)
  if (length(cohort_names) < 2) {
    stop("Meta-analysis requires at least two cohort fgsea results.")
  }

  all_pathways <- unique(unlist(lapply(fgsea_results_list, function(df) df$pathway)))
  pval_mat <- vapply(
    fgsea_results_list,
    function(df) df$pval[match(all_pathways, df$pathway)],
    numeric(length(all_pathways))
  )
  if (is.vector(pval_mat)) {
    pval_mat <- matrix(pval_mat, ncol = 1)
  }
  nes_mat <- vapply(
    fgsea_results_list,
    function(df) df$NES[match(all_pathways, df$pathway)],
    numeric(length(all_pathways))
  )
  if (is.vector(nes_mat)) {
    nes_mat <- matrix(nes_mat, ncol = 1)
  }
  colnames(pval_mat) <- cohort_names
  colnames(nes_mat) <- cohort_names

  fisher_p <- apply(pval_mat, 1, function(row) {
    row <- row[is.finite(row)]
    if (length(row) < 2) {
      return(NA_real_)
    }
    stats::pchisq(-2 * sum(log(row)), df = 2 * length(row), lower.tail = FALSE)
  })

  meta <- data.frame(
    pathway = all_pathways,
    fisher_p = fisher_p,
    fisher_padj = stats::p.adjust(fisher_p, method = "BH"),
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
  meta <- meta[order(meta$fisher_padj, meta$fisher_p), ]

  out_path <- file.path(output_dir, paste0(file_prefix, "_meta_", collection, "_fisher.csv"))
  readr::write_csv(meta, out_path)
  message("Wrote Fisher meta: ", out_path)

  sig <- meta %>% dplyr::filter(fisher_padj < fdr_cutoff)
  message("Meta pathways FDR < ", fdr_cutoff, ": ", nrow(sig))

  export_fgsea_meta_dotplot(
    meta, output_dir, file_prefix, collection, fdr_cutoff
  )
  export_fgsea_meta_heatmap(
    meta, output_dir, file_prefix, collection, fdr_cutoff
  )

  invisible(meta)
}

export_fgsea_meta_dotplot <- function(
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
    dplyr::filter(.data$fisher_padj < fdr_cutoff) %>%
    dplyr::arrange(.data$fisher_padj, dplyr::desc(abs(.data$mean_NES))) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(plot_df) == 0) {
    plot_df <- meta %>%
      dplyr::arrange(.data$fisher_padj, dplyr::desc(abs(.data$mean_NES))) %>%
      dplyr::slice_head(n = top_n)
  }

  plot_df$pathway_label <- factor(
    plot_df$pathway_label,
    levels = plot_df$pathway_label[order(plot_df$mean_NES, decreasing = TRUE)]
  )

  title <- paste0(
    file_prefix, " meta — Fisher ", collection, " (High vs Low risk)"
  )
  subtitle <- if (any(meta$fisher_padj < fdr_cutoff, na.rm = TRUE)) {
    paste0(
      "FDR < ", fdr_cutoff, "; top ", nrow(plot_df),
      " by |mean NES|; 4 cohorts (fgsea preranked)"
    )
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
      colour = -log10(.data$fisher_padj),
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
      x = "Mean NES across cohorts (negative = down in High risk)",
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
    paste0(file_prefix, "_meta_", collection, "_fisher_dotplot.png")
  )
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = 150)
  message("Wrote ", png_path)
  invisible(png_path)
}

export_fgsea_meta_heatmap <- function(
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
    dplyr::filter(.data$fisher_padj < fdr_cutoff) %>%
    dplyr::arrange(.data$fisher_padj) %>%
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
      title = paste0(
        file_prefix, " meta — per-cohort NES (", collection, ")"
      ),
      subtitle = paste0(
        "Fisher FDR < ", fdr_cutoff, "; top ", nrow(plot_df),
        " pathways; red = up in High, blue = down in High"
      ),
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
    paste0(file_prefix, "_meta_", collection, "_fisher_heatmap.png")
  )
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = 150)
  message("Wrote ", png_path)
  invisible(png_path)
}
