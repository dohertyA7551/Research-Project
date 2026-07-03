# Shared helpers: always analyse working copies, never modify originals.
# Source from cohort Rmds after base_dir is set.

setup_analysis_dirs <- function(base_dir, cohort) {
  work_dir <- file.path(base_dir, "analysis_working", cohort)
  output_dir <- file.path(base_dir, "analysis_output", cohort)
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  list(work_dir = work_dir, output_dir = output_dir)
}

ensure_working_copy <- function(src_path, work_dir) {
  if (!file.exists(src_path)) {
    stop("Original dataset not found: ", src_path)
  }
  dest_path <- file.path(work_dir, basename(src_path))
  copied <- file.copy(src_path, dest_path, overwrite = TRUE)
  if (!copied) {
    stop("Failed to copy original to working directory: ", src_path)
  }
  message("Working copy: ", dest_path)
  dest_path
}

ensure_working_copies <- function(src_paths, work_dir) {
  dest_paths <- vapply(unname(src_paths), ensure_working_copy, character(1), work_dir = work_dir)
  names(dest_paths) <- names(src_paths)
  dest_paths
}

# Wilcoxon High vs Low enrichment per cell type (raw p-values).
wilcox_enrichment_diff <- function(scores_long) {
  diff <- scores_long %>%
    dplyr::group_by(.data$cell_type) %>%
    dplyr::summarise(
      mean_high = mean(.data$enrichment[.data$risk_grp == "High"], na.rm = TRUE),
      mean_low  = mean(.data$enrichment[.data$risk_grp == "Low"], na.rm = TRUE),
      log2FC    = log2((mean_high + 1e-10) / (mean_low + 1e-10)),
      pval      = stats::wilcox.test(.data$enrichment ~ .data$risk_grp)$p.value,
      .groups   = "drop"
    )
  diff$neg_log10_p <- -log10(diff$pval)
  diff
}

# Volcano axis limits: tighter than +/-2, but always show p <= p_cutoff points.
# Extreme log2FC artefacts (near-zero denominators) are excluded from scaling only.
volcano_axis_limits <- function(
  diff,
  fc_cap = 5,
  xmax_ceiling = 2,
  p_cutoff = 0.05001
) {
  p_use <- diff$pval

  fc <- diff$log2FC[is.finite(diff$log2FC)]
  p_pos <- p_use[is.finite(p_use) & p_use > 0]

  fc_body <- fc[abs(fc) <= fc_cap]
  if (length(fc_body) == 0) {
    fc_body <- fc
  }

  xmax <- min(
    xmax_ceiling,
    max(1.25, as.numeric(stats::quantile(abs(fc_body), 0.95, na.rm = TRUE)) * 1.1)
  )

  sig_mask <- is.finite(p_use) & is.finite(diff$log2FC) & p_use <= p_cutoff
  sig_fc <- abs(diff$log2FC[sig_mask])
  if (length(sig_fc) > 0) {
    xmax <- min(xmax_ceiling, max(xmax, max(sig_fc, na.rm = TRUE) * 1.08))
  }

  ymax <- max(1.5, max(-log10(p_pos), na.rm = TRUE) + 0.75)

  list(xlim = c(-xmax, xmax), ylim = c(0, ymax))
}

# Top n cell types up (High > Low) and down (High < Low).
# When prioritize_significant = TRUE, rows meeting volcano criteria
# (pval <= p_cutoff and |log2FC| > fc_cutoff) are listed first, ordered by p-value;
# remaining slots are filled by |log2FC| among non-significant hits.
top_log2fc_table <- function(
  diff,
  n = 5,
  p_cutoff = 0.05001,
  fc_cutoff = 0.5,
  prioritize_significant = TRUE
) {
  clean <- diff[is.finite(diff$log2FC) & !is.na(diff$pval), , drop = FALSE]
  if (!"neg_log10_p" %in% names(clean)) {
    clean$neg_log10_p <- -log10(clean$pval)
  }
  clean$significant <- clean$pval <= p_cutoff & abs(clean$log2FC) > fc_cutoff

  cols <- intersect(
    c("cell_type", "mean_high", "mean_low", "log2FC", "pval", "neg_log10_p", "significant"),
    names(clean)
  )

  rank_one_direction <- function(df, positive = TRUE) {
    if (positive) {
      df <- df[df$log2FC > 0, , drop = FALSE]
    } else {
      df <- df[df$log2FC < 0, , drop = FALSE]
    }
    if (nrow(df) == 0) {
      return(df)
    }

    if (prioritize_significant) {
      sig <- df[df$significant, , drop = FALSE]
      non <- df[!df$significant, , drop = FALSE]

      if (nrow(sig) > 0) {
        sig <- sig[order(sig$pval, -abs(sig$log2FC)), , drop = FALSE]
      }
      if (nrow(non) > 0) {
        if (positive) {
          non <- non[order(-non$log2FC), , drop = FALSE]
        } else {
          non <- non[order(non$log2FC), , drop = FALSE]
        }
      }

      df <- rbind(sig, non)
    } else if (positive) {
      df <- df[order(-df$log2FC), , drop = FALSE]
    } else {
      df <- df[order(df$log2FC), , drop = FALSE]
    }

    df <- utils::head(df, n)
    df$direction <- if (positive) "up" else "down"
    df$rank <- seq_len(nrow(df))
    df
  }

  up <- rank_one_direction(clean, positive = TRUE)
  down <- rank_one_direction(clean, positive = FALSE)

  dplyr::bind_rows(up, down) %>%
    dplyr::select("direction", "rank", dplyr::all_of(setdiff(cols, c("direction", "rank"))))
}

export_top_log2fc_table <- function(diff, out_path, n = 5, ...) {
  tbl <- top_log2fc_table(diff, n = n, ...)
  readr::write_csv(tbl, out_path)
  message("Saved table: ", out_path)
  tbl
}

# Expression matrices: samples in columns, genes in rows.
prepare_expr_samples_x_genes <- function(expr_df, sample_col) {
  sample_ids <- expr_df[[sample_col]]
  gene_mat <- as.matrix(expr_df[, setdiff(names(expr_df), sample_col), drop = FALSE])
  mat <- t(gene_mat)
  colnames(mat) <- sample_ids
  mat
}

prepare_expr_genes_x_samples <- function(expr_df, gene_col) {
  genes <- expr_df[[gene_col]]
  mat <- as.matrix(expr_df[, setdiff(names(expr_df), gene_col), drop = FALSE])
  rownames(mat) <- genes
  mat
}

filter_impute_expr <- function(expr_mat, max_na_frac = 0.2) {
  na_frac <- rowMeans(is.na(expr_mat))
  expr_mat <- expr_mat[na_frac <= max_na_frac, , drop = FALSE]
  for (i in seq_len(nrow(expr_mat))) {
    row <- expr_mat[i, ]
    if (any(is.na(row))) {
      row[is.na(row)] <- stats::median(row, na.rm = TRUE)
      expr_mat[i, ] <- row
    }
  }
  expr_mat
}

# Average rows sharing the same gene symbol (e.g. readr ... suffix duplicates).
collapse_duplicate_rownames <- function(expr_mat) {
  base <- sub("\\.\\.\\.[0-9]+$", "", rownames(expr_mat))
  if (!any(duplicated(base))) {
    rownames(expr_mat) <- base
    return(expr_mat)
  }
  split_idx <- split(seq_len(nrow(expr_mat)), base)
  out <- vapply(
    names(split_idx),
    function(g) {
      rows <- expr_mat[split_idx[[g]], , drop = FALSE]
      if (nrow(rows) == 1) rows[1, ] else colMeans(rows)
    },
    numeric(ncol(expr_mat))
  )
  t(out)
}

load_taxonomy_expression_matrix <- function(rlog_path) {
  expr_df <- readr::read_tsv(rlog_path, show_col_types = FALSE)
  names(expr_df)[1] <- "sample_id"
  expr_mat <- prepare_expr_samples_x_genes(expr_df, "sample_id")
  expr_mat <- filter_impute_expr(expr_mat)
  collapse_duplicate_rownames(expr_mat)
}

build_taxonomy_sample_meta <- function(taxonomy_map_path, survival_path) {
  taxonomy_map <- readr::read_tsv(taxonomy_map_path, show_col_types = FALSE)
  survival_df <- readRDS(survival_path)
  taxonomy_map %>%
    dplyr::transmute(
      sample_id = .data$Patient,
      patient_id = .data$Code
    ) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      survival_df %>% dplyr::select(.data$patient_id, .data$risk_grp, .data$rsf_risk),
      by = "patient_id"
    )
}

# Limma: High vs Low risk on log-normalised expression (genes x samples).
# Uses Bioconductor limma when available; otherwise a linear-model + BH fallback.
run_limma_high_vs_low <- function(expr_mat, sample_meta, output_dir, file_prefix, min_group = 3) {
  meta <- sample_meta %>%
    dplyr::filter(.data$risk_grp %in% c("Low", "High")) %>%
    dplyr::filter(.data$sample_id %in% colnames(expr_mat))

  n_low <- sum(meta$risk_grp == "Low")
  n_high <- sum(meta$risk_grp == "High")
  if (n_low < min_group || n_high < min_group) {
    stop(
      "Insufficient samples for ", file_prefix, ": Low=", n_low, ", High=", n_high,
      " (need >= ", min_group, " per group)"
    )
  }

  expr <- expr_mat[, meta$sample_id, drop = FALSE]
  risk_grp <- factor(
    meta$risk_grp[match(colnames(expr), meta$sample_id)],
    levels = c("Low", "High")
  )

  if (requireNamespace("limma", quietly = TRUE)) {
    design <- stats::model.matrix(~ 0 + risk_grp)
    colnames(design) <- make.names(colnames(design))
    contrast <- limma::makeContrasts(
      HighvsLow = risk_grpHigh - risk_grpLow,
      levels = design
    )
    fit <- limma::lmFit(expr, design)
    fit2 <- limma::contrasts.fit(fit, contrast)
    fit2 <- limma::eBayes(fit2)
    results <- limma::topTable(
      fit2,
      coef = "HighvsLow",
      number = Inf,
      adjust.method = "BH",
      sort.by = "P"
    )
    results$gene <- rownames(results)
  } else {
    message("Package 'limma' not installed — using linear-model t-test fallback with BH adjustment.")
    results <- de_two_group_fallback(expr, risk_grp)
  }

  all_path <- file.path(output_dir, paste0(file_prefix, "_limma_high_vs_low_all_genes.csv"))
  readr::write_csv(results, all_path)
  message("Saved: ", all_path, " (n=", n_low, " Low, ", n_high, " High)")

  up <- results %>%
    dplyr::filter(.data$logFC > 0) %>%
    dplyr::arrange(dplyr::desc(.data$logFC)) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::mutate(direction = "up", rank = dplyr::row_number())

  down <- results %>%
    dplyr::filter(.data$logFC < 0) %>%
    dplyr::arrange(.data$logFC) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::mutate(direction = "down", rank = dplyr::row_number())

  top10 <- dplyr::bind_rows(up, down) %>%
    dplyr::select(.data$direction, .data$rank, dplyr::everything())

  top_path <- file.path(output_dir, paste0(file_prefix, "_limma_top5_logfc_up_down.csv"))
  readr::write_csv(top10, top_path)
  message("Saved: ", top_path)

  list(results = results, top = top10, n_low = n_low, n_high = n_high)
}

de_two_group_fallback <- function(expr_mat, risk_grp) {
  low_idx  <- which(risk_grp == "Low")
  high_idx <- which(risk_grp == "High")
  low_mat  <- expr_mat[, low_idx, drop = FALSE]
  high_mat <- expr_mat[, high_idx, drop = FALSE]

  logFC <- rowMeans(high_mat) - rowMeans(low_mat)
  AveExpr <- rowMeans(expr_mat)

  pvals <- vapply(seq_len(nrow(expr_mat)), function(i) {
    stats::t.test(high_mat[i, ], low_mat[i, ], var.equal = FALSE)$p.value
  }, numeric(1))

  t_vals <- vapply(seq_len(nrow(expr_mat)), function(i) {
    unname(stats::t.test(high_mat[i, ], low_mat[i, ], var.equal = FALSE)$statistic)
  }, numeric(1))

  adj_p <- stats::p.adjust(pvals, method = "BH")
  results <- data.frame(
    logFC = logFC,
    AveExpr = AveExpr,
    t = t_vals,
    P.Value = pvals,
    adj.P.Val = adj_p,
    B = NA_real_,
    gene = rownames(expr_mat),
    row.names = rownames(expr_mat),
    stringsAsFactors = FALSE
  )
  results[order(results$P.Value), , drop = FALSE]
}

build_colossus_id_map <- function(master_df, sample_ids) {
  strip_fp <- function(x) stringr::str_replace(x, "-FP.*$", "")
  master_col_map <- dplyr::bind_rows(
    master_df %>%
      dplyr::filter(!is.na(.data$patient_id_g), .data$patient_id_g != "") %>%
      dplyr::transmute(patient_id = .data$patient_id_g, col_id = strip_fp(.data$COLOSSUS_ID)),
    master_df %>%
      dplyr::filter(!is.na(.data$patient_id_g), .data$patient_id_g != "") %>%
      dplyr::transmute(patient_id = .data$patient_id_g, col_id = strip_fp(.data$Old_COLOSSUS_ID)),
    master_df %>%
      dplyr::filter(!is.na(.data$patient_id_g), .data$patient_id_g != "") %>%
      dplyr::transmute(patient_id = .data$patient_id_g, col_id = strip_fp(.data$Sample_ID)),
    master_df %>%
      dplyr::filter(!is.na(.data$patient_id_g), .data$patient_id_g != "") %>%
      dplyr::transmute(patient_id = .data$patient_id_g, col_id = strip_fp(.data$Alternative_ID))
  ) %>%
    dplyr::filter(!is.na(.data$col_id), .data$col_id != "", .data$col_id != "NA") %>%
    dplyr::distinct(.data$col_id, .keep_all = TRUE)

  data.frame(sample = sample_ids, stringsAsFactors = FALSE) %>%
    dplyr::left_join(master_col_map, by = c("sample" = "col_id"))
}

# Y-axis limit for limma FDR volcanoes: fit the bulk of points; only stretch to the
# FDR cutoff when genes actually pass FDR (avoids empty space for null cohorts).
limma_volcano_ymax <- function(y, adj_p, label_y = numeric(0), fdr_cutoff = 0.05) {
  y_finite <- y[is.finite(y) & y >= 0]
  if (length(y_finite) == 0) {
    return(-log10(fdr_cutoff) * 1.05)
  }

  y_body <- stats::quantile(y_finite, 0.99, na.rm = TRUE) * 1.2
  y_labels <- if (length(label_y)) {
    max(label_y[is.finite(label_y)], na.rm = TRUE) * 1.12
  } else {
    0
  }
  has_fdr <- !is.null(adj_p) && any(is.finite(adj_p) & adj_p <= fdr_cutoff, na.rm = TRUE)
  fdr_y <- -log10(fdr_cutoff)

  ymax <- max(y_body, y_labels, 0.35, na.rm = TRUE)
  if (has_fdr) {
    sig_y <- y[is.finite(y) & is.finite(adj_p) & adj_p <= fdr_cutoff]
    ymax <- max(
      ymax,
      fdr_y * 1.08,
      if (length(sig_y)) max(sig_y, na.rm = TRUE) * 1.08 else fdr_y * 1.08,
      na.rm = TRUE
    )
  }
  ymax
}

prepare_limma_volcano_df <- function(
  results,
  fc_cutoff = 1,
  top_labels = 12,
  fc_display_cap = 5,
  fdr_cutoff = 0.05
) {
  df <- results %>%
    dplyr::filter(is.finite(.data$logFC), !is.na(.data$P.Value)) %>%
    dplyr::mutate(
      adj.P.Val = if ("adj.P.Val" %in% names(results)) {
        .data$adj.P.Val
      } else {
        stats::p.adjust(.data$P.Value, method = "BH")
      },
      neg_log10_adjP = -log10(pmax(.data$adj.P.Val, .Machine$double.xmin)),
      fc_class = dplyr::case_when(
        .data$logFC >= fc_cutoff ~ "Higher in High risk",
        .data$logFC <= -fc_cutoff ~ "Higher in Low risk",
        TRUE ~ "Other"
      )
    )

  label_df <- df %>%
    dplyr::filter(abs(.data$logFC) <= fc_display_cap) %>%
    dplyr::arrange(dplyr::desc(abs(.data$logFC))) %>%
    dplyr::slice_head(n = top_labels)

  has_fdr_hits <- any(is.finite(df$adj.P.Val) & df$adj.P.Val <= fdr_cutoff, na.rm = TRUE)

  list(
    df = df,
    label_df = label_df,
    has_fdr_hits = has_fdr_hits,
    fdr_y = -log10(fdr_cutoff)
  )
}

compute_limma_volcano_limits <- function(
  results,
  fc_cutoff = 1,
  top_labels = 12,
  fc_display_cap = 5,
  fdr_cutoff = 0.05
) {
  prep <- prepare_limma_volcano_df(
    results,
    fc_cutoff = fc_cutoff,
    top_labels = top_labels,
    fc_display_cap = fc_display_cap,
    fdr_cutoff = fdr_cutoff
  )
  df <- prep$df
  label_df <- prep$label_df

  fc_body <- df$logFC[is.finite(df$logFC) & abs(df$logFC) <= fc_display_cap]
  label_fc <- label_df$logFC
  sig_fc <- df$logFC[!is.na(df$adj.P.Val) & df$adj.P.Val <= fdr_cutoff]
  sig_body <- sig_fc[is.finite(sig_fc) & abs(sig_fc) <= fc_display_cap]
  sig_q <- if (length(sig_body) >= 1) {
    stats::quantile(abs(sig_body), 0.95, na.rm = TRUE) * 1.1
  } else {
    0
  }
  xmax <- max(
    fc_cutoff * 1.5,
    2,
    stats::quantile(abs(fc_body), 0.95, na.rm = TRUE) * 1.15,
    sig_q,
    if (length(label_fc)) max(abs(label_fc), na.rm = TRUE) * 1.12 else 0,
    na.rm = TRUE
  )

  list(
    xmax = min(fc_display_cap, xmax),
    ymax = limma_volcano_ymax(
      df$neg_log10_adjP,
      df$adj.P.Val,
      label_y = label_df$neg_log10_adjP,
      fdr_cutoff = fdr_cutoff
    ),
    has_fdr_hits = prep$has_fdr_hits
  )
}

# Volcano + heatmap for limma DE — fold-change first (suited to small High-risk groups).
# Y-axis uses BH FDR-adjusted p-values (adj.P.Val from limma).
export_limma_figures <- function(
  results,
  expr_mat,
  sample_meta,
  output_dir,
  file_prefix,
  cohort_title,
  fc_cutoff = 1,
  top_labels = 12,
  heatmap_genes = 25,
  fc_display_cap = 5,
  fdr_cutoff = 0.05,
  xmax_override = NULL,
  ymax_override = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for limma figures.")
  }

  prep <- prepare_limma_volcano_df(
    results,
    fc_cutoff = fc_cutoff,
    top_labels = top_labels,
    fc_display_cap = fc_display_cap,
    fdr_cutoff = fdr_cutoff
  )
  df <- prep$df
  label_df <- prep$label_df
  has_fdr_hits <- prep$has_fdr_hits
  fdr_y <- prep$fdr_y

  limits <- compute_limma_volcano_limits(
    results,
    fc_cutoff = fc_cutoff,
    top_labels = top_labels,
    fc_display_cap = fc_display_cap,
    fdr_cutoff = fdr_cutoff
  )
  xmax <- if (is.null(xmax_override)) limits$xmax else xmax_override
  ymax <- if (is.null(ymax_override)) limits$ymax else ymax_override

  line_caption <- if (has_fdr_hits) {
    paste0("Dashed line: BH FDR = ", fdr_cutoff)
  } else {
    paste0("No genes at BH FDR ≤ ", fdr_cutoff, " (FDR line omitted; y-axis scaled to data)")
  }

  meta_n <- sample_meta %>%
    dplyr::filter(.data$risk_grp %in% c("Low", "High"), .data$sample_id %in% colnames(expr_mat))
  n_low <- sum(meta_n$risk_grp == "Low")
  n_high <- sum(meta_n$risk_grp == "High")

  vol <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$logFC, y = .data$neg_log10_adjP, colour = .data$fc_class, size = abs(.data$logFC))
  ) +
    ggplot2::geom_point(alpha = 0.55) +
    ggplot2::scale_size_continuous(range = c(0.8, 4), guide = "none") +
    ggplot2::scale_colour_manual(
      values = c(
        "Higher in High risk" = "#B2182B",
        "Higher in Low risk"  = "#2166AC",
        "Other"               = "grey75"
      ),
      name = NULL
    ) +
    ggplot2::geom_vline(
      xintercept = c(-fc_cutoff, fc_cutoff),
      linetype = "dashed",
      colour = "black",
      linewidth = 0.4
    )

  if (has_fdr_hits) {
    vol <- vol + ggplot2::geom_hline(
      yintercept = fdr_y,
      linetype = "dashed",
      colour = "black",
      linewidth = 0.4
    )
  }

  vol <- vol +
    ggplot2::coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax), clip = "off", expand = FALSE) +
    ggplot2::labs(
      title = paste0("Differential expression: High vs Low DFS risk"),
      subtitle = paste0(
        cohort_title, " | n = ", n_low + n_high, " (Low ", n_low, ", High ", n_high,
        ") | Labels: top ", top_labels, " genes by |logFC|"
      ),
      caption = line_caption,
      x = "log2 fold change (High / Low)",
      y = expression("-log"[10] * "(FDR-adjusted p-value)")
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      plot.caption = ggplot2::element_text(size = 9, colour = "grey30", hjust = 0),
      plot.margin = ggplot2::margin(15, 50, 10, 15),
      legend.position = "bottom"
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    vol <- vol + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = .data$gene),
      size = 3,
      max.overlaps = Inf,
      segment.size = 0.2,
      colour = "black",
      show.legend = FALSE
    )
  } else {
    vol <- vol + ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(label = .data$gene),
      size = 3,
      hjust = -0.1,
      colour = "black",
      show.legend = FALSE
    )
  }

  vol_png <- file.path(output_dir, paste0(file_prefix, "_limma_volcano_logfc.png"))
  vol_pdf <- file.path(output_dir, paste0(file_prefix, "_limma_volcano_logfc.pdf"))
  ggplot2::ggsave(vol_png, vol, width = 12, height = 7, dpi = 300, bg = "white")
  ggplot2::ggsave(vol_pdf, vol, width = 12, height = 7, bg = "white")
  message("Saved: ", vol_png)

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    warning("Package 'pheatmap' not installed — skipping heatmap.")
    return(invisible(list(volcano = vol)))
  }

  meta <- sample_meta %>%
    dplyr::filter(.data$risk_grp %in% c("Low", "High"), .data$sample_id %in% colnames(expr_mat))
  if ("rsf_risk" %in% names(meta)) {
    meta <- meta %>% dplyr::arrange(.data$risk_grp, .data$rsf_risk)
  } else {
    meta <- meta %>% dplyr::arrange(.data$risk_grp, .data$sample_id)
  }

  top_genes <- results %>%
    dplyr::filter(is.finite(.data$logFC), abs(.data$logFC) <= fc_display_cap) %>%
    dplyr::arrange(dplyr::desc(abs(.data$logFC))) %>%
    dplyr::slice_head(n = heatmap_genes) %>%
    dplyr::pull(.data$gene)

  if (length(top_genes) < 2) {
    warning("Fewer than 2 genes for heatmap in ", file_prefix)
    return(invisible(list(volcano = vol)))
  }

  hm_mat <- expr_mat[top_genes, meta$sample_id, drop = FALSE]
  hm_z <- t(scale(t(hm_mat)))
  hm_z[!is.finite(hm_z)] <- 0

  ann_col <- data.frame(
    risk_grp = factor(meta$risk_grp, levels = c("Low", "High")),
    row.names = meta$sample_id
  )
  ann_colors <- list(risk_grp = c(Low = "#004B23", High = "#800000"))
  hm_colors <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  n_low_hm <- sum(meta$risk_grp == "Low")
  gaps_col <- if (n_low_hm > 0L && n_low_hm < ncol(hm_z)) n_low_hm else NULL

  hm_title <- paste0(
    cohort_title, ": top ", length(top_genes), " genes by |logFC|\n",
    "n = ", nrow(meta), " (Low ", n_low_hm, ", High ", sum(meta$risk_grp == "High"), ")"
  )

  draw_hm <- function() {
    pheatmap::pheatmap(
      hm_z,
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      annotation_col = ann_col,
      annotation_colors = ann_colors,
      color = hm_colors,
      main = hm_title,
      fontsize_row = 9,
      fontsize_col = 7,
      border_color = NA,
      gaps_col = gaps_col,
      show_colnames = FALSE
    )
  }

  hm_png <- file.path(output_dir, paste0(file_prefix, "_limma_heatmap_top_logfc.png"))
  hm_pdf <- file.path(output_dir, paste0(file_prefix, "_limma_heatmap_top_logfc.pdf"))
  grDevices::png(hm_png, width = 10, height = max(5, length(top_genes) * 0.22), units = "in", res = 300)
  draw_hm()
  grDevices::dev.off()
  grDevices::pdf(hm_pdf, width = 10, height = max(5, length(top_genes) * 0.22))
  draw_hm()
  grDevices::dev.off()
  message("Saved: ", hm_png)

  invisible(list(volcano = vol, heatmap_genes = top_genes))
}
