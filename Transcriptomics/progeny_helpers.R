# PROGENy pathway activity via decoupleR, with optional class-specific (risk group) normalization.
source("analysis_helpers.R")

resolve_progeny_paths <- function() {
  if (file.exists("stage2/normalized_genes.csv")) {
    list(
      base = ".",
      proteomics = "../proteomics",
      spatial = "../spatial",
      integrated = "../Integrated "
    )
  } else if (file.exists("Transcriptomics/stage2/normalized_genes.csv")) {
    list(
      base = "Transcriptomics",
      proteomics = "proteomics",
      spatial = "spatial",
      integrated = "Integrated "
    )
  } else {
    stop("Run from Transcriptomics/ or its parent directory.")
  }
}

load_stage2_expression_meta <- function(p, survival_df) {
  expr_df <- readr::read_csv(
    file.path(p$base, "stage2/normalized_genes.csv"),
    show_col_types = FALSE
  )
  clinical <- readr::read_csv(
    file.path(p$base, "stage2/clinical.csv"),
    show_col_types = FALSE
  )
  mapping <- readr::read_csv(
    file.path(p$spatial, "master_patient_mapping.csv"),
    show_col_types = FALSE
  )
  expr_mat <- filter_impute_expr(
    prepare_expr_samples_x_genes(expr_df, "patient_id")
  )
  sample_ids <- colnames(expr_mat)
  sample_meta <- data.frame(
    sample_id = sample_ids,
    patient_id = mapping$patient_id_g[
      match(clinical$r_code[match(sample_ids, clinical$patient_id)], mapping$r_code)
    ],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(
      survival_df %>% dplyr::select(patient_id, risk_grp, rsf_risk),
      by = "patient_id"
    ) %>%
    dplyr::mutate(cohort = "Stage2")
  list(expr = expr_mat, meta = sample_meta)
}

load_retrospective_expression_meta <- function(p, survival_df) {
  expr_df <- readr::read_tsv(
    file.path(p$base, "retrospective/Leuven_rlog_values.txt"),
    show_col_types = FALSE
  )
  expr_mat <- filter_impute_expr(
    prepare_expr_genes_x_samples(expr_df, "Geneid")
  )
  sample_ids <- colnames(expr_mat)
  sample_meta <- data.frame(
    sample_id = sample_ids,
    patient_id = sample_ids,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(
      survival_df %>% dplyr::select(patient_id, risk_grp, rsf_risk),
      by = "patient_id"
    ) %>%
    dplyr::mutate(cohort = "Retrospective")
  list(expr = expr_mat, meta = sample_meta)
}

load_colossus_expression_meta <- function(p, survival_df) {
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required for Colossus sample mapping.")
  }
  expr_df <- readr::read_csv(
    file.path(p$base, "colossus/rna.csv"),
    show_col_types = FALSE
  )
  master_df <- readr::read_csv(
    file.path(p$integrated, "All_patinet_IDS_concatonated_masterdoc.csv"),
    show_col_types = FALSE
  )
  expr_mat <- filter_impute_expr(
    prepare_expr_samples_x_genes(expr_df, "patient_id")
  )
  id_map <- build_colossus_id_map(master_df, colnames(expr_mat)) %>%
    dplyr::rename(sample_id = sample)
  sample_meta <- id_map %>%
    dplyr::left_join(
      survival_df %>% dplyr::select(patient_id, risk_grp, rsf_risk),
      by = "patient_id"
    ) %>%
    dplyr::mutate(cohort = "Colossus")
  list(expr = expr_mat, meta = sample_meta)
}

load_taxonomy_expression_meta <- function(p, survival_path) {
  rlog_path <- file.path(p$base, "Taxonomy_calls/Taxonomy_manuela_with_rna_rlog.txt")
  map_path <- file.path(
    p$base, "Taxonomy_calls",
    "Manuela_and_Belfast_RNA_classifications_CMS_CRIS.txt"
  )
  expr_mat <- load_taxonomy_expression_matrix(rlog_path)
  sample_meta <- build_taxonomy_sample_meta(map_path, survival_path) %>%
    dplyr::mutate(cohort = "Taxonomy")
  list(expr = expr_mat, meta = sample_meta)
}

combine_legacy_expression <- function(cohort_data) {
  genes <- Reduce(intersect, lapply(cohort_data, function(x) rownames(x$expr)))
  if (length(genes) == 0) {
    stop("No shared genes across legacy cohort expression matrices.")
  }

  expr_parts <- list()
  meta_parts <- list()
  for (nm in names(cohort_data)) {
    ex <- cohort_data[[nm]]$expr[genes, , drop = FALSE]
    meta <- cohort_data[[nm]]$meta %>%
      dplyr::filter(.data$sample_id %in% colnames(ex)) %>%
      dplyr::mutate(
        sample_uid = paste0(.data$cohort, "::", .data$sample_id)
      )
    ex <- ex[, meta$sample_id, drop = FALSE]
    colnames(ex) <- meta$sample_uid
    meta$sample_id <- meta$sample_uid
    expr_parts[[nm]] <- ex
    meta_parts[[nm]] <- meta
  }

  meta <- dplyr::bind_rows(meta_parts) %>%
    dplyr::filter(.data$risk_grp %in% c("Low", "High"))
  expr <- do.call(cbind, expr_parts)
  expr <- expr[, meta$sample_id, drop = FALSE]

  list(
    expr = expr,
    meta = meta,
    shared_genes = genes
  )
}

class_specific_normalize_expr <- function(
  expr_mat,
  sample_meta,
  class_col = "risk_grp",
  classes = c("Low", "High"),
  method = c("zscore", "center")
) {
  method <- match.arg(method)
  sample_meta <- sample_meta %>%
    dplyr::filter(
      .data$sample_id %in% colnames(expr_mat),
      .data[[class_col]] %in% classes
    )
  if (nrow(sample_meta) == 0) {
    stop("No samples with valid risk groups for class-specific normalization.")
  }

  out <- expr_mat[, sample_meta$sample_id, drop = FALSE]
  for (cl in classes) {
    ids <- sample_meta$sample_id[sample_meta[[class_col]] == cl]
    if (length(ids) < 2) {
      warning("Fewer than 2 samples in class ", cl, "; skipping normalization for that class.")
      next
    }
    block <- out[, ids, drop = FALSE]
    if (method == "zscore") {
      scaled <- t(scale(t(block)))
      scaled[!is.finite(scaled)] <- 0
      out[, ids] <- scaled
    } else {
      row_means <- rowMeans(block, na.rm = TRUE)
      out[, ids] <- block - row_means
    }
  }

  list(expr = out, sample_meta = sample_meta)
}

run_progeny_mlm <- function(
  expr_mat,
  organism = "human",
  top = 500L,
  minsize = 5L
) {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Package 'decoupleR' is required for PROGENy analysis.")
  }
  net <- decoupleR::get_progeny(organism = organism, top = top)
  acts <- decoupleR::run_mlm(
    mat = expr_mat,
    net = net,
    .source = "source",
    .target = "target",
    .mor = "weight",
    minsize = minsize
  )
  list(activities = acts, network = net)
}

progeny_activities_to_long <- function(activities, sample_meta) {
  acts <- activities %>%
    dplyr::rename(
      sample_id = condition,
      pathway = source,
      score = score
    ) %>%
    dplyr::inner_join(
      sample_meta %>%
        dplyr::select(sample_id, risk_grp, cohort, patient_id),
      by = "sample_id"
    )
  acts
}

progeny_pathway_wilcox <- function(acts_long) {
  out <- acts_long %>%
    dplyr::group_by(pathway) %>%
    dplyr::summarise(
      mean_high = mean(.data$score[.data$risk_grp == "High"], na.rm = TRUE),
      mean_low = mean(.data$score[.data$risk_grp == "Low"], na.rm = TRUE),
      log2FC = .data$mean_high - .data$mean_low,
      pval = tryCatch(
        stats::wilcox.test(.data$score ~ .data$risk_grp)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    )
  out$neg_log10_p <- -log10(out$pval)
  out
}

export_progeny_pathway_volcano <- function(
  diff,
  title,
  subtitle,
  out_path
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for PROGENy volcano plots.")
  }
  plot_df <- diff[is.finite(diff$log2FC) & !is.na(diff$pval), , drop = FALSE]
  plot_df$sig <- plot_df$pval <= 0.05 & abs(plot_df$log2FC) > 0.25
  xmax <- min(2, max(0.5, stats::quantile(abs(plot_df$log2FC), 0.98, na.rm = TRUE) * 1.1))
  ymax <- max(1.5, max(-log10(plot_df$pval), na.rm = TRUE) + 0.5)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = log2FC, y = -log10(pval), colour = sig, label = pathway)
  ) +
    ggplot2::geom_point(size = 2.8, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#d73027", `FALSE` = "grey55")) +
    ggplot2::coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax)) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Pathway score difference (High − Low)",
      y = expression(-log[10] * "(Wilcoxon p-value)"),
      colour = "p <= 0.05 & |diff| > 0.25"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = pathway),
      size = 3,
      max.overlaps = 20,
      show.legend = FALSE
    )
  }

  ggplot2::ggsave(out_path, p, width = 10, height = 7, dpi = 200, bg = "white")
  message("Saved: ", out_path)
  invisible(out_path)
}

export_progeny_heatmap <- function(
  acts_long,
  out_path,
  title = "PROGENy pathway activities"
) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("pheatmap is required for PROGENy heatmaps.")
  }
  mat <- acts_long %>%
    dplyr::select(sample_id, pathway, score) %>%
    tidyr::pivot_wider(names_from = pathway, values_from = score) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()
  mat <- t(scale(t(mat)))
  mat[!is.finite(mat)] <- 0

  ann <- acts_long %>%
    dplyr::distinct(sample_id, risk_grp, cohort) %>%
    tibble::column_to_rownames("sample_id")

  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors.use <- grDevices::colorRampPalette(colors = colors)(100)
  breaks <- c(
    seq(-2, 0, length.out = ceiling(100 / 2) + 1),
    seq(0.05, 2, length.out = floor(100 / 2))
  )

  grDevices::png(out_path, width = 1200, height = 900, res = 120)
  pheatmap::pheatmap(
    mat,
    annotation_row = ann[, c("risk_grp", "cohort"), drop = FALSE],
    color = colors.use,
    breaks = breaks,
    border_color = NA,
    main = title,
    fontsize_row = 6,
    fontsize_col = 9,
    cluster_cols = TRUE
  )
  grDevices::dev.off()
  message("Saved: ", out_path)
  invisible(out_path)
}

run_progeny_risk_analysis <- function(
  expr_mat,
  sample_meta,
  output_dir,
  file_prefix,
  title,
  apply_class_norm = TRUE,
  class_norm_method = "zscore"
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  expr_use <- expr_mat
  norm_note <- "No class-specific normalization (raw log-expression)."
  if (apply_class_norm) {
    norm_out <- class_specific_normalize_expr(
      expr_mat,
      sample_meta,
      method = class_norm_method
    )
    expr_use <- norm_out$expr
    sample_meta <- norm_out$sample_meta
    norm_note <- paste0(
      "Class-specific ", class_norm_method,
      " normalization within Low and High risk groups (per gene)."
    )
  } else {
    sample_meta <- sample_meta %>%
      dplyr::filter(
        .data$sample_id %in% colnames(expr_mat),
        .data$risk_grp %in% c("Low", "High")
      )
    expr_use <- expr_mat[, sample_meta$sample_id, drop = FALSE]
  }

  expr_df_out <- as.data.frame(expr_use, check.names = FALSE)
  expr_df_out$gene <- rownames(expr_use)
  expr_df_out <- expr_df_out[, c("gene", setdiff(names(expr_df_out), "gene")), drop = FALSE]
  if (nrow(expr_df_out) > 5000L) {
    message("Skipping full expression export (>5000 genes).")
  } else {
    readr::write_csv(
      expr_df_out,
      file.path(output_dir, paste0(file_prefix, "_expression_classNorm.csv"))
    )
  }

  progeny <- run_progeny_mlm(expr_use)
  acts_long <- progeny_activities_to_long(progeny$activities, sample_meta)

  readr::write_csv(
    acts_long,
    file.path(output_dir, paste0(file_prefix, "_progeny_activities_long.csv"))
  )

  acts_wide <- acts_long %>%
    dplyr::select(sample_id, risk_grp, cohort, pathway, score) %>%
    tidyr::pivot_wider(names_from = pathway, values_from = score)
  readr::write_csv(
    acts_wide,
    file.path(output_dir, paste0(file_prefix, "_progeny_activities_wide.csv"))
  )

  diff <- progeny_pathway_wilcox(acts_long)
  readr::write_csv(
    diff,
    file.path(output_dir, paste0(file_prefix, "_progeny_wilcox_pathways.csv"))
  )

  n_low <- sum(sample_meta$risk_grp == "Low")
  n_high <- sum(sample_meta$risk_grp == "High")
  export_progeny_pathway_volcano(
    diff,
    paste0("PROGENy: High vs Low risk (", title, ")"),
    paste0(
      norm_note, " | n = ", nrow(sample_meta),
      " (Low ", n_low, ", High ", n_high, ")"
    ),
    file.path(output_dir, paste0(file_prefix, "_progeny_volcano.png"))
  )
  export_progeny_heatmap(
    acts_long,
    file.path(output_dir, paste0(file_prefix, "_progeny_heatmap.png")),
    title = paste0("PROGENy — ", title)
  )

  summary_lines <- c(
    paste0("PROGENy analysis: ", title),
    norm_note,
    paste0("Samples: ", nrow(sample_meta), " (Low=", n_low, ", High=", n_high, ")"),
    paste0("Genes: ", nrow(expr_use)),
    paste0("Pathways: ", dplyr::n_distinct(acts_long$pathway)),
    paste0(
      "Significant pathways (p <= 0.05, |score diff| > 0.25): ",
      sum(diff$pval <= 0.05 & abs(diff$log2FC) > 0.25, na.rm = TRUE),
      " / ", nrow(diff)
    ),
    "",
    "Top pathways by p-value:",
    paste0(
      "  ",
      head(diff$pathway[order(diff$pval)], 5),
      ": p=",
      signif(head(diff$pval[order(diff$pval)], 5), 3),
      collapse = "\n"
    )
  )
  writeLines(summary_lines, file.path(output_dir, paste0(file_prefix, "_summary.txt")))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(
    expr = expr_use,
    sample_meta = sample_meta,
    acts_long = acts_long,
    diff = diff
  ))
}

run_progeny_legacy_and_taxonomy <- function(
  legacy_output_dir = NULL,
  taxonomy_output_dir = NULL
) {
  p <- resolve_progeny_paths()
  if (is.null(legacy_output_dir)) {
    legacy_output_dir <- file.path(p$base, "analysis_output", "legacy_combined")
  }
  if (is.null(taxonomy_output_dir)) {
    taxonomy_output_dir <- file.path(p$base, "analysis_output", "taxonomy")
  }

  survival_path <- file.path(p$proteomics, "survival_df_with_risk.rds")
  survival_df <- readRDS(survival_path)

  legacy_parts <- list(
    Stage2 = load_stage2_expression_meta(p, survival_df),
    Retrospective = load_retrospective_expression_meta(p, survival_df),
    Colossus = load_colossus_expression_meta(p, survival_df)
  )
  legacy <- combine_legacy_expression(legacy_parts)
  message(
    "Legacy combined expression: ",
    length(legacy$shared_genes), " genes, ",
    ncol(legacy$expr), " samples"
  )

  legacy_out <- run_progeny_risk_analysis(
    legacy$expr,
    legacy$meta,
    legacy_output_dir,
    file_prefix = "Legacy_classNorm",
    title = "Stage 2 + Retrospective + Colossus",
    apply_class_norm = TRUE
  )

  tax <- load_taxonomy_expression_meta(p, survival_path)
  tax_meta <- tax$meta %>%
    dplyr::filter(.data$sample_id %in% colnames(tax$expr))
  tax_out <- run_progeny_risk_analysis(
    tax$expr[, tax_meta$sample_id, drop = FALSE],
    tax_meta,
    taxonomy_output_dir,
    file_prefix = "Taxonomy_classNorm",
    title = "Taxonomy RNA-Seq",
    apply_class_norm = TRUE
  )

  invisible(list(legacy = legacy_out, taxonomy = tax_out))
}
