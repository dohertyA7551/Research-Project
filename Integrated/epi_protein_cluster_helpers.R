# Epithelial proteomics clustering from limma DE (High vs Low RSF; G1 vs G4 TMA groups).
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

source("panel_predicted_observed_network_helpers.R")

DEFAULT_N_CLUSTERS <- 6L
DEFAULT_TOP_MARKERS <- 10L
DEFAULT_TMA_GROUPS_RDS <- file.path(PROTEOMICS_DIR, "tma_composition_patient_groups.rds")
G1_GROUP <- "G1_all_low"
G4_GROUP <- "G4_high_3plus"

IMMUNE_COMPARTMENTS <- data.frame(
  suffix = c("helperT", "cytT", "regT"),
  label = c("Helper T cells", "Cytotoxic T cells", "Regulatory T cells"),
  stringsAsFactors = FALSE
)

compartment_label <- function(compartment_suffix) {
  idx <- match(compartment_suffix, IMMUNE_COMPARTMENTS$suffix)
  if (!is.na(idx)) {
    return(IMMUNE_COMPARTMENTS$label[idx])
  }
  compartment_suffix
}

cluster_markers_by_profile <- function(mat, n_clusters = DEFAULT_N_CLUSTERS) {
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("stats package required.")
  }
  mat <- mat[rowSums(is.finite(mat)) > 0, , drop = FALSE]
  if (nrow(mat) < 3L || ncol(mat) < 4L) {
    stop("Need >=3 markers and >=4 patients for clustering.")
  }
  mat_z <- t(scale(t(mat)))
  mat_z[!is.finite(mat_z)] <- 0
  cor_mat <- stats::cor(t(mat_z), method = "spearman", use = "pairwise.complete.obs")
  cor_mat[!is.finite(cor_mat)] <- 0
  diag(cor_mat) <- 1
  d <- stats::as.dist(1 - cor_mat)
  hc <- stats::hclust(d, method = "average")
  k <- min(as.integer(n_clusters), nrow(mat))
  clusters <- stats::cutree(hc, k = k)
  list(hclust = hc, clusters = clusters, cor_mat = cor_mat, mat_z = mat_z)
}

load_epi_protein_with_tma_groups <- function(
  groups_rds = DEFAULT_TMA_GROUPS_RDS,
  compartment_suffix = DEFAULT_COMPARTMENT_SUFFIX
) {
  prot <- load_epi_protein_matrix(compartment_suffix = compartment_suffix)
  if (!file.exists(groups_rds)) {
    stop("Missing TMA composition groups: ", groups_rds)
  }
  patient_grp <- readRDS(groups_rds)$patient_groups %>%
    dplyr::select(.data$patient_id, .data$tma_comp_grp)
  prot$meta <- prot$meta %>%
    dplyr::left_join(patient_grp, by = "patient_id")
  prot
}

subset_protein_two_groups <- function(
  prot,
  group_col,
  group_low,
  group_high,
  min_patients = DEFAULT_MIN_PATIENTS
) {
  meta <- prot$meta %>%
    dplyr::filter(.data[[group_col]] %in% c(group_low, group_high)) %>%
    dplyr::mutate(
      contrast_grp = dplyr::if_else(
        .data[[group_col]] == group_high,
        "High",
        "Low"
      )
    )
  pts <- meta$patient_id
  if (length(pts) < min_patients) {
    stop(
      "Too few patients for ", group_low, " vs ", group_high, ": ", length(pts)
    )
  }
  list(
    matrix = prot$matrix[, pts, drop = FALSE],
    meta = meta,
    markers = rownames(prot$matrix),
    patients = pts,
    group_col = group_col,
    group_low = group_low,
    group_high = group_high
  )
}

run_limma_two_group_matrix <- function(
  expr_mat,
  meta,
  group_col,
  group_low,
  group_high,
  feature_col = "marker"
) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required.")
  }
  meta <- meta %>%
    dplyr::filter(.data[[group_col]] %in% c(group_low, group_high)) %>%
    dplyr::filter(.data$patient_id %in% colnames(expr_mat))
  n_low <- sum(meta[[group_col]] == group_low)
  n_high <- sum(meta[[group_col]] == group_high)
  if (n_low < 3L || n_high < 3L) {
    stop(
      "Need >=3 patients per group; ",
      group_low, "=", n_low, ", ", group_high, "=", n_high
    )
  }
  expr <- expr_mat[, meta$patient_id, drop = FALSE]
  grp <- factor(
    meta[[group_col]][match(colnames(expr), meta$patient_id)],
    levels = c(group_low, group_high)
  )
  design <- stats::model.matrix(~ grp)
  fit <- limma::lmFit(expr, design)
  fit2 <- limma::eBayes(fit)
  coef_name <- paste0("grp", group_high)
  results <- limma::topTable(
    fit2,
    coef = coef_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  results[[feature_col]] <- rownames(results)
  list(
    results = results,
    n_low = n_low,
    n_high = n_high,
    group_low = group_low,
    group_high = group_high
  )
}

summarize_marker_clusters <- function(
  clusters,
  de_df,
  top_n = DEFAULT_TOP_MARKERS,
  up_direction = "up_in_high",
  down_direction = "down_in_high"
) {
  cl_df <- tibble::tibble(
    marker = names(clusters),
    cluster = unname(clusters)
  ) %>%
    dplyr::left_join(de_df, by = "marker") %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::summarise(
      n_markers = dplyr::n(),
      mean_logFC = mean(.data$logFC, na.rm = TRUE),
      median_logFC = median(.data$logFC, na.rm = TRUE),
      min_p = min(.data$P.Value, na.rm = TRUE),
      n_up = sum(.data$logFC > 0, na.rm = TRUE),
      n_down = sum(.data$logFC < 0, na.rm = TRUE),
      markers = paste(.data$marker, collapse = ";"),
      top_by_logFC = paste(
        .data$marker[order(-abs(.data$logFC))][seq_len(min(top_n, dplyr::n()))],
        collapse = ";"
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        .data$mean_logFC > 0.05 & .data$n_up >= .data$n_down ~ up_direction,
        .data$mean_logFC < -0.05 & .data$n_down > .data$n_up ~ down_direction,
        TRUE ~ "mixed"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(.data$mean_logFC))
  cl_df
}

export_cluster_heatmap <- function(
  mat_z,
  hc,
  de_df,
  clusters,
  out_path,
  title = "Epithelial proteomics marker clusters (High vs Low risk DE)"
) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("Skipping heatmap (pheatmap not installed).")
    return(invisible(NULL))
  }
  ord <- hc$order
  mat_ord <- mat_z[ord, , drop = FALSE]
  ann_row <- data.frame(
    logFC = de_df$logFC[match(rownames(mat_ord), de_df$marker)],
    cluster = factor(clusters[rownames(mat_ord)]),
    row.names = rownames(mat_ord)
  )
  risk_ann <- data.frame(
    risk_grp = factor(
      c("Low", "High")[match(colnames(mat_z), colnames(mat_z))],
      levels = c("Low", "High")
    ),
    row.names = colnames(mat_z)
  )
  # rebuild column annotation from meta passed separately via colnames order - fix below in caller
  invisible(NULL)
}

export_cluster_heatmap_with_meta <- function(
  mat_z,
  hc,
  de_df,
  clusters,
  sample_meta,
  out_path,
  title = NULL,
  group_col = "risk_grp",
  group_levels = c("Low", "High")
) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("Skipping heatmap (pheatmap not installed).")
    return(invisible(NULL))
  }
  ord <- hc$order
  mat_ord <- mat_z[ord, , drop = FALSE]
  ann_row <- data.frame(
    logFC = de_df$logFC[match(rownames(mat_ord), de_df$marker)],
    cluster = factor(clusters[rownames(mat_ord)]),
    row.names = rownames(mat_ord)
  )
  col_meta <- sample_meta[match(colnames(mat_ord), sample_meta$patient_id), , drop = FALSE]
  ann_col <- data.frame(
    group = factor(col_meta[[group_col]], levels = group_levels),
    row.names = colnames(mat_ord)
  )
  pheatmap::pheatmap(
    mat_ord,
    cluster_rows = hc,
    cluster_cols = TRUE,
    annotation_row = ann_row,
    annotation_col = ann_col,
    show_colnames = FALSE,
    main = title %||% "Epithelial markers: z-scored abundance",
    filename = out_path,
    width = 10,
    height = 8
  )
  message("Wrote ", out_path)
  invisible(out_path)
}

pick_top_clusters <- function(cluster_summary, up_direction, down_direction) {
  top_up_cluster <- cluster_summary %>%
    dplyr::filter(.data$direction == up_direction) %>%
    dplyr::slice_max(.data$mean_logFC, n = 1, with_ties = FALSE)
  if (nrow(top_up_cluster) == 0L) {
    top_up_cluster <- cluster_summary %>% dplyr::slice_max(.data$mean_logFC, n = 1)
  }
  top_down_cluster <- cluster_summary %>%
    dplyr::filter(.data$direction == down_direction) %>%
    dplyr::slice_min(.data$mean_logFC, n = 1, with_ties = FALSE)
  if (nrow(top_down_cluster) == 0L) {
    top_down_cluster <- cluster_summary %>% dplyr::slice_min(.data$mean_logFC, n = 1)
  }
  list(up = top_up_cluster, down = top_down_cluster)
}

run_protein_de_clusters_internal <- function(
  aligned,
  de,
  output_dir,
  contrast_label,
  group_col,
  group_levels,
  group_low_label,
  group_high_label,
  up_direction,
  down_direction,
  y_axis_label,
  n_clusters = DEFAULT_N_CLUSTERS,
  top_n = DEFAULT_TOP_MARKERS,
  compartment_suffix = DEFAULT_COMPARTMENT_SUFFIX
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  plots_dir <- file.path(output_dir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  de_df <- de$results %>%
    dplyr::arrange(dplyr::desc(.data$logFC))

  cl <- cluster_markers_by_profile(aligned$matrix, n_clusters = n_clusters)
  cluster_assign <- tibble::tibble(
    marker = names(cl$clusters),
    cluster = unname(cl$clusters)
  ) %>%
    dplyr::left_join(de_df, by = "marker") %>%
    dplyr::arrange(.data$cluster, dplyr::desc(.data$logFC))

  cluster_summary <- summarize_marker_clusters(
    cl$clusters,
    de_df,
    top_n = top_n,
    up_direction = up_direction,
    down_direction = down_direction
  )
  tops <- pick_top_clusters(cluster_summary, up_direction, down_direction)
  top_up_cluster <- tops$up
  top_down_cluster <- tops$down

  top_up_markers <- cluster_assign %>%
    dplyr::filter(.data$cluster == top_up_cluster$cluster[1]) %>%
    dplyr::arrange(dplyr::desc(.data$logFC))
  top_down_markers <- cluster_assign %>%
    dplyr::filter(.data$cluster == top_down_cluster$cluster[1]) %>%
    dplyr::arrange(.data$logFC)

  readr::write_csv(de_df, file.path(output_dir, "protein_epi_limma.csv"))
  readr::write_csv(cluster_assign, file.path(output_dir, "marker_cluster_assignments.csv"))
  readr::write_csv(cluster_summary, file.path(output_dir, "cluster_summary.csv"))
  readr::write_csv(top_up_markers, file.path(output_dir, "top_upregulated_cluster_markers.csv"))
  readr::write_csv(top_down_markers, file.path(output_dir, "top_downregulated_cluster_markers.csv"))

  n_low <- sum(aligned$meta[[group_col]] == group_levels[1])
  n_high <- sum(aligned$meta[[group_col]] == group_levels[2])
  export_cluster_heatmap_with_meta(
    cl$mat_z,
    cl$hclust,
    de_df,
    cl$clusters,
    aligned$meta,
    file.path(plots_dir, "marker_cluster_heatmap.png"),
    title = paste0(
      compartment_label(compartment_suffix),
      " markers (n=",
      length(aligned$patients),
      " patients: ",
      n_low, " ", group_low_label, " / ",
      n_high, " ", group_high_label, ")"
    ),
    group_col = group_col,
    group_levels = group_levels
  )

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    bar_df <- cluster_summary %>%
      dplyr::mutate(
        cluster_label = paste0("Cluster ", .data$cluster, " (", .data$direction, ")")
      )
    p <- ggplot2::ggplot(
      bar_df,
      ggplot2::aes(
        x = reorder(.data$cluster_label, .data$mean_logFC),
        y = .data$mean_logFC,
        fill = .data$direction
      )
    ) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(
        values = stats::setNames(
          c("#B2182B", "#2166AC", "grey60"),
          c(up_direction, down_direction, "mixed")
        )
      ) +
      ggplot2::labs(
        x = NULL,
        y = y_axis_label,
        title = paste0("Epithelial protein clusters — ", contrast_label),
        fill = "Direction"
      ) +
      ggplot2::theme_bw()
    ggplot2::ggsave(
      file.path(plots_dir, "cluster_mean_logFC_barplot.png"),
      p,
      width = 8,
      height = 5,
      dpi = 150
    )
    message("Wrote cluster mean logFC barplot.")
  }

  top_de_up <- de_df %>% dplyr::filter(.data$logFC > 0) %>% dplyr::slice_head(n = top_n)
  top_de_down <- de_df %>% dplyr::filter(.data$logFC < 0) %>% dplyr::slice_tail(n = top_n)
  readr::write_csv(top_de_up, file.path(output_dir, "top_individual_up.csv"))
  readr::write_csv(top_de_down, file.path(output_dir, "top_individual_down.csv"))

  log_lines <- c(
    paste0("Epithelial proteomics clusters — ", contrast_label),
    paste0("Compartment: ", compartment_suffix),
    paste0(
      "Patients: ", length(aligned$patients),
      " (", group_low_label, "=", n_low,
      ", ", group_high_label, "=", n_high, ")"
    ),
    paste0("Markers: ", nrow(de_df)),
    paste0("Clusters (hclust on Spearman profile): ", n_clusters),
    "",
    paste0(
      "Top UP-regulated cluster: ", top_up_cluster$cluster[1],
      " (mean logFC=", signif(top_up_cluster$mean_logFC[1], 3),
      ", n=", top_up_cluster$n_markers[1], ")"
    ),
    paste0("  Markers: ", top_up_cluster$markers[1]),
    "",
    paste0(
      "Top DOWN-regulated cluster: ", top_down_cluster$cluster[1],
      " (mean logFC=", signif(top_down_cluster$mean_logFC[1], 3),
      ", n=", top_down_cluster$n_markers[1], ")"
    ),
    paste0("  Markers: ", top_down_cluster$markers[1]),
    "",
    "Individual top markers (limma, not clustered):",
    paste0("  Up in ", group_high_label, ": ", paste(top_de_up$marker, collapse = ", ")),
    paste0("  Down in ", group_high_label, ": ", paste(top_de_down$marker, collapse = ", ")),
    "",
    paste0("Output: ", output_dir)
  )
  writeLines(log_lines, file.path(output_dir, "summary.txt"))
  cat(paste(log_lines, collapse = "\n"), "\n")

  invisible(list(
    de = de_df,
    clusters = cluster_assign,
    cluster_summary = cluster_summary,
    top_up = top_up_markers,
    top_down = top_down_markers,
    top_up_cluster = top_up_cluster,
    top_down_cluster = top_down_cluster,
    contrast_label = contrast_label
  ))
}

run_epi_protein_de_clusters <- function(
  output_dir = NULL,
  n_clusters = DEFAULT_N_CLUSTERS,
  top_n = DEFAULT_TOP_MARKERS,
  compartment_suffix = DEFAULT_COMPARTMENT_SUFFIX,
  seed = 42L
) {
  set.seed(seed)
  if (is.null(output_dir)) {
    output_dir <- file.path(INTEGRATED_DIR, "results", "epi_protein_high_vs_low_clusters")
  }

  prot <- load_epi_protein_matrix(compartment_suffix = compartment_suffix)
  aligned <- subset_protein_high_low(prot)
  de <- run_limma_high_vs_low_matrix(
    aligned$matrix,
    aligned$meta,
    feature_col = "marker"
  )

  run_protein_de_clusters_internal(
    aligned = aligned,
    de = de,
    output_dir = output_dir,
    contrast_label = "High vs Low RSF risk",
    group_col = "risk_grp",
    group_levels = c("Low", "High"),
    group_low_label = "Low",
    group_high_label = "High",
    up_direction = "up_in_high",
    down_direction = "down_in_high",
    y_axis_label = "Mean logFC (High vs Low)",
    n_clusters = n_clusters,
    top_n = top_n,
    compartment_suffix = compartment_suffix
  )
}

run_epi_protein_de_clusters_g1_vs_g4 <- function(
  output_dir = NULL,
  groups_rds = DEFAULT_TMA_GROUPS_RDS,
  group_low = G1_GROUP,
  group_high = G4_GROUP,
  n_clusters = DEFAULT_N_CLUSTERS,
  top_n = DEFAULT_TOP_MARKERS,
  compartment_suffix = DEFAULT_COMPARTMENT_SUFFIX,
  min_patients = 10L,
  seed = 42L
) {
  set.seed(seed)
  if (is.null(output_dir)) {
    output_dir <- file.path(
      INTEGRATED_DIR,
      "results",
      "epi_protein_G1_vs_G4_clusters"
    )
  }

  prot <- load_epi_protein_with_tma_groups(
    groups_rds = groups_rds,
    compartment_suffix = compartment_suffix
  )
  aligned <- subset_protein_two_groups(
    prot,
    group_col = "tma_comp_grp",
    group_low = group_low,
    group_high = group_high,
    min_patients = min_patients
  )
  de <- run_limma_two_group_matrix(
    aligned$matrix,
    aligned$meta,
    group_col = "tma_comp_grp",
    group_low = group_low,
    group_high = group_high,
    feature_col = "marker"
  )

  run_protein_de_clusters_internal(
    aligned = aligned,
    de = de,
    output_dir = output_dir,
    contrast_label = paste0(group_low, " vs ", group_high),
    group_col = "tma_comp_grp",
    group_levels = c(group_low, group_high),
    group_low_label = group_low,
    group_high_label = group_high,
    up_direction = "up_in_g4",
    down_direction = "down_in_g4",
    y_axis_label = paste0("Mean logFC (", group_high, " vs ", group_low, ")"),
    n_clusters = n_clusters,
    top_n = top_n,
    compartment_suffix = compartment_suffix
  )
}

run_protein_de_clusters_both_contrasts <- function(
  output_root = NULL,
  n_clusters = DEFAULT_N_CLUSTERS,
  top_n = DEFAULT_TOP_MARKERS,
  compartment_suffix = DEFAULT_COMPARTMENT_SUFFIX,
  seed = 42L
) {
  if (is.null(output_root)) {
    output_root <- file.path(INTEGRATED_DIR, "results", "protein_de_clusters")
  }
  high_low <- run_epi_protein_de_clusters(
    output_dir = file.path(output_root, "high_vs_low"),
    n_clusters = n_clusters,
    top_n = top_n,
    compartment_suffix = compartment_suffix,
    seed = seed
  )
  g1_g4 <- run_epi_protein_de_clusters_g1_vs_g4(
    output_dir = file.path(output_root, "G1_vs_G4"),
    n_clusters = n_clusters,
    top_n = top_n,
    compartment_suffix = compartment_suffix,
    seed = seed
  )

  combined_lines <- c(
    "Protein DE clusters — both contrasts",
    paste0("Compartment: ", compartment_label(compartment_suffix), " (", compartment_suffix, ")"),
    "",
    "=== High risk vs Low risk (RSF) ===",
    paste0(
      "Top UP cluster ", high_low$top_up_cluster$cluster[1],
      ": mean logFC=", signif(high_low$top_up_cluster$mean_logFC[1], 3),
      " | ", high_low$top_up_cluster$markers[1]
    ),
    paste0(
      "Top DOWN cluster ", high_low$top_down_cluster$cluster[1],
      ": mean logFC=", signif(high_low$top_down_cluster$mean_logFC[1], 3),
      " | ", high_low$top_down_cluster$markers[1]
    ),
    "",
    "=== G1_all_low vs G4_high_3plus (TMA composition) ===",
    paste0(
      "Top UP cluster ", g1_g4$top_up_cluster$cluster[1],
      ": mean logFC=", signif(g1_g4$top_up_cluster$mean_logFC[1], 3),
      " | ", g1_g4$top_up_cluster$markers[1]
    ),
    paste0(
      "Top DOWN cluster ", g1_g4$top_down_cluster$cluster[1],
      ": mean logFC=", signif(g1_g4$top_down_cluster$mean_logFC[1], 3),
      " | ", g1_g4$top_down_cluster$markers[1]
    ),
    "",
    paste0("Outputs: ", output_root)
  )
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  writeLines(combined_lines, file.path(output_root, "combined_summary.txt"))
  cat("\n", paste(combined_lines, collapse = "\n"), "\n", sep = "")

  invisible(list(high_vs_low = high_low, g1_vs_g4 = g1_g4))
}

run_protein_de_clusters_immune_compartments <- function(
  output_root = NULL,
  compartments = IMMUNE_COMPARTMENTS,
  n_clusters = DEFAULT_N_CLUSTERS,
  top_n = DEFAULT_TOP_MARKERS,
  seed = 42L
) {
  if (is.null(output_root)) {
    output_root <- file.path(INTEGRATED_DIR, "results", "protein_de_clusters_immune")
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  results <- list()
  summary_rows <- list()

  for (i in seq_len(nrow(compartments))) {
    suffix <- compartments$suffix[i]
    label <- compartments$label[i]
    comp_dir <- file.path(output_root, suffix)
    message("\n=== ", label, " (", suffix, ") ===")

    res <- tryCatch(
      run_protein_de_clusters_both_contrasts(
        output_root = comp_dir,
        n_clusters = n_clusters,
        top_n = top_n,
        compartment_suffix = suffix,
        seed = seed
      ),
      error = function(e) {
        message("FAILED ", suffix, ": ", conditionMessage(e))
        NULL
      }
    )
    results[[suffix]] <- res
    if (is.null(res)) next

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      compartment = suffix,
      compartment_label = label,
      contrast = "high_vs_low",
      top_up_cluster = res$high_vs_low$top_up_cluster$cluster[1],
      top_up_mean_logFC = res$high_vs_low$top_up_cluster$mean_logFC[1],
      top_up_markers = res$high_vs_low$top_up_cluster$markers[1],
      top_down_cluster = res$high_vs_low$top_down_cluster$cluster[1],
      top_down_mean_logFC = res$high_vs_low$top_down_cluster$mean_logFC[1],
      top_down_markers = res$high_vs_low$top_down_cluster$markers[1],
      stringsAsFactors = FALSE
    )
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      compartment = suffix,
      compartment_label = label,
      contrast = "G1_vs_G4",
      top_up_cluster = res$g1_vs_g4$top_up_cluster$cluster[1],
      top_up_mean_logFC = res$g1_vs_g4$top_up_cluster$mean_logFC[1],
      top_up_markers = res$g1_vs_g4$top_up_cluster$markers[1],
      top_down_cluster = res$g1_vs_g4$top_down_cluster$cluster[1],
      top_down_mean_logFC = res$g1_vs_g4$top_down_cluster$mean_logFC[1],
      top_down_markers = res$g1_vs_g4$top_down_cluster$markers[1],
      stringsAsFactors = FALSE
    )
  }

  summary_df <- if (length(summary_rows) > 0L) {
    dplyr::bind_rows(summary_rows)
  } else {
    data.frame()
  }
  if (nrow(summary_df) > 0L) {
    readr::write_csv(summary_df, file.path(output_root, "immune_compartments_summary.csv"))
  }

  master_lines <- c(
    "Protein DE clusters — immune compartments (helperT, cytT, regT)",
    paste0("Output root: ", output_root),
    ""
  )
  for (i in seq_len(nrow(compartments))) {
    suffix <- compartments$suffix[i]
    label <- compartments$label[i]
    res <- results[[suffix]]
    master_lines <- c(master_lines, paste0("--- ", label, " (", suffix, ") ---"))
    if (is.null(res)) {
      master_lines <- c(master_lines, "  (analysis failed)", "")
      next
    }
    master_lines <- c(
      master_lines,
      "  High vs Low:",
      paste0(
        "    UP cluster ", res$high_vs_low$top_up_cluster$cluster[1],
        " (logFC=", signif(res$high_vs_low$top_up_cluster$mean_logFC[1], 3),
        "): ", res$high_vs_low$top_up_cluster$markers[1]
      ),
      paste0(
        "    DOWN cluster ", res$high_vs_low$top_down_cluster$cluster[1],
        " (logFC=", signif(res$high_vs_low$top_down_cluster$mean_logFC[1], 3),
        "): ", res$high_vs_low$top_down_cluster$markers[1]
      ),
      "  G1 vs G4:",
      paste0(
        "    UP cluster ", res$g1_vs_g4$top_up_cluster$cluster[1],
        " (logFC=", signif(res$g1_vs_g4$top_up_cluster$mean_logFC[1], 3),
        "): ", res$g1_vs_g4$top_up_cluster$markers[1]
      ),
      paste0(
        "    DOWN cluster ", res$g1_vs_g4$top_down_cluster$cluster[1],
        " (logFC=", signif(res$g1_vs_g4$top_down_cluster$mean_logFC[1], 3),
        "): ", res$g1_vs_g4$top_down_cluster$markers[1]
      ),
      ""
    )
  }
  writeLines(master_lines, file.path(output_root, "master_summary.txt"))
  cat("\n", paste(master_lines, collapse = "\n"), "\n", sep = "")

  invisible(list(by_compartment = results, summary = summary_df))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
