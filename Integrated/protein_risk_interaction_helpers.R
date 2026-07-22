# Compare protein co-abundance interactions: Low vs High RSF risk.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

source("protein_de_string_helpers.R")
source("cms_stratified_proteomics_helpers.R")

DEFAULT_CORR_EDGE <- 0.35
DEFAULT_CORR_WEAK <- 0.15
DEFAULT_DIFF_P <- 0.05

HIGH_HIGHLIGHT_COLORS <- c(
  "Present in High only" = "#F46D43",
  "Stronger in High" = "#D73027",
  "Other" = "#C8C8C8"
)
LOW_HIGHLIGHT_COLORS <- c(
  "Present in Low only" = "#1B7837",
  "Stronger in Low" = "#5AAE61",
  "Other" = "#C8C8C8"
)

network_layout_coords <- function(edge_df) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph is required.")
  }
  if (nrow(edge_df) == 0L) {
    return(list(node_df = data.frame(), layout = NULL))
  }
  nodes <- unique(c(edge_df$marker_a, edge_df$marker_b))
  g <- igraph::graph_from_data_frame(
    edge_df[, c("marker_a", "marker_b")],
    directed = FALSE,
    vertices = nodes
  )
  lay <- igraph::layout_with_fr(g)
  list(
    node_df = data.frame(
      marker = nodes,
      x = lay[, 1],
      y = lay[, 2],
      stringsAsFactors = FALSE
    ),
    layout = lay
  )
}

pair_key <- function(a, b) {
  paste(pmin(a, b), pmax(a, b), sep = "|")
}

cor_diff_test <- function(x1, y1, x2, y2) {
  ok1 <- is.finite(x1) & is.finite(y1)
  ok2 <- is.finite(x2) & is.finite(y2)
  n1 <- sum(ok1)
  n2 <- sum(ok2)
  if (n1 < 6L || n2 < 6L) {
    return(list(r_low = NA_real_, r_high = NA_real_, p_diff = NA_real_, n_low = n1, n_high = n2))
  }
  r1 <- stats::cor(x1[ok1], y1[ok1], method = "spearman")
  r2 <- stats::cor(x2[ok2], y2[ok2], method = "spearman")
  if (!is.finite(r1) || !is.finite(r2) || abs(r1) >= 1 || abs(r2) >= 1) {
    return(list(r_low = r1, r_high = r2, p_diff = NA_real_, n_low = n1, n_high = n2))
  }
  z1 <- atanh(r1)
  z2 <- atanh(r2)
  se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  z <- (z1 - z2) / se
  p <- 2 * stats::pnorm(-abs(z))
  list(r_low = r1, r_high = r2, p_diff = p, n_low = n1, n_high = n2)
}

compare_risk_protein_interactions <- function(
  mat,
  meta,
  corr_edge = DEFAULT_CORR_EDGE,
  corr_weak = DEFAULT_CORR_WEAK,
  diff_p = DEFAULT_DIFF_P
) {
  low_pts <- meta$patient_id[meta$risk_grp == "Low"]
  high_pts <- meta$patient_id[meta$risk_grp == "High"]
  mat_low <- mat[, low_pts, drop = FALSE]
  mat_high <- mat[, high_pts, drop = FALSE]
  markers <- rownames(mat)
  n <- length(markers)
  if (n < 2L) {
    return(data.frame())
  }

  rows <- list()
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      a <- markers[i]
      b <- markers[j]
      tst <- cor_diff_test(
        mat_low[a, ], mat_low[b, ],
        mat_high[a, ], mat_high[b, ]
      )
      r_low <- tst$r_low
      r_high <- tst$r_high
      category <- dplyr::case_when(
        is.finite(r_low) && is.finite(r_high) &&
          abs(r_low) >= corr_edge && abs(r_high) < corr_weak ~ "present_in_Low_only",
        is.finite(r_low) && is.finite(r_high) &&
          abs(r_high) >= corr_edge && abs(r_low) < corr_weak ~ "present_in_High_only",
        is.finite(r_low) && is.finite(r_high) &&
          r_low * r_high < 0 && abs(r_low) >= corr_weak && abs(r_high) >= corr_weak ~ "sign_flip",
        is.finite(tst$p_diff) && tst$p_diff <= diff_p && is.finite(r_low) && is.finite(r_high) &&
          abs(r_low) > abs(r_high) + 0.1 ~ "stronger_in_Low",
        is.finite(tst$p_diff) && tst$p_diff <= diff_p && is.finite(r_low) && is.finite(r_high) &&
          abs(r_high) > abs(r_low) + 0.1 ~ "stronger_in_High",
        TRUE ~ "not_different"
      )
      if (category == "not_different") next
      rows[[length(rows) + 1L]] <- data.frame(
        marker_a = a,
        marker_b = b,
        r_low = r_low,
        r_high = r_high,
        delta_r = r_low - r_high,
        p_diff = tst$p_diff,
        n_low = tst$n_low,
        n_high = tst$n_high,
        category = category,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    return(data.frame())
  }
  dplyr::bind_rows(rows) %>%
    dplyr::arrange(.data$category, dplyr::desc(abs(.data$delta_r)))
}

plot_differential_interaction_network <- function(
  edge_df,
  title,
  subtitle = NULL,
  highlight = c("all", "high", "low"),
  node_df = NULL
) {
  highlight <- match.arg(highlight)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }
  if (nrow(edge_df) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No differential interactions", size = 4) +
        ggplot2::labs(title = title, subtitle = subtitle) +
        ggplot2::theme_void()
    )
  }

  if (is.null(node_df)) {
    node_df <- network_layout_coords(edge_df)$node_df
  }
  edge_plot <- edge_df %>%
    dplyr::left_join(node_df, by = c("marker_a" = "marker")) %>%
    dplyr::rename(x1 = x, y1 = y) %>%
    dplyr::left_join(node_df, by = c("marker_b" = "marker")) %>%
    dplyr::rename(x2 = x, y2 = y)

  if (highlight == "high") {
    edge_plot <- edge_plot %>%
      dplyr::mutate(
        plot_category = dplyr::case_when(
          .data$category == "present_in_High_only" ~ "Present in High only",
          .data$category == "stronger_in_High" ~ "Stronger in High",
          TRUE ~ "Other"
        )
      )
    cat_cols <- HIGH_HIGHLIGHT_COLORS
    edge_other <- edge_plot %>% dplyr::filter(.data$plot_category == "Other")
    edge_hi <- edge_plot %>% dplyr::filter(.data$plot_category != "Other")
  } else if (highlight == "low") {
    edge_plot <- edge_plot %>%
      dplyr::mutate(
        plot_category = dplyr::case_when(
          .data$category == "present_in_Low_only" ~ "Present in Low only",
          .data$category == "stronger_in_Low" ~ "Stronger in Low",
          TRUE ~ "Other"
        )
      )
    cat_cols <- LOW_HIGHLIGHT_COLORS
    edge_other <- edge_plot %>% dplyr::filter(.data$plot_category == "Other")
    edge_hi <- edge_plot %>% dplyr::filter(.data$plot_category != "Other")
  } else {
    edge_plot <- edge_plot %>%
      dplyr::mutate(plot_category = .data$category)
    cat_cols <- c(
      present_in_Low_only = "#2166AC",
      present_in_High_only = "#B2182B",
      stronger_in_Low = "#4393C3",
      stronger_in_High = "#D6604D",
      sign_flip = "#9970AB"
    )
    edge_other <- edge_plot[0, ]
    edge_hi <- edge_plot
  }

  edge_plot <- edge_plot %>%
    dplyr::mutate(
      plot_category = factor(
        .data$plot_category,
        levels = names(cat_cols)
      )
    )
  if (nrow(edge_hi) > 0L) {
    edge_hi <- edge_hi %>%
      dplyr::mutate(
        plot_category = factor(.data$plot_category, levels = names(cat_cols))
      )
  }
  if (nrow(edge_other) > 0L) {
    edge_other <- edge_other %>%
      dplyr::mutate(
        plot_category = factor("Other", levels = names(cat_cols))
      )
  }

  p <- ggplot2::ggplot()
  if (nrow(edge_other) > 0L) {
    p <- p +
      ggplot2::geom_segment(
        data = edge_other,
        ggplot2::aes(
          x = x1, y = y1, xend = x2, yend = y2,
          linewidth = abs(delta_r)
        ),
        colour = cat_cols[["Other"]],
        alpha = 0.45
      )
  }
  if (nrow(edge_hi) > 0L) {
    p <- p +
      ggplot2::geom_segment(
        data = edge_hi,
        ggplot2::aes(
          x = x1, y = y1, xend = x2, yend = y2,
          colour = plot_category,
          linewidth = abs(delta_r)
        ),
        alpha = 0.95
      )
  }
  p +
    ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(x = x, y = y),
      size = 4,
      shape = 21,
      fill = "grey92",
      colour = "grey25"
    ) +
    ggplot2::geom_text(
      data = node_df,
      ggplot2::aes(x = x, y = y, label = marker),
      size = 2.6,
      vjust = -1.1
    ) +
    ggplot2::scale_colour_manual(values = cat_cols, name = NULL, drop = FALSE) +
    ggplot2::scale_linewidth_continuous(range = c(0.4, 1.8), guide = "none") +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8, colour = "grey35"),
      legend.position = "bottom"
    )
}

export_risk_interaction_comparison <- function(
  compartment_suffix,
  output_dir = NULL,
  corr_edge = DEFAULT_CORR_EDGE,
  corr_weak = DEFAULT_CORR_WEAK,
  diff_p = DEFAULT_DIFF_P
) {
  label <- compartment_display_label(compartment_suffix)
  if (is.null(output_dir)) {
    output_dir <- file.path(
      INTEGRATED_DIR,
      "results",
      "protein_risk_interaction_diff",
      compartment_suffix
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  prot <- load_epi_protein_matrix(compartment_suffix = compartment_suffix)
  aligned <- subset_protein_high_low(prot)
  diff_edges <- compare_risk_protein_interactions(
    aligned$matrix,
    aligned$meta,
    corr_edge = corr_edge,
    corr_weak = corr_weak,
    diff_p = diff_p
  )
  readr::write_csv(diff_edges, file.path(output_dir, "differential_interactions_low_vs_high.csv"))

  n_low <- sum(aligned$meta$risk_grp == "Low")
  n_high <- sum(aligned$meta$risk_grp == "High")
  layout <- network_layout_coords(diff_edges)
  node_df <- layout$node_df
  sub_base <- paste0(
    "Co-abundance (Spearman); Low n=", n_low, ", High n=", n_high,
    "; ", nrow(diff_edges), " differential pair(s)"
  )

  p_all <- plot_differential_interaction_network(
    diff_edges,
    title = paste0(label, ": interactions differing by risk"),
    subtitle = sub_base,
    highlight = "all",
    node_df = node_df
  )
  p_high <- plot_differential_interaction_network(
    diff_edges,
    title = paste0(label, ": High-risk interactions"),
    subtitle = paste0(sub_base, " | orange = High only, red = stronger in High"),
    highlight = "high",
    node_df = node_df
  )
  p_low <- plot_differential_interaction_network(
    diff_edges,
    title = paste0(label, ": Low-risk interactions"),
    subtitle = paste0(sub_base, " | dark green = Low only, light green = stronger in Low"),
    highlight = "low",
    node_df = node_df
  )

  ggplot2::ggsave(
    file.path(output_dir, "differential_interactions_network.png"),
    p_all, width = 10, height = 8, dpi = 180, bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, "differential_interactions_network_high.png"),
    p_high, width = 10, height = 8, dpi = 180, bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, "differential_interactions_network_low.png"),
    p_low, width = 10, height = 8, dpi = 180, bg = "white"
  )
  message("Saved high/low highlight networks in ", output_dir)

  summarize_edges <- function(df, cat) {
    sub <- df %>% dplyr::filter(.data$category == cat)
    if (nrow(sub) == 0L) return("")
    pairs <- paste0(sub$marker_a, "-", sub$marker_b, " (r_Low=", signif(sub$r_low, 2),
                    ", r_High=", signif(sub$r_high, 2), ")", collapse = "; ")
    paste0(cat, ": ", pairs)
  }

  summary_lines <- c(
    paste0("Differential protein interactions: ", label),
    paste0("Low=", n_low, " patients | High=", n_high, " patients"),
    paste0("Total differential pairs: ", nrow(diff_edges)),
    "",
    summarize_edges(diff_edges, "present_in_Low_only"),
    summarize_edges(diff_edges, "present_in_High_only"),
    summarize_edges(diff_edges, "stronger_in_Low"),
    summarize_edges(diff_edges, "stronger_in_High"),
    summarize_edges(diff_edges, "sign_flip"),
    "",
    paste0("Output: ", output_dir),
    "  differential_interactions_network_high.png",
    "  differential_interactions_network_low.png"
  )
  summary_lines <- summary_lines[nzchar(summary_lines)]
  writeLines(summary_lines, file.path(output_dir, "summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(
    edges = diff_edges,
    plot_all = p_all,
    plot_high = p_high,
    plot_low = p_low,
    node_df = node_df,
    output_dir = output_dir
  ))
}

export_all_compartment_risk_interactions <- function(
  output_root = NULL,
  compartments = ALL_PROTEIN_COMPARTMENTS
) {
  if (is.null(output_root)) {
    output_root <- file.path(INTEGRATED_DIR, "results", "protein_risk_interaction_diff")
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  results <- list()
  summary_rows <- list()
  plots_high <- list()
  plots_low <- list()

  for (i in seq_len(nrow(compartments))) {
    suffix <- compartments$suffix[i]
    message("=== ", compartments$label[i], " ===")
    res <- export_risk_interaction_comparison(
      compartment_suffix = suffix,
      output_dir = file.path(output_root, suffix)
    )
    results[[suffix]] <- res
    if (nrow(res$edges) > 0L) {
      summary_rows[[length(summary_rows) + 1L]] <- res$edges %>%
        dplyr::mutate(
          compartment = suffix,
          compartment_label = compartments$label[i]
        )
    }
    plots_high[[suffix]] <- res$plot_high
    plots_low[[suffix]] <- res$plot_low
  }

  if (length(summary_rows) > 0L) {
    readr::write_csv(
      dplyr::bind_rows(summary_rows),
      file.path(output_root, "all_compartments_differential_interactions.csv")
    )
  }

  if (requireNamespace("cowplot", quietly = TRUE) && length(plots_high) > 0L) {
    save_master_panel <- function(plot_list, title_text, filename) {
      master <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
      master_title <- cowplot::plot_grid(
        cowplot::ggdraw() +
          cowplot::draw_label(title_text, fontface = "bold", size = 14, x = 0.01, hjust = 0),
        master,
        ncol = 1,
        rel_heights = c(0.04, 1)
      )
      out <- file.path(output_root, filename)
      ggplot2::ggsave(out, master_title, width = 14, height = 18, dpi = 160, bg = "white")
      message("Saved: ", out)
    }
    save_master_panel(
      plots_high,
      "High-risk protein interactions (orange = present in High only, red = stronger in High)",
      "differential_interactions_all_compartments_high.png"
    )
    save_master_panel(
      plots_low,
      "Low-risk protein interactions (dark green = present in Low only, light green = stronger in Low)",
      "differential_interactions_all_compartments_low.png"
    )
  }

  invisible(results)
}

export_cms_risk_interaction_comparison <- function(
  cms_group,
  compartment_suffix,
  output_dir = NULL,
  corr_edge = DEFAULT_CORR_EDGE,
  corr_weak = DEFAULT_CORR_WEAK,
  diff_p = DEFAULT_DIFF_P,
  min_low = 6L,
  min_high = 6L
) {
  label <- compartment_display_label(compartment_suffix)
  if (is.null(output_dir)) {
    output_dir <- file.path(
      INTEGRATED_DIR,
      "results",
      "protein_risk_interaction_diff_cms",
      cms_group,
      compartment_suffix
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  prot <- load_proteomics_with_cms(compartment_suffix = compartment_suffix)
  aligned <- subset_protein_cms_risk(
    prot,
    cms_group = cms_group,
    min_patients = min_low + min_high,
    min_low = min_low,
    min_high = min_high
  )
  diff_edges <- compare_risk_protein_interactions(
    aligned$matrix,
    aligned$meta,
    corr_edge = corr_edge,
    corr_weak = corr_weak,
    diff_p = diff_p
  )
  readr::write_csv(
    diff_edges,
    file.path(output_dir, "differential_interactions_low_vs_high.csv")
  )

  n_low <- aligned$n_low
  n_high <- aligned$n_high
  layout <- network_layout_coords(diff_edges)
  node_df <- layout$node_df
  sub_base <- paste0(
    cms_group, " | Co-abundance (Spearman); Low n=", n_low, ", High n=", n_high,
    "; ", nrow(diff_edges), " differential pair(s)"
  )

  p_all <- plot_differential_interaction_network(
    diff_edges,
    title = paste0(label, " (", cms_group, "): interactions differing by risk"),
    subtitle = sub_base,
    highlight = "all",
    node_df = node_df
  )
  p_high <- plot_differential_interaction_network(
    diff_edges,
    title = paste0(label, " (", cms_group, "): High-risk interactions"),
    subtitle = paste0(sub_base, " | orange = High only, red = stronger in High"),
    highlight = "high",
    node_df = node_df
  )
  p_low <- plot_differential_interaction_network(
    diff_edges,
    title = paste0(label, " (", cms_group, "): Low-risk interactions"),
    subtitle = paste0(sub_base, " | dark green = Low only, light green = stronger in Low"),
    highlight = "low",
    node_df = node_df
  )

  ggplot2::ggsave(
    file.path(output_dir, "differential_interactions_network.png"),
    p_all, width = 10, height = 8, dpi = 180, bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, "differential_interactions_network_high.png"),
    p_high, width = 10, height = 8, dpi = 180, bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, "differential_interactions_network_low.png"),
    p_low, width = 10, height = 8, dpi = 180, bg = "white"
  )

  summarize_edges <- function(df, cat) {
    sub <- df %>% dplyr::filter(.data$category == cat)
    if (nrow(sub) == 0L) return("")
    pairs <- paste0(
      sub$marker_a, "-", sub$marker_b,
      " (r_Low=", signif(sub$r_low, 2),
      ", r_High=", signif(sub$r_high, 2), ")",
      collapse = "; "
    )
    paste0(cat, ": ", pairs)
  }

  summary_lines <- c(
    paste0("Differential protein interactions: ", label, " (", cms_group, ")"),
    paste0("Low=", n_low, " patients | High=", n_high, " patients"),
    paste0("Total differential pairs: ", nrow(diff_edges)),
    "",
    summarize_edges(diff_edges, "present_in_Low_only"),
    summarize_edges(diff_edges, "present_in_High_only"),
    summarize_edges(diff_edges, "stronger_in_Low"),
    summarize_edges(diff_edges, "stronger_in_High"),
    summarize_edges(diff_edges, "sign_flip"),
    "",
    paste0("Output: ", output_dir)
  )
  summary_lines <- summary_lines[nzchar(summary_lines)]
  writeLines(summary_lines, file.path(output_dir, "summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(
    cms_group = cms_group,
    edges = diff_edges,
    plot_all = p_all,
    plot_high = p_high,
    plot_low = p_low,
    node_df = node_df,
    output_dir = output_dir,
    n_low = n_low,
    n_high = n_high
  ))
}

export_all_cms_compartment_risk_interactions <- function(
  cms_groups = CMS_STRATIFIED_GROUPS,
  output_root = NULL,
  compartments = ALL_PROTEIN_COMPARTMENTS,
  min_low = 6L,
  min_high = 6L
) {
  if (is.null(output_root)) {
    output_root <- file.path(INTEGRATED_DIR, "results", "protein_risk_interaction_diff_cms")
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  all_results <- list()
  all_summary_rows <- list()
  cms_counts <- list()

  for (cms_group in cms_groups) {
    message("\n========== ", cms_group, " ==========")
    cms_out <- file.path(output_root, cms_group)
    dir.create(cms_out, recursive = TRUE, showWarnings = FALSE)

    results <- list()
    summary_rows <- list()
    plots_high <- list()
    plots_low <- list()

    for (i in seq_len(nrow(compartments))) {
      suffix <- compartments$suffix[i]
      message("=== ", cms_group, " | ", compartments$label[i], " ===")
      res <- tryCatch(
        export_cms_risk_interaction_comparison(
          cms_group = cms_group,
          compartment_suffix = suffix,
          output_dir = file.path(cms_out, suffix),
          min_low = min_low,
          min_high = min_high
        ),
        error = function(e) {
          message("Skipped ", cms_group, " / ", suffix, ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(res)) next
      results[[suffix]] <- res
      cms_counts[[paste(cms_group, suffix, sep = "_")]] <- data.frame(
        cms = cms_group,
        compartment = suffix,
        compartment_label = compartments$label[i],
        n_low = res$n_low,
        n_high = res$n_high,
        n_pairs = nrow(res$edges),
        stringsAsFactors = FALSE
      )
      if (nrow(res$edges) > 0L) {
        summary_rows[[length(summary_rows) + 1L]] <- res$edges %>%
          dplyr::mutate(
            cms_subtype = cms_group,
            compartment = suffix,
            compartment_label = compartments$label[i]
          )
      }
      plots_high[[suffix]] <- res$plot_high
      plots_low[[suffix]] <- res$plot_low
    }

    if (length(summary_rows) > 0L) {
      cms_edges <- dplyr::bind_rows(summary_rows)
      readr::write_csv(
        cms_edges,
        file.path(cms_out, paste0(cms_group, "_all_compartments_differential_interactions.csv"))
      )
      all_summary_rows[[cms_group]] <- cms_edges
    }

    if (requireNamespace("cowplot", quietly = TRUE) && length(plots_high) > 0L) {
      save_master_panel <- function(plot_list, title_text, filename) {
        master <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
        master_title <- cowplot::plot_grid(
          cowplot::ggdraw() +
            cowplot::draw_label(title_text, fontface = "bold", size = 14, x = 0.01, hjust = 0),
          master,
          ncol = 1,
          rel_heights = c(0.04, 1)
        )
        out <- file.path(cms_out, filename)
        ggplot2::ggsave(out, master_title, width = 14, height = 18, dpi = 160, bg = "white")
        message("Saved: ", out)
      }
      save_master_panel(
        plots_high,
        paste0(
          cms_group,
          " — High-risk protein interactions (orange = present in High only, red = stronger in High)"
        ),
        paste0(cms_group, "_differential_interactions_all_compartments_high.png")
      )
      save_master_panel(
        plots_low,
        paste0(
          cms_group,
          " — Low-risk protein interactions (dark green = present in Low only, light green = stronger in Low)"
        ),
        paste0(cms_group, "_differential_interactions_all_compartments_low.png")
      )
    }

    all_results[[cms_group]] <- results
  }

  if (length(all_summary_rows) > 0L) {
    readr::write_csv(
      dplyr::bind_rows(all_summary_rows),
      file.path(output_root, "all_cms_compartments_differential_interactions.csv")
    )
  }
  if (length(cms_counts) > 0L) {
    readr::write_csv(
      dplyr::bind_rows(cms_counts),
      file.path(output_root, "cms_compartment_pair_counts.csv")
    )
  }

  invisible(all_results)
}
