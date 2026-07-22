# Publication-style figures for metacluster vs PROGENy/xCell correlation research.
# Works from existing CSV outputs or fresh run objects.

MC_COHORT_COLS <- c(
  Stage2 = "#377eb8",
  Retrospective = "#ff7f00",
  Colossus = "#4daf4a",
  Taxonomy = "#e41a1c"
)

mc_plot_dir <- function(base = NULL) {
  if (is.null(base)) {
    p <- resolve_progeny_paths()
    base <- file.path(p$base, "analysis_output", "metacluster_correlation")
  }
  fig_dir <- file.path(base, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  fig_dir
}

load_mc_correlation_outputs <- function(output_dir = NULL) {
  if (is.null(output_dir)) {
    p <- resolve_progeny_paths()
    output_dir <- file.path(p$base, "analysis_output", "metacluster_correlation")
  }
  list(
    output_dir = output_dir,
    progeny_meta = readr::read_csv(
      file.path(output_dir, "MC_progeny_fisher_meta.csv"),
      show_col_types = FALSE
    ),
    progeny_sig = readr::read_csv(
      file.path(output_dir, "MC_progeny_fisher_meta_sig.csv"),
      show_col_types = FALSE
    ),
    xcell_meta = readr::read_csv(
      file.path(output_dir, "MC_xcell_fisher_meta.csv"),
      show_col_types = FALSE
    )
  )
}

load_mc_rsf_direction <- function(rsf_path = NULL) {
  if (is.null(rsf_path)) {
    p <- resolve_progeny_paths()
    rsf_path <- file.path(
      p$proteomics,
      "results",
      "tma_composition_groups",
      "rsf_metacluster_direction.csv"
    )
  }
  if (!file.exists(rsf_path)) {
    stop("Missing RSF direction file: ", rsf_path)
  }
  readr::read_csv(rsf_path, show_col_types = FALSE) %>%
    dplyr::mutate(
      rsf_group = dplyr::if_else(
        .data$spearman_rho_rsf > 0,
        "Higher in High RSF risk",
        "Higher in Low RSF risk"
      ),
      rsf_group = factor(
        .data$rsf_group,
        levels = c("Higher in Low RSF risk", "Higher in High RSF risk")
      )
    )
}

export_mc_progeny_meta_heatmap <- function(
  progeny_meta,
  out_path,
  fdr_label = 0.05,
  title = "Metacluster vs PROGENy (Fisher meta mean Spearman rho)",
  metaclusters = NULL,
  mc_order = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  plot_df <- progeny_meta %>%
    dplyr::filter(is.finite(.data$mean_rho)) %>%
    {if (!is.null(metaclusters)) dplyr::filter(., .data$metacluster %in% metaclusters) else .} %>%
    dplyr::mutate(
      neg_log10_fdr = -log10(pmax(.data$fisher_padj, .Machine$double.xmin)),
      sig = .data$fisher_padj <= fdr_label & .data$sign_consistent %in% TRUE
    )
  if (nrow(plot_df) == 0) return(invisible(NULL))

  if (!is.null(mc_order)) {
    plot_df$metacluster <- factor(
      plot_df$metacluster,
      levels = rev(intersect(mc_order, unique(plot_df$metacluster)))
    )
  }

  n_mc <- length(unique(plot_df$metacluster))
  n_path <- length(unique(plot_df$feature))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$feature,
      y = .data$metacluster,
      fill = .data$mean_rho
    )
  ) +
    ggplot2::geom_tile(colour = "grey92", linewidth = 0.3) +
    ggplot2::geom_point(
      data = plot_df %>% dplyr::filter(.data$sig),
      ggplot2::aes(size = .data$neg_log10_fdr),
      shape = 21,
      fill = NA,
      colour = "black",
      stroke = 0.8
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Meta rho"
    ) +
    ggplot2::scale_size_continuous(range = c(1.5, 4), name = expression(-log[10] * "FDR")) +
    ggplot2::labs(
      title = title,
      subtitle = "Rings = sign-consistent Fisher FDR <= 0.05",
      x = NULL,
      y = "Metacluster"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  fig_h <- max(4.5, 0.35 * n_mc + 2.5)
  fig_w <- max(9, 0.55 * n_path + 4)
  ggplot2::ggsave(out_path, p, width = fig_w, height = fig_h, dpi = 180, bg = "white")
  message("Wrote ", out_path)
  invisible(out_path)
}

export_mc_progeny_meta_heatmaps_by_rsf <- function(
  progeny_meta,
  rsf_path,
  fig_dir,
  fdr_label = 0.05
) {
  rsf_dir <- load_mc_rsf_direction(rsf_path)
  mc_order <- rsf_dir %>%
    dplyr::arrange(dplyr::desc(abs(.data$spearman_rho_rsf))) %>%
    dplyr::pull(.data$metacluster)

  groups <- list(
    low_risk = rsf_dir %>%
      dplyr::filter(.data$spearman_rho_rsf < 0) %>%
      dplyr::pull(.data$metacluster),
    high_risk = rsf_dir %>%
      dplyr::filter(.data$spearman_rho_rsf > 0) %>%
      dplyr::pull(.data$metacluster)
  )

  export_mc_progeny_meta_heatmap(
    progeny_meta,
    file.path(fig_dir, "MC_progeny_meta_heatmap_low_rsf_risk_mcs.png"),
    fdr_label = fdr_label,
    title = "PROGENy meta correlations — MCs higher in Low RSF risk",
    metaclusters = groups$low_risk,
    mc_order = mc_order
  )
  export_mc_progeny_meta_heatmap(
    progeny_meta,
    file.path(fig_dir, "MC_progeny_meta_heatmap_high_rsf_risk_mcs.png"),
    fdr_label = fdr_label,
    title = "PROGENy meta correlations — MCs higher in High RSF risk",
    metaclusters = groups$high_risk,
    mc_order = mc_order
  )
  invisible(fig_dir)
}

export_mc_meta_bubble <- function(
  meta_df,
  out_path,
  title,
  omics_label = "PROGENy",
  top_n = 30L,
  fdr_label = 0.05
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  plot_df <- meta_df %>%
    dplyr::filter(is.finite(.data$mean_rho), is.finite(.data$fisher_padj)) %>%
    dplyr::arrange(.data$fisher_padj) %>%
    utils::head(top_n) %>%
    dplyr::mutate(
      label = paste0(.data$metacluster, " ~ ", .data$feature),
      neg_log10_fdr = -log10(pmax(.data$fisher_padj, .Machine$double.xmin)),
      sig = .data$fisher_padj <= fdr_label
    )
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$mean_rho,
      y = .data$neg_log10_fdr,
      size = .data$n_cohorts,
      colour = .data$sig
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    ggplot2::geom_hline(
      yintercept = -log10(fdr_label),
      linetype = "dashed",
      colour = "grey70"
    ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "#B2182B", `FALSE` = "#636363"),
      labels = c(`TRUE` = paste0("FDR <= ", fdr_label), `FALSE` = "NS"),
      name = NULL
    ) +
    ggplot2::scale_size_continuous(range = c(2, 7), name = "Cohorts") +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Top ", nrow(plot_df), " ", omics_label, " pairs by Fisher meta FDR"),
      x = "Mean Spearman rho (meta)",
      y = expression(-log[10] * "(Fisher FDR)")
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    sig_df <- plot_df %>% dplyr::filter(.data$sig | .data$neg_log10_fdr > 1)
    if (nrow(sig_df) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = sig_df,
        ggplot2::aes(label = .data$label),
        size = 2.8,
        max.overlaps = 15,
        show.legend = FALSE
      )
    }
  }

  ggplot2::ggsave(out_path, p, width = 9, height = 6.5, dpi = 180, bg = "white")
  message("Wrote ", out_path)
  invisible(out_path)
}

export_mc_forest_plot <- function(
  meta_rows,
  out_path,
  title = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  cohort_cols <- grep("^rho_", names(meta_rows), value = TRUE)
  if (length(cohort_cols) == 0) return(invisible(NULL))

  plot_df <- meta_rows %>%
    dplyr::mutate(pair = paste0(.data$metacluster, " ~ ", .data$feature)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(cohort_cols),
      names_to = "cohort",
      values_to = "rho"
    ) %>%
    dplyr::mutate(
      cohort = sub("^rho_", "", .data$cohort),
      cohort = factor(.data$cohort, levels = names(MC_COHORT_COLS))
    ) %>%
    dplyr::filter(is.finite(.data$rho))

  if (nrow(plot_df) == 0) return(invisible(NULL))
  if (is.null(title)) {
    title <- "Cohort Spearman rho for meta-significant MC ~ PROGENy pairs"
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$rho, y = .data$pair, colour = .data$cohort)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::geom_line(
      ggplot2::aes(group = .data$pair),
      colour = "grey75",
      linewidth = 0.4
    ) +
    ggplot2::scale_colour_manual(values = MC_COHORT_COLS, name = "Cohort") +
    ggplot2::labs(
      title = title,
      subtitle = "Patient-level Spearman rho per cohort",
      x = "Spearman rho",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  h <- max(4, 0.55 * length(unique(plot_df$pair)) + 2)
  ggplot2::ggsave(out_path, p, width = 9, height = h, dpi = 180, bg = "white")
  message("Wrote ", out_path)
  invisible(out_path)
}

export_mc_rsf_lollipop <- function(
  rsf_path,
  out_path,
  top_n = 12L
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (!file.exists(rsf_path)) {
    message("Missing RSF direction file: ", rsf_path)
    return(invisible(NULL))
  }
  df <- readr::read_csv(rsf_path, show_col_types = FALSE) %>%
    dplyr::arrange(dplyr::desc(abs(.data$spearman_rho_rsf))) %>%
    utils::head(top_n) %>%
    dplyr::mutate(
      metacluster = factor(.data$metacluster, levels = rev(.data$metacluster))
    )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$spearman_rho_rsf, y = .data$metacluster, colour = .data$spearman_rho_rsf > 0)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data$spearman_rho_rsf, yend = .data$metacluster),
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "#B2182B", `FALSE` = "#1B4332"),
      labels = c(`TRUE` = "Higher in High RSF risk", `FALSE` = "Higher in Low RSF risk"),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Metacluster abundance vs RSF clinical risk score",
      subtitle = "Spearman rho across all patients with spatial + risk data",
      x = "Spearman rho with rsf_risk",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggplot2::ggsave(out_path, p, width = 9, height = max(5, 0.35 * nrow(df) + 2), dpi = 180, bg = "white")
  message("Wrote ", out_path)
  invisible(out_path)
}

build_mc_progeny_scatter_data <- function(
  mc_df,
  cohort_loaders,
  pairs
) {
  rows <- list()
  for (cn in names(cohort_loaders)) {
    dat <- cohort_loaders[[cn]]()
    omics_wide <- run_progeny_cohort_patient_wide(dat$expr, dat$meta)
    if (is.null(omics_wide)) next
    combined <- join_mc_with_omics(mc_df, omics_wide)
    for (i in seq_len(nrow(pairs))) {
      mc <- pairs$metacluster[i]
      feat <- pairs$feature[i]
      mc_col <- paste0(MC_PREFIX, sub("^MC_", "", mc))
      if (!mc_col %in% names(combined) || !feat %in% names(combined)) next
      df <- data.frame(
        cohort = cn,
        metacluster = mc,
        feature = feat,
        pair = paste0(mc, " ~ ", feat),
        mc_abundance = combined[[mc_col]],
        omics_score = combined[[feat]],
        risk_grp = mc_df$risk_grp[match(rownames(combined), mc_df$patient_id)],
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1L]] <- df
    }
  }
  out <- dplyr::bind_rows(rows)
  if (nrow(out) > 0) {
    out$cohort <- factor(out$cohort, levels = names(MC_COHORT_COLS))
  }
  out
}

export_mc_progeny_scatter_facets <- function(
  scatter_df,
  out_path,
  title = "Metacluster abundance vs PROGENy (patient-level)"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(scatter_df) == 0) {
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(
    scatter_df,
    ggplot2::aes(
      x = .data$mc_abundance,
      y = .data$omics_score,
      colour = .data$cohort
    )
  ) +
    ggplot2::geom_point(alpha = 0.65, size = 1.8) +
    ggplot2::geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.7,
      alpha = 0.15,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~pair, scales = "free", ncol = min(3, length(unique(scatter_df$pair)))) +
    ggplot2::scale_colour_manual(values = MC_COHORT_COLS, name = "Cohort") +
    ggplot2::labs(
      title = title,
      subtitle = "Points = patients; line = cohort-specific linear trend",
      x = "Metacluster abundance",
      y = "PROGENy pathway score"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      legend.position = "bottom"
    )

  n_panels <- length(unique(scatter_df$pair))
  ncol <- min(3, n_panels)
  ggplot2::ggsave(
    out_path,
    p,
    width = 4.2 * ncol,
    height = max(5, 3.8 * ceiling(n_panels / ncol)),
    dpi = 180,
    bg = "white"
  )
  message("Wrote ", out_path)
  invisible(out_path)
}

export_mc_cohort_progeny_heatmap <- function(
  cor_df,
  cohort,
  out_path,
  title = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  sub <- cor_df %>%
    dplyr::filter(.data$cohort == cohort, is.finite(.data$rho))
  if (nrow(sub) == 0) return(invisible(NULL))
  if (is.null(title)) {
    title <- paste0(cohort, ": MC vs PROGENy Spearman rho")
  }

  p <- ggplot2::ggplot(
    sub,
    ggplot2::aes(x = .data$feature, y = .data$metacluster, fill = .data$rho)
  ) +
    ggplot2::geom_tile(colour = "grey92", linewidth = 0.25) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-max(abs(sub$rho), na.rm = TRUE), max(abs(sub$rho), na.rm = TRUE)),
      name = "Spearman rho"
    ) +
    ggplot2::labs(title = title, x = NULL, y = "Metacluster") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      panel.grid = ggplot2::element_blank()
    )

  ggplot2::ggsave(out_path, p, width = 11, height = 7, dpi = 180, bg = "white")
  message("Wrote ", out_path)
  invisible(out_path)
}

export_mc_clinical_risk_boxplot <- function(
  mc_df,
  mc_ids = c("MC_13", "MC_8", "MC_17", "MC_0"),
  out_path
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (!"risk_grp" %in% names(mc_df)) return(invisible(NULL))

  long <- mc_df %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with(MC_PREFIX),
      names_to = "mc_column",
      values_to = "abundance"
    ) %>%
    dplyr::mutate(
      metacluster = mc_short_name(.data$mc_column),
      risk_grp = factor(.data$risk_grp, levels = c("Low", "High"))
    ) %>%
    dplyr::filter(.data$metacluster %in% mc_ids, is.finite(.data$abundance))

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$risk_grp, y = .data$abundance, fill = .data$risk_grp)
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0.35, width = 0.55) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.25, size = 0.8, colour = "grey20") +
    ggplot2::facet_wrap(~metacluster, scales = "free_y", ncol = 2) +
    ggplot2::scale_fill_manual(values = c(Low = "#4393C3", High = "#D6604D"), name = "RSF risk") +
    ggplot2::labs(
      title = "Metacluster abundance by clinical RSF risk group",
      subtitle = "Spatial metaclusters linked to PROGENy meta hits + MC_0 (high in High risk)",
      x = NULL,
      y = "Metacluster abundance"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )

  ggplot2::ggsave(out_path, p, width = 9, height = 7, dpi = 180, bg = "white")
  message("Wrote ", out_path)
  invisible(out_path)
}

run_metacluster_correlation_figures <- function(
  output_dir = NULL,
  include_patient_scatter = FALSE
) {
  if (!requireNamespace("readr", quietly = TRUE)) library(readr)
  if (!requireNamespace("dplyr", quietly = TRUE)) library(dplyr)

  p <- resolve_progeny_paths()
  if (is.null(output_dir)) {
    output_dir <- file.path(p$base, "analysis_output", "metacluster_correlation")
  }
  fig_dir <- mc_plot_dir(output_dir)
  dat <- load_mc_correlation_outputs(output_dir)

  export_mc_progeny_meta_heatmap(
    dat$progeny_meta,
    file.path(fig_dir, "MC_progeny_meta_heatmap.png")
  )
  rsf_path <- file.path(p$proteomics, "results", "tma_composition_groups", "rsf_metacluster_direction.csv")
  export_mc_progeny_meta_heatmaps_by_rsf(
    dat$progeny_meta,
    rsf_path,
    fig_dir
  )
  export_mc_meta_bubble(
    dat$progeny_meta,
    file.path(fig_dir, "MC_progeny_meta_bubble.png"),
    title = "Metacluster vs PROGENy — Fisher meta",
    omics_label = "PROGENy"
  )
  export_mc_meta_bubble(
    dat$xcell_meta %>% dplyr::arrange(.data$fisher_padj) %>% utils::head(40),
    file.path(fig_dir, "MC_xcell_meta_bubble_top40.png"),
    title = "Metacluster vs xCell — top Fisher meta pairs",
    omics_label = "xCell",
    top_n = 40L
  )

  if (nrow(dat$progeny_sig) > 0) {
    export_mc_forest_plot(
      dat$progeny_sig,
      file.path(fig_dir, "MC_progeny_sig_forest.png")
    )
  }

  rsf_path <- file.path(p$proteomics, "results", "tma_composition_groups", "rsf_metacluster_direction.csv")
  export_mc_rsf_lollipop(rsf_path, file.path(fig_dir, "MC_rsf_risk_lollipop.png"))

  survival_df <- readRDS(file.path(p$proteomics, "survival_df_with_risk.rds"))
  mc_df <- load_patient_metaclusters(survival_df)
  export_mc_clinical_risk_boxplot(
    mc_df,
    out_path = file.path(fig_dir, "MC_abundance_by_risk_grp.png")
  )

  if (include_patient_scatter && nrow(dat$progeny_sig) > 0) {
    message("Building patient scatter panels (re-runs PROGENy per cohort; may take several minutes)...")
    cohort_loaders <- list(
      Stage2 = function() load_stage2_expression_meta(p, survival_df),
      Retrospective = function() load_retrospective_expression_meta(p, survival_df),
      Colossus = function() load_colossus_expression_meta(p, survival_df),
      Taxonomy = function() load_taxonomy_expression_meta(
        p, file.path(p$proteomics, "survival_df_with_risk.rds")
      )
    )
    scatter_df <- build_mc_progeny_scatter_data(mc_df, cohort_loaders, dat$progeny_sig)
    export_mc_progeny_scatter_facets(
      scatter_df,
      file.path(fig_dir, "MC_progeny_sig_scatter_by_cohort.png")
    )
  }

  for (cn in c("Stage2", "Retrospective", "Colossus", "Taxonomy")) {
    f <- file.path(output_dir, paste0(cn, "_mc_progeny_spearman.csv"))
    if (!file.exists(f)) next
    cor_df <- readr::read_csv(f, show_col_types = FALSE)
    export_mc_cohort_progeny_heatmap(
      cor_df,
      cn,
      file.path(fig_dir, paste0(cn, "_MC_progeny_cohort_heatmap.png"))
    )
  }

  summary_lines <- c(
    "Metacluster correlation figures",
    paste0("Output: ", fig_dir),
    "",
    "MC_progeny_meta_heatmap.png — MC x PROGENy meta rho tile",
    "MC_progeny_meta_bubble.png — meta bubble volcano",
    "MC_progeny_sig_forest.png — cohort rhos for FDR-significant pairs",
    "MC_progeny_sig_scatter_by_cohort.png — patient scatter for sig pairs",
    "MC_rsf_risk_lollipop.png — MC vs rsf_risk Spearman",
    "MC_abundance_by_risk_grp.png — boxplots for key MCs",
    "{Cohort}_MC_progeny_cohort_heatmap.png — per-cohort rho heatmaps",
    "MC_xcell_meta_bubble_top40.png — top xCell meta pairs"
  )
  writeLines(summary_lines, file.path(fig_dir, "figures_index.txt"))
  message("Figures written to ", fig_dir)
  invisible(fig_dir)
}
