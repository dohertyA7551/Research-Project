# Mergeability checks for High vs Low Hallmark GSEA across cohorts.
# Analogous to xCell AD + density: per-pathway gene logFC distributions from limma.

DEFAULT_GSEA_COHORTS <- c("Stage2", "Retrospective", "Colossus", "Taxonomy")

gsea_cohort_colors <- function() {
  c(
    Stage2 = "#377eb8",
    Retrospective = "#ff7f00",
    Colossus = "#4daf4a",
    Taxonomy = "#e41a1c"
  )
}

resolve_gsea_base_dir <- function() {
  if (file.exists("analysis_output/stage2/Stage2_limma_high_vs_low_all_genes.csv")) {
    "."
  } else if (file.exists("Transcriptomics/analysis_output/stage2/Stage2_limma_high_vs_low_all_genes.csv")) {
    "Transcriptomics"
  } else {
    stop("Run from Transcriptomics/ or its parent directory.")
  }
}

gsea_cohort_configs <- function(base_dir) {
  list(
    Stage2 = list(
      prefix = "Stage2",
      dir = "stage2",
      limma = "Stage2_limma_high_vs_low_all_genes.csv",
      fgsea = "Stage2_fgsea_hallmark_all.csv"
    ),
    Retrospective = list(
      prefix = "Retrospective",
      dir = "retrospective",
      limma = "Retrospective_limma_high_vs_low_all_genes.csv",
      fgsea = "Retrospective_fgsea_hallmark_all.csv"
    ),
    Colossus = list(
      prefix = "Colossus",
      dir = "colossus",
      limma = "Colossus_limma_high_vs_low_all_genes.csv",
      fgsea = "Colossus_fgsea_hallmark_all.csv"
    ),
    Taxonomy = list(
      prefix = "Taxonomy",
      dir = "taxonomy",
      limma = "Taxonomy_limma_high_vs_low_all_genes.csv",
      fgsea = "Taxonomy_fgsea_hallmark_all.csv"
    )
  )
}

ad_test_k_samples <- function(values_by_cohort) {
  if (!requireNamespace("kSamples", quietly = TRUE)) {
    stop("Package 'kSamples' is required for Anderson-Darling tests.")
  }
  values_by_cohort <- lapply(values_by_cohort, function(v) v[is.finite(v)])
  values_by_cohort <- values_by_cohort[vapply(values_by_cohort, length, 0) > 0]
  if (length(values_by_cohort) < 2) {
    return(list(ad = NA_real_, t_ad = NA_real_, p = NA_real_, k = length(values_by_cohort)))
  }
  k <- length(values_by_cohort)
  ns <- vapply(values_by_cohort, length, integer(1))
  if (sum(ns) < 4 || min(ns) < 2) {
    return(list(ad = NA_real_, t_ad = NA_real_, p = NA_real_, k = k, ns = paste(ns, collapse = ",")))
  }
  res <- kSamples::ad.test(values_by_cohort, method = "asymptotic")
  list(
    ad = as.numeric(res$ad[1, 1]),
    t_ad = as.numeric(res$ad[1, 2]),
    p = as.numeric(res$ad[1, 3]),
    k = k,
    ns = paste(ns, collapse = ",")
  )
}

load_hallmark_high_vs_low_data <- function(
  base_dir = resolve_gsea_base_dir(),
  cohorts = DEFAULT_GSEA_COHORTS
) {
  if (!requireNamespace("readr", quietly = TRUE)) library(readr)
  cfgs <- gsea_cohort_configs(base_dir)
  pathways <- hallmark_pathways()

  limma_list <- list()
  fgsea_list <- list()
  for (cn in cohorts) {
    cfg <- cfgs[[cn]]
    out_dir <- file.path(base_dir, "analysis_output", cfg$dir)
    limma_path <- file.path(out_dir, cfg$limma)
    fgsea_path <- file.path(out_dir, cfg$fgsea)
    if (!file.exists(limma_path)) {
      stop("Missing limma CSV: ", limma_path)
    }
    if (!file.exists(fgsea_path)) {
      stop("Missing fgsea CSV: ", fgsea_path, " — run run_gsea_all.R first.")
    }
    limma_list[[cn]] <- readr::read_csv(limma_path, show_col_types = FALSE)
    fgsea_list[[cn]] <- readr::read_csv(fgsea_path, show_col_types = FALSE)
  }

  list(
    base_dir = base_dir,
    cohorts = cohorts,
    pathways = pathways,
    limma_list = limma_list,
    fgsea_list = fgsea_list
  )
}

pathways_in_all_cohorts <- function(
  pathways,
  limma_list,
  cohorts = names(limma_list),
  min_genes = 15L
) {
  gene_sets <- lapply(cohorts, function(cn) {
    unique(limma_list[[cn]]$gene[!is.na(limma_list[[cn]]$gene)])
  })
  names(gene_sets) <- cohorts

  counts <- vapply(names(pathways), function(pw) {
    genes <- pathways[[pw]]
    min(vapply(cohorts, function(cn) sum(genes %in% gene_sets[[cn]]), integer(1)))
  }, integer(1))

  names(pathways)[counts >= min_genes]
}

pathway_gene_logfc_vectors <- function(
  pathway_genes,
  limma_list,
  cohorts = names(limma_list),
  stat_col = "logFC"
) {
  stats::setNames(
    lapply(cohorts, function(cn) {
      df <- limma_list[[cn]]
      sub <- df[df$gene %in% pathway_genes, , drop = FALSE]
      v <- sub[[stat_col]]
      v <- v[is.finite(v)]
      v
    }),
    cohorts
  )
}

zscore_within_cohort <- function(vecs) {
  lapply(vecs, function(v) {
    if (length(v) < 2 || stats::sd(v) == 0) {
      return(rep(0, length(v)))
    }
    as.numeric(scale(v))
  })
}

run_hallmark_logfc_ad_by_pathway <- function(
  data,
  min_genes = 15L,
  ad_alpha = 0.05,
  zscore_within_pathway = FALSE
) {
  pathways <- data$pathways
  limma_list <- data$limma_list
  fgsea_list <- data$fgsea_list
  cohorts <- data$cohorts

  pathway_names <- pathways_in_all_cohorts(
    pathways, limma_list, cohorts = cohorts, min_genes = min_genes
  )

  rows <- lapply(pathway_names, function(pw) {
    genes <- pathways[[pw]]
    vecs <- pathway_gene_logfc_vectors(genes, limma_list, cohorts = cohorts)
    if (zscore_within_pathway) {
      vecs <- zscore_within_cohort(vecs)
    }
    tst <- ad_test_k_samples(vecs)

    nes <- vapply(cohorts, function(cn) {
      fg <- fgsea_list[[cn]]
      val <- fg$NES[match(pw, fg$pathway)]
      if (length(val) == 0) NA_real_ else val
    }, numeric(1))

    out <- data.frame(
      pathway = pw,
      pathway_label = gsub("^HALLMARK_", "", pw),
      ad_stat = tst$ad,
      t_ad = tst$t_ad,
      ad_p = tst$p,
      n_cohorts = tst$k,
      n_genes_min = min(vapply(vecs, length, integer(1))),
      pass_merge = is.finite(tst$p) && tst$p >= ad_alpha,
      mean_NES = mean(nes, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    out$pathway_label <- gsub("_", " ", out$pathway_label)
    for (cn in cohorts) {
      out[[paste0("n_genes_", cn)]] <- length(vecs[[cn]])
      out[[paste0("NES_", cn)]] <- nes[cn]
    }
    out
  })

  dplyr::bind_rows(rows)
}

run_hallmark_nes_profile_ad <- function(
  data,
  ad_alpha = 0.05
) {
  cohorts <- data$cohorts
  fgsea_list <- data$fgsea_list
  all_pw <- unique(unlist(lapply(fgsea_list, function(df) df$pathway)))

  vecs <- stats::setNames(
    lapply(cohorts, function(cn) {
      fg <- fgsea_list[[cn]]
      fg$NES[match(all_pw, fg$pathway)]
    }),
    cohorts
  )
  vecs <- lapply(vecs, function(v) v[is.finite(v)])

  tst <- ad_test_k_samples(vecs)
  data.frame(
    test = "Hallmark NES profile (all pathways)",
    n_pathways = length(all_pw),
    ad_stat = tst$ad,
    ad_p = tst$p,
    n_cohorts = tst$k,
    pass_merge = is.finite(tst$p) && tst$p >= ad_alpha,
    stringsAsFactors = FALSE
  )
}

build_pathway_logfc_long <- function(
  data,
  pathways = NULL,
  zscore_within_pathway = FALSE
) {
  if (is.null(pathways)) {
    pathways <- data$pathways
  }
  limma_list <- data$limma_list
  cohorts <- data$cohorts

  rows <- lapply(names(pathways), function(pw) {
    genes <- pathways[[pw]]
    vecs <- pathway_gene_logfc_vectors(genes, limma_list, cohorts = cohorts)
    if (zscore_within_pathway) {
      vecs <- zscore_within_cohort(vecs)
    }
    dplyr::bind_rows(lapply(cohorts, function(cn) {
      data.frame(
        pathway = pw,
        pathway_label = gsub("_", " ", gsub("^HALLMARK_", "", pw)),
        cohort = cn,
        logFC = vecs[[cn]],
        stringsAsFactors = FALSE
      )
    }))
  })
  dplyr::bind_rows(rows)
}

pick_hallmark_density_pathways <- function(
  ad_results,
  meta_path = NULL,
  n_pathways = 12L
) {
  ad_ok <- ad_results %>%
    dplyr::filter(is.finite(.data$ad_p)) %>%
    dplyr::arrange(.data$ad_p)

  pass <- ad_ok %>% dplyr::filter(.data$pass_merge)
  fail <- ad_ok %>% dplyr::filter(!.data$pass_merge)

  picks <- character(0)
  if (nrow(fail) > 0) {
    picks <- c(picks, fail$pathway[seq_len(min(4, nrow(fail)))])
  }
  if (nrow(pass) > 0) {
    picks <- c(picks, pass$pathway[seq_len(min(4, nrow(pass)))])
  }

  if (!is.null(meta_path) && file.exists(meta_path)) {
    meta <- readr::read_csv(meta_path, show_col_types = FALSE)
    meta <- meta[order(meta$fisher_padj, meta$fisher_p), ]
    picks <- c(picks, meta$pathway[seq_len(min(6, nrow(meta)))])
  } else {
    ad_fail_top <- ad_ok %>% dplyr::filter(!.data$pass_merge)
    if (nrow(ad_fail_top) > 0) {
      by_nes <- ad_fail_top[order(-abs(ad_fail_top$mean_NES)), ]
      picks <- c(picks, by_nes$pathway[seq_len(min(4, nrow(by_nes)))])
    }
  }

  picks <- unique(picks)
  picks[seq_len(min(n_pathways, length(picks)))]
}

export_hallmark_ad_summary_plot <- function(
  ad_results,
  output_dir,
  ad_alpha = 0.05
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  plot_df <- ad_results %>%
    dplyr::filter(is.finite(.data$ad_p)) %>%
    dplyr::mutate(
      neg_log10_p = -log10(pmax(.data$ad_p, .Machine$double.xmin)),
      pass = ifelse(.data$pass_merge, "Pass", "Fail")
    )
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$neg_log10_p,
      y = reorder(.data$pathway_label, .data$ad_p),
      fill = .data$pass
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(
      xintercept = -log10(ad_alpha),
      linetype = "dashed",
      colour = "grey40"
    ) +
    ggplot2::scale_fill_manual(values = c(Pass = "#4daf4a", Fail = "#e41a1c")) +
    ggplot2::labs(
      title = "Anderson-Darling: Hallmark gene logFC by pathway",
      subtitle = paste0(
        "Per pathway: do limma logFC distributions (pathway genes) match across cohorts? ",
        "Pass if AD p \u2265 ", ad_alpha
      ),
      x = expression(-log[10] * "(AD p-value)"),
      y = NULL,
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  h <- max(6, min(20, 0.28 * nrow(plot_df) + 2))
  out <- file.path(output_dir, "hallmark_ad_by_pathway_summary.png")
  ggplot2::ggsave(out, p, width = 9, height = h, dpi = 150, bg = "white")
  message("Wrote ", out)
  invisible(out)
}

export_hallmark_logfc_density_plots <- function(
  logfc_long,
  pathways,
  output_dir,
  zscore_within_pathway = FALSE,
  file_prefix = "hallmark_logfc_density"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for density plots.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  plot_df <- logfc_long %>%
    dplyr::filter(.data$pathway %in% pathways) %>%
    dplyr::mutate(
      cohort = factor(.data$cohort, levels = DEFAULT_GSEA_COHORTS),
      pathway_label = factor(
        .data$pathway_label,
        levels = unique(.data$pathway_label[match(pathways, .data$pathway)])
      )
    )

  cohort_cols <- gsea_cohort_colors()
  n_ct <- length(pathways)
  ncol <- if (n_ct <= 4) 2L else if (n_ct <= 9) 3L else 4L

  if (zscore_within_pathway) {
    x_lab <- "Z-scored logFC (within pathway & cohort)"
    suffix <- "zscore"
    scale_note <- "logFC z-scored within each pathway in each cohort (shape comparison)."
  } else {
    x_lab <- "limma logFC (High vs Low)"
    suffix <- "raw"
    scale_note <- "Raw limma logFC for genes in each Hallmark set; facet axes are free."
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$logFC, colour = .data$cohort, fill = .data$cohort)
  ) +
    ggplot2::geom_density(alpha = 0.18, linewidth = 0.9) +
    ggplot2::geom_rug(
      ggplot2::aes(colour = .data$cohort),
      alpha = 0.25,
      linewidth = 0.2,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~pathway_label, scales = "free", ncol = ncol) +
    ggplot2::scale_colour_manual(values = cohort_cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = cohort_cols, drop = FALSE) +
    ggplot2::labs(
      title = "Hallmark pathway gene logFC distributions by cohort",
      subtitle = paste0(scale_note, " High vs Low DFS risk limma contrast."),
      x = x_lab,
      y = "Density",
      colour = "Cohort",
      fill = "Cohort"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(size = 8, face = "bold"),
      legend.position = "bottom"
    )

  panel_path <- file.path(output_dir, paste0(file_prefix, "_", suffix, "_panel.png"))
  panel_h <- max(8, 2.2 * ceiling(n_ct / ncol))
  ggplot2::ggsave(panel_path, p, width = 4.2 * ncol, height = panel_h, dpi = 180)
  message("Wrote ", panel_path)
  invisible(panel_path)
}

export_hallmark_nes_profile_density <- function(
  data,
  output_dir
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  cohorts <- data$cohorts
  fgsea_list <- data$fgsea_list

  plot_df <- dplyr::bind_rows(lapply(cohorts, function(cn) {
    fg <- fgsea_list[[cn]]
    data.frame(
      cohort = cn,
      NES = fg$NES,
      pathway = fg$pathway,
      stringsAsFactors = FALSE
    )
  })) %>%
    dplyr::filter(is.finite(.data$NES)) %>%
    dplyr::mutate(cohort = factor(.data$cohort, levels = cohorts))

  cohort_cols <- gsea_cohort_colors()
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$NES, colour = .data$cohort, fill = .data$cohort)
  ) +
    ggplot2::geom_density(alpha = 0.2, linewidth = 1) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
    ggplot2::scale_colour_manual(values = cohort_cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = cohort_cols, drop = FALSE) +
    ggplot2::labs(
      title = "Hallmark NES profile across cohorts",
      subtitle = paste0(
        "One NES per pathway per cohort (", nrow(plot_df) / length(cohorts),
        " pathways). Compares overall enrichment fingerprint, not per-pathway merge."
      ),
      x = "fgsea NES (High vs Low)",
      y = "Density",
      colour = "Cohort",
      fill = "Cohort"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  out <- file.path(output_dir, "hallmark_nes_profile_density.png")
  ggplot2::ggsave(out, p, width = 9, height = 5.5, dpi = 180)
  message("Wrote ", out)
  invisible(out)
}

run_gsea_hallmark_ad_and_merge_check <- function(
  output_dir = NULL,
  ad_alpha = 0.05,
  min_genes = 15L,
  n_density_pathways = 12L
) {
  if (!requireNamespace("readr", quietly = TRUE)) library(readr)
  if (!requireNamespace("dplyr", quietly = TRUE)) library(dplyr)

  base_dir <- resolve_gsea_base_dir()
  if (is.null(output_dir)) {
    output_dir <- file.path(base_dir, "analysis_output", "gsea_hallmark_merge")
  }
  dist_dir <- file.path(output_dir, "distribution_checks")
  dir.create(dist_dir, recursive = TRUE, showWarnings = FALSE)

  data <- load_hallmark_high_vs_low_data(base_dir = base_dir)
  meta_path <- file.path(base_dir, "analysis_output", "Transcriptomics_meta_hallmark_fisher.csv")

  ad_raw <- run_hallmark_logfc_ad_by_pathway(
    data,
    min_genes = min_genes,
    ad_alpha = ad_alpha,
    zscore_within_pathway = FALSE
  )
  ad_z <- run_hallmark_logfc_ad_by_pathway(
    data,
    min_genes = min_genes,
    ad_alpha = ad_alpha,
    zscore_within_pathway = TRUE
  )
  ad_raw$pass_merge_zscore = ad_z$pass_merge[match(ad_raw$pathway, ad_z$pathway)]
  ad_raw$ad_p_zscore = ad_z$ad_p[match(ad_raw$pathway, ad_z$pathway)]

  nes_ad <- run_hallmark_nes_profile_ad(data, ad_alpha = ad_alpha)

  readr::write_csv(ad_raw, file.path(output_dir, "hallmark_ad_by_pathway.csv"))
  readr::write_csv(nes_ad, file.path(output_dir, "hallmark_nes_profile_ad.csv"))

  export_hallmark_ad_summary_plot(ad_raw, output_dir, ad_alpha = ad_alpha)
  export_hallmark_nes_profile_density(data, dist_dir)

  density_pws <- pick_hallmark_density_pathways(
    ad_raw,
    meta_path = meta_path,
    n_pathways = n_density_pathways
  )
  readr::write_csv(
    data.frame(pathway = density_pws, stringsAsFactors = FALSE),
    file.path(dist_dir, "hallmark_density_pathways.csv")
  )

  path_sub <- stats::setNames(data$pathways[density_pws], density_pws)
  for (zscore in c(FALSE, TRUE)) {
    long <- build_pathway_logfc_long(
      data,
      pathways = path_sub,
      zscore_within_pathway = zscore
    )
    export_hallmark_logfc_density_plots(
      long,
      pathways = density_pws,
      output_dir = dist_dir,
      zscore_within_pathway = zscore
    )
  }

  n_pass_raw <- sum(ad_raw$pass_merge, na.rm = TRUE)
  n_pass_z <- sum(ad_raw$pass_merge_zscore, na.rm = TRUE)
  n_tested <- nrow(ad_raw)

  summary_lines <- c(
    "Hallmark High vs Low GSEA — mergeability (Anderson-Darling)",
    "",
    "Contrast: limma High vs Low DFS risk (same as run_gsea_all.R).",
    "Per-pathway AD: k-sample test on limma logFC for genes in each Hallmark set",
    "  (>= min_genes present in all four cohorts). Analogous to xCell AD on sample scores.",
    "NES profile AD: single test on the vector of all Hallmark NES values per cohort.",
    "",
    paste0("Pathways tested: ", n_tested),
    paste0("Pass merge (raw logFC AD p >= ", ad_alpha, "): ", n_pass_raw),
    paste0("Pass merge (z-scored within pathway AD p >= ", ad_alpha, "): ", n_pass_z),
    paste0(
      "NES profile pass: ",
      ifelse(isTRUE(nes_ad$pass_merge), "YES", "NO"),
      " (AD p = ",
      signif(nes_ad$ad_p, 4),
      ")"
    ),
    "",
    "Interpretation:",
    "  Pass = gene-level logFC distributions within the pathway are not significantly",
    "  different across cohorts (safer to meta-combine at pathway level).",
    "  Fail = cohort heterogeneity; inspect density plots before trusting Fisher meta.",
    "",
    "Outputs:",
    "  hallmark_ad_by_pathway.csv",
    "  hallmark_nes_profile_ad.csv",
    "  hallmark_ad_by_pathway_summary.png",
    "  distribution_checks/hallmark_logfc_density_{raw,zscore}_panel.png",
    "  distribution_checks/hallmark_nes_profile_density.png"
  )
  writeLines(summary_lines, file.path(output_dir, "hallmark_merge_summary.txt"))

  list(
    ad_by_pathway = ad_raw,
    nes_profile_ad = nes_ad,
    n_pass_raw = n_pass_raw,
    n_pass_z = n_pass_z,
    n_tested = n_tested,
    can_merge = n_pass_raw > 0,
    output_dir = output_dir
  )
}
