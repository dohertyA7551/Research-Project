# Transcriptomics (xCell, PROGENy, limma/GSEA) by TMA composition group.
source("analysis_helpers.R")
source("progeny_helpers.R")
source("xcell_combine_helpers.R")
source("gsea_helpers.R")
source("pbine_meta_helpers.R")

load_tma_comp_patient_groups <- function(
  groups_rds = "/work_space/files/proteomics/tma_composition_patient_groups.rds",
  exclude_other = TRUE
) {
  obj <- readRDS(groups_rds)
  pg <- obj$patient_groups
  if (exclude_other) {
    pg <- pg %>% dplyr::filter(.data$tma_comp_grp != "Other_mixed")
  }
  pg %>%
    dplyr::mutate(
      tma_de_extreme = dplyr::case_when(
        .data$n_tma_high == 0L ~ "all_low",
        .data$n_tma_high == .data$n_tma ~ "all_high",
        TRUE ~ NA_character_
      ),
      tma_de_pole = dplyr::case_when(
        .data$tma_comp_grp %in% c("G1_all_low", "G2_low_2plus") ~ "low_enriched",
        .data$tma_comp_grp %in% c("G3_high_2plus", "G4_high_3plus") ~ "high_enriched",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(
      .data$patient_id,
      .data$tma_comp_grp,
      .data$tma_de_extreme,
      .data$tma_de_pole,
      .data$n_tma,
      .data$n_tma_low,
      .data$n_tma_high,
      .data$rsf_risk,
      .data$risk_grp
    )
}

join_comp_groups_meta <- function(sample_meta, comp_groups) {
  cg <- comp_groups %>%
    dplyr::select(
      .data$patient_id,
      .data$tma_comp_grp,
      .data$tma_de_extreme,
      .data$tma_de_pole
    )
  sample_meta %>%
    dplyr::left_join(cg, by = "patient_id") %>%
    dplyr::filter(!is.na(.data$tma_comp_grp))
}

kruskal_by_group <- function(values, groups) {
  ok <- is.finite(values) & !is.na(groups)
  if (sum(ok) < 10L || dplyr::n_distinct(groups[ok]) < 2L) {
    return(NA_real_)
  }
  tryCatch(
    stats::kruskal.test(values[ok] ~ groups[ok])$p.value,
    error = function(e) NA_real_
  )
}

run_tma_comp_xcell_cohort <- function(scores_long, comp_groups, cohort, output_dir) {
  cg <- comp_groups %>%
    dplyr::select(.data$patient_id, .data$tma_comp_grp)
  df <- scores_long %>%
    dplyr::inner_join(cg, by = "patient_id") %>%
    dplyr::filter(.data$cohort == !!cohort | is.na(.data$cohort))
  if (!"cohort" %in% names(df) || all(is.na(df$cohort))) {
    df$cohort <- cohort
  }
  df <- df %>% dplyr::filter(.data$cohort == cohort)

  out <- df %>%
    dplyr::group_by(.data$cell_type) %>%
    dplyr::summarise(
      kruskal_p = kruskal_by_group(.data$enrichment, .data$tma_comp_grp),
      .groups = "drop"
    ) %>%
    dplyr::mutate(neg_log10_p = -log10(.data$kruskal_p), cohort = cohort)

  readr::write_csv(out, file.path(output_dir, paste0(cohort, "_xcell_kruskal_4groups.csv")))
  out
}

run_tma_comp_progeny_cohort <- function(expr_mat, sample_meta, comp_groups, cohort, output_dir) {
  meta <- join_comp_groups_meta(sample_meta, comp_groups)
  meta <- meta %>% dplyr::filter(.data$sample_id %in% colnames(expr_mat))
  expr <- expr_mat[, meta$sample_id, drop = FALSE]
  progeny <- run_progeny_mlm(expr)
  acts <- progeny_activities_to_long(progeny$activities, meta) %>%
    dplyr::left_join(
      meta %>% dplyr::select(.data$sample_id, .data$tma_comp_grp),
      by = "sample_id"
    )

  out <- acts %>%
    dplyr::group_by(.data$pathway) %>%
    dplyr::summarise(
      kruskal_p = kruskal_by_group(.data$score, .data$tma_comp_grp),
      .groups = "drop"
    ) %>%
    dplyr::mutate(neg_log10_p = -log10(.data$kruskal_p), cohort = cohort)

  readr::write_csv(out, file.path(output_dir, paste0(cohort, "_progeny_kruskal_4groups.csv")))
  readr::write_csv(
    acts,
    file.path(output_dir, paste0(cohort, "_progeny_activities_long.csv"))
  )
  out
}

run_tma_comp_limma_gsea_cohort <- function(
  expr_mat,
  sample_meta,
  comp_groups,
  cohort,
  output_dir,
  file_prefix,
  min_group = 4L
) {
  meta <- join_comp_groups_meta(sample_meta, comp_groups)
  meta <- meta %>% dplyr::filter(.data$sample_id %in% colnames(expr_mat))
  grp_tab <- table(meta$tma_comp_grp)
  if (length(grp_tab) < 2L || min(grp_tab) < min_group) {
    message("Skipping limma/GSEA ", cohort, ": insufficient samples per group.")
    return(invisible(NULL))
  }

  expr <- expr_mat[, meta$sample_id, drop = FALSE]
  comp_grp <- factor(meta$tma_comp_grp, levels = names(grp_tab)[order(names(grp_tab))])

  if (requireNamespace("limma", quietly = TRUE)) {
    design <- stats::model.matrix(~ 0 + comp_grp)
    colnames(design) <- make.names(colnames(design))
    fit <- limma::lmFit(expr, design)
    fit2 <- limma::eBayes(fit)
    ref <- colnames(design)[1L]
    contrast_specs <- setdiff(colnames(design), ref)
    if (length(contrast_specs) > 0L) {
      cn <- limma::makeContrasts(
        contrasts = paste0(contrast_specs, "-", ref),
        levels = design
      )
      fit_c <- limma::contrasts.fit(fit, cn)
      fit_c <- limma::eBayes(fit_c)
      res <- limma::topTable(fit_c, number = Inf, sort.by = "F") %>%
        tibble::rownames_to_column("gene")
      readr::write_csv(
        res,
        file.path(output_dir, paste0(file_prefix, "_limma_4groups_vs_G1.csv"))
      )
      rank_col <- if ("t" %in% names(res)) "t" else "logFC"
      ranks <- res[[rank_col]]
      names(ranks) <- res$gene
      ranks <- sort(ranks, decreasing = TRUE)
      gsea_h <- tryCatch({
        run_fgsea_preranked(ranks, hallmark_pathways())
      }, error = function(e) {
        message(e$message)
        NULL
      })
      if (!is.null(gsea_h)) {
        readr::write_csv(
          gsea_h,
          file.path(output_dir, paste0(file_prefix, "_fgsea_hallmark_4groups.csv"))
        )
      }
    }
  } else {
    message("limma not available; skipping DE/GSEA for ", cohort)
  }
  invisible(meta)
}

export_tma_comp_de_volcano <- function(
  results,
  output_dir,
  file_prefix,
  cohort_title,
  group_high_label,
  group_low_label,
  n_low,
  n_high,
  fc_cutoff = 1,
  top_labels = 12,
  fdr_cutoff = 0.05
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for volcano plots.")
  }

  prep <- prepare_limma_volcano_df(
    results,
    fc_cutoff = fc_cutoff,
    top_labels = top_labels,
    fdr_cutoff = fdr_cutoff
  )
  df <- prep$df %>%
    dplyr::mutate(
      fc_class = dplyr::case_when(
        .data$logFC >= fc_cutoff ~ paste0("Higher in ", group_high_label),
        .data$logFC <= -fc_cutoff ~ paste0("Higher in ", group_low_label),
        TRUE ~ "Other"
      )
    )
  label_df <- prep$label_df
  has_fdr_hits <- prep$has_fdr_hits
  fdr_y <- prep$fdr_y

  limits <- compute_limma_volcano_limits(
    results,
    fc_cutoff = fc_cutoff,
    top_labels = top_labels,
    fdr_cutoff = fdr_cutoff
  )

  line_caption <- if (has_fdr_hits) {
    paste0("Dashed line: BH FDR = ", fdr_cutoff)
  } else {
    paste0("No genes at BH FDR <= ", fdr_cutoff, " (FDR line omitted; y-axis scaled to data)")
  }

  vol <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$logFC,
      y = .data$neg_log10_adjP,
      colour = .data$fc_class,
      size = abs(.data$logFC)
    )
  ) +
    ggplot2::geom_point(alpha = 0.55) +
    ggplot2::scale_size_continuous(range = c(0.8, 4), guide = "none") +
    ggplot2::scale_colour_manual(
      values = stats::setNames(
        c("#B2182B", "#2166AC", "grey75"),
        c(
          paste0("Higher in ", group_high_label),
          paste0("Higher in ", group_low_label),
          "Other"
        )
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
    ggplot2::coord_cartesian(
      xlim = c(-limits$xmax, limits$xmax),
      ylim = c(0, limits$ymax),
      clip = "off",
      expand = FALSE
    ) +
    ggplot2::labs(
      title = paste0("TMA composition: ", group_high_label, " vs ", group_low_label),
      subtitle = paste0(
        cohort_title, " | n = ", n_low + n_high,
        " (", group_low_label, " ", n_low, ", ", group_high_label, " ", n_high,
        ") | Labels: top ", top_labels, " genes by |logFC|"
      ),
      caption = line_caption,
      x = paste0("log2 fold change (", group_high_label, " / ", group_low_label, ")"),
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
  }

  vol_png <- file.path(output_dir, paste0(file_prefix, "_limma_volcano.png"))
  vol_pdf <- file.path(output_dir, paste0(file_prefix, "_limma_volcano.pdf"))
  ggplot2::ggsave(vol_png, vol, width = 10, height = 7, dpi = 200, bg = "white")
  ggplot2::ggsave(vol_pdf, vol, width = 10, height = 7, bg = "white")
  message("Saved: ", vol_png)
  invisible(vol)
}

run_tma_comp_limma_two_group_cohort <- function(
  expr_mat,
  sample_meta,
  comp_groups,
  cohort,
  output_dir,
  group_col,
  level_low,
  level_high,
  file_suffix,
  group_low_label,
  group_high_label,
  min_group = 3L,
  fc_cutoff = 1,
  fdr_cutoff = 0.05,
  run_gsea = TRUE
) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    message("limma not available; skipping ", file_suffix, " for ", cohort)
    return(invisible(NULL))
  }

  meta <- join_comp_groups_meta(sample_meta, comp_groups) %>%
    dplyr::filter(.data[[group_col]] %in% c(level_low, level_high)) %>%
    dplyr::filter(.data$sample_id %in% colnames(expr_mat))

  n_low <- sum(meta[[group_col]] == level_low)
  n_high <- sum(meta[[group_col]] == level_high)
  if (n_low < min_group || n_high < min_group) {
    message(
      "Skipping ", cohort, " ", file_suffix,
      ": n_", level_low, "=", n_low, ", n_", level_high, "=", n_high
    )
    return(invisible(NULL))
  }

  expr <- expr_mat[, meta$sample_id, drop = FALSE]
  grp <- factor(
    meta[[group_col]][match(colnames(expr), meta$sample_id)],
    levels = c(level_low, level_high)
  )

  design <- stats::model.matrix(~ 0 + grp)
  colnames(design) <- c("Low", "High")
  contrast <- limma::makeContrasts(HighVsLow = High - Low, levels = design)
  fit <- limma::lmFit(expr, design)
  fit2 <- limma::contrasts.fit(fit, contrast)
  fit2 <- limma::eBayes(fit2)
  results <- limma::topTable(
    fit2,
    coef = "HighVsLow",
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  results$gene <- rownames(results)

  file_prefix <- paste0(cohort, "_", file_suffix)
  all_path <- file.path(output_dir, paste0(file_prefix, "_all_genes.csv"))
  readr::write_csv(results, all_path)
  message("Saved: ", all_path, " (", group_low_label, "=", n_low, ", ", group_high_label, "=", n_high, ")")

  sig <- results %>%
    dplyr::filter(.data$adj.P.Val <= fdr_cutoff, abs(.data$logFC) >= fc_cutoff)
  readr::write_csv(
    sig,
    file.path(output_dir, paste0(file_prefix, "_sig.csv"))
  )

  export_tma_comp_de_volcano(
    results,
    output_dir,
    file_prefix = file_prefix,
    cohort_title = cohort,
    group_high_label = group_high_label,
    group_low_label = group_low_label,
    n_low = n_low,
    n_high = n_high,
    fc_cutoff = fc_cutoff,
    fdr_cutoff = fdr_cutoff
  )

  if (isTRUE(run_gsea)) {
    ranks <- results$t
    names(ranks) <- results$gene
    ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
    gsea_h <- tryCatch(
      run_fgsea_preranked(ranks, hallmark_pathways()),
      error = function(e) {
        message(e$message)
        NULL
      }
    )
    if (!is.null(gsea_h)) {
      readr::write_csv(
        gsea_h,
        file.path(output_dir, paste0(file_prefix, "_fgsea_hallmark.csv"))
      )
    }
  }

  invisible(list(results = results, n_low = n_low, n_high = n_high))
}

run_tma_comp_de_contrasts_cohort <- function(
  expr_mat,
  sample_meta,
  comp_groups,
  cohort,
  output_dir,
  min_group = 3L
) {
  extreme <- run_tma_comp_limma_two_group_cohort(
    expr_mat, sample_meta, comp_groups, cohort, output_dir,
    group_col = "tma_de_extreme",
    level_low = "all_low",
    level_high = "all_high",
    file_suffix = "limma_all_low_vs_all_high",
    group_low_label = "all TMAs Low",
    group_high_label = "all TMAs High",
    min_group = min_group
  )
  pole <- run_tma_comp_limma_two_group_cohort(
    expr_mat, sample_meta, comp_groups, cohort, output_dir,
    group_col = "tma_de_pole",
    level_low = "low_enriched",
    level_high = "high_enriched",
    file_suffix = "limma_low_enriched_vs_high_enriched",
    group_low_label = "low-enriched (G1+G2)",
    group_high_label = "high-enriched (G3+G4)",
    min_group = min_group
  )
  invisible(list(all_low_vs_all_high = extreme, low_vs_high_enriched = pole))
}

export_tma_comp_de_meta_volcano <- function(
  meta,
  output_dir,
  file_prefix,
  group_high_label,
  group_low_label,
  n_cohorts,
  fc_cutoff = 1,
  fdr_cutoff = 0.05,
  top_labels = 15L
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for meta volcano plots.")
  }

  df <- meta %>%
    dplyr::filter(is.finite(.data$mean_logFC), is.finite(.data$meta_p)) %>%
    dplyr::mutate(
      neg_log10_fdr = -log10(pmax(.data$meta_padj, .Machine$double.xmin)),
      fc_class = dplyr::case_when(
        .data$mean_logFC >= fc_cutoff & .data$sign_consistent %in% TRUE ~
          paste0("Higher in ", group_high_label),
        .data$mean_logFC <= -fc_cutoff & .data$sign_consistent %in% TRUE ~
          paste0("Higher in ", group_low_label),
        TRUE ~ "Other"
      )
    )
  has_fdr <- any(df$meta_padj <= fdr_cutoff & df$sign_consistent %in% TRUE, na.rm = TRUE)
  label_df <- df %>%
    dplyr::filter(.data$sign_consistent %in% TRUE) %>%
    dplyr::arrange(dplyr::desc(abs(.data$mean_logFC))) %>%
    dplyr::slice_head(n = top_labels)

  limits <- compute_limma_volcano_limits(
    data.frame(
      logFC = df$mean_logFC,
      P.Value = df$meta_p,
      adj.P.Val = df$meta_padj,
      gene = df$gene
    ),
    fc_cutoff = fc_cutoff,
    top_labels = top_labels,
    fdr_cutoff = fdr_cutoff
  )

  vol <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$mean_logFC,
      y = .data$neg_log10_fdr,
      colour = .data$fc_class,
      size = abs(.data$mean_logFC)
    )
  ) +
    ggplot2::geom_point(alpha = 0.55) +
    ggplot2::scale_size_continuous(range = c(0.8, 4), guide = "none") +
    ggplot2::scale_colour_manual(
      values = stats::setNames(
        c("#B2182B", "#2166AC", "grey75"),
        c(
          paste0("Higher in ", group_high_label),
          paste0("Higher in ", group_low_label),
          "Other"
        )
      ),
      name = NULL
    ) +
    ggplot2::geom_vline(
      xintercept = c(-fc_cutoff, fc_cutoff),
      linetype = "dashed",
      colour = "black",
      linewidth = 0.4
    )

  if (has_fdr) {
    vol <- vol + ggplot2::geom_hline(
      yintercept = -log10(fdr_cutoff),
      linetype = "dashed",
      colour = "black",
      linewidth = 0.4
    )
  }

  vol <- vol +
    ggplot2::coord_cartesian(
      xlim = c(-limits$xmax, limits$xmax),
      ylim = c(0, limits$ymax),
      clip = "off",
      expand = FALSE
    ) +
    ggplot2::labs(
      title = paste0(
        "Meta-analysis (Pbine): ", group_high_label, " vs ", group_low_label
      ),
      subtitle = paste0(
        "Internal limma per cohort, then Pbine meta (method=Int) | ",
        n_cohorts, " cohorts | sign-consistent genes highlighted"
      ),
      caption = paste0(
        "x = mean log2 FC across cohorts; y = Pbine combined FDR; ",
        "requires consistent direction in contributing cohorts for coloured hits"
      ),
      x = paste0("Mean log2 FC (", group_high_label, " / ", group_low_label, ")"),
      y = expression("-log"[10] * "(Pbine meta FDR)")
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      plot.caption = ggplot2::element_text(size = 9, colour = "grey30", hjust = 0),
      legend.position = "bottom"
    )

  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(label_df) > 0L) {
    vol <- vol + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = .data$gene),
      size = 3,
      max.overlaps = 25,
      segment.size = 0.2,
      colour = "black",
      show.legend = FALSE
    )
  }

  vol_png <- file.path(output_dir, paste0(file_prefix, "_meta_volcano.png"))
  vol_pdf <- file.path(output_dir, paste0(file_prefix, "_meta_volcano.pdf"))
  ggplot2::ggsave(vol_png, vol, width = 10, height = 7, dpi = 200, bg = "white")
  ggplot2::ggsave(vol_pdf, vol, width = 10, height = 7, bg = "white")
  message("Saved: ", vol_png)
  invisible(vol)
}

compile_tma_comp_de_meta <- function(
  output_dir,
  contrast_suffix,
  meta_file_prefix,
  group_low_label,
  group_high_label,
  cohort_names = c("Stage2", "Retrospective", "Colossus", "Taxonomy"),
  min_cohorts = 2L,
  fc_cutoff = 1,
  fdr_cutoff = 0.05,
  pbine_method = "Int",
  run_gene_pbine_meta = FALSE,
  gene_rowwise_max_rows = 500L
) {
  cohort_results <- list()
  fgsea_results <- list()
  for (cn in cohort_names) {
    gene_path <- file.path(
      output_dir,
      paste0(cn, "_", contrast_suffix, "_all_genes.csv")
    )
    if (file.exists(gene_path)) {
      cohort_results[[cn]] <- readr::read_csv(gene_path, show_col_types = FALSE)
    }
    gsea_path <- file.path(
      output_dir,
      paste0(cn, "_", contrast_suffix, "_fgsea_hallmark.csv")
    )
    if (file.exists(gsea_path)) {
      fgsea_results[[cn]] <- readr::read_csv(gsea_path, show_col_types = FALSE)
    }
  }

  gsea_meta <- NULL
  if (length(fgsea_results) >= min_cohorts) {
    message("Pbine GSEA meta: ", contrast_suffix, " (", length(fgsea_results), " cohorts)")
    gsea_meta <- run_fgsea_meta_pbine(
      fgsea_results,
      output_dir = output_dir,
      file_prefix = meta_file_prefix,
      collection = "hallmark",
      fdr_cutoff = fdr_cutoff,
      method = pbine_method
    )
  } else {
    message(
      "Skipping Pbine GSEA meta for ", contrast_suffix,
      ": only ", length(fgsea_results), " cohort(s) with internal GSEA."
    )
  }

  gene_meta <- NULL
  if (isTRUE(run_gene_pbine_meta) && length(cohort_results) >= min_cohorts) {
    gene_meta <- pbine_meta_limma_genes(
      cohort_results,
      min_cohorts = min_cohorts,
      method = pbine_method,
      rowwise_max_rows = gene_rowwise_max_rows
    )
    if (!is.null(gene_meta)) {
      readr::write_csv(
        gene_meta,
        file.path(output_dir, paste0(meta_file_prefix, "_meta_pbine_all_genes.csv"))
      )
      sig <- gene_meta %>%
        dplyr::filter(
          .data$meta_padj <= fdr_cutoff,
          abs(.data$mean_logFC) >= fc_cutoff,
          .data$sign_consistent %in% TRUE
        )
      readr::write_csv(
        sig,
        file.path(output_dir, paste0(meta_file_prefix, "_meta_pbine_sig.csv"))
      )
      export_tma_comp_de_meta_volcano(
        gene_meta,
        output_dir,
        file_prefix = meta_file_prefix,
        group_high_label = group_high_label,
        group_low_label = group_low_label,
        n_cohorts = length(cohort_results),
        fc_cutoff = fc_cutoff,
        fdr_cutoff = fdr_cutoff
      )
      message(
        contrast_suffix, " gene meta: ",
        nrow(sig), " sign-consistent genes at Pbine FDR <= ", fdr_cutoff
      )
    }
  } else if (!isTRUE(run_gene_pbine_meta)) {
    message(
      "Gene-level Pbine meta skipped for ", contrast_suffix,
      " (set run_gene_pbine_meta=TRUE to enable; very slow for full transcriptomes)."
    )
  } else {
    message(
      "Skipping Pbine gene meta for ", contrast_suffix,
      ": only ", length(cohort_results), " cohort(s) with internal DE."
    )
  }

  invisible(list(genes = gene_meta, gsea = gsea_meta, n_cohorts = length(cohort_results)))
}

run_tma_comp_de_meta_all <- function(
  output_dir = "/work_space/files/Transcriptomics/analysis_output/tma_composition",
  cohort_names = c("Stage2", "Retrospective", "Colossus", "Taxonomy"),
  min_cohorts = 2L,
  fc_cutoff = 1,
  fdr_cutoff = 0.05,
  pbine_method = "Int",
  run_gene_pbine_meta = FALSE,
  gene_rowwise_max_rows = 500L
) {
  extreme <- compile_tma_comp_de_meta(
    output_dir,
    contrast_suffix = "limma_all_low_vs_all_high",
    meta_file_prefix = "meta_all_low_vs_all_high",
    group_low_label = "all TMAs Low",
    group_high_label = "all TMAs High",
    cohort_names = cohort_names,
    min_cohorts = min_cohorts,
    fc_cutoff = fc_cutoff,
    fdr_cutoff = fdr_cutoff,
    pbine_method = pbine_method,
    run_gene_pbine_meta = run_gene_pbine_meta,
    gene_rowwise_max_rows = gene_rowwise_max_rows
  )
  pole <- compile_tma_comp_de_meta(
    output_dir,
    contrast_suffix = "limma_low_enriched_vs_high_enriched",
    meta_file_prefix = "meta_low_enriched_vs_high_enriched",
    group_low_label = "low-enriched (G1+G2)",
    group_high_label = "high-enriched (G3+G4)",
    cohort_names = cohort_names,
    min_cohorts = min_cohorts,
    fc_cutoff = fc_cutoff,
    fdr_cutoff = fdr_cutoff,
    pbine_method = pbine_method,
    run_gene_pbine_meta = run_gene_pbine_meta,
    gene_rowwise_max_rows = gene_rowwise_max_rows
  )
  invisible(list(all_low_vs_all_high = extreme, low_vs_high_enriched = pole))
}

run_tma_comp_transcriptomics_all <- function(
  output_dir = "/work_space/files/Transcriptomics/analysis_output/tma_composition",
  exclude_other = TRUE
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  comp_groups <- load_tma_comp_patient_groups(exclude_other = exclude_other)
  p <- resolve_progeny_paths()
  survival_path <- file.path(p$proteomics, "survival_df_with_risk.rds")
  survival_df <- readRDS(survival_path)

  cohort_loaders <- list(
    Stage2 = function() load_stage2_expression_meta(p, survival_df),
    Retrospective = function() load_retrospective_expression_meta(p, survival_df),
    Colossus = function() load_colossus_expression_meta(p, survival_df),
    Taxonomy = function() load_taxonomy_expression_meta(p, survival_path)
  )

  progeny_results <- list()
  de_results <- list()
  for (cn in names(cohort_loaders)) {
    message("TMA comp PROGENy + GSEA: ", cn)
    dat <- cohort_loaders[[cn]]()
    progeny_results[[cn]] <- run_tma_comp_progeny_cohort(
      dat$expr, dat$meta, comp_groups, cn, output_dir
    )
    run_tma_comp_limma_gsea_cohort(
      dat$expr, dat$meta, comp_groups, cn, output_dir, cn
    )
    message("TMA comp DE contrasts: ", cn)
    de_results[[cn]] <- run_tma_comp_de_contrasts_cohort(
      dat$expr, dat$meta, comp_groups, cn, output_dir
    )
  }

  xcell_results <- list()
  xcell_loaders <- list(
    Stage2 = function() load_xcell_scores_stage2(p, survival_df),
    Retrospective = function() load_xcell_scores_retrospective(p, survival_df),
    Colossus = function() load_xcell_scores_colossus(p, survival_df),
    Taxonomy = function() load_xcell_scores_taxonomy(p, survival_df)
  )
  for (cn in names(xcell_loaders)) {
    message("TMA comp xCell: ", cn)
    scores_long <- xcell_loaders[[cn]]()
    xcell_results[[cn]] <- run_tma_comp_xcell_cohort(
      scores_long, comp_groups, cn, output_dir
    )
  }

  readr::write_csv(comp_groups, file.path(output_dir, "patient_groups_used.csv"))

  message("TMA comp DE Pbine meta (internal cohort limma -> combined)")
  meta_results <- run_tma_comp_de_meta_all(output_dir = output_dir)

  summary_lines <- c(
    "Transcriptomics by TMA composition group",
    "",
    "Method: DE is run internally within each cohort (no cross-dataset expression pooling).",
    "Cross-cohort synthesis uses Pbine (Lin et al. 2022, method=Int) on per-cohort limma P-values.",
    "For 4 cohorts: hierarchical Pbine (Stage2+Retro, Colossus+Taxonomy, then combined).",
    "Pbine source: Transcriptomics/vendor/Pbine/ (workspace-local, not system-installed).",
    "Hallmark GSEA meta uses the same Pbine approach on per-cohort fgsea P-values.",
    "Gene-level Pbine meta is off by default (very slow at full transcriptome scale).",
    "",
    paste0("Patients with composition labels (excl Other): ", nrow(comp_groups)),
    paste(capture.output(print(table(comp_groups$tma_comp_grp))), collapse = "\n"),
    "",
    "DE contrast patient counts (assignment rules):",
    paste0("  all_low vs all_high: ", sum(comp_groups$tma_de_extreme == "all_low", na.rm = TRUE),
           " vs ", sum(comp_groups$tma_de_extreme == "all_high", na.rm = TRUE)),
    paste0("  low_enriched vs high_enriched: ",
           sum(comp_groups$tma_de_pole == "low_enriched", na.rm = TRUE),
           " vs ", sum(comp_groups$tma_de_pole == "high_enriched", na.rm = TRUE)),
    "",
    "Per-cohort (internal): {cohort}_progeny_kruskal_4groups.csv, {cohort}_xcell_kruskal_4groups.csv",
    "Per-cohort DE: {cohort}_limma_all_low_vs_all_high_* and {cohort}_limma_low_enriched_vs_high_enriched_*",
    "Meta (compiled, Pbine): meta_all_low_vs_all_high_* and meta_low_enriched_vs_high_enriched_*"
  )
  writeLines(summary_lines, file.path(output_dir, "transcriptomics_summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(
    progeny = progeny_results,
    xcell = xcell_results,
    de = de_results,
    meta = meta_results
  ))
}
