# G1 (all TMAs low) vs G4 (>=3 TMAs high) transcriptomics at cohort level,
# then merge xCell / PROGENy scores after distribution checks; KEGG GSEA + Pbine meta.
source("analysis_helpers.R")
source("progeny_helpers.R")
source("gsea_helpers.R")
source("pbine_meta_helpers.R")

G1_GRP <- "G1_all_low"
G4_GRP <- "G4_high_3plus"
G1_LABEL <- "All TMAs Low"
G4_LABEL <- ">=3 TMAs High"
COHORTS <- c("Stage2", "Retrospective", "Colossus", "Taxonomy")

load_g1_g4_comp_groups <- function(
  groups_rds = "/work_space/files/proteomics/tma_composition_patient_groups.rds"
) {
  pg <- readRDS(groups_rds)$patient_groups
  pg %>%
    dplyr::filter(.data$tma_comp_grp %in% c(G1_GRP, G4_GRP)) %>%
    dplyr::select(.data$patient_id, .data$tma_comp_grp, dplyr::any_of(c("n_tma", "n_tma_low", "n_tma_high")))
}

join_g1_g4_meta <- function(sample_meta, comp_groups) {
  sample_meta %>%
    dplyr::left_join(
      comp_groups %>% dplyr::select(.data$patient_id, .data$tma_comp_grp),
      by = "patient_id"
    ) %>%
    dplyr::filter(.data$tma_comp_grp %in% c(G1_GRP, G4_GRP))
}

cohort_output_dir <- function(base_dir, cohort) {
  out <- file.path(base_dir, cohort)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

count_g1_g4 <- function(meta) {
  c(
    G1 = sum(meta$tma_comp_grp == G1_GRP, na.rm = TRUE),
    G4 = sum(meta$tma_comp_grp == G4_GRP, na.rm = TRUE)
  )
}

wilcox_g1_vs_g4 <- function(
  df,
  feature_col,
  value_col = "score",
  group_col = "tma_comp_grp",
  min_per_group = 1L
) {
  df <- df %>% dplyr::filter(.data[[group_col]] %in% c(G1_GRP, G4_GRP))
  df %>%
    dplyr::group_by(.data[[feature_col]]) %>%
    dplyr::summarise(
      n_g1 = sum(.data[[group_col]] == G1_GRP),
      n_g4 = sum(.data[[group_col]] == G4_GRP),
      mean_g1 = mean(.data[[value_col]][.data[[group_col]] == G1_GRP], na.rm = TRUE),
      mean_g4 = mean(.data[[value_col]][.data[[group_col]] == G4_GRP], na.rm = TRUE),
      diff_g4_minus_g1 = .data$mean_g4 - .data$mean_g1,
      log2FC = log2((.data$mean_g4 + 1e-10) / (.data$mean_g1 + 1e-10)),
      pval = if (n_g1 >= min_per_group && n_g4 >= min_per_group) {
        tryCatch(
          stats::wilcox.test(
            .data[[value_col]] ~ .data[[group_col]],
            exact = FALSE
          )$p.value,
          error = function(e) NA_real_
        )
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(padj = stats::p.adjust(.data$pval, method = "BH"))
}

export_g1_g4_volcano <- function(
  diff,
  feature_col,
  title,
  subtitle,
  out_path,
  x_col = "diff_g4_minus_g1",
  fc_cutoff = 0.25,
  label_top = 12L
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required.")
  }
  feats <- diff[[feature_col]]
  plot_df <- diff %>%
    dplyr::filter(is.finite(.data[[x_col]]), !is.na(.data$pval)) %>%
    dplyr::mutate(
      neg_log10_p = -log10(pmax(.data$pval, .Machine$double.xmin)),
      feature = .data[[feature_col]],
      sig = .data$pval <= 0.05 & abs(.data[[x_col]]) >= fc_cutoff
    )
  if (nrow(plot_df) == 0) {
    message("No data for volcano: ", out_path)
    return(invisible(NULL))
  }
  xmax <- min(3, max(0.5, stats::quantile(abs(plot_df[[x_col]]), 0.98, na.rm = TRUE) * 1.15))
  ymax <- max(1.5, max(plot_df$neg_log10_p, na.rm = TRUE) * 1.08)
  label_df <- plot_df %>%
    dplyr::arrange(dplyr::desc(abs(.data[[x_col]]))) %>%
    dplyr::slice_head(n = label_top)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[[x_col]], y = .data$neg_log10_p, colour = .data$sig)
  ) +
    ggplot2::geom_point(size = 2, alpha = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60")) +
    ggplot2::coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax)) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = paste0(G4_LABEL, " - ", G1_LABEL),
      y = expression(-log[10] * "(p-value)"),
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
      ggplot2::aes(label = .data$feature),
      size = 2.8,
      max.overlaps = 20,
      show.legend = FALSE
    )
  }
  ggplot2::ggsave(out_path, p, width = 8, height = 6, dpi = 150, bg = "white")
  message("Saved: ", out_path)
  invisible(p)
}

ad_test_k_samples <- function(values_by_cohort) {
  if (!requireNamespace("kSamples", quietly = TRUE)) {
    stop("Package kSamples required for Anderson-Darling tests.")
  }
  values_by_cohort <- lapply(values_by_cohort, function(v) v[is.finite(v)])
  values_by_cohort <- values_by_cohort[vapply(values_by_cohort, length, 0) > 0]
  if (length(values_by_cohort) < 2L) {
    return(list(ad = NA_real_, p = NA_real_, k = length(values_by_cohort)))
  }
  ns <- vapply(values_by_cohort, length, integer(1))
  if (sum(ns) < 4L || min(ns) < 2L) {
    return(list(ad = NA_real_, p = NA_real_, k = length(values_by_cohort)))
  }
  res <- kSamples::ad.test(values_by_cohort, method = "asymptotic")
  list(ad = as.numeric(res$ad[1, 1]), p = as.numeric(res$ad[1, 3]), k = length(values_by_cohort))
}

run_per_feature_ad <- function(
  long_df,
  feature_col,
  value_col = "score",
  cohort_col = "cohort",
  cohorts = COHORTS,
  ad_alpha = 0.05
) {
  feats <- sort(unique(long_df[[feature_col]]))
  out <- lapply(feats, function(f) {
    sub <- long_df[long_df[[feature_col]] == f, , drop = FALSE]
    vecs <- stats::setNames(
      lapply(cohorts, function(cn) sub[[value_col]][sub[[cohort_col]] == cn]),
      cohorts
    )
    vecs <- vecs[vapply(vecs, function(v) sum(is.finite(v)) > 0, logical(1))]
    tst <- ad_test_k_samples(vecs)
    data.frame(
      feature = f,
      ad_stat = tst$ad,
      ad_p = tst$p,
      n_cohorts = tst$k,
      pass_merge = is.finite(tst$p) && tst$p >= ad_alpha,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

normalize_by_cohort_quantiles <- function(x, batch) {
  x <- as.numeric(x)
  batch <- as.character(batch)
  out <- x
  ok <- is.finite(x)
  if (sum(ok) < 4L) return(out)
  idx_ok <- which(ok)
  batches <- split(idx_ok, batch[idx_ok])
  n_per <- vapply(batches, length, integer(1))
  if (length(batches) < 2L || any(n_per < 2L)) return(out)
  max_n <- max(n_per)
  probs_ref <- (seq_len(max_n) - 0.5) / max_n
  qmat <- matrix(NA_real_, nrow = max_n, ncol = length(batches))
  for (j in seq_along(batches)) {
    v <- sort(x[batches[[j]]])
    n <- length(v)
    probs_j <- (seq_len(n) - 0.5) / n
    qmat[, j] <- stats::approx(probs_j, v, xout = probs_ref, rule = 2)$y
  }
  ref <- rowMeans(qmat, na.rm = TRUE)
  for (idx in batches) {
    v <- x[idx]
    n <- length(v)
    ranks <- rank(v, ties.method = "average")
    probs_j <- (ranks - 0.5) / n
    out[idx] <- stats::approx(probs_ref, ref, xout = probs_j, rule = 2)$y
  }
  out
}

export_distribution_density <- function(
  long_df,
  feature_col,
  value_col = "score",
  cohort_col = "cohort",
  group_col = "tma_comp_grp",
  out_path,
  n_panel = 12L,
  title = "Score distributions by cohort"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  top_feats <- long_df %>%
    dplyr::count(.data[[feature_col]], name = "n") %>%
    dplyr::slice_max(.data$n, n = n_panel, with_ties = FALSE) %>%
    dplyr::pull(.data[[feature_col]])
  plot_df <- long_df %>%
    dplyr::filter(.data[[feature_col]] %in% top_feats) %>%
    dplyr::mutate(
      feature = .data[[feature_col]],
      cohort = .data[[cohort_col]],
      group = .data[[group_col]]
    )
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[[value_col]], colour = .data$cohort, fill = .data$cohort)
  ) +
    ggplot2::geom_density(alpha = 0.15, linewidth = 0.7) +
    ggplot2::facet_wrap(~feature, scales = "free", ncol = 3) +
    ggplot2::labs(title = title, x = NULL, y = "Density", colour = NULL, fill = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(size = 8),
      legend.position = "bottom"
    )
  ggplot2::ggsave(out_path, p, width = 10, height = 8, dpi = 150, bg = "white")
  message("Saved: ", out_path)
  invisible(out_path)
}

export_ad_summary_plot <- function(ad_df, out_path, title, ad_alpha = 0.05) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  plot_df <- ad_df %>%
    dplyr::filter(is.finite(.data$ad_p)) %>%
    dplyr::mutate(
      neg_log10_p = -log10(pmax(.data$ad_p, .Machine$double.xmin)),
      pass = ifelse(.data$pass_merge, "Pass", "Fail")
    )
  if (nrow(plot_df) == 0) return(invisible(NULL))
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$neg_log10_p, y = reorder(.data$feature, .data$ad_p), fill = .data$pass)
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(
      xintercept = -log10(ad_alpha),
      linetype = "dashed",
      colour = "grey40"
    ) +
    ggplot2::scale_fill_manual(values = c(Pass = "#4daf4a", Fail = "#e41a1c")) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("AD p >= ", ad_alpha, " → OK to merge across cohorts"),
      x = expression(-log[10] * "(AD p-value)"),
      y = NULL,
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  h <- max(5, min(14, 0.22 * nrow(plot_df) + 2))
  ggplot2::ggsave(out_path, p, width = 8, height = h, dpi = 150, bg = "white")
  message("Saved: ", out_path)
  invisible(out_path)
}

run_xcell2_on_expr <- function(expr_mat, min_shared = 0.5) {
  if (!requireNamespace("xCell2", quietly = TRUE)) {
    stop("Package xCell2 required.")
  }
  data("PanCancer.xCell2Ref", package = "xCell2", envir = environment())
  ref_genes <- xCell2::getGenesUsed(PanCancer.xCell2Ref)
  overlap <- length(intersect(rownames(expr_mat), ref_genes)) / length(ref_genes)
  if (overlap < min_shared) {
    warning("Gene overlap with xCell2 ref is low: ", round(overlap, 3))
  }
  xCell2::xCell2Analysis(mix = expr_mat, xcell2object = PanCancer.xCell2Ref)
}

xcell_wide_to_long <- function(xcell_wide, meta) {
  # xcell_wide: cell types x samples (matrix or df)
  if (is.matrix(xcell_wide)) {
    df <- as.data.frame(xcell_wide, check.names = FALSE)
    df$cell_type <- rownames(xcell_wide)
    long <- tidyr::pivot_longer(df, -cell_type, names_to = "sample_id", values_to = "score")
  } else {
    stop("Expected matrix from xCell2.")
  }
  long %>%
    dplyr::inner_join(
      meta %>% dplyr::select(.data$sample_id, .data$patient_id, .data$tma_comp_grp, .data$cohort),
      by = "sample_id"
    )
}

run_cohort_xcell2 <- function(expr_mat, meta, cohort, out_dir, force = FALSE) {
  csv_path <- file.path(out_dir, paste0(cohort, "_xcell2_scores.csv"))
  if (file.exists(csv_path) && !force) {
    message("Loading cached xCell2: ", csv_path)
    wide <- readr::read_csv(csv_path, show_col_types = FALSE)
    mat <- as.matrix(wide[, setdiff(names(wide), "cell_type"), drop = FALSE])
    rownames(mat) <- wide$cell_type
    return(xcell_wide_to_long(mat, meta))
  }
  message("Running xCell2: ", cohort, " (", ncol(expr_mat), " samples)")
  meta_f <- meta %>% dplyr::filter(.data$sample_id %in% colnames(expr_mat))
  expr <- expr_mat[, meta_f$sample_id, drop = FALSE]
  res <- run_xcell2_on_expr(expr)
  wide <- as.data.frame(res, check.names = FALSE)
  wide$cell_type <- rownames(res)
  readr::write_csv(wide, csv_path)
  xcell_wide_to_long(res, meta_f)
}

run_cohort_progeny <- function(expr_mat, meta, cohort, out_dir) {
  meta_f <- meta %>% dplyr::filter(.data$sample_id %in% colnames(expr_mat))
  expr <- expr_mat[, meta_f$sample_id, drop = FALSE]
  pg <- run_progeny_mlm(expr)
  long <- pg$activities %>%
    dplyr::rename(
      sample_id = .data$condition,
      pathway = .data$source,
      score = .data$score
    ) %>%
    dplyr::inner_join(
      meta_f %>% dplyr::select(.data$sample_id, .data$patient_id, .data$tma_comp_grp, .data$cohort),
      by = "sample_id"
    )
  readr::write_csv(
    long,
    file.path(out_dir, paste0(cohort, "_progeny_activities.csv"))
  )
  long
}

run_cohort_limma_g1_g4 <- function(expr_mat, meta, cohort, out_dir, min_group = 3L) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma required.")
  }
  meta_f <- meta %>% dplyr::filter(.data$sample_id %in% colnames(expr_mat))
  n <- count_g1_g4(meta_f)
  if (n["G1"] < min_group || n["G4"] < min_group) {
    message("Skip limma ", cohort, ": G1=", n["G1"], " G4=", n["G4"])
    return(NULL)
  }
  expr <- expr_mat[, meta_f$sample_id, drop = FALSE]
  grp <- factor(meta_f$tma_comp_grp, levels = c(G1_GRP, G4_GRP))
  design <- stats::model.matrix(~ 0 + grp)
  colnames(design) <- c("G1", "G4")
  contrast <- limma::makeContrasts(G4vsG1 = G4 - G1, levels = design)
  fit <- limma::lmFit(expr, design)
  fit2 <- limma::contrasts.fit(fit, contrast)
  fit2 <- limma::eBayes(fit2)
  results <- limma::topTable(fit2, coef = "G4vsG1", number = Inf, adjust.method = "BH", sort.by = "P")
  results$gene <- rownames(results)
  out_csv <- file.path(out_dir, paste0(cohort, "_limma_g1_vs_g4_all_genes.csv"))
  readr::write_csv(results, out_csv)
  message("Saved: ", out_csv, " (G1=", n["G1"], ", G4=", n["G4"], ")")

  export_g1_g4_volcano(
    results %>% dplyr::mutate(
      diff_g4_minus_g1 = .data$logFC,
      pval = .data$P.Value,
      feature = .data$gene
    ),
    "feature",
    paste0(cohort, ": G1 vs G4 DE"),
    paste0("G1=", n["G1"], ", G4=", n["G4"]),
    file.path(out_dir, paste0(cohort, "_limma_g1_vs_g4_volcano.png")),
    x_col = "diff_g4_minus_g1",
    fc_cutoff = 1
  )
  list(results = results, n_g1 = n["G1"], n_g4 = n["G4"])
}

run_cohort_fgsea_kegg <- function(limma_results, cohort, out_dir) {
  ranks <- limma_results$t
  names(ranks) <- limma_results$gene
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  pw <- kegg_pathways()
  fg <- run_fgsea_preranked(ranks, pw)
  export_fgsea_tables(fg, out_dir, cohort, collection = "kegg", fdr_cutoff = 0.05)
  export_fgsea_dotplot(
    fg,
    out_dir,
    cohort,
    collection = "kegg",
    cohort_title = paste0(cohort, " | G1 vs G4"),
    fdr_cutoff = 0.05
  )
  fg
}

merge_scores_quantile <- function(long_df, feature_col, pass_features, value_col = "score") {
  long_df %>%
    dplyr::filter(.data[[feature_col]] %in% pass_features) %>%
    dplyr::group_by(.data[[feature_col]]) %>%
    dplyr::mutate(
      score_qn = normalize_by_cohort_quantiles(.data[[value_col]], .data$cohort)
    ) %>%
    dplyr::ungroup()
}

run_merged_analysis <- function(
  long_qn,
  feature_col,
  value_col = "score_qn",
  out_dir,
  prefix,
  title_prefix
) {
  diff <- wilcox_g1_vs_g4(long_qn, feature_col, value_col = value_col)
  readr::write_csv(diff, file.path(out_dir, paste0(prefix, "_g1_vs_g4_wilcoxon.csv")))
  n <- count_g1_g4(long_qn)
  export_g1_g4_volcano(
    diff %>% dplyr::mutate(feature = .data[[feature_col]]),
    "feature",
    paste0(title_prefix, " (merged cohorts)"),
    paste0("G1=", n["G1"], ", G4=", n["G4"], " | quantile-normalized"),
    file.path(out_dir, paste0(prefix, "_g1_vs_g4_volcano.png"))
  )
  diff
}

run_tma_g1_vs_g4_transcriptomics_all <- function(
  output_dir = "/work_space/files/Transcriptomics/analysis_output/tma_g1_vs_g4",
  ad_alpha = 0.05,
  min_group = 2L,
  force_xcell2 = FALSE,
  cohorts = COHORTS
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dist_dir <- file.path(output_dir, "distribution_checks")
  merged_dir <- file.path(output_dir, "merged")
  meta_dir <- file.path(output_dir, "meta")
  dir.create(dist_dir, recursive = FALSE, showWarnings = FALSE)
  dir.create(merged_dir, recursive = FALSE, showWarnings = FALSE)
  dir.create(meta_dir, recursive = FALSE, showWarnings = FALSE)

  p <- resolve_progeny_paths()
  surv <- readRDS(file.path(p$proteomics, "survival_df_with_risk.rds"))
  comp_groups <- load_g1_g4_comp_groups()
  readr::write_csv(comp_groups, file.path(output_dir, "g1_g4_patients_used.csv"))

  loaders <- list(
    Stage2 = function() load_stage2_expression_meta(p, surv),
    Retrospective = function() load_retrospective_expression_meta(p, surv),
    Colossus = function() load_colossus_expression_meta(p, surv),
    Taxonomy = function() load_taxonomy_expression_meta(p, file.path(p$proteomics, "survival_df_with_risk.rds"))
  )

  xcell_long <- list()
  progeny_long <- list()
  fgsea_list <- list()
  cohort_n <- list()

  for (cn in cohorts) {
    message("\n========== ", cn, " ==========")
    out_c <- cohort_output_dir(output_dir, cn)
    dat <- loaders[[cn]]()
    meta <- join_g1_g4_meta(dat$meta, comp_groups)
    n <- count_g1_g4(meta)
    cohort_n[[cn]] <- n
    message(cn, " G1=", n["G1"], " G4=", n["G4"])

    if (n["G1"] < 1L || n["G4"] < 1L) {
      message("Skipping ", cn, " (no G1 or G4 RNA samples).")
      next
    }

    if (n["G4"] < min_group) {
      message(
        "Including ", cn, " for xCell/PROGENy (contributes to merged analysis); ",
        "per-cohort limma/GSEA skipped (G4 n=", n["G4"], " < ", min_group, ")."
      )
    }

    expr <- dat$expr[, meta$sample_id, drop = FALSE]

    xc <- run_cohort_xcell2(expr, meta, cn, out_c, force = force_xcell2 || cn == "Taxonomy")
    xcell_long[[cn]] <- xc
    xdiff <- wilcox_g1_vs_g4(xc, "cell_type", value_col = "score")
    readr::write_csv(xdiff, file.path(out_c, paste0(cn, "_xcell_g1_vs_g4.csv")))
    export_g1_g4_volcano(
      xdiff %>% dplyr::mutate(feature = .data$cell_type),
      "feature",
      paste0(cn, ": xCell2"),
      paste0("G1=", n["G1"], ", G4=", n["G4"]),
      file.path(out_c, paste0(cn, "_xcell_g1_vs_g4_volcano.png"))
    )

    pg <- run_cohort_progeny(expr, meta, cn, out_c)
    progeny_long[[cn]] <- pg
    pdiff <- wilcox_g1_vs_g4(pg, "pathway", value_col = "score")
    readr::write_csv(pdiff, file.path(out_c, paste0(cn, "_progeny_g1_vs_g4.csv")))
    export_g1_g4_volcano(
      pdiff %>% dplyr::mutate(feature = .data$pathway),
      "feature",
      paste0(cn, ": PROGENy"),
      paste0("G1=", n["G1"], ", G4=", n["G4"]),
      file.path(out_c, paste0(cn, "_progeny_g1_vs_g4_volcano.png"))
    )

    lim <- NULL
    if (n["G1"] >= min_group && n["G4"] >= min_group) {
      lim <- run_cohort_limma_g1_g4(expr, meta, cn, out_c, min_group = min_group)
      if (!is.null(lim)) {
        fg <- run_cohort_fgsea_kegg(lim$results, cn, out_c)
        fgsea_list[[cn]] <- fg
      }
    }
  }

  if (length(xcell_long) >= 2L) {
    xcell_all <- dplyr::bind_rows(xcell_long)
    readr::write_csv(xcell_all, file.path(dist_dir, "xcell_all_cohorts_long.csv"))
    x_ad <- run_per_feature_ad(xcell_all, "cell_type", value_col = "score", ad_alpha = ad_alpha)
    readr::write_csv(x_ad, file.path(dist_dir, "xcell_ad_by_cell_type.csv"))
    export_ad_summary_plot(
      x_ad,
      file.path(dist_dir, "xcell_ad_summary.png"),
      "xCell2: cohort distribution similarity (Anderson-Darling)",
      ad_alpha = ad_alpha
    )
    export_distribution_density(
      xcell_all,
      "cell_type",
      out_path = file.path(dist_dir, "xcell_density_panel.png"),
      title = "xCell2 scores by cohort (top cell types)"
    )
    pass_x <- x_ad$feature[x_ad$pass_merge]
    if (length(pass_x) >= 1L) {
      x_qn <- merge_scores_quantile(xcell_all, "cell_type", pass_x)
      readr::write_csv(x_qn, file.path(merged_dir, "xcell_quantile_normalized.csv"))
      run_merged_analysis(
        x_qn, "cell_type", out_dir = merged_dir, prefix = "xcell",
        title_prefix = "xCell2"
      )
    } else {
      message("Too few xCell types passed AD (", length(pass_x), "); skip merged xCell.")
    }
  }

  if (length(progeny_long) >= 2L) {
    pg_all <- dplyr::bind_rows(progeny_long)
    readr::write_csv(pg_all, file.path(dist_dir, "progeny_all_cohorts_long.csv"))
    p_ad <- run_per_feature_ad(pg_all, "pathway", value_col = "score", ad_alpha = ad_alpha)
    readr::write_csv(p_ad, file.path(dist_dir, "progeny_ad_by_pathway.csv"))
    export_ad_summary_plot(
      p_ad,
      file.path(dist_dir, "progeny_ad_summary.png"),
      "PROGENy: cohort distribution similarity (Anderson-Darling)",
      ad_alpha = ad_alpha
    )
    export_distribution_density(
      pg_all,
      "pathway",
      out_path = file.path(dist_dir, "progeny_density_panel.png"),
      title = "PROGENy pathway scores by cohort"
    )
    pass_p <- p_ad$feature[p_ad$pass_merge]
    if (length(pass_p) >= 1L) {
      p_qn <- merge_scores_quantile(pg_all, "pathway", pass_p)
      readr::write_csv(p_qn, file.path(merged_dir, "progeny_quantile_normalized.csv"))
      run_merged_analysis(
        p_qn, "pathway", out_dir = merged_dir, prefix = "progeny",
        title_prefix = "PROGENy"
      )
    } else {
      message("Too few PROGENy pathways passed AD (", length(pass_p), "); skip merged PROGENy.")
    }
  }

  if (length(fgsea_list) >= 2L) {
    message("\nPbine KEGG GSEA meta (", length(fgsea_list), " cohorts)")
    meta_fg <- run_fgsea_meta_pbine(
      fgsea_list,
      output_dir = meta_dir,
      file_prefix = "g1_vs_g4",
      collection = "kegg",
      fdr_cutoff = 0.05
    )
  } else {
    meta_fg <- NULL
    message("Skip KEGG Pbine meta: need >=2 cohorts with internal GSEA.")
  }

  summary_lines <- c(
    "TMA composition transcriptomics: G1 (all TMAs low) vs G4 (>=3 TMAs high)",
    "",
    paste0(
      "Design: per-cohort xCell2, PROGENy; limma+KEGG fgsea when >=",
      min_group,
      " samples per group; merge scores after AD checks; Pbine on KEGG p-values."
    ),
    paste0("Output root: ", output_dir),
    "",
    "RNA samples per cohort (G1 / G4):",
    paste0("  ", names(cohort_n), ": ", vapply(cohort_n, function(x) paste(x, collapse = "/"), character(1)), collapse = "\n"),
    "",
    "Folders:",
    "  {cohort}/           per-cohort xCell, PROGENy, limma, KEGG fgsea",
    "  distribution_checks/  AD tests + density panels",
    "  merged/             quantile-normalized pooled Wilcoxon (features passing AD)",
    "  meta/               Pbine KEGG GSEA meta",
    "",
    paste0("AD merge threshold: p >= ", ad_alpha)
  )
  writeLines(summary_lines, file.path(output_dir, "summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(
    xcell = xcell_long,
    progeny = progeny_long,
    fgsea = fgsea_list,
    meta = meta_fg,
    cohort_n = cohort_n
  ))
}
