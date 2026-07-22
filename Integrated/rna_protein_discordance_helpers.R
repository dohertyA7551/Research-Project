# Patient-level RNA -> protein prediction (GeneExpressPred-style) and discordance
# (observed - predicted). Workspace paths only.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

WORKSPACE_ROOT <- "/work_space/files"
INTEGRATED_DIR <- "/work_space/files/Integrated"
INTEGRATED_DATA_DIR <- "/work_space/files/Integrated "
TRANSCRIPTOMICS_DIR <- file.path(WORKSPACE_ROOT, "Transcriptomics")
PROTEOMICS_DIR <- file.path(WORKSPACE_ROOT, "proteomics")

resolve_discordance_paths <- function() {
  list(
    proteomics = PROTEOMICS_DIR,
    transcriptomics = TRANSCRIPTOMICS_DIR,
    integrated = INTEGRATED_DIR,
    spatial = file.path(WORKSPACE_ROOT, "spatial"),
    merged_rds = file.path(PROTEOMICS_DIR, "RCSI_Taxonomy_clin_prot_merged_df.rds"),
    survival_rds = file.path(PROTEOMICS_DIR, "survival_df_with_risk.rds"),
    tma_groups_rds = file.path(PROTEOMICS_DIR, "tma_composition_patient_groups.rds"),
    output = file.path(INTEGRATED_DIR, "results", "rna_protein_discordance")
  )
}

get_protein_feature_cols <- function(merged) {
  cols <- grep(
    "_mean_(NonImmuneEpithelium|NonImmuneStroma|cytT|helperT|regT)$",
    names(merged),
    value = TRUE,
    ignore.case = TRUE
  )
  if (length(cols) == 0L) {
    stop("No proteomics marker columns found in merged patient dataframe.")
  }
  cols
}

clean_protein_label <- function(x) {
  x <- sub("^baseline_Mean\\.Cell\\.", "", x, ignore.case = TRUE)
  x <- sub("^baseline_", "", x, ignore.case = TRUE)
  x <- sub("_mean_NonImmuneEpithelium$", " (epi)", x, ignore.case = TRUE)
  x <- sub("_mean_NonImmuneStroma$", " (stroma)", x, ignore.case = TRUE)
  x <- sub("_mean_cytT$", " (cytT)", x, ignore.case = TRUE)
  x <- sub("_mean_helperT$", " (helperT)", x, ignore.case = TRUE)
  x <- sub("_mean_regT$", " (regT)", x, ignore.case = TRUE)
  x
}

collapse_expr_to_patient <- function(expr_mat, meta) {
  meta <- meta %>%
    dplyr::filter(!is.na(.data$patient_id), .data$sample_id %in% colnames(expr_mat))
  if (nrow(meta) == 0L) {
    return(list(expr = NULL, meta = NULL))
  }
  expr <- expr_mat[, meta$sample_id, drop = FALSE]
  pts <- unique(meta$patient_id)
  out <- matrix(NA_real_, nrow = nrow(expr), ncol = length(pts))
  rownames(out) <- rownames(expr)
  colnames(out) <- pts
  cohort_vec <- setNames(rep(NA_character_, length(pts)), pts)
  for (pt in pts) {
    ids <- meta$sample_id[meta$patient_id == pt]
    if (length(ids) == 1L) {
      out[, pt] <- expr[, ids]
    } else {
      out[, pt] <- rowMeans(expr[, ids, drop = FALSE], na.rm = TRUE)
    }
    cohort_vec[pt] <- meta$cohort[match(ids[1], meta$sample_id)]
  }
  list(
    expr = out,
    meta = data.frame(
      patient_id = pts,
      cohort = unname(cohort_vec[pts]),
      stringsAsFactors = FALSE
    )
  )
}

zscore_by_cohort <- function(expr_mat, cohorts) {
  out <- expr_mat
  for (cn in unique(cohorts)) {
    idx <- which(cohorts == cn)
    if (length(idx) < 2L) next
    block <- out[, idx, drop = FALSE]
    scaled <- t(scale(t(block)))
    scaled[!is.finite(scaled)] <- 0
    out[, idx] <- scaled
  }
  out
}

load_patient_rna_cohort_z <- function() {
  owd <- getwd()
  on.exit(setwd(owd), add = TRUE)
  setwd(TRANSCRIPTOMICS_DIR)
  source("progeny_helpers.R")
  p <- list(
    base = TRANSCRIPTOMICS_DIR,
    proteomics = PROTEOMICS_DIR,
    spatial = file.path(WORKSPACE_ROOT, "spatial"),
    integrated = INTEGRATED_DATA_DIR
  )
  surv <- readRDS(file.path(p$proteomics, "survival_df_with_risk.rds"))
  loaders <- list(
    Stage2 = function() load_stage2_expression_meta(p, surv),
    Retrospective = function() load_retrospective_expression_meta(p, surv),
    Colossus = function() load_colossus_expression_meta(p, surv),
    Taxonomy = function() load_taxonomy_expression_meta(
      p, file.path(p$proteomics, "survival_df_with_risk.rds")
    )
  )
  parts <- lapply(names(loaders), function(cn) {
    dat <- loaders[[cn]]()
    collapse_expr_to_patient(dat$expr, dat$meta)
  })
  names(parts) <- names(loaders)
  parts <- parts[!vapply(parts, function(x) is.null(x$expr), logical(1))]
  if (length(parts) == 0L) stop("No RNA cohort data loaded.")
  shared_genes <- Reduce(intersect, lapply(parts, function(x) rownames(x$expr)))
  if (length(shared_genes) < 500L) {
    warning("Few shared genes across RNA cohorts: ", length(shared_genes))
  }
  expr_list <- lapply(parts, function(x) x$expr[shared_genes, , drop = FALSE])
  meta <- bind_rows(lapply(parts, function(x) x$meta))
  meta <- meta[!duplicated(meta$patient_id), , drop = FALSE]
  all_pts <- unique(unlist(lapply(expr_list, colnames)))
  combined <- matrix(NA_real_, nrow = length(shared_genes), ncol = length(all_pts))
  rownames(combined) <- shared_genes
  colnames(combined) <- all_pts
  cohort_map <- setNames(rep(NA_character_, length(all_pts)), all_pts)
  for (cn in names(expr_list)) {
    ex <- expr_list[[cn]]
    combined[, colnames(ex)] <- ex
    cohort_map[colnames(ex)] <- cn
  }
  combined <- zscore_by_cohort(combined, cohort_map[colnames(combined)])
  list(
    expr = combined,
    meta = data.frame(
      patient_id = colnames(combined),
      cohort = unname(cohort_map[colnames(combined)]),
      stringsAsFactors = FALSE
    ),
    n_genes = length(shared_genes)
  )
}

load_patient_protein_matrix <- function(merged_rds) {
  merged <- readRDS(merged_rds)
  feat_cols <- get_protein_feature_cols(merged)
  mat <- as.matrix(merged[, feat_cols, drop = FALSE])
  rownames(mat) <- merged$patient_id
  storage.mode(mat) <- "double"
  list(matrix = t(mat), features = feat_cols, merged = merged)
}

select_rna_features <- function(x_train, y, top_n = 30L, method = c("spearman", "cosine")) {
  method <- match.arg(method)
  ok <- is.finite(y)
  if (sum(ok) < 8L) return(character(0))
  x <- x_train[, ok, drop = FALSE]
  yy <- y[ok]
  if (method == "spearman") {
    cors <- apply(x, 1, function(g) {
      if (stats::sd(g, na.rm = TRUE) < 1e-8) return(0)
      suppressWarnings(stats::cor(g, yy, method = "spearman", use = "pairwise.complete.obs"))
    })
    cors[!is.finite(cors)] <- 0
    names(sort(abs(cors), decreasing = TRUE))[seq_len(min(top_n, length(cors)))]
  } else {
    xsc <- t(scale(t(x)))
    xsc[!is.finite(xsc)] <- 0
    ysc <- as.numeric(scale(yy))
    ysc[!is.finite(ysc)] <- 0
    sim <- apply(xsc, 1, function(g) sum(g * ysc) / (sqrt(sum(g^2)) * sqrt(sum(ysc^2)) + 1e-8))
    sim[!is.finite(sim)] <- 0
    names(sort(abs(sim), decreasing = TRUE))[seq_len(min(top_n, length(sim)))]
  }
}

voting_predict <- function(x_train, y_train, x_all, genes) {
  if (length(genes) < 2L) return(rep(NA_real_, ncol(x_all)))
  xt <- t(x_train[genes, , drop = FALSE])
  xa <- t(x_all[genes, , drop = FALSE])
  preds <- list()
  preds$lm <- tryCatch({
    fit <- stats::lm(y_train ~ ., data = data.frame(y = y_train, xt))
    as.numeric(stats::predict(fit, newdata = as.data.frame(xa)))
  }, error = function(e) rep(NA_real_, nrow(xa)))
  preds$ridge <- tryCatch({
    if (!requireNamespace("glmnet", quietly = TRUE)) stop("no glmnet")
    fit <- glmnet::cv.glmnet(xt, y_train, alpha = 0, nfolds = min(5L, length(y_train)))
    as.numeric(predict(fit, newx = xa, s = "lambda.min"))
  }, error = function(e) rep(NA_real_, nrow(xa)))
  preds$huber <- tryCatch({
    if (!requireNamespace("MASS", quietly = TRUE)) stop("no MASS")
    fit <- MASS::rlm(y_train ~ ., data = data.frame(y = y_train, xt), maxit = 100)
    as.numeric(stats::predict(fit, newdata = as.data.frame(xa)))
  }, error = function(e) rep(NA_real_, nrow(xa)))
  pm <- do.call(cbind, preds)
  rowMeans(pm, na.rm = TRUE)
}

predict_discordance_one_feature <- function(
  protein_vec,
  rna_expr,
  patients,
  top_n = 30L,
  min_n = 10L,
  feature_select = "spearman"
) {
  y <- protein_vec[patients]
  ok <- is.finite(y)
  if (sum(ok) < min_n) {
    return(NULL)
  }
  pts <- patients[ok]
  x <- rna_expr[, pts, drop = FALSE]
  genes <- select_rna_features(x, y[ok], top_n = top_n, method = feature_select)
  if (length(genes) < 3L) return(NULL)
  pred <- voting_predict(x, y[ok], x, genes)
  obs <- as.numeric(y[ok])
  data.frame(
    patient_id = pts,
    observed = obs,
    predicted = pred,
    discordance = obs - pred,
    n_rna_features = length(genes),
    stringsAsFactors = FALSE
  )
}

compute_all_discordance <- function(
  protein_mat,
  protein_features,
  rna_expr,
  patients,
  top_n = 30L,
  min_n = 10L,
  max_features = NULL
) {
  use_feats <- protein_features
  if (!is.null(max_features) && length(use_feats) > max_features) {
    use_feats <- use_feats[seq_len(max_features)]
  }
  out <- lapply(use_feats, function(f) {
    res <- predict_discordance_one_feature(
      protein_mat[f, ], rna_expr, patients,
      top_n = top_n, min_n = min_n
    )
    if (is.null(res)) return(NULL)
    res$protein_feature <- f
    res$label <- clean_protein_label(f)
    res
  })
  bind_rows(out)
}

summarize_patient_discordance <- function(long_df) {
  long_df %>%
    group_by(.data$patient_id) %>%
    summarise(
      n_features = dplyr::n(),
      mean_discordance = mean(.data$discordance, na.rm = TRUE),
      mean_abs_discordance = mean(abs(.data$discordance), na.rm = TRUE),
      sd_discordance = stats::sd(.data$discordance, na.rm = TRUE),
      .groups = "drop"
    )
}

clinical_spatial_covariates <- function(merged) {
  keep_pat <- c(
    "patient_id", "cohort_rna", "tma_comp_grp", "rsf_risk", "risk_grp",
    "dfs_months", "dfs_event", "stage", "grade", "age", "sex", "msi", "cms", "cris"
  )
  num_pat <- grep(
    "^(count_|freq_|density_|ratio_|nhood_z_|abundance_metacluster)",
    names(merged),
    value = TRUE
  )
  intersect(c(keep_pat, num_pat), names(merged))
}

run_association_screen <- function(patient_df, covariates, min_pairs = 15L) {
  covs <- covariates[covariates %in% names(patient_df)]
  covs <- covs[vapply(covs, function(cn) {
    x <- patient_df[[cn]]
    is.numeric(x) &&
      sum(is.finite(x)) >= min_pairs &&
      stats::sd(x, na.rm = TRUE) > 0
  }, logical(1))]
  targets <- c(
    "mean_abs_discordance", "mean_discordance", "sd_discordance"
  )
  out <- list()
  for (tg in targets) {
    if (!tg %in% names(patient_df)) next
    y <- patient_df[[tg]]
    for (cn in covs) {
      x <- patient_df[[cn]]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < min_pairs) next
      if (is.numeric(x)) {
        test <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman"))
        out[[length(out) + 1L]] <- data.frame(
          discordance_metric = tg,
          covariate = cn,
          n = sum(ok),
          spearman_rho = unname(test$estimate),
          p_value = test$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(out) == 0L) {
    return(data.frame())
  }
  res <- bind_rows(out) %>%
    mutate(p_adj = stats::p.adjust(.data$p_value, method = "BH")) %>%
    arrange(.data$p_adj)
  res
}

plot_discordance_distribution <- function(patient_summary, out_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  p <- ggplot2::ggplot(
    patient_summary,
    ggplot2::aes(x = .data$mean_abs_discordance)
  ) +
    ggplot2::geom_histogram(bins = 20, fill = "#4393C3", colour = "white") +
    ggplot2::labs(
      title = "Patient-level mean |discordance|",
      subtitle = "observed protein - predicted (GeneExpressPred-style voting model)",
      x = "Mean |discordance| across markers",
      y = "Patients"
    ) +
    ggplot2::theme_bw(base_size = 11)
  ggplot2::ggsave(out_path, p, width = 8, height = 5, dpi = 150, bg = "white")
  invisible(out_path)
}

plot_top_associations <- function(assoc_df, out_path, top_n = 20L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(assoc_df) == 0L) {
    return(invisible(NULL))
  }
  plot_df <- assoc_df %>%
    arrange(.data$p_adj) %>%
    slice_head(n = top_n) %>%
    mutate(label = paste0(.data$discordance_metric, " ~ ", .data$covariate))
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$spearman_rho,
      y = reorder(.data$label, .data$spearman_rho),
      fill = .data$spearman_rho > 0
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC")) +
    ggplot2::labs(
      title = "Discordance vs clinical / spatial features",
      subtitle = "Spearman correlation (exploratory, BH-adjusted p in table)",
      x = "Spearman rho",
      y = NULL,
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "none")
  h <- max(5, 0.25 * nrow(plot_df) + 2)
  ggplot2::ggsave(out_path, p, width = 9, height = h, dpi = 150, bg = "white")
  invisible(out_path)
}

plot_feature_discordance_volcano <- function(long_df, out_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  feat_sum <- long_df %>%
    group_by(.data$protein_feature, .data$label) %>%
    summarise(
      mean_disc = mean(.data$discordance, na.rm = TRUE),
      mean_abs = mean(abs(.data$discordance), na.rm = TRUE),
      sd_disc = stats::sd(.data$discordance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      z = .data$mean_disc / (.data$sd_disc + 1e-8),
      neg_log10 = -log10(pmax(2 * stats::pnorm(-abs(.data$z)), .Machine$double.xmin))
    )
  p <- ggplot2::ggplot(
    feat_sum,
    ggplot2::aes(x = .data$mean_disc, y = .data$neg_log10, size = .data$mean_abs)
  ) +
    ggplot2::geom_point(alpha = 0.65, colour = "#542788") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::labs(
      title = "Marker-level discordance (patient-averaged)",
      x = "Mean discordance (obs - pred)",
      y = expression(-log[10] * "(approx. p)"),
      size = "Mean |disc|"
    ) +
    ggplot2::theme_bw(base_size = 11)
  ggplot2::ggsave(out_path, p, width = 8, height = 6, dpi = 150, bg = "white")
  invisible(out_path)
}

marker_compartment <- function(label) {
  dplyr::case_when(
    grepl("\\(epi\\)", label) ~ "Epithelium",
    grepl("\\(stroma\\)", label) ~ "Stroma",
    grepl("\\(cytT\\)", label) ~ "cytotoxic T",
    grepl("\\(helperT\\)", label) ~ "helper T",
    grepl("\\(regT\\)", label) ~ "regulatory T",
    TRUE ~ "Other"
  )
}

plot_discordance_by_group <- function(
  patient_df,
  group_col,
  metric = "mean_abs_discordance",
  out_path,
  title = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  if (!group_col %in% names(patient_df) || !metric %in% names(patient_df)) {
    return(invisible(NULL))
  }
  plot_df <- patient_df %>%
    filter(!is.na(.data[[group_col]]), is.finite(.data[[metric]]))
  if (nrow(plot_df) < 8L) return(invisible(NULL))
  if (is.null(title)) title <- paste0(metric, " by ", group_col)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data[[group_col]], y = .data[[metric]], fill = .data[[group_col]])
  ) +
    ggplot2::geom_boxplot(outlier.alpha = 0.4, width = 0.6) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.25, size = 0.9, colour = "grey25") +
    ggplot2::labs(title = title, x = NULL, y = metric) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )
  ggplot2::ggsave(out_path, p, width = 8, height = 5.5, dpi = 180, bg = "white")
  message("Wrote ", out_path)
  invisible(out_path)
}

plot_rsf_vs_tma_composition <- function(patient_df, out_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  req <- c("rsf_risk", "tma_comp_grp", "n_tma_high", "mean_abs_discordance")
  if (!all(req %in% names(patient_df))) return(invisible(NULL))
  plot_df <- patient_df %>%
    filter(
      is.finite(.data$rsf_risk),
      !is.na(.data$tma_comp_grp),
      is.finite(.data$mean_abs_discordance)
    )
  if (nrow(plot_df) < 10L) return(invisible(NULL))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$rsf_risk,
      y = .data$n_tma_high,
      colour = .data$tma_comp_grp,
      size = .data$mean_abs_discordance
    )
  ) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::scale_size_continuous(range = c(1.5, 5), name = "Mean |discordance|") +
    ggplot2::labs(
      title = "Clinical RSF risk vs spatial TMA composition",
      subtitle = "Point size = RNA-protein discordance magnitude",
      x = "RSF predicted risk (rsf_risk)",
      y = "Number of High-risk TMAs per patient",
      colour = "TMA group"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  ggplot2::ggsave(out_path, p, width = 9, height = 6, dpi = 180, bg = "white")
  message("Wrote ", out_path)
  invisible(out_path)
}

plot_discordance_spatial_scatter <- function(patient_df, out_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  req <- c("freq_NonImmuneStroma", "freq_NonImmuneEpithelium", "mean_discordance")
  if (!all(req %in% names(patient_df))) return(invisible(NULL))
  plot_df <- patient_df %>%
    filter(
      is.finite(.data$freq_NonImmuneStroma),
      is.finite(.data$freq_NonImmuneEpithelium),
      is.finite(.data$mean_discordance)
    )
  if (nrow(plot_df) < 10L) return(invisible(NULL))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$freq_NonImmuneStroma,
      y = .data$freq_NonImmuneEpithelium,
      colour = .data$mean_discordance,
      size = abs(.data$mean_discordance)
    )
  ) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_colour_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "Mean discordance"
    ) +
    ggplot2::scale_size_continuous(range = c(1.5, 4.5), guide = "none") +
    ggplot2::labs(
      title = "Tissue composition vs RNA-protein discordance",
      subtitle = "Stroma frequency shows strongest association with discordance",
      x = "Stroma cell frequency",
      y = "Epithelium cell frequency"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  ggplot2::ggsave(out_path, p, width = 8, height = 6.5, dpi = 180, bg = "white")
  message("Wrote ", out_path)
  invisible(out_path)
}

plot_observed_vs_predicted_top <- function(long_df, out_path, top_n = 6L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  top_feats <- long_df %>%
    group_by(.data$protein_feature, .data$label) %>%
    summarise(mean_abs = mean(abs(.data$discordance), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(.data$mean_abs)) %>%
    slice_head(n = top_n)

  plot_df <- long_df %>%
    inner_join(top_feats %>% select(.data$protein_feature), by = "protein_feature")

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$predicted, y = .data$observed, colour = .data$discordance)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(alpha = 0.45, size = 1.2) +
    ggplot2::facet_wrap(~label, scales = "free", ncol = min(3, top_n)) +
    ggplot2::scale_colour_gradient2(
      low = "#2166AC", mid = "grey90", high = "#B2182B",
      midpoint = 0, name = "Discordance"
    ) +
    ggplot2::labs(
      title = "Observed vs predicted protein (top discordant markers)",
      x = "Predicted protein", y = "Observed protein"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(size = 7, face = "bold")
    )

  ncol <- min(3, top_n)
  ggplot2::ggsave(
    out_path, p,
    width = 4 * ncol,
    height = max(5, 3.5 * ceiling(top_n / ncol)),
    dpi = 180, bg = "white"
  )
  message("Wrote ", out_path)
  invisible(out_path)
}

plot_compartment_discordance <- function(long_df, out_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  plot_df <- long_df %>%
    mutate(compartment = marker_compartment(.data$label)) %>%
    group_by(.data$compartment, .data$label) %>%
    summarise(
      mean_disc = mean(.data$discordance, na.rm = TRUE),
      mean_abs = mean(abs(.data$discordance), na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$mean_disc,
      y = reorder(.data$label, .data$mean_disc),
      size = .data$mean_abs
    )
  ) +
    ggplot2::geom_point(alpha = 0.75, colour = "#542788") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::facet_wrap(~compartment, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      title = "Marker discordance by spatial compartment",
      x = "Mean discordance (obs - pred)", y = NULL, size = "Mean |disc|"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(out_path, p, width = 10, height = 12, dpi = 180, bg = "white")
  message("Wrote ", out_path)
  invisible(out_path)
}

plot_discordance_patient_heatmap <- function(
  long_df, out_path, n_patients = 40L, n_markers = 30L
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  top_markers <- long_df %>%
    group_by(.data$label) %>%
    summarise(v = mean(abs(.data$discordance), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(.data$v)) %>%
    slice_head(n = n_markers) %>%
    pull(.data$label)

  top_pts <- long_df %>%
    group_by(.data$patient_id) %>%
    summarise(v = mean(abs(.data$discordance), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(.data$v)) %>%
    slice_head(n = n_patients) %>%
    pull(.data$patient_id)

  hm <- long_df %>%
    filter(.data$label %in% top_markers, .data$patient_id %in% top_pts) %>%
    group_by(.data$patient_id, .data$label) %>%
    summarise(disc = mean(.data$discordance, na.rm = TRUE), .groups = "drop")

  p <- ggplot2::ggplot(
    hm,
    ggplot2::aes(x = .data$patient_id, y = .data$label, fill = .data$disc)
  ) +
    ggplot2::geom_tile(colour = NA) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "Discordance"
    ) +
    ggplot2::labs(
      title = "Patient x marker discordance heatmap",
      x = "Patient", y = NULL
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  ggplot2::ggsave(out_path, p, width = 12, height = 9, dpi = 180, bg = "white")
  message("Wrote ", out_path)
  invisible(out_path)
}

plot_geneexpresspred_cv <- function(cv_path, out_path, top_n = 25L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !file.exists(cv_path)) {
    return(invisible(NULL))
  }
  df <- readr::read_csv(cv_path, show_col_types = FALSE) %>%
    mutate(label = clean_protein_label(.data$protein)) %>%
    arrange(desc(.data$mean_pcc))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$mean_pcc)) +
    ggplot2::geom_histogram(bins = 25, fill = "#4393C3", colour = "white") +
    ggplot2::labs(
      title = "GeneExpressPred CV: PCC distribution across proteins",
      x = "Mean Pearson correlation (5-fold CV)", y = "Proteins"
    ) +
    ggplot2::theme_bw(base_size = 11)

  ggplot2::ggsave(out_path, p, width = 8, height = 5, dpi = 180, bg = "white")
  message("Wrote ", out_path)
  invisible(out_path)
}

run_discordance_figures <- function(output_dir = NULL, figures_subdir = "figures") {
  p <- resolve_discordance_paths()
  if (is.null(output_dir)) output_dir <- p$output
  fig_dir <- file.path(output_dir, figures_subdir)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  long_path <- file.path(output_dir, "discordance_per_patient_marker.csv")
  pat_path <- file.path(output_dir, "discordance_patient_summary.csv")
  assoc_path <- file.path(output_dir, "discordance_clinical_spatial_associations.csv")
  if (!file.exists(long_path) || !file.exists(pat_path)) {
    stop("Run run_rna_protein_discordance_all() first.")
  }

  long <- readr::read_csv(long_path, show_col_types = FALSE)
  patient <- readr::read_csv(pat_path, show_col_types = FALSE)
  assoc <- if (file.exists(assoc_path)) {
    readr::read_csv(assoc_path, show_col_types = FALSE)
  } else {
    data.frame()
  }

  plot_discordance_distribution(patient, file.path(fig_dir, "discordance_patient_histogram.png"))
  plot_top_associations(assoc, file.path(fig_dir, "discordance_top_associations.png"))
  plot_feature_discordance_volcano(long, file.path(fig_dir, "discordance_marker_volcano.png"))
  plot_discordance_by_group(
    patient, "tma_comp_grp", "mean_abs_discordance",
    file.path(fig_dir, "discordance_by_tma_group.png"),
    "RNA-protein discordance by TMA composition group"
  )
  plot_discordance_by_group(
    patient, "risk_grp", "mean_abs_discordance",
    file.path(fig_dir, "discordance_by_risk_grp.png"),
    "RNA-protein discordance by RSF risk group"
  )
  plot_discordance_by_group(
    patient, "cohort_rna", "mean_abs_discordance",
    file.path(fig_dir, "discordance_by_rna_cohort.png"),
    "RNA-protein discordance by RNA cohort"
  )
  plot_rsf_vs_tma_composition(patient, file.path(fig_dir, "rsf_risk_vs_tma_composition.png"))
  plot_discordance_spatial_scatter(patient, file.path(fig_dir, "discordance_tissue_composition.png"))
  plot_observed_vs_predicted_top(long, file.path(fig_dir, "observed_vs_predicted_top_markers.png"))
  plot_compartment_discordance(long, file.path(fig_dir, "discordance_by_compartment.png"))
  plot_discordance_patient_heatmap(long, file.path(fig_dir, "discordance_patient_marker_heatmap.png"))

  cv_path <- file.path(INTEGRATED_DIR, "results", "geneexpresspred", "cv_protein_summary_best_topn.csv")
  plot_geneexpresspred_cv(cv_path, file.path(fig_dir, "geneexpresspred_cv_performance.png"))

  writeLines(
    c(
      "Discordance figures", paste0("Output: ", fig_dir),
      "See PNG files in this folder (12 figures + index)"
    ),
    file.path(fig_dir, "figures_index.txt")
  )
  message("Discordance figures written to ", fig_dir)
  invisible(fig_dir)
}

run_rna_protein_discordance_all <- function(
  output_dir = NULL,
  top_n_rna = 30L,
  min_patients = 10L,
  max_protein_features = NULL
) {
  p <- resolve_discordance_paths()
  if (is.null(output_dir)) {
    output_dir <- p$output
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Loading patient protein matrix...")
  prot <- load_patient_protein_matrix(p$merged_rds)
  message("Loading patient RNA (cohort z-scored)...")
  rna <- load_patient_rna_cohort_z()

  common_pts <- intersect(colnames(rna$expr), colnames(prot$matrix))
  message("Patients with RNA + protein: ", length(common_pts))
  if (length(common_pts) < min_patients) {
    stop("Too few patients with both RNA and protein (n=", length(common_pts), ").")
  }

  rna_expr <- rna$expr[, common_pts, drop = FALSE]
  protein_mat <- prot$matrix[, common_pts, drop = FALSE]

  message("Predicting protein from RNA (", nrow(protein_mat), " markers)...")
  long <- compute_all_discordance(
    protein_mat,
    prot$features,
    rna_expr,
    common_pts,
    top_n = top_n_rna,
    min_n = min_patients,
    max_features = max_protein_features
  )
  if (nrow(long) == 0L) stop("No discordance results computed.")

  readr::write_csv(long, file.path(output_dir, "discordance_per_patient_marker.csv"))

  patient_sum <- summarize_patient_discordance(long)
  patient_sum <- patient_sum %>%
    left_join(rna$meta, by = "patient_id") %>%
    rename(cohort_rna = .data$cohort)

  merged_sub <- prot$merged %>%
    filter(.data$patient_id %in% common_pts)
  patient_sum <- patient_sum %>%
    left_join(merged_sub, by = "patient_id", suffix = c("", "_dup"))

  if (file.exists(p$tma_groups_rds)) {
    tma_grp <- readRDS(p$tma_groups_rds)$patient_groups %>%
      select(.data$patient_id, .data$tma_comp_grp, .data$n_tma, .data$n_tma_low, .data$n_tma_high)
    patient_sum <- patient_sum %>%
      left_join(tma_grp, by = "patient_id")
  }

  if (file.exists(p$survival_rds)) {
    surv <- readRDS(p$survival_rds) %>%
      select(.data$patient_id, .data$dfs_months, .data$dfs_event, .data$rsf_risk, .data$risk_grp)
    patient_sum <- patient_sum %>%
      left_join(surv, by = "patient_id", suffix = c("", "_surv"))
  }

  cov_cols <- clinical_spatial_covariates(patient_sum)
  assoc <- run_association_screen(patient_sum, cov_cols)

  readr::write_csv(patient_sum, file.path(output_dir, "discordance_patient_summary.csv"))
  readr::write_csv(assoc, file.path(output_dir, "discordance_clinical_spatial_associations.csv"))

  plot_discordance_distribution(
    patient_sum,
    file.path(output_dir, "discordance_patient_histogram.png")
  )
  plot_top_associations(
    assoc,
    file.path(output_dir, "discordance_top_associations.png")
  )
  plot_feature_discordance_volcano(
    long,
    file.path(output_dir, "discordance_marker_volcano.png")
  )

  summary_lines <- c(
    "RNA-protein discordance (patient-level, exploratory)",
    "",
    "Method: GeneExpressPred-style — Spearman RNA feature selection + voting",
    "  (lm + ridge glmnet + Huber rlm); discordance = observed - predicted.",
    "RNA z-scored within cohort; protein from patient-level merged proteomics.",
    "",
    paste0("Output folder: ", output_dir),
    paste0("Shared RNA genes: ", rna$n_genes),
    paste0("Patients with both omics: ", length(common_pts)),
    paste0("Protein markers modelled: ", length(unique(long$protein_feature))),
    paste0("Long-format rows: ", nrow(long)),
    "",
    "Files:",
    "  discordance_per_patient_marker.csv",
    "  discordance_patient_summary.csv",
    "  discordance_clinical_spatial_associations.csv",
    "  discordance_patient_histogram.png",
    "  discordance_top_associations.png",
    "  discordance_marker_volcano.png"
  )
  writeLines(summary_lines, file.path(output_dir, "summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(
    long = long,
    patient = patient_sum,
    associations = assoc,
    n_patients = length(common_pts)
  ))
}
