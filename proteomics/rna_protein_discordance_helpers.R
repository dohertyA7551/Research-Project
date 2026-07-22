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

  readr::write_csv(long, file.path(output_dir, "discordance_per_patient_marker.csv"))
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
