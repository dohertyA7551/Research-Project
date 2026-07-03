# Spearman correlations: spatial metacluster abundance vs PROGENy / xCell (per cohort + Fisher meta).
source("analysis_helpers.R")
source("progeny_helpers.R")
source("xcell_combine_helpers.R")

MC_PREFIX <- "abundance_metacluster_leiden_1.0_"

mc_short_name <- function(col) {
  sub(paste0("^", MC_PREFIX), "", col)
}

load_patient_metaclusters <- function(survival_df) {
  mc_cols <- grep(paste0("^", MC_PREFIX), names(survival_df), value = TRUE)
  if (length(mc_cols) == 0) {
    stop("No metacluster abundance columns found in survival_df.")
  }
  out <- survival_df %>%
    dplyr::select(
      dplyr::any_of(c("patient_id", "rsf_risk", "risk_grp")),
      dplyr::all_of(mc_cols)
    ) %>%
    dplyr::distinct(.data$patient_id, .keep_all = TRUE)
  out
}

aggregate_omics_to_patient <- function(
  scores_long,
  patient_col = "patient_id",
  feature_col,
  value_col = "score",
  sample_col = "sample_id"
) {
  if (!patient_col %in% names(scores_long)) {
    stop("Column not found: ", patient_col)
  }
  scores_long %>%
    dplyr::filter(!is.na(.data[[patient_col]])) %>%
    dplyr::group_by(.data[[patient_col]], .data[[feature_col]]) %>%
    dplyr::summarise(
      value = mean(.data[[value_col]], na.rm = TRUE),
      n_samples = dplyr::n(),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = .data[[feature_col]],
      values_from = value
    ) %>%
    tibble::column_to_rownames(patient_col) %>%
    as.data.frame(check.names = FALSE)
}

run_progeny_cohort_patient_wide <- function(expr_mat, sample_meta) {
  sample_meta <- sample_meta %>%
    dplyr::filter(
      !is.na(.data$patient_id),
      .data$sample_id %in% colnames(expr_mat)
    )
  if (nrow(sample_meta) < 3L) {
    return(NULL)
  }
  expr <- expr_mat[, sample_meta$sample_id, drop = FALSE]
  progeny <- run_progeny_mlm(expr)
  acts_long <- progeny_activities_to_long(progeny$activities, sample_meta)
  aggregate_omics_to_patient(
    acts_long,
    patient_col = "patient_id",
    feature_col = "pathway",
    value_col = "score",
    sample_col = "sample_id"
  )
}

xcell_long_to_patient_wide <- function(scores_long) {
  aggregate_omics_to_patient(
    scores_long,
    patient_col = "patient_id",
    feature_col = "cell_type",
    value_col = "enrichment",
    sample_col = "sample"
  )
}

join_mc_with_omics <- function(mc_df, omics_patient_wide) {
  common <- intersect(mc_df$patient_id, rownames(omics_patient_wide))
  if (length(common) == 0L) {
    stop("No overlapping patients between metaclusters and omics scores.")
  }
  mc_sub <- mc_df %>%
    dplyr::filter(.data$patient_id %in% common) %>%
    tibble::column_to_rownames("patient_id") %>%
    as.data.frame(check.names = FALSE)
  omics_sub <- omics_patient_wide[common, , drop = FALSE]
  cbind(mc_sub, omics_sub)
}

spearman_mc_omics <- function(
  combined_df,
  mc_cols,
  feature_cols,
  cohort,
  omics_type,
  min_n = 10L
) {
  rows <- list()
  for (mc in mc_cols) {
    mc_vals <- combined_df[[mc]]
    for (feat in feature_cols) {
      feat_vals <- combined_df[[feat]]
      ok <- is.finite(mc_vals) & is.finite(feat_vals)
      n <- sum(ok)
      if (n < min_n) {
        next
      }
      ct <- tryCatch(
        stats::cor.test(mc_vals[ok], feat_vals[ok], method = "spearman", exact = FALSE),
        error = function(e) NULL
      )
      if (is.null(ct)) {
        next
      }
      rows[[length(rows) + 1L]] <- data.frame(
        cohort = cohort,
        omics_type = omics_type,
        metacluster = mc_short_name(mc),
        mc_column = mc,
        feature = feat,
        rho = unname(ct$estimate),
        pval = ct$p.value,
        n = n,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) {
    return(tibble::tibble())
  }
  dplyr::bind_rows(rows)
}

fisher_meta_mc_omics <- function(cohort_results, min_cohorts = 2L) {
  if (length(cohort_results) < min_cohorts) {
    stop("Fisher meta requires at least ", min_cohorts, " cohort result tables.")
  }
  all <- dplyr::bind_rows(cohort_results)
  keys <- all %>%
    dplyr::distinct(.data$omics_type, .data$metacluster, .data$feature)

  meta_rows <- lapply(seq_len(nrow(keys)), function(i) {
    row <- keys[i, , drop = FALSE]
    sub <- all %>%
      dplyr::filter(
        .data$omics_type == row$omics_type,
        .data$metacluster == row$metacluster,
        .data$feature == row$feature
      )
    pvals <- sub$pval[is.finite(sub$pval)]
    rhos <- sub$rho[is.finite(sub$rho)]
    n_cohorts <- length(pvals)
    fisher_p <- if (n_cohorts >= min_cohorts) {
      stats::pchisq(-2 * sum(log(pvals)), df = 2 * n_cohorts, lower.tail = FALSE)
    } else {
      NA_real_
    }
    mean_rho <- mean(rhos, na.rm = TRUE)
    sign_agree <- if (length(rhos) >= 2L) {
      sum(rhos > 0) == length(rhos) || sum(rhos < 0) == length(rhos)
    } else {
      NA
    }
    out <- data.frame(
      omics_type = row$omics_type,
      metacluster = row$metacluster,
      feature = row$feature,
      fisher_p = fisher_p,
      fisher_padj = NA_real_,
      n_cohorts = n_cohorts,
      mean_rho = mean_rho,
      sign_consistent = sign_agree,
      stringsAsFactors = FALSE
    )
    for (cn in unique(sub$cohort)) {
      cs <- sub[sub$cohort == cn, , drop = FALSE]
      out[[paste0("rho_", cn)]] <- cs$rho[1]
      out[[paste0("pval_", cn)]] <- cs$pval[1]
      out[[paste0("n_", cn)]] <- cs$n[1]
    }
    out
  })
  meta <- dplyr::bind_rows(meta_rows)
  ok_fisher <- is.finite(meta$fisher_p)
  meta$fisher_padj[ok_fisher] <- stats::p.adjust(meta$fisher_p[ok_fisher], method = "BH")
  meta[order(meta$fisher_padj, meta$fisher_p, na.last = TRUE), , drop = FALSE]
}

export_mc_omics_heatmap <- function(
  meta_sig,
  out_path,
  title = "Metacluster vs omics (Fisher meta)"
) {
  if (nrow(meta_sig) == 0) {
    message("Skipping heatmap (no significant pairs): ", out_path)
    return(invisible(NULL))
  }
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    warning("pheatmap not installed; skipping heatmap.")
    return(invisible(NULL))
  }
  hm <- meta_sig %>%
    dplyr::mutate(label = paste0(.data$metacluster, " | ", .data$feature)) %>%
    dplyr::select(label, mean_rho) %>%
    tibble::column_to_rownames("label") %>%
    as.matrix()
  grDevices::png(out_path, width = 1200, height = max(600, 20 * nrow(hm)), res = 120)
  pheatmap::pheatmap(
    hm,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = title,
    border_color = NA
  )
  grDevices::dev.off()
  message("Saved: ", out_path)
  invisible(out_path)
}

run_metacluster_progeny_xcell_correlation <- function(
  output_dir = NULL,
  min_n = 10L,
  fisher_fdr = 0.05,
  require_sign_consistency = TRUE
) {
  p <- resolve_progeny_paths()
  if (is.null(output_dir)) {
    output_dir <- file.path(p$base, "analysis_output", "metacluster_correlation")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  survival_path <- file.path(p$proteomics, "survival_df_with_risk.rds")
  survival_df <- readRDS(survival_path)
  mc_df <- load_patient_metaclusters(survival_df)
  mc_cols <- grep(paste0("^", MC_PREFIX), names(mc_df), value = TRUE)

  cohort_loaders <- list(
    Stage2 = function() load_stage2_expression_meta(p, survival_df),
    Retrospective = function() load_retrospective_expression_meta(p, survival_df),
    Colossus = function() load_colossus_expression_meta(p, survival_df),
    Taxonomy = function() load_taxonomy_expression_meta(p, survival_path)
  )

  progeny_cohort_results <- list()
  for (cn in names(cohort_loaders)) {
    message("PROGENy + MC correlations: ", cn)
    dat <- cohort_loaders[[cn]]()
    omics_wide <- run_progeny_cohort_patient_wide(dat$expr, dat$meta)
    if (is.null(omics_wide) || ncol(omics_wide) == 0) {
      warning("Skipping ", cn, " PROGENy (no patient-level scores).")
      next
    }
    combined <- join_mc_with_omics(mc_df, omics_wide)
    feat_cols <- colnames(omics_wide)
    cor_df <- spearman_mc_omics(
      combined,
      mc_cols = intersect(mc_cols, colnames(combined)),
      feature_cols = feat_cols,
      cohort = cn,
      omics_type = "PROGENy",
      min_n = min_n
    )
    if (nrow(cor_df) > 0L) {
      readr::write_csv(
        cor_df,
        file.path(output_dir, paste0(cn, "_mc_progeny_spearman.csv"))
      )
      progeny_cohort_results[[cn]] <- cor_df
      message("  pairs tested: ", nrow(cor_df))
    }
  }

  progeny_meta <- fisher_meta_mc_omics(progeny_cohort_results)
  readr::write_csv(
    progeny_meta,
    file.path(output_dir, "MC_progeny_fisher_meta.csv")
  )

  progeny_sig <- progeny_meta %>%
    dplyr::filter(.data$fisher_padj <= fisher_fdr)
  if (require_sign_consistency) {
    progeny_sig <- progeny_sig %>%
      dplyr::filter(.data$sign_consistent %in% TRUE | .data$n_cohorts < 2L)
  }
  readr::write_csv(
    progeny_sig,
    file.path(output_dir, "MC_progeny_fisher_meta_sig.csv")
  )
  export_mc_omics_heatmap(
    utils::head(progeny_sig, 40),
    file.path(output_dir, "MC_progeny_fisher_meta_sig_heatmap.png"),
    title = "Significant MC vs PROGENy (Fisher meta)"
  )

  xcell_loaders <- list(
    Stage2 = function() load_xcell_scores_stage2(p, survival_df),
    Retrospective = function() load_xcell_scores_retrospective(p, survival_df),
    Colossus = function() load_xcell_scores_colossus(p, survival_df),
    Taxonomy = function() load_xcell_scores_taxonomy(p, survival_df)
  )

  xcell_cohort_results <- list()
  for (cn in names(xcell_loaders)) {
    message("xCell + MC correlations: ", cn)
    scores_long <- xcell_loaders[[cn]]()
    scores_long <- scores_long %>%
      dplyr::filter(!is.na(.data$patient_id))
    if (nrow(scores_long) == 0L) {
      warning("Skipping ", cn, " xCell (no mapped patients).")
      next
    }
    omics_wide <- xcell_long_to_patient_wide(scores_long)
    combined <- join_mc_with_omics(mc_df, omics_wide)
    feat_cols <- intersect(colnames(omics_wide), colnames(combined))
    cor_df <- spearman_mc_omics(
      combined,
      mc_cols = intersect(mc_cols, colnames(combined)),
      feature_cols = feat_cols,
      cohort = cn,
      omics_type = "xCell",
      min_n = min_n
    )
    if (nrow(cor_df) > 0L) {
      readr::write_csv(
        cor_df,
        file.path(output_dir, paste0(cn, "_mc_xcell_spearman.csv"))
      )
      xcell_cohort_results[[cn]] <- cor_df
      message("  pairs tested: ", nrow(cor_df))
    }
  }

  xcell_meta <- fisher_meta_mc_omics(xcell_cohort_results)
  readr::write_csv(
    xcell_meta,
    file.path(output_dir, "MC_xcell_fisher_meta.csv")
  )

  xcell_sig <- xcell_meta %>%
    dplyr::filter(.data$fisher_padj <= fisher_fdr)
  if (require_sign_consistency) {
    xcell_sig <- xcell_sig %>%
      dplyr::filter(.data$sign_consistent %in% TRUE | .data$n_cohorts < 2L)
  }
  readr::write_csv(
    xcell_sig,
    file.path(output_dir, "MC_xcell_fisher_meta_sig.csv")
  )
  export_mc_omics_heatmap(
    utils::head(xcell_sig, 40),
    file.path(output_dir, "MC_xcell_fisher_meta_sig_heatmap.png"),
    title = "Significant MC vs xCell (Fisher meta)"
  )

  legacy_xcell <- xcell_cohort_results[names(xcell_cohort_results) %in% LEGACY_XCELL_COHORTS]
  if (length(legacy_xcell) >= 2L) {
    legacy_meta <- fisher_meta_mc_omics(legacy_xcell)
    readr::write_csv(
      legacy_meta,
      file.path(output_dir, "MC_xcell_legacy3_fisher_meta.csv")
    )
  }

  summary_lines <- c(
    "Metacluster abundance vs PROGENy / xCell (Spearman, patient-level)",
    paste0("Metaclusters: ", length(mc_cols)),
    paste0("Min patients per test: ", min_n),
    paste0("Fisher FDR threshold: ", fisher_fdr),
    "",
    paste0("PROGENy cohorts: ", paste(names(progeny_cohort_results), collapse = ", ")),
    paste0(
      "PROGENy Fisher sig (padj <= ", fisher_fdr, "): ",
      nrow(progeny_sig), " / ", nrow(progeny_meta)
    ),
    "",
    paste0("xCell cohorts: ", paste(names(xcell_cohort_results), collapse = ", ")),
    paste0(
      "xCell Fisher sig (padj <= ", fisher_fdr, "): ",
      nrow(xcell_sig), " / ", nrow(xcell_meta)
    ),
    ""
  )

  if (nrow(progeny_sig) > 0L) {
    top_p <- progeny_sig %>%
      dplyr::arrange(.data$fisher_padj) %>%
      utils::head(8)
    summary_lines <- c(
      summary_lines,
      "Top PROGENy associations:",
      paste0(
        "  ",
        top_p$metacluster,
        " ~ ",
        top_p$feature,
        ": meta rho=",
        signif(top_p$mean_rho, 3),
        ", Fisher padj=",
        signif(top_p$fisher_padj, 3)
      )
    )
  } else {
    summary_lines <- c(summary_lines, "Top PROGENy associations: none at FDR threshold")
  }

  if (nrow(xcell_sig) > 0L) {
    top_x <- xcell_sig %>%
      dplyr::arrange(.data$fisher_padj) %>%
      utils::head(8)
    summary_lines <- c(
      summary_lines,
      "",
      "Top xCell associations:",
      paste0(
        "  ",
        top_x$metacluster,
        " ~ ",
        top_x$feature,
        ": meta rho=",
        signif(top_x$mean_rho, 3),
        ", Fisher padj=",
        signif(top_x$fisher_padj, 3)
      )
    )
  } else {
    summary_lines <- c(summary_lines, "", "Top xCell associations: none at FDR threshold")
  }

  writeLines(summary_lines, file.path(output_dir, "summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(
    progeny_cohort = progeny_cohort_results,
    progeny_meta = progeny_meta,
    progeny_sig = progeny_sig,
    xcell_cohort = xcell_cohort_results,
    xcell_meta = xcell_meta,
    xcell_sig = xcell_sig
  ))
}
