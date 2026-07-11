# xCell enrichment: Anderson-Darling k-sample tests across cohorts and optional pooling.
# Requires kSamples (ad.test). Install once if missing:
#   install.packages(c("SuppDists", "kSamples"))

resolve_transcriptomics_paths <- function() {
  if (file.exists("Stage2_transcriptomics_xcell.csv")) {
    list(
      base = ".",
      proteomics = "../proteomics",
      spatial = "../spatial",
      integrated = "../Integrated "
    )
  } else if (file.exists("Transcriptomics/Stage2_transcriptomics_xcell.csv")) {
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

load_xcell_scores_stage2 <- function(p, survival_df) {
  stage2_res <- readr::read_csv(
    file.path(p$base, "Stage2_transcriptomics_xcell.csv"),
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
  r_codes <- clinical$r_code[match(colnames(stage2_res)[-1], clinical$patient_id)]
  patient_id_g <- mapping$patient_id_g[match(r_codes, mapping$r_code)]
  id_map <- data.frame(
    sample = colnames(stage2_res)[-1],
    patient_id = patient_id_g,
    stringsAsFactors = FALSE
  )
  stage2_res %>%
    dplyr::rename(cell_type = 1) %>%
    tidyr::pivot_longer(-cell_type, names_to = "sample", values_to = "enrichment") %>%
    dplyr::left_join(id_map, by = "sample") %>%
    dplyr::left_join(
      survival_df %>% dplyr::select(patient_id, risk_grp, rsf_risk),
      by = "patient_id"
    ) %>%
    dplyr::filter(risk_grp %in% c("Low", "High")) %>%
    dplyr::mutate(cohort = "Stage2")
}

load_xcell_scores_retrospective <- function(p, survival_df) {
  xcell_r <- readr::read_csv(
    file.path(p$base, "reterospective_transcriptomics_xcell_results.csv"),
    show_col_types = FALSE
  )
  xcell_r %>%
    dplyr::rename(cell_type = 1) %>%
    tidyr::pivot_longer(-cell_type, names_to = "sample", values_to = "enrichment") %>%
    dplyr::mutate(patient_id = sample) %>%
    dplyr::left_join(
      survival_df %>% dplyr::select(patient_id, risk_grp, rsf_risk),
      by = "patient_id"
    ) %>%
    dplyr::filter(risk_grp %in% c("Low", "High")) %>%
    dplyr::mutate(cohort = "Retrospective")
}

load_xcell_scores_colossus <- function(p, survival_df) {
  xcell_c <- readr::read_csv(
    file.path(p$base, "collosus_transcriptomics_xcell_results.csv"),
    show_col_types = FALSE
  )
  master_df <- readr::read_csv(
    file.path(p$integrated, "All_patinet_IDS_concatonated_masterdoc.csv"),
    show_col_types = FALSE
  )
  strip_fp <- function(x) stringr::str_replace(x, "-FP.*$", "")
  master_col_map <- dplyr::bind_rows(
    master_df %>%
      dplyr::filter(!is.na(patient_id_g), patient_id_g != "") %>%
      dplyr::transmute(patient_id = patient_id_g, col_id = strip_fp(COLOSSUS_ID)),
    master_df %>%
      dplyr::filter(!is.na(patient_id_g), patient_id_g != "") %>%
      dplyr::transmute(patient_id = patient_id_g, col_id = strip_fp(Old_COLOSSUS_ID)),
    master_df %>%
      dplyr::filter(!is.na(patient_id_g), patient_id_g != "") %>%
      dplyr::transmute(patient_id = patient_id_g, col_id = strip_fp(Sample_ID)),
    master_df %>%
      dplyr::filter(!is.na(patient_id_g), patient_id_g != "") %>%
      dplyr::transmute(patient_id = patient_id_g, col_id = strip_fp(Alternative_ID))
  ) %>%
    dplyr::filter(!is.na(col_id), col_id != "", col_id != "NA") %>%
    dplyr::distinct(col_id, .keep_all = TRUE)
  samples_c <- colnames(xcell_c)[-1]
  id_map_c <- data.frame(sample = samples_c, stringsAsFactors = FALSE) %>%
    dplyr::left_join(master_col_map, by = c("sample" = "col_id"))
  xcell_c %>%
    dplyr::rename(cell_type = 1) %>%
    tidyr::pivot_longer(-cell_type, names_to = "sample", values_to = "enrichment") %>%
    dplyr::left_join(id_map_c, by = "sample") %>%
    dplyr::left_join(
      survival_df %>% dplyr::select(patient_id, risk_grp, rsf_risk),
      by = "patient_id"
    ) %>%
    dplyr::filter(risk_grp %in% c("Low", "High")) %>%
    dplyr::mutate(cohort = "Colossus")
}

load_xcell_scores_taxonomy <- function(p, survival_df) {
  taxonomy_xcell <- file.path(
    p$base, "Taxonomy_calls",
    "QUB-Taxonomy-classifieR (RNA-Seq)",
    "classifieR results", "xCell2022-03-14.csv"
  )
  taxonomy_map <- file.path(
    p$base, "Taxonomy_calls",
    "Manuela_and_Belfast_RNA_classifications_CMS_CRIS.txt"
  )
  xcell_wide <- readr::read_csv(taxonomy_xcell, show_col_types = FALSE)
  map_df <- readr::read_tsv(taxonomy_map, show_col_types = FALSE)
  names(xcell_wide)[1] <- "sample"
  xcell_wide <- xcell_wide %>%
    dplyr::mutate(sample = stringr::str_replace(sample, stringr::fixed("."), "-"))
  id_map <- map_df %>%
    dplyr::transmute(sample_id = Patient, patient_id = Code) %>%
    dplyr::distinct()
  xcell_wide %>%
    dplyr::rename(sample_id = sample) %>%
    tidyr::pivot_longer(
      -sample_id,
      names_to = "cell_type",
      values_to = "enrichment"
    ) %>%
    dplyr::left_join(id_map, by = "sample_id") %>%
    dplyr::left_join(
      survival_df %>% dplyr::select(patient_id, risk_grp, rsf_risk),
      by = "patient_id"
    ) %>%
    dplyr::filter(risk_grp %in% c("Low", "High")) %>%
    dplyr::transmute(
      cell_type, sample = sample_id, enrichment,
      patient_id, risk_grp, rsf_risk, cohort = "Taxonomy"
    )
}

# Legacy cohorts share xCell v1-style labels; Taxonomy uses xCell2022 separately.
LEGACY_XCELL_COHORTS <- c("Stage2", "Retrospective", "Colossus")

load_all_xcell_scores_long <- function(...) {
  if (!requireNamespace("kSamples", quietly = TRUE)) {
    stop("Package 'kSamples' is required for Anderson-Darling tests.")
  }
  load_xcell_scores_long(...)
}

load_xcell_scores_long <- function(
  harmonize_cell_types = TRUE,
  cohorts = NULL
) {
  if (!is.null(cohorts) && "Taxonomy" %in% cohorts &&
      !requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required when loading Taxonomy xCell scores.")
  }
  p <- resolve_transcriptomics_paths()
  survival_df <- readRDS(file.path(p$proteomics, "survival_df_with_risk.rds"))
  all_cohorts <- c(LEGACY_XCELL_COHORTS, "Taxonomy")
  if (is.null(cohorts)) {
    cohorts <- all_cohorts
  }
  cohorts <- intersect(cohorts, all_cohorts)
  if (length(cohorts) == 0) {
    stop("No valid cohorts requested.")
  }

  loaders <- list(
    Stage2 = function() load_xcell_scores_stage2(p, survival_df),
    Retrospective = function() load_xcell_scores_retrospective(p, survival_df),
    Colossus = function() load_xcell_scores_colossus(p, survival_df),
    Taxonomy = function() load_xcell_scores_taxonomy(p, survival_df)
  )
  out <- dplyr::bind_rows(lapply(cohorts, function(cn) loaders[[cn]]())) %>%
    dplyr::mutate(
      cell_type_raw = .data$cell_type,
      cohort = factor(
        .data$cohort,
        levels = intersect(all_cohorts, cohorts)
      )
    )
  if (harmonize_cell_types) {
    out <- harmonize_xcell_cell_types(out, p)
  }
  out
}

load_xcell_cell_type_mapping <- function(p = resolve_transcriptomics_paths()) {
  map_path <- file.path(p$base, "xcell_cell_type_unified.csv")
  if (!file.exists(map_path)) {
    stop("Missing cell type mapping: ", map_path)
  }
  readr::read_csv(map_path, show_col_types = FALSE)
}

harmonize_xcell_cell_types <- function(scores_long, p = resolve_transcriptomics_paths()) {
  mapping <- load_xcell_cell_type_mapping(p)
  scores_long %>%
    dplyr::left_join(
      mapping %>% dplyr::select(source_name, unified_cell_type),
      by = c("cell_type" = "source_name")
    ) %>%
    dplyr::mutate(
      cell_type = dplyr::coalesce(.data$unified_cell_type, .data$cell_type)
    ) %>%
    dplyr::select(-unified_cell_type)
}

xcell_cell_types_in_all_cohorts <- function(
  scores_long,
  cohorts = c("Stage2", "Retrospective", "Colossus", "Taxonomy")
) {
  scores_long %>%
    dplyr::filter(.data$cohort %in% cohorts) %>%
    dplyr::count(.data$cell_type, .data$cohort) %>%
    dplyr::count(.data$cell_type, name = "n_cohorts") %>%
    dplyr::filter(.data$n_cohorts == length(cohorts)) %>%
    dplyr::pull(.data$cell_type) %>%
    sort()
}

ad_test_k_samples <- function(values_by_cohort) {
  values_by_cohort <- lapply(values_by_cohort, function(v) v[is.finite(v)])
  values_by_cohort <- values_by_cohort[vapply(values_by_cohort, length, 0) > 0]
  if (length(values_by_cohort) < 2) {
    return(list(ad = NA_real_, t_ad = NA_real_, p = NA_real_, k = length(values_by_cohort)))
  }
  k <- length(values_by_cohort)
  ns <- vapply(values_by_cohort, length, integer(1))
  if (sum(ns) < 4 || min(ns) < 2) {
    return(list(ad = NA_real_, t_ad = NA_real_, p = NA_real_, k = k))
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

run_xcell_ad_by_cell_type <- function(
  scores_long,
  cohorts = c("Stage2", "Retrospective", "Colossus", "Taxonomy"),
  ad_alpha = 0.05,
  require_all_cohorts = TRUE
) {
  cohorts <- intersect(cohorts, unique(as.character(scores_long$cohort)))
  n_per_cohort <- scores_long %>%
    dplyr::filter(.data$cohort %in% cohorts) %>%
    dplyr::distinct(.data$cohort, .data$sample) %>%
    dplyr::count(.data$cohort, name = "n_samples")

  cell_types <- if (require_all_cohorts) {
    xcell_cell_types_in_all_cohorts(scores_long, cohorts)
  } else {
    sort(unique(scores_long$cell_type))
  }

  rows <- lapply(cell_types, function(ct) {
    sub <- scores_long %>% dplyr::filter(.data$cell_type == ct)
    vecs <- stats::setNames(
      lapply(cohorts, function(cn) sub$enrichment[sub$cohort == cn]),
      cohorts
    )
    vecs <- vecs[vapply(vecs, function(v) sum(is.finite(v)) > 0, logical(1))]
    tst <- ad_test_k_samples(vecs)
    ns <- vapply(
      cohorts,
      function(cn) sum(sub$cohort == cn & is.finite(sub$enrichment)),
      integer(1)
    )
    out <- data.frame(
      cell_type = ct,
      ad_stat = tst$ad,
      t_ad = tst$t_ad,
      ad_p = tst$p,
      n_cohorts = tst$k,
      pass_merge = is.finite(tst$p) && tst$p >= ad_alpha,
      stringsAsFactors = FALSE
    )
    for (cn in cohorts) {
      out[[paste0("n_", cn)]] <- ns[cn]
    }
    out
  })

  results <- dplyr::bind_rows(rows)
  list(
    results = results,
    n_samples = n_per_cohort,
    n_cell_types = nrow(results),
    n_pass = sum(results$pass_merge, na.rm = TRUE),
    passing_cell_types = results$cell_type[results$pass_merge],
    cohorts = cohorts
  )
}

export_xcell_ad_cell_type_plot <- function(ad_out, output_dir, ad_alpha = 0.05) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  plot_df <- ad_out$results %>%
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
      y = reorder(.data$cell_type, .data$ad_p),
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
      title = "Anderson-Darling k-sample: xCell2 scores by cell type",
      subtitle = paste0(
        "Per cell type: do score distributions match across cohorts? ",
        "Pass if AD p \u2265 ", ad_alpha
      ),
      x = expression(-log[10] * "(AD p-value)"),
      y = NULL,
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  h <- max(5, min(18, 0.22 * nrow(plot_df) + 2))
  out <- file.path(output_dir, "xcell_ad_by_cell_type_summary.png")
  ggplot2::ggsave(out, p, width = 8, height = h, dpi = 150, bg = "white")
  message("Wrote ", out)
  invisible(out)
}

prepare_pooled_xcell_scores <- function(scores_long) {
  scores_long %>%
    dplyr::mutate(
      sample_uid = paste(.data$cohort, .data$sample, sep = "::")
    )
}

build_xcell_combined_scores <- function(
  scores_long,
  cell_types = NULL,
  cohorts = NULL,
  require_all_cohorts = TRUE
) {
  if (!is.null(cohorts)) {
    scores_long <- scores_long %>%
      dplyr::filter(.data$cohort %in% cohorts)
  }
  cohorts_present <- sort(unique(as.character(scores_long$cohort)))
  if (length(cohorts_present) == 0) {
    stop("No samples left after cohort filtering.")
  }
  if (is.null(cell_types)) {
    cell_types <- if (require_all_cohorts) {
      xcell_cell_types_in_all_cohorts(scores_long, cohorts = cohorts_present)
    } else {
      sort(unique(scores_long$cell_type))
    }
  }
  if (length(cell_types) == 0) {
    stop("No cell types available for the combined xCell dataset.")
  }
  scores_long %>%
    dplyr::filter(.data$cell_type %in% cell_types) %>%
    prepare_pooled_xcell_scores()
}

normalize_values_by_cohort_quantiles <- function(x, batch) {
  x <- as.numeric(x)
  batch <- as.character(batch)
  out <- x
  ok <- is.finite(x)
  if (sum(ok) < 4) {
    return(out)
  }

  idx_ok <- which(ok)
  batches <- split(idx_ok, batch[idx_ok])
  n_per <- vapply(batches, length, integer(1))
  if (length(batches) < 2 || any(n_per < 2)) {
    return(out)
  }

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

quantile_normalize_xcell_cross_cohort <- function(scores_long) {
  if (!"sample_uid" %in% names(scores_long)) {
    scores_long <- prepare_pooled_xcell_scores(scores_long)
  }
  scores_long %>%
    dplyr::group_by(.data$cell_type) %>%
    dplyr::mutate(
      enrichment_raw = .data$enrichment,
      enrichment = normalize_values_by_cohort_quantiles(
        .data$enrichment,
        .data$cohort
      )
    ) %>%
    dplyr::ungroup()
}

export_combined_xcell_differential <- function(
  scores_qn,
  output_dir,
  file_prefix = "Combined_QN",
  title_suffix = "quantile-normalized pooled cohorts"
) {
  diff_combined <- wilcox_enrichment_diff(scores_qn)
  meta <- scores_qn %>%
    dplyr::distinct(.data$sample_uid, .data$risk_grp, .data$cohort)
  n_samples <- nrow(meta)
  n_low <- sum(meta$risk_grp == "Low")
  n_high <- sum(meta$risk_grp == "High")

  export_xcell_volcano(
    diff_combined,
    paste0("Cell Type Enrichment: High vs Low Risk (", title_suffix, ")"),
    paste0(
      "Stage 2 + Retrospective + Colossus + Taxonomy | xCell | n = ", n_samples,
      " (Low ", n_low, ", High ", n_high, ")"
    ),
    file.path(output_dir, paste0(file_prefix, "_risk_enrichment_volcano.png"))
  )
  export_top_log2fc_table(
    diff_combined,
    file.path(output_dir, paste0(file_prefix, "_top5_log2fc_up_down.csv"))
  )
  readr::write_csv(
    diff_combined,
    file.path(output_dir, paste0(file_prefix, "_wilcox_all_cell_types.csv"))
  )
  diff_combined
}

run_xcell_quantile_combined_analysis <- function(
  output_dir = NULL,
  cell_types = NULL
) {
  if (!requireNamespace("readr", quietly = TRUE)) library(readr)
  p <- resolve_transcriptomics_paths()
  if (is.null(output_dir)) {
    output_dir <- file.path(p$base, "analysis_output", "combined")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  scores_long <- load_all_xcell_scores_long()
  scores <- build_xcell_combined_scores(scores_long, cell_types = cell_types)
  cell_types_used <- sort(unique(scores$cell_type))
  scores_qn <- quantile_normalize_xcell_cross_cohort(scores)

  long_path <- file.path(output_dir, "Combined_QN_xcell_scores_long.csv")
  readr::write_csv(scores_qn, long_path)
  message("Wrote ", long_path)

  wide_path <- file.path(output_dir, "Combined_QN_xcell_scores_wide.csv")
  scores_qn %>%
    dplyr::select(
      .data$sample_uid, .data$cohort, .data$sample, .data$patient_id,
      .data$risk_grp, .data$rsf_risk, .data$cell_type, .data$enrichment
    ) %>%
    tidyr::pivot_wider(
      names_from = .data$cell_type,
      values_from = .data$enrichment
    ) %>%
    readr::write_csv(wide_path)
  message("Wrote ", wide_path)

  readr::write_csv(
    data.frame(cell_type = cell_types_used, stringsAsFactors = FALSE),
    file.path(output_dir, "Combined_QN_cell_types.csv")
  )

  sample_meta <- scores_qn %>%
    dplyr::distinct(
      .data$sample_uid, .data$cohort, .data$sample,
      .data$patient_id, .data$risk_grp, .data$rsf_risk
    )
  readr::write_csv(
    sample_meta,
    file.path(output_dir, "Combined_QN_sample_metadata.csv")
  )

  diff_combined <- export_combined_xcell_differential(
    scores_qn,
    output_dir,
    file_prefix = "Combined_QN"
  )

  summary_lines <- c(
    "Combined xCell analysis with cross-cohort quantile normalization",
    "",
    paste0("Cell types (harmonized, present in all 4 cohorts): ", length(cell_types_used)),
    paste0("Samples: ", nrow(sample_meta)),
    paste0(
      "Per cohort: ",
      paste(
        sample_meta %>% dplyr::count(.data$cohort, name = "n") %>%
          dplyr::transmute(label = paste0(.data$cohort, "=", .data$n)) %>%
          dplyr::pull(.data$label),
        collapse = "; "
      )
    ),
    paste0(
      "Risk groups: Low=", sum(sample_meta$risk_grp == "Low"),
      ", High=", sum(sample_meta$risk_grp == "High")
    ),
    "",
    "Normalization:",
    "  For each unified cell type, cohort-specific score distributions were aligned",
    "  to a common reference quantile function (mean quantiles across Stage 2,",
    "  Retrospective, Colossus, and Taxonomy). Original harmonized scores are kept",
    "  in enrichment_raw; enrichment holds quantile-normalized values used for testing.",
    "",
    "Differential testing:",
    "  Wilcoxon rank-sum High vs Low per cell type on quantile-normalized scores.",
    "",
    paste0(
      "Significant at p <= 0.05 and |log2FC| > 0.5: ",
      sum(diff_combined$pval <= 0.05 & abs(diff_combined$log2FC) > 0.5, na.rm = TRUE),
      " / ", nrow(diff_combined)
    )
  )
  summary_path <- file.path(output_dir, "Combined_QN_summary.txt")
  writeLines(summary_lines, summary_path)
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(
    scores_harmonized = scores,
    scores_qn = scores_qn,
    diff_combined = diff_combined,
    cell_types = cell_types_used,
    sample_meta = sample_meta
  ))
}

run_xcell_quantile_legacy_combined_analysis <- function(output_dir = NULL) {
  if (!requireNamespace("readr", quietly = TRUE)) library(readr)
  p <- resolve_transcriptomics_paths()
  if (is.null(output_dir)) {
    output_dir <- file.path(p$base, "analysis_output", "legacy_combined")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  scores_long <- load_xcell_scores_long(
    harmonize_cell_types = FALSE,
    cohorts = LEGACY_XCELL_COHORTS
  )
  scores <- build_xcell_combined_scores(
    scores_long,
    cohorts = LEGACY_XCELL_COHORTS
  )
  cell_types_used <- sort(unique(scores$cell_type))
  scores_qn <- quantile_normalize_xcell_cross_cohort(scores)

  long_path <- file.path(output_dir, "Legacy_QN_xcell_scores_long.csv")
  readr::write_csv(scores_qn, long_path)
  message("Wrote ", long_path)

  wide_path <- file.path(output_dir, "Legacy_QN_xcell_scores_wide.csv")
  scores_qn %>%
    dplyr::select(
      sample_uid, cohort, sample, patient_id,
      risk_grp, rsf_risk, cell_type, enrichment
    ) %>%
    tidyr::pivot_wider(
      names_from = cell_type,
      values_from = enrichment
    ) %>%
    readr::write_csv(wide_path)
  message("Wrote ", wide_path)

  readr::write_csv(
    data.frame(cell_type = cell_types_used, stringsAsFactors = FALSE),
    file.path(output_dir, "Legacy_QN_cell_types.csv")
  )

  sample_meta <- scores_qn %>%
    dplyr::distinct(
      sample_uid, cohort, sample, patient_id, risk_grp, rsf_risk
    )
  readr::write_csv(
    sample_meta,
    file.path(output_dir, "Legacy_QN_sample_metadata.csv")
  )

  diff_combined <- export_combined_xcell_differential(
    scores_qn,
    output_dir,
    file_prefix = "Legacy_QN",
    title_suffix = "Stage 2 + Retrospective + Colossus (quantile-normalized)"
  )

  summary_lines <- c(
    "Legacy combined xCell: Stage 2 + Retrospective + Colossus only",
    "Taxonomy is analysed separately (different xCell version / cell-type labels).",
    "",
    paste0("Cell types (shared across 3 legacy cohorts): ", length(cell_types_used)),
    paste0("Samples: ", nrow(sample_meta)),
    paste0(
      "Per cohort: ",
      paste(
        sample_meta %>%
          dplyr::count(cohort, name = "n") %>%
          dplyr::transmute(label = paste0(cohort, "=", n)) %>%
          dplyr::pull(label),
        collapse = "; "
      )
    ),
    paste0(
      "Risk groups: Low=", sum(sample_meta$risk_grp == "Low"),
      ", High=", sum(sample_meta$risk_grp == "High")
    ),
    "",
    "Normalization:",
    "  Per cell type, quantile alignment across Stage 2, Retrospective, and Colossus",
    "  only (enrichment_raw retains pre-normalization scores).",
    "",
    "Differential testing:",
    "  Wilcoxon rank-sum High vs Low on quantile-normalized scores.",
    "",
    paste0(
      "Significant at p <= 0.05 and |log2FC| > 0.5: ",
      sum(diff_combined$pval <= 0.05 & abs(diff_combined$log2FC) > 0.5, na.rm = TRUE),
      " / ", nrow(diff_combined)
    ),
    "",
    "Taxonomy outputs: analysis_output/taxonomy/"
  )
  summary_path <- file.path(output_dir, "Legacy_QN_summary.txt")
  writeLines(summary_lines, summary_path)
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(
    scores = scores,
    scores_qn = scores_qn,
    diff_combined = diff_combined,
    cell_types = cell_types_used,
    sample_meta = sample_meta
  ))
}

export_xcell_volcano <- function(
  differential_results,
  title,
  subtitle,
  out_path
) {
  if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
    stop("Package 'EnhancedVolcano' is required for volcano plots.")
  }
  diff <- differential_results[
    is.finite(differential_results$log2FC) & !is.na(differential_results$pval),
  ]
  axes <- volcano_axis_limits(diff)
  p <- suppressWarnings(
    EnhancedVolcano::EnhancedVolcano(
      diff,
      lab = diff$cell_type,
      x = "log2FC",
      y = "pval",
      pCutoff = 0.05001,
      FCcutoff = 0.5,
      title = title,
      subtitle = subtitle,
      xlab = "Log2 Fold Change (High / Low)",
      ylab = "-Log10(p-value)",
      pointSize = 2.0,
      labSize = 3.5,
      drawConnectors = TRUE,
      widthConnectors = 0.4,
      maxoverlapsConnectors = Inf,
      xlim = axes$xlim,
      ylim = axes$ylim
    ) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(20, 60, 20, 60),
        legend.position = "bottom"
      )
  )
  ggplot2::ggsave(out_path, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  message("Saved: ", out_path)
  invisible(out_path)
}

export_xcell_ad_summary_plot <- function(ad_out, output_dir, ad_alpha = 0.05) {
  export_xcell_ad_cell_type_plot(ad_out, output_dir, ad_alpha)
}

xcell_cohort_colors <- function() {
  c(
    Stage2 = "#377eb8",
    Retrospective = "#ff7f00",
    Colossus = "#4daf4a",
    Taxonomy = "#e41a1c"
  )
}

default_xcell_density_cell_types <- function() {
  c(
    "B cells",
    "Exhausted CD8 T cells",
    "Macrophages",
    "CD8 T cells",
    "Tregs",
    "NK cells",
    "Fibroblasts",
    "Endothelial cells",
    "Monocytes",
    "Plasma cells",
    "Mast cells",
    "DC"
  )
}

pick_xcell_density_cell_types <- function(
  scores_long,
  cell_types = NULL,
  n_cell_types = 12L,
  require_all_cohorts = TRUE
) {
  available_all <- if (require_all_cohorts) {
    xcell_cell_types_in_all_cohorts(scores_long)
  } else {
    sort(unique(scores_long$cell_type))
  }

  if (!is.null(cell_types)) {
    available <- intersect(cell_types, available_all)
    missing <- setdiff(cell_types, available)
    if (length(missing) > 0) {
      message(
        "Cell types skipped (missing in one or more cohorts): ",
        paste(missing, collapse = ", ")
      )
    }
    if (length(available) == 0) {
      stop("None of the requested cell types are present in all four cohorts.")
    }
    return(available)
  }

  defaults <- intersect(default_xcell_density_cell_types(), available_all)
  if (length(defaults) >= n_cell_types) {
    return(defaults[seq_len(n_cell_types)])
  }

  extra <- setdiff(available_all, defaults)
  c(defaults, extra)[seq_len(min(n_cell_types, length(defaults) + length(extra)))]
}

export_xcell_cohort_density_plots <- function(
  scores_long,
  output_dir,
  cell_types = NULL,
  n_cell_types = 12L,
  z_score_within_cell_type = FALSE,
  file_prefix = "xcell_cohort_density"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for density plots.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cell_types <- pick_xcell_density_cell_types(
    scores_long,
    cell_types = cell_types,
    n_cell_types = n_cell_types
  )
  plot_df <- scores_long %>%
    dplyr::filter(.data$cell_type %in% cell_types) %>%
    dplyr::mutate(
      cohort = factor(
        .data$cohort,
        levels = c("Stage2", "Retrospective", "Colossus", "Taxonomy")
      ),
      cell_type = factor(.data$cell_type, levels = cell_types)
    )

  cohort_n <- plot_df %>%
    dplyr::distinct(.data$cohort, .data$sample) %>%
    dplyr::count(.data$cohort, name = "n_samples")
  cohort_labels <- stats::setNames(
    paste0(cohort_n$cohort, " (n=", cohort_n$n_samples, ")"),
    as.character(cohort_n$cohort)
  )

  if (z_score_within_cell_type) {
    plot_df <- plot_df %>%
      dplyr::group_by(.data$cell_type) %>%
      dplyr::mutate(
        enrichment = as.numeric(scale(.data$enrichment))
      ) %>%
      dplyr::ungroup()
    x_lab <- "Z-scored xCell enrichment (within cell type)"
    suffix <- "zscore"
    scale_note <- "Scores z-scored within each cell type to compare shape across cohorts."
  } else {
    x_lab <- "xCell enrichment score"
    suffix <- "raw"
    scale_note <- "Raw xCell scores; facet axes are free. Legacy cohorts use xCell v1-style names mapped to Taxonomy xCell2022 labels."
  }

  cohort_cols <- xcell_cohort_colors()
  n_ct <- length(cell_types)
  ncol <- if (n_ct <= 4) 2L else if (n_ct <= 9) 3L else 4L

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$enrichment, colour = .data$cohort, fill = .data$cohort)
  ) +
    ggplot2::geom_density(alpha = 0.18, linewidth = 0.9) +
    ggplot2::geom_rug(
      ggplot2::aes(colour = .data$cohort),
      alpha = 0.35,
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~cell_type, scales = "free", ncol = ncol) +
    ggplot2::scale_colour_manual(values = cohort_cols, labels = cohort_labels, drop = FALSE) +
    ggplot2::scale_fill_manual(values = cohort_cols, labels = cohort_labels, drop = FALSE) +
    ggplot2::labs(
      title = "xCell enrichment distributions by cohort",
      subtitle = paste0(
        scale_note, " Showing ", n_ct,
        " cell types present in all four cohorts."
      ),
      x = x_lab,
      y = "Density",
      colour = "Cohort",
      fill = "Cohort"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      legend.position = "bottom"
    )

  panel_path <- file.path(
    output_dir,
    paste0(file_prefix, "_", suffix, "_panel.png")
  )
  panel_h <- max(8, 2.2 * ceiling(n_ct / ncol))
  ggplot2::ggsave(panel_path, p, width = 4.2 * ncol, height = panel_h, dpi = 180)
  message("Wrote ", panel_path)

  for (ct in cell_types) {
    ct_df <- plot_df %>% dplyr::filter(.data$cell_type == ct)
    ct_slug <- gsub("[^A-Za-z0-9]+", "_", ct)
    ct_slug <- gsub("_+", "_", ct_slug)
    ct_slug <- gsub("^_|_$", "", ct_slug)
    p_one <- ggplot2::ggplot(
      ct_df,
      ggplot2::aes(x = .data$enrichment, colour = .data$cohort, fill = .data$cohort)
    ) +
      ggplot2::geom_density(alpha = 0.22, linewidth = 1) +
      ggplot2::geom_rug(
        ggplot2::aes(colour = .data$cohort),
        alpha = 0.4,
        linewidth = 0.3,
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_manual(values = cohort_cols, labels = cohort_labels, drop = FALSE) +
      ggplot2::scale_fill_manual(values = cohort_cols, labels = cohort_labels, drop = FALSE) +
      ggplot2::labs(
        title = ct,
        subtitle = paste0(scale_note, " All four cohorts shown."),
        x = x_lab,
        y = "Density",
        colour = "Cohort",
        fill = "Cohort"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom"
      )
    one_path <- file.path(
      output_dir,
      paste0(file_prefix, "_", suffix, "_", ct_slug, ".png")
    )
    ggplot2::ggsave(one_path, p_one, width = 8, height = 5, dpi = 180)
  }
  message(
    "Wrote ", length(cell_types), " individual density plots to ", output_dir
  )

  readr::write_csv(
    data.frame(cell_type = cell_types, stringsAsFactors = FALSE),
    file.path(output_dir, paste0(file_prefix, "_cell_types.csv"))
  )

  invisible(list(
    cell_types = cell_types,
    panel_path = panel_path,
    plot_df = plot_df
  ))
}

run_xcell_ad_and_combined <- function(
  output_dir = NULL,
  ad_alpha = 0.05,
  require_all_cohorts = TRUE
) {
  if (!requireNamespace("readr", quietly = TRUE)) library(readr)
  p <- resolve_transcriptomics_paths()
  if (is.null(output_dir)) {
    output_dir <- file.path(p$base, "analysis_output", "combined")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  scores_long <- load_all_xcell_scores_long()
  ad_out <- run_xcell_ad_by_cell_type(
    scores_long,
    ad_alpha = ad_alpha,
    require_all_cohorts = require_all_cohorts
  )
  ad_path <- file.path(output_dir, "xcell_ad_by_cell_type.csv")
  readr::write_csv(ad_out$results, ad_path)
  message("Wrote ", ad_path)
  readr::write_csv(
    ad_out$n_samples,
    file.path(output_dir, "xcell_ad_cohort_sample_counts.csv")
  )
  if (length(ad_out$passing_cell_types) > 0L) {
    readr::write_csv(
      data.frame(cell_type = ad_out$passing_cell_types, stringsAsFactors = FALSE),
      file.path(output_dir, "xcell_ad_passing_cell_types.csv")
    )
  }

  can_combine <- ad_out$n_pass > 0L

  summary_lines <- c(
    paste0("Cohorts: ", paste(ad_out$cohorts, collapse = ", ")),
    paste0("Cell types tested (present in all cohorts): ", ad_out$n_cell_types),
    paste0("Pass AD (p >= ", ad_alpha, "): ", ad_out$n_pass),
    paste0("Samples per cohort: ", paste(
      ad_out$n_samples$cohort, ad_out$n_samples$n_samples, sep = "=", collapse = "; "
    )),
    "",
    "Anderson-Darling k-sample tests (k = number of cohorts):",
    "  One test per cell type — compares score distributions across datasets.",
    "  Not stratified by High/Low risk.",
    "",
    paste0("Cell types OK to merge: ", length(ad_out$passing_cell_types)),
    if (length(ad_out$passing_cell_types) > 0L) {
      paste0("  ", paste(ad_out$passing_cell_types, collapse = ", "))
    } else {
      "  (none)"
    },
    "",
    paste0("Alpha: ", ad_alpha)
  )
  writeLines(summary_lines, file.path(output_dir, "xcell_ad_summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")
  export_xcell_ad_summary_plot(ad_out, output_dir, ad_alpha)

  result <- list(
    scores_long = scores_long,
    ad_out = ad_out,
    passing_cell_types = ad_out$passing_cell_types,
    can_combine = can_combine
  )

  run_combined_volcano <- function(scores, label_suffix, file_prefix) {
    diff_combined <- wilcox_enrichment_diff(scores)
    n_samples <- dplyr::n_distinct(scores$sample_uid)
    n_low <- sum(scores$risk_grp == "Low") / dplyr::n_distinct(scores$cell_type)
    n_high <- sum(scores$risk_grp == "High") / dplyr::n_distinct(scores$cell_type)
    export_xcell_volcano(
      diff_combined,
      paste0("Cell Type Enrichment: High vs Low Risk (", label_suffix, ")"),
      paste0(
        "Pooled Stage 2 + Retrospective + Colossus + Taxonomy | xCell | n = ", n_samples,
        " (Low ", n_low, ", High ", n_high, ")"
      ),
      file.path(output_dir, paste0(file_prefix, "_risk_enrichment_volcano.png"))
    )
    export_top_log2fc_table(
      diff_combined,
      file.path(output_dir, paste0(file_prefix, "_top5_log2fc_up_down.csv"))
    )
    readr::write_csv(
      diff_combined,
      file.path(output_dir, paste0(file_prefix, "_wilcox_all_cell_types.csv"))
    )
    diff_combined
  }

  if (can_combine) {
    pooled <- build_xcell_combined_scores(
      scores_long,
      cell_types = ad_out$passing_cell_types
    )
    result$diff_combined <- run_combined_volcano(
      pooled,
      paste0(
        "AD-passing cell types pooled (n=",
        length(ad_out$passing_cell_types), ")"
      ),
      "Combined"
    )
    return(result)
  }

  message(
    "Combined volcano skipped: no cell types passed AD (p >= ", ad_alpha, ")."
  )
  result
}
