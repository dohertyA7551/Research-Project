# Proteomics differential abundance by TMA composition group (4 groups + optional Other).
source("tma_composition_groups_helpers.R")

PROTEOMICS_COMPARTMENTS <- data.frame(
  id = c("helperT", "cytT", "regT", "NonImmuneEpithelium", "NonImmuneStroma"),
  suffix = c("helperT", "cytT", "regT", "NonImmuneEpithelium", "NonImmuneStroma"),
  label = c(
    "Helper T cells",
    "Cytotoxic T cells",
    "Regulatory T cells",
    "Non-immune epithelium (tumour)",
    "Stromal cells"
  ),
  stringsAsFactors = FALSE
)

get_proteomics_cell_cols <- function(merged) {
  cell_cols <- grep(
    "_mean_(NonImmuneEpithelium|NonImmuneStroma|cytT|helperT|regT)$",
    colnames(merged),
    value = TRUE,
    ignore.case = TRUE
  )
  if (length(cell_cols) == 0L) {
    cell_cols <- grep("_mean_", colnames(merged), value = TRUE)
    cell_cols <- cell_cols[grepl("cytT|helperT|regT|Epithelium|Stroma", cell_cols, ignore.case = TRUE)]
  }
  cell_cols
}

get_proteomics_cols_for_compartment <- function(merged, suffix) {
  pat <- paste0("_mean_", suffix, "$")
  cols <- grep(pat, colnames(merged), value = TRUE, ignore.case = TRUE)
  if (length(cols) == 0L) {
    stop("No proteomics columns for compartment suffix: ", suffix)
  }
  cols
}

clean_feature_label <- function(x) {
  x <- sub("^baseline_Mean\\.Cell\\.", "", x, ignore.case = TRUE)
  x <- sub("^baseline_", "", x, ignore.case = TRUE)
  x <- sub("_mean_NonImmuneEpithelium$", " (epi)", x, ignore.case = TRUE)
  x <- sub("_mean_NonImmuneStroma$", " (stroma)", x, ignore.case = TRUE)
  x <- sub("_mean_cytT$", " (cytT)", x, ignore.case = TRUE)
  x <- sub("_mean_helperT$", " (helperT)", x, ignore.case = TRUE)
  x <- sub("_mean_regT$", " (regT)", x, ignore.case = TRUE)
  x
}

clean_marker_only <- function(x) {
  x <- sub("^baseline_Mean\\.Cell\\.", "", x, ignore.case = TRUE)
  x <- sub("^baseline_", "", x, ignore.case = TRUE)
  x <- sub("_mean_(NonImmuneEpithelium|NonImmuneStroma|cytT|helperT|regT)$", "", x, ignore.case = TRUE)
  x
}

build_tma_extreme_groups <- function(patient_grp) {
  patient_grp %>%
    dplyr::mutate(
      tma_extreme = dplyr::case_when(
        .data$n_tma_high == 0L ~ "all_low",
        .data$n_tma_high == .data$n_tma ~ "all_high",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(.data$tma_extreme))
}

export_proteomics_limma_volcano <- function(
  results,
  output_dir,
  file_prefix,
  title,
  subtitle,
  group_high = "G4_high_3plus",
  group_low = "G1_all_low",
  group_high_label = NULL,
  group_low_label = NULL,
  fc_cutoff = 0.5,
  p_cutoff = 0.05,
  fdr_cutoff = 0.05,
  top_labels = 15L,
  label_fn = clean_feature_label
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for volcano plots.")
  }
  high_lab <- if (is.null(group_high_label)) group_high else group_high_label
  low_lab <- if (is.null(group_low_label)) group_low else group_low_label

  df <- results %>%
    dplyr::filter(is.finite(.data$logFC), !is.na(.data$P.Value)) %>%
    dplyr::mutate(
      gene = .data$feature,
      label = label_fn(.data$feature),
      neg_log10_p = -log10(pmax(.data$P.Value, .Machine$double.xmin)),
      fc_class = dplyr::case_when(
        .data$logFC >= fc_cutoff ~ paste0("Higher in ", high_lab),
        .data$logFC <= -fc_cutoff ~ paste0("Higher in ", low_lab),
        TRUE ~ "Other"
      )
    )
  label_df <- df %>%
    dplyr::arrange(dplyr::desc(abs(.data$logFC))) %>%
    dplyr::slice_head(n = top_labels)

  xmax <- min(3, max(1.5, stats::quantile(abs(df$logFC), 0.98, na.rm = TRUE) * 1.1))
  ymax <- max(2, max(df$neg_log10_p, na.rm = TRUE) * 1.08, -log10(p_cutoff) * 1.15)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = logFC, y = neg_log10_p, colour = fc_class)
  ) +
    ggplot2::geom_point(size = 1.8, alpha = 0.65) +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", colour = "grey40") +
    ggplot2::geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed",
      colour = "grey40"
    ) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(
        c("#B2182B", "#2166AC", "grey75"),
        c(paste0("Higher in ", high_lab), paste0("Higher in ", low_lab), "Other")
      )
    ) +
    ggplot2::coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax)) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = paste0("log2 fold change (", high_lab, " / ", low_lab, ")"),
      y = expression("-log"[10] * "(p-value)"),
      colour = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = label),
      size = 2.8,
      max.overlaps = 25,
      show.legend = FALSE
    )
  }

  out_png <- file.path(output_dir, paste0(file_prefix, "_limma_volcano.png"))
  out_pdf <- file.path(output_dir, paste0(file_prefix, "_limma_volcano.pdf"))
  ggplot2::ggsave(out_png, p, width = 10, height = 7, dpi = 200, bg = "white")
  ggplot2::ggsave(out_pdf, p, width = 10, height = 7, bg = "white")
  message("Saved: ", out_png)
  invisible(list(plot = p, results = df))
}

run_tma_comp_proteomics_g1_vs_g4 <- function(
  groups_rds = "/work_space/files/proteomics/tma_composition_patient_groups.rds",
  merged_rds = "/work_space/files/proteomics/RCSI_Taxonomy_clin_prot_merged_df.rds",
  output_dir = "/work_space/files/proteomics/results/tma_composition_groups/proteomics",
  group_low = "G1_all_low",
  group_high = "G4_high_3plus",
  min_group = 5L
) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for proteomics DE.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  patient_grp <- readRDS(groups_rds)$patient_groups %>%
    dplyr::filter(.data$tma_comp_grp %in% c(group_low, group_high))

  n_low <- sum(patient_grp$tma_comp_grp == group_low)
  n_high <- sum(patient_grp$tma_comp_grp == group_high)
  if (n_low < min_group || n_high < min_group) {
    stop(
      "Insufficient samples for ", group_high, " vs ", group_low,
      ": n_low=", n_low, ", n_high=", n_high
    )
  }

  merged <- readRDS(merged_rds)
  cell_cols <- get_proteomics_cell_cols(merged)

  analysis_df <- merged %>%
    dplyr::inner_join(
      patient_grp %>% dplyr::select(patient_id, tma_comp_grp),
      by = "patient_id"
    )

  expr_mat <- t(as.matrix(analysis_df[, cell_cols, drop = FALSE]))
  colnames(expr_mat) <- analysis_df$patient_id

  comp_grp <- factor(
    analysis_df$tma_comp_grp[match(colnames(expr_mat), analysis_df$patient_id)],
    levels = c(group_low, group_high)
  )

  design <- stats::model.matrix(~ comp_grp)
  fit <- limma::lmFit(expr_mat, design)
  fit2 <- limma::eBayes(fit)

  results <- limma::topTable(
    fit2,
    coef = paste0("comp_grp", group_high),
    number = Inf,
    sort.by = "P"
  ) %>%
    tibble::rownames_to_column("feature")

  readr::write_csv(
    results,
    file.path(output_dir, "proteomics_limma_G1_all_low_vs_G4_high_3plus.csv")
  )

  sig <- results %>% dplyr::filter(.data$adj.P.Val <= 0.05, abs(.data$logFC) >= 0.5)
  readr::write_csv(
    sig,
    file.path(output_dir, "proteomics_limma_G1_all_low_vs_G4_high_3plus_sig.csv")
  )

  export_proteomics_limma_volcano(
    results,
    output_dir,
    file_prefix = "proteomics_G1_all_low_vs_G4_high_3plus",
    title = "Proteomics: TMA composition extremes",
    subtitle = "All TMAs Low vs >=3 TMAs High",
    group_high = group_high,
    group_low = group_low,
    group_high_label = ">=3 TMAs High",
    group_low_label = "All TMAs Low"
  )

  summary_lines <- c(
    paste0("Limma: ", group_low, " vs ", group_high),
    paste0("Features: ", length(cell_cols)),
    paste0("n ", group_low, ": ", n_low, " | n ", group_high, ": ", n_high),
    paste0("Significant (FDR<=0.05, |logFC|>=0.5): ", nrow(sig)),
    "",
    "Outputs:",
    "  proteomics_limma_G1_all_low_vs_G4_high_3plus.csv",
    "  proteomics_G1_all_low_vs_G4_high_3plus_limma_volcano.png"
  )
  writeLines(
    summary_lines,
    file.path(output_dir, "proteomics_G1_vs_G4_summary.txt")
  )
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(results)
}

run_tma_comp_proteomics_all_low_vs_all_high_by_compartment <- function(
  groups_rds = "/work_space/files/proteomics/tma_composition_patient_groups.rds",
  merged_rds = "/work_space/files/proteomics/RCSI_Taxonomy_clin_prot_merged_df.rds",
  output_dir = "/work_space/files/proteomics/results/tma_composition_groups/proteomics_all_low_vs_all_high_by_compartment",
  group_low = "all_low",
  group_high = "all_high",
  min_group = 5L
) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for proteomics DE.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  patient_grp <- readRDS(groups_rds)$patient_groups %>%
    build_tma_extreme_groups()

  n_low <- sum(patient_grp$tma_extreme == group_low)
  n_high <- sum(patient_grp$tma_extreme == group_high)
  if (n_low < min_group || n_high < min_group) {
    stop(
      "Insufficient samples for ", group_high, " vs ", group_low,
      ": n_low=", n_low, ", n_high=", n_high
    )
  }

  merged <- readRDS(merged_rds)
  analysis_df <- merged %>%
    dplyr::inner_join(
      patient_grp %>%
        dplyr::select(.data$patient_id, .data$tma_extreme, .data$n_tma, .data$n_tma_high),
      by = "patient_id"
    )

  all_results <- list()
  summary_rows <- list()

  for (i in seq_len(nrow(PROTEOMICS_COMPARTMENTS))) {
    comp <- PROTEOMICS_COMPARTMENTS[i, ]
    cell_cols <- get_proteomics_cols_for_compartment(merged, comp$suffix)
    message("\n--- ", comp$label, " (", length(cell_cols), " markers) ---")

    expr_mat <- t(as.matrix(analysis_df[, cell_cols, drop = FALSE]))
    colnames(expr_mat) <- analysis_df$patient_id
    comp_grp <- factor(
      analysis_df$tma_extreme[match(colnames(expr_mat), analysis_df$patient_id)],
      levels = c(group_low, group_high)
    )

    design <- stats::model.matrix(~ comp_grp)
    fit <- limma::lmFit(expr_mat, design)
    fit2 <- limma::eBayes(fit)
    results <- limma::topTable(
      fit2,
      coef = paste0("comp_grp", group_high),
      number = Inf,
      sort.by = "P"
    ) %>%
      tibble::rownames_to_column("feature") %>%
      dplyr::mutate(
        marker = clean_marker_only(.data$feature),
        compartment = comp$id,
        compartment_label = comp$label
      )

    sig <- results %>%
      dplyr::filter(.data$adj.P.Val <= 0.05, abs(.data$logFC) >= 0.5)
    prefix <- paste0("proteomics_", comp$id, "_all_low_vs_all_high")
    readr::write_csv(
      results,
      file.path(output_dir, paste0(prefix, "_limma.csv"))
    )
    readr::write_csv(
      sig,
      file.path(output_dir, paste0(prefix, "_limma_sig.csv"))
    )

    export_proteomics_limma_volcano(
      results,
      output_dir,
      file_prefix = prefix,
      title = paste0("Proteomics: ", comp$label),
      subtitle = "All TMAs Low vs All TMAs High",
      group_high = group_high,
      group_low = group_low,
      group_high_label = "All TMAs High",
      group_low_label = "All TMAs Low",
      label_fn = clean_marker_only
    )

    all_results[[comp$id]] <- results
    summary_rows[[comp$id]] <- data.frame(
      compartment = comp$id,
      compartment_label = comp$label,
      n_markers = length(cell_cols),
      n_sig_fdr05_fc05 = nrow(sig),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- dplyr::bind_rows(summary_rows)
  readr::write_csv(summary_df, file.path(output_dir, "compartment_summary.csv"))

  summary_lines <- c(
    "Proteomics: all TMAs Low vs all TMAs High — by cell compartment",
    paste0("n all_low: ", n_low, " | n all_high: ", n_high),
    "",
    "Each limma run uses markers measured in ONE compartment only",
    "(helper T, cytotoxic T, regulatory T, epithelium, stroma).",
    "",
    paste0("Output folder: ", output_dir),
    "",
    capture.output(print(summary_df)),
    "",
    "Per compartment:",
    "  proteomics_{compartment}_all_low_vs_all_high_limma.csv",
    "  proteomics_{compartment}_all_low_vs_all_high_limma_sig.csv",
    "  proteomics_{compartment}_all_low_vs_all_high_limma_volcano.png"
  )
  writeLines(summary_lines, file.path(output_dir, "summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(by_compartment = all_results, summary = summary_df))
}

run_tma_comp_proteomics_limma <- function(
  groups_rds = "/work_space/files/proteomics/tma_composition_patient_groups.rds",
  merged_rds = "/work_space/files/proteomics/RCSI_Taxonomy_clin_prot_merged_df.rds",
  output_dir = "/work_space/files/proteomics/results/tma_composition_groups/proteomics",
  exclude_other = TRUE,
  min_group = 5L
) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for proteomics DE.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  grp_obj <- readRDS(groups_rds)
  patient_grp <- grp_obj$patient_groups
  if (exclude_other) {
    patient_grp <- patient_grp %>%
      dplyr::filter(.data$tma_comp_grp != "Other_mixed")
  }

  merged <- readRDS(merged_rds)
  cell_cols <- get_proteomics_cell_cols(merged)

  analysis_df <- merged %>%
    dplyr::inner_join(
      patient_grp %>% dplyr::select(patient_id, tma_comp_grp),
      by = "patient_id"
    ) %>%
    dplyr::filter(!is.na(.data$tma_comp_grp))

  grp_tab <- table(analysis_df$tma_comp_grp)
  if (any(grp_tab < min_group)) {
    warning("Some groups below min_group=", min_group, ": ", paste(names(grp_tab), grp_tab, collapse = ", "))
  }

  expr_mat <- t(as.matrix(analysis_df[, cell_cols, drop = FALSE]))
  colnames(expr_mat) <- analysis_df$patient_id
  comp_grp <- factor(
    analysis_df$tma_comp_grp[match(colnames(expr_mat), analysis_df$patient_id)],
    levels = levels(patient_grp$tma_comp_grp)
  )
  comp_grp <- droplevels(comp_grp)

  design <- stats::model.matrix(~ 0 + comp_grp)
  colnames(design) <- make.names(colnames(design))

  fit <- limma::lmFit(expr_mat, design)
  fit2 <- limma::eBayes(fit)

  coef_names <- colnames(design)
  all_results <- list()
  for (i in seq_along(coef_names)) {
    if (i == 1L) {
      next
    }
    ref <- coef_names[1L]
    cn <- limma::makeContrasts(
      contrasts = paste0(coef_names[i], "-", ref),
      levels = design
    )
    fit_c <- limma::contrasts.fit(fit, cn)
    fit_c <- limma::eBayes(fit_c)
    res <- limma::topTable(fit_c, number = Inf, sort.by = "P") %>%
      tibble::rownames_to_column("feature")
    grp_label <- sub("^comp_grp", "", coef_names[i])
    res$contrast <- paste0(grp_label, "_vs_", sub("^comp_grp", "", ref))
    all_results[[length(all_results) + 1L]] <- res
  }

  results <- dplyr::bind_rows(all_results)
  readr::write_csv(results, file.path(output_dir, "proteomics_limma_vs_G1_all_low.csv"))

  kw_rows <- lapply(cell_cols, function(feat) {
    vals <- analysis_df[[feat]]
    kw <- tryCatch(
      stats::kruskal.test(vals ~ analysis_df$tma_comp_grp)$p.value,
      error = function(e) NA_real_
    )
    tibble(feature = feat, kruskal_p = kw)
  })
  kw_df <- dplyr::bind_rows(kw_rows) %>%
    mutate(neg_log10_p = -log10(.data$kruskal_p)) %>%
    arrange(.data$kruskal_p)
  readr::write_csv(kw_df, file.path(output_dir, "proteomics_kruskal_4groups.csv"))

  sig <- kw_df %>% filter(.data$kruskal_p <= 0.05)
  summary_lines <- c(
    "Proteomics by TMA composition group",
    paste0("Cell-type features: ", length(cell_cols)),
    paste0("Patients: ", nrow(analysis_df)),
    paste(capture.output(print(grp_tab)), collapse = "\n"),
    paste0("Kruskal FDR<=0.05 features: ", nrow(sig)),
    "",
    "Contrasts vs G1_all_low written to proteomics_limma_vs_G1_all_low.csv"
  )
  writeLines(summary_lines, file.path(output_dir, "proteomics_summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(limma = results, kruskal = kw_df))
}
