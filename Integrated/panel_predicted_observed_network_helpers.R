# Predicted (bulk RNA + STRING PPI) vs observed (epithelial proteomics co-abundance)
# for the TMA CODEX marker panel, High vs Low RSF risk.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

WORKSPACE_ROOT <- "/work_space/files"
TRANSCRIPTOMICS_DIR <- file.path(WORKSPACE_ROOT, "Transcriptomics")
PROTEOMICS_DIR <- file.path(WORKSPACE_ROOT, "proteomics")
INTEGRATED_DIR <- file.path(WORKSPACE_ROOT, "Integrated")
INTEGRATED_DATA_DIR <- file.path(WORKSPACE_ROOT, "Integrated ")

DEFAULT_COMPARTMENT_SUFFIX <- "NonImmuneEpithelium"
DEFAULT_STRING_SCORE <- 400L
DEFAULT_CORR_THRESHOLD <- 0.3
DEFAULT_MIN_PATIENTS <- 20L
DEFAULT_HUB_TOP_N <- 10L

PANEL_MANUAL_GENE_MAP <- c(
  CD3 = "CD3D",
  COLIV = "COL4A1",
  BCATENIN = "CTNNB1",
  CMYC = "MYC",
  KI67 = "MKI67",
  GRANZYMEB = "GZMB",
  HLA1 = "HLA-A",
  BCLXL = "BCL2L1",
  FLIP_CS = "CFLAR",
  PD1 = "PDCD1",
  CIAP1 = "BIRC2",
  S6 = "RPS6",
  SMA = "ACTA2",
  CASP3CLEAVED = "CASP3",
  PCK26 = "KRT26",
  SR2B = "HTR2B",
  MUC5 = "MUC5AC",
  AE1 = "SLC4A1",
  PCAD = "CDH17",
  NAK = "NAGK"
)

resolve_panel_network_paths <- function() {
  list(
    proteomics = PROTEOMICS_DIR,
    transcriptomics = TRANSCRIPTOMICS_DIR,
    integrated = INTEGRATED_DIR,
    integrated_data = INTEGRATED_DATA_DIR,
    spatial = file.path(WORKSPACE_ROOT, "spatial"),
    merged_rds = file.path(PROTEOMICS_DIR, "RCSI_Taxonomy_clin_prot_merged_df.rds"),
    survival_rds = file.path(PROTEOMICS_DIR, "survival_df_with_risk.rds"),
    panel_de_csv = file.path(WORKSPACE_ROOT, "Features", "marker_de_patient_level_limma.csv"),
    output = file.path(INTEGRATED_DIR, "results", "panel_predicted_vs_observed")
  )
}

marker_from_epi_col <- function(col) {
  marker_from_compartment_col(col, DEFAULT_COMPARTMENT_SUFFIX)
}

marker_from_compartment_col <- function(col, compartment_suffix) {
  pat <- paste0(
    "^baseline_Mean\\.Cell\\.(.+)_mean_",
    compartment_suffix,
    "$"
  )
  sub(pat, "\\1", col, ignore.case = TRUE)
}

marker_primary_gene <- function(marker) {
  marker <- as.character(marker)
  if (marker %in% names(PANEL_MANUAL_GENE_MAP)) {
    return(unname(PANEL_MANUAL_GENE_MAP[[marker]]))
  }
  toupper(marker)
}

load_panel_marker_table <- function(panel_csv = NULL) {
  paths <- resolve_panel_network_paths()
  if (is.null(panel_csv)) {
    panel_csv <- paths$panel_de_csv
  }
  if (!file.exists(panel_csv)) {
    stop("Missing panel marker list: ", panel_csv)
  }
  markers <- readr::read_csv(panel_csv, show_col_types = FALSE)$marker
  markers <- unique(as.character(markers))
  tibble::tibble(
    marker = markers,
    gene_symbol = vapply(markers, marker_primary_gene, character(1))
  )
}

load_epi_protein_matrix <- function(
  merged_rds = NULL,
  survival_rds = NULL,
  compartment_suffix = DEFAULT_COMPARTMENT_SUFFIX
) {
  paths <- resolve_panel_network_paths()
  if (is.null(merged_rds)) merged_rds <- paths$merged_rds
  if (is.null(survival_rds)) survival_rds <- paths$survival_rds
  merged <- readRDS(merged_rds)
  surv <- readRDS(survival_rds)
  pat <- intersect(merged$patient_id, surv$patient_id)
  merged <- merged[match(pat, merged$patient_id), , drop = FALSE]
  surv <- surv[match(pat, surv$patient_id), , drop = FALSE]

  pat_suffix <- paste0("_mean_", compartment_suffix, "$")
  epi_cols <- grep(pat_suffix, names(merged), value = TRUE, ignore.case = TRUE)
  if (length(epi_cols) == 0L) {
    stop("No proteomics columns for compartment suffix: ", compartment_suffix)
  }

  mat <- as.matrix(merged[, epi_cols, drop = FALSE])
  rownames(mat) <- merged$patient_id
  storage.mode(mat) <- "double"
  markers <- vapply(epi_cols, marker_from_compartment_col, character(1), compartment_suffix)

  mat <- t(mat)
  rownames(mat) <- markers
  meta <- data.frame(
    patient_id = merged$patient_id,
    risk_grp = surv$risk_grp,
    rsf_risk = surv$rsf_risk,
    stringsAsFactors = FALSE
  )
  list(matrix = mat, meta = meta, markers = markers, feature_cols = epi_cols)
}

collapse_expr_to_patient_raw <- function(expr_mat, sample_meta) {
  sample_meta <- sample_meta %>%
    dplyr::filter(
      !is.na(.data$patient_id),
      .data$sample_id %in% colnames(expr_mat)
    )
  if (nrow(sample_meta) == 0L) {
    return(list(expr = NULL, meta = NULL))
  }
  expr <- expr_mat[, sample_meta$sample_id, drop = FALSE]
  pts <- unique(sample_meta$patient_id)
  out <- matrix(NA_real_, nrow = nrow(expr), ncol = length(pts))
  rownames(out) <- rownames(expr)
  colnames(out) <- pts
  cohort_vec <- setNames(rep(NA_character_, length(pts)), pts)
  for (pt in pts) {
    ids <- sample_meta$sample_id[sample_meta$patient_id == pt]
    if (length(ids) == 1L) {
      out[, pt] <- expr[, ids]
    } else {
      out[, pt] <- rowMeans(expr[, ids, drop = FALSE], na.rm = TRUE)
    }
    cohort_vec[pt] <- sample_meta$cohort[match(ids[1], sample_meta$sample_id)]
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

limma_cohort_paths <- function() {
  base <- file.path(TRANSCRIPTOMICS_DIR, "analysis_output")
  list(
    Stage2 = file.path(base, "stage2", "Stage2_limma_high_vs_low_all_genes.csv"),
    Retrospective = file.path(base, "retrospective", "Retrospective_limma_high_vs_low_all_genes.csv"),
    Colossus = file.path(base, "colossus", "Colossus_limma_high_vs_low_all_genes.csv"),
    Taxonomy = file.path(base, "taxonomy", "Taxonomy_limma_high_vs_low_all_genes.csv")
  )
}

fisher_combine_pvals <- function(pvals) {
  pvals <- as.numeric(pvals)
  pvals <- pvals[is.finite(pvals) & pvals > 0 & pvals <= 1]
  if (length(pvals) == 0L) {
    return(NA_real_)
  }
  if (length(pvals) == 1L) {
    return(pvals[1])
  }
  stats::pchisq(-2 * sum(log(pvals)), df = 2 * length(pvals), lower.tail = FALSE)
}

load_rna_panel_de_meta <- function(panel_df) {
  paths <- limma_cohort_paths()
  genes <- unique(panel_df$gene_symbol)
  pieces <- list()
  for (cn in names(paths)) {
    if (!file.exists(paths[[cn]])) {
      message("Skipping missing limma file: ", paths[[cn]])
      next
    }
    df <- readr::read_csv(paths[[cn]], show_col_types = FALSE)
    df <- df %>% dplyr::filter(.data$gene %in% genes)
    if (nrow(df) == 0L) next
    pieces[[cn]] <- df %>%
      dplyr::select(.data$gene, .data$logFC, .data$P.Value, .data$adj.P.Val) %>%
      dplyr::mutate(cohort = cn)
  }
  if (length(pieces) == 0L) {
    stop("No cohort limma results found for panel genes.")
  }
  long <- dplyr::bind_rows(pieces)
  meta_df <- long %>%
    dplyr::group_by(.data$gene) %>%
    dplyr::summarise(
      logFC = mean(.data$logFC, na.rm = TRUE),
      P.Value = fisher_combine_pvals(.data$P.Value),
      adj.P.Val = fisher_combine_pvals(.data$adj.P.Val),
      n_cohorts = dplyr::n(),
      cohorts = paste(sort(unique(.data$cohort)), collapse = ";"),
      .groups = "drop"
    )
  list(results = meta_df, by_cohort = long)
}

subset_protein_high_low <- function(prot, min_patients = DEFAULT_MIN_PATIENTS) {
  meta <- prot$meta %>%
    dplyr::filter(.data$risk_grp %in% c("Low", "High"))
  pts <- meta$patient_id
  if (length(pts) < min_patients) {
    stop("Too few protein patients with High/Low risk: ", length(pts))
  }
  list(
    matrix = prot$matrix[, pts, drop = FALSE],
    meta = meta,
    markers = rownames(prot$matrix),
    patients = pts
  )
}

run_limma_high_vs_low_matrix <- function(expr_mat, meta, feature_col = "feature") {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required.")
  }
  meta <- meta %>%
    dplyr::filter(.data$risk_grp %in% c("Low", "High")) %>%
    dplyr::filter(.data$patient_id %in% colnames(expr_mat))
  n_low <- sum(meta$risk_grp == "Low")
  n_high <- sum(meta$risk_grp == "High")
  if (n_low < 3L || n_high < 3L) {
    stop("Need >=3 patients per risk group; Low=", n_low, ", High=", n_high)
  }
  expr <- expr_mat[, meta$patient_id, drop = FALSE]
  risk_grp <- factor(
    meta$risk_grp[match(colnames(expr), meta$patient_id)],
    levels = c("Low", "High")
  )
  has_cohort <- "cohort" %in% names(meta) &&
    length(unique(meta$cohort[!is.na(meta$cohort)])) > 1L
  if (has_cohort) {
    cohort <- factor(meta$cohort[match(colnames(expr), meta$patient_id)])
    design <- stats::model.matrix(~ 0 + risk_grp + cohort)
  } else {
    design <- stats::model.matrix(~ 0 + risk_grp)
  }
  colnames(design) <- make.names(colnames(design))
  contrast <- limma::makeContrasts(
    HighvsLow = risk_grpHigh - risk_grpLow,
    levels = design
  )
  fit <- limma::lmFit(expr, design)
  fit2 <- limma::contrasts.fit(fit, contrast)
  fit2 <- limma::eBayes(fit2)
  results <- limma::topTable(
    fit2,
    coef = "HighvsLow",
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  results[[feature_col]] <- rownames(results)
  list(
    results = results,
    n_low = n_low,
    n_high = n_high,
    cohort_adjusted = has_cohort
  )
}

build_spearman_network <- function(mat, corr_threshold = DEFAULT_CORR_THRESHOLD) {
  markers <- rownames(mat)
  n <- ncol(mat)
  if (n < 4L) {
    return(list(edges = data.frame(), corr = matrix(NA, 0, 0)))
  }
  corr <- stats::cor(t(mat), method = "spearman", use = "pairwise.complete.obs")
  diag(corr) <- NA
  edges <- which(
    upper.tri(corr) & is.finite(corr) & abs(corr) >= corr_threshold,
    arr.ind = TRUE
  )
  if (nrow(edges) == 0L) {
    return(list(edges = data.frame(), corr = corr))
  }
  edge_df <- data.frame(
    from = markers[edges[, 1]],
    to = markers[edges[, 2]],
    spearman_r = corr[edges],
    abs_r = abs(corr[edges]),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(dplyr::desc(.data$abs_r))
  list(edges = edge_df, corr = corr)
}

build_string_panel_network <- function(
  panel_df,
  de_df,
  score_threshold = DEFAULT_STRING_SCORE,
  string_version = "12.0"
) {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("Package 'STRINGdb' is required.")
  }
  genes <- unique(panel_df$gene_symbol)
  gene_de <- de_df %>%
    dplyr::filter(.data$gene %in% genes) %>%
    dplyr::distinct(.data$gene, .keep_all = TRUE)
  map_in <- panel_df %>%
    dplyr::distinct(.data$marker, .data$gene_symbol) %>%
    dplyr::left_join(
      gene_de %>% dplyr::rename(gene_symbol = gene),
      by = "gene_symbol"
    )
  string_db <- STRINGdb::STRINGdb$new(
    version = string_version,
    species = 9606L,
    score_threshold = score_threshold
  )
  gene_tbl <- map_in %>%
    dplyr::distinct(.data$gene_symbol) %>%
    dplyr::filter(!is.na(.data$gene_symbol), nzchar(.data$gene_symbol))
  mapped <- string_db$map(
    as.data.frame(gene_tbl),
    "gene_symbol",
    removeUnmappedRows = TRUE
  )
  if (nrow(mapped) < 2L) {
    return(list(
      mapped = mapped,
      interactions = data.frame(),
      string_ids = character(0)
    ))
  }
  string_ids <- unique(mapped$STRING_id)
  interactions <- tryCatch(
    string_db$get_interactions(string_ids),
    error = function(e) {
      message("STRING get_interactions failed: ", conditionMessage(e))
      data.frame()
    }
  )
  list(
    string_db = string_db,
    mapped = mapped,
    interactions = interactions,
    string_ids = string_ids,
    map_in = map_in
  )
}

string_interactions_to_gene_edges <- function(interactions, mapped) {
  if (is.null(interactions) || nrow(interactions) == 0L) {
    return(data.frame())
  }
  id2gene <- mapped %>%
    dplyr::distinct(.data$STRING_id, .data$gene_symbol) %>%
    dplyr::group_by(.data$STRING_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  lookup <- setNames(id2gene$gene_symbol, id2gene$STRING_id)
  from_gene <- lookup[interactions$from]
  to_gene <- lookup[interactions$to]
  ok <- !is.na(from_gene) & !is.na(to_gene) & from_gene != to_gene
  out <- data.frame(
    from_gene = from_gene[ok],
    to_gene = to_gene[ok],
    combined_score = interactions$combined_score[ok],
    stringsAsFactors = FALSE
  )
  out <- out %>%
    dplyr::mutate(
      gene_a = pmin(.data$from_gene, .data$to_gene),
      gene_b = pmax(.data$from_gene, .data$to_gene)
    ) %>%
    dplyr::group_by(.data$gene_a, .data$gene_b) %>%
    dplyr::slice_max(.data$combined_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$gene_a, -.data$gene_b)
  out
}

gene_edge_to_marker_edge <- function(gene_edges, panel_df) {
  if (nrow(gene_edges) == 0L) {
    return(data.frame())
  }
  gene_to_marker <- panel_df %>%
    dplyr::group_by(.data$gene_symbol) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  lookup <- setNames(gene_to_marker$marker, gene_to_marker$gene_symbol)
  data.frame(
    from = lookup[gene_edges$from_gene],
    to = lookup[gene_edges$to_gene],
    from_gene = gene_edges$from_gene,
    to_gene = gene_edges$to_gene,
    combined_score = gene_edges$combined_score,
    edge_type = "predicted_string",
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(.data$from), !is.na(.data$to))
}

compute_network_centrality <- function(edges, nodes, id_col = "node") {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.")
  }
  if (length(nodes) == 0L) {
    return(data.frame())
  }
  node_df <- data.frame(node = nodes, stringsAsFactors = FALSE)
  names(node_df)[1] <- id_col
  if (is.null(edges) || nrow(edges) == 0L) {
    node_df$degree <- 0L
    node_df$betweenness <- 0
    return(node_df)
  }
  g <- igraph::graph_from_data_frame(
    edges[, c("from", "to"), drop = FALSE],
    directed = FALSE,
    vertices = nodes
  )
  node_df$degree <- igraph::degree(g)
  node_df$betweenness <- igraph::betweenness(g, directed = FALSE)
  node_df
}

annotate_nodes_with_de <- function(
  nodes_df,
  id_col,
  de_df,
  id_in_de,
  panel_df = NULL
) {
  de_sub <- de_df %>%
    dplyr::distinct(.data[[id_in_de]], .keep_all = TRUE)
  out <- nodes_df %>%
    dplyr::left_join(
      de_sub,
      by = setNames(id_in_de, id_col)
    )
  if (!is.null(panel_df)) {
    out <- out %>%
      dplyr::left_join(
        panel_df %>% dplyr::select(.data$marker, .data$gene_symbol),
        by = setNames("marker", id_col)
      )
  }
  out
}

compare_predicted_observed_edges <- function(
  predicted_edges,
  observed_edges,
  corr_lookup = NULL
) {
  pred_key <- if (nrow(predicted_edges) > 0L) {
    paste(
      pmin(predicted_edges$from, predicted_edges$to),
      pmax(predicted_edges$from, predicted_edges$to),
      sep = "|"
    )
  } else {
    character(0)
  }
  obs_key <- if (nrow(observed_edges) > 0L) {
    paste(
      pmin(observed_edges$from, observed_edges$to),
      pmax(observed_edges$from, observed_edges$to),
      sep = "|"
    )
  } else {
    character(0)
  }
  all_keys <- unique(c(pred_key, obs_key))
  if (length(all_keys) == 0L) {
    return(data.frame())
  }
  split_key <- strsplit(all_keys, "|", fixed = TRUE)
  out <- data.frame(
    marker_a = vapply(split_key, `[`, character(1), 1),
    marker_b = vapply(split_key, `[`, character(1), 2),
    in_predicted = all_keys %in% pred_key,
    in_observed = all_keys %in% obs_key,
    stringsAsFactors = FALSE
  )
  if (nrow(predicted_edges) > 0L) {
    pred_map <- setNames(predicted_edges$combined_score, pred_key)
    out$string_score <- pred_map[all_keys]
    out$from_gene <- NA_character_
    out$to_gene <- NA_character_
    for (i in seq_len(nrow(predicted_edges))) {
      k <- pred_key[i]
      idx <- which(all_keys == k)
      out$from_gene[idx] <- predicted_edges$from_gene[i]
      out$to_gene[idx] <- predicted_edges$to_gene[i]
    }
  } else {
    out$string_score <- NA_real_
    out$from_gene <- NA_character_
    out$to_gene <- NA_character_
  }
  if (!is.null(corr_lookup) && length(corr_lookup) > 0L) {
    out$spearman_r <- corr_lookup[all_keys]
  } else if (nrow(observed_edges) > 0L) {
    obs_map <- setNames(observed_edges$spearman_r, obs_key)
    out$spearman_r <- obs_map[all_keys]
  } else {
    out$spearman_r <- NA_real_
  }
  out$concordance <- dplyr::case_when(
    out$in_predicted & out$in_observed ~ "supported",
    out$in_predicted & !out$in_observed ~ "predicted_only",
    !out$in_predicted & out$in_observed ~ "observed_only",
    TRUE ~ "neither"
  )
  out
}

jaccard_top_terms <- function(a, b) {
  a <- unique(na.omit(a))
  b <- unique(na.omit(b))
  if (length(a) == 0L && length(b) == 0L) return(NA_real_)
  if (length(a) == 0L || length(b) == 0L) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

run_string_enrichment <- function(string_db, string_ids) {
  if (length(string_ids) < 3L) {
    return(data.frame())
  }
  tryCatch(
    string_db$get_enrichment(string_ids),
    error = function(e) {
      message("STRING enrichment failed: ", conditionMessage(e))
      data.frame()
    }
  )
}

export_de_scatter <- function(
  panel_df,
  rna_de,
  prot_de,
  out_path,
  title = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Skipping DE scatter (ggplot2 not available)")
    return(invisible(NULL))
  }
  df <- panel_df %>%
    dplyr::left_join(
      rna_de %>%
        dplyr::select(gene, rna_logFC = logFC, rna_adj.P.Val = adj.P.Val),
      by = c("gene_symbol" = "gene")
    ) %>%
    dplyr::left_join(
      prot_de %>%
        dplyr::select(marker, protein_logFC = logFC, protein_adj.P.Val = adj.P.Val),
      by = "marker"
    )
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = rna_logFC, y = protein_logFC, label = marker)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::geom_text(check_overlap = TRUE, hjust = -0.1, vjust = 0, size = 2.8) +
    ggplot2::labs(
      x = "Bulk RNA logFC (High vs Low)",
      y = "Epithelial protein logFC (High vs Low)",
      title = title %||% "Panel marker DE concordance"
    ) +
    ggplot2::theme_bw()
  ggplot2::ggsave(out_path, p, width = 8, height = 6, dpi = 150)
  message("Wrote ", out_path)
  invisible(out_path)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

run_panel_predicted_vs_observed <- function(
  output_dir = NULL,
  corr_threshold = DEFAULT_CORR_THRESHOLD,
  string_score = DEFAULT_STRING_SCORE,
  hub_top_n = DEFAULT_HUB_TOP_N,
  fdr_report = 0.05,
  seed = 42L
) {
  set.seed(seed)
  paths <- resolve_panel_network_paths()
  if (is.null(output_dir)) {
    output_dir <- paths$output
  }
  dirs <- list(
    root = output_dir,
    gene_lists = file.path(output_dir, "gene_lists"),
    networks = file.path(output_dir, "networks"),
    enrichment = file.path(output_dir, "enrichment"),
    comparison = file.path(output_dir, "comparison"),
    plots = file.path(output_dir, "plots")
  )
  lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  log_lines <- c(
    "Panel predicted (RNA + STRING) vs observed (epithelial proteomics)",
    paste0("Compartment: ", DEFAULT_COMPARTMENT_SUFFIX),
    "RNA DE: Fisher meta across 4 transcriptomic cohort limma tables (panel genes)",
    "Protein DE: patient-level limma on epithelial abundances (High vs Low RSF risk)",
    paste0("STRING score >= ", string_score),
    paste0("Observed edge threshold: |Spearman r| >= ", corr_threshold),
    paste0("Seed: ", seed),
    ""
  )

  panel_df <- load_panel_marker_table()
  log_lines <- c(log_lines, paste0("Panel markers: ", nrow(panel_df)))

  prot <- load_epi_protein_matrix()
  aligned <- subset_protein_high_low(prot)
  rna_meta <- load_rna_panel_de_meta(panel_df)
  rna_de <- list(results = rna_meta$results)
  log_lines <- c(
    log_lines,
    paste0(
      "Protein patients (High/Low): ", length(aligned$patients),
      " (Low=", sum(aligned$meta$risk_grp == "Low"),
      ", High=", sum(aligned$meta$risk_grp == "High"), ")"
    ),
    paste0("RNA panel genes in meta-analysis: ", nrow(rna_de$results)),
    ""
  )

  prot_de <- run_limma_high_vs_low_matrix(
    aligned$matrix,
    aligned$meta,
    feature_col = "marker"
  )
  log_lines <- c(
    log_lines,
    paste0("Protein limma cohort-adjusted: ", prot_de$cohort_adjusted),
    ""
  )

  readr::write_csv(prot_de$results, file.path(dirs$gene_lists, "protein_epi_limma_all_markers.csv"))
  readr::write_csv(rna_de$results, file.path(dirs$gene_lists, "rna_panel_limma_meta_all_genes.csv"))
  readr::write_csv(rna_meta$by_cohort, file.path(dirs$gene_lists, "rna_panel_limma_by_cohort.csv"))

  de_compare <- panel_df %>%
    dplyr::left_join(
      prot_de$results %>%
        dplyr::select(marker, protein_logFC = logFC, protein_adj.P.Val = adj.P.Val, protein_P.Value = P.Value),
      by = "marker"
    ) %>%
    dplyr::left_join(
      rna_de$results %>%
        dplyr::select(gene, rna_logFC = logFC, rna_adj.P.Val = adj.P.Val, rna_P.Value = P.Value),
      by = c("gene_symbol" = "gene")
    )
  readr::write_csv(de_compare, file.path(dirs$gene_lists, "panel_de_rna_vs_protein.csv"))

  readr::write_csv(
    de_compare %>% dplyr::filter(.data$rna_logFC > 0 | .data$protein_logFC > 0),
    file.path(dirs$gene_lists, "panel_up_in_high_rna_or_protein.csv")
  )
  readr::write_csv(
    de_compare %>% dplyr::filter(.data$rna_logFC < 0 | .data$protein_logFC < 0),
    file.path(dirs$gene_lists, "panel_down_in_high_rna_or_protein.csv")
  )

  string_net <- build_string_panel_network(panel_df, rna_de$results, score_threshold = string_score)
  gene_edges <- string_interactions_to_gene_edges(string_net$interactions, string_net$mapped)
  pred_marker_edges <- gene_edge_to_marker_edge(gene_edges, panel_df)
  readr::write_csv(gene_edges, file.path(dirs$networks, "predicted_string_gene_edges.csv"))
  readr::write_csv(pred_marker_edges, file.path(dirs$networks, "predicted_string_marker_edges.csv"))

  obs_net <- build_spearman_network(aligned$matrix, corr_threshold = corr_threshold)
  obs_edges <- obs_net$edges
  if (nrow(obs_edges) > 0L) {
    obs_edges$edge_type <- "observed_spearman"
  }
  readr::write_csv(obs_edges, file.path(dirs$networks, "observed_epi_spearman_edges.csv"))

  pred_nodes <- compute_network_centrality(pred_marker_edges, panel_df$marker, id_col = "marker")
  obs_nodes <- compute_network_centrality(obs_edges, panel_df$marker, id_col = "marker")
  pred_nodes <- pred_nodes %>%
    dplyr::left_join(panel_df, by = "marker") %>%
    dplyr::left_join(
      rna_de$results %>%
        dplyr::select(gene, rna_logFC = logFC, rna_adj.P.Val = adj.P.Val, rna_P.Value = P.Value),
      by = c("gene_symbol" = "gene")
    ) %>%
    dplyr::left_join(
      prot_de$results %>%
        dplyr::select(marker, protein_logFC = logFC, protein_adj.P.Val = adj.P.Val, protein_P.Value = P.Value),
      by = "marker"
    )
  obs_nodes <- obs_nodes %>%
    dplyr::left_join(panel_df, by = "marker") %>%
    dplyr::left_join(
      prot_de$results %>%
        dplyr::select(marker, protein_logFC = logFC, protein_adj.P.Val = adj.P.Val, protein_P.Value = P.Value),
      by = "marker"
    ) %>%
    dplyr::left_join(
      rna_de$results %>%
        dplyr::select(gene, rna_logFC = logFC, rna_adj.P.Val = adj.P.Val, rna_P.Value = P.Value),
      by = c("gene_symbol" = "gene")
    )
  readr::write_csv(pred_nodes, file.path(dirs$networks, "predicted_string_nodes.csv"))
  readr::write_csv(obs_nodes, file.path(dirs$networks, "observed_epi_spearman_nodes.csv"))

  pred_hubs <- pred_nodes %>%
    dplyr::arrange(dplyr::desc(.data$degree), dplyr::desc(.data$betweenness)) %>%
    dplyr::slice_head(n = hub_top_n)
  obs_hubs <- obs_nodes %>%
    dplyr::arrange(dplyr::desc(.data$degree), dplyr::desc(.data$betweenness)) %>%
    dplyr::slice_head(n = hub_top_n)
  readr::write_csv(pred_hubs, file.path(dirs$comparison, "predicted_top_hubs.csv"))
  readr::write_csv(obs_hubs, file.path(dirs$comparison, "observed_top_hubs.csv"))

  hub_concordance <- pred_nodes %>%
    dplyr::select(marker, gene_symbol, predicted_degree = degree, predicted_betweenness = betweenness) %>%
    dplyr::left_join(
      obs_nodes %>%
        dplyr::select(marker, observed_degree = degree, observed_betweenness = betweenness),
      by = "marker"
    ) %>%
    dplyr::mutate(
      in_predicted_top = .data$marker %in% pred_hubs$marker,
      in_observed_top = .data$marker %in% obs_hubs$marker,
      hub_concordance = .data$in_predicted_top & .data$in_observed_top
    ) %>%
    dplyr::arrange(dplyr::desc(.data$hub_concordance), dplyr::desc(.data$predicted_degree))
  readr::write_csv(hub_concordance, file.path(dirs$comparison, "hub_gene_concordance.csv"))

  corr_keys <- if (nrow(obs_edges) > 0L) {
    setNames(
      obs_edges$spearman_r,
      paste(
        pmin(obs_edges$from, obs_edges$to),
        pmax(obs_edges$from, obs_edges$to),
        sep = "|"
      )
    )
  } else {
    NULL
  }
  edge_disc <- compare_predicted_observed_edges(
    pred_marker_edges,
    obs_edges,
    corr_lookup = corr_keys
  )
  readr::write_csv(edge_disc, file.path(dirs$comparison, "edge_discordance.csv"))

  enrich_dirs <- list(
    rna_up = de_compare %>% dplyr::filter(.data$rna_logFC > 0) %>% dplyr::pull(gene_symbol) %>% unique(),
    rna_down = de_compare %>% dplyr::filter(.data$rna_logFC < 0) %>% dplyr::pull(gene_symbol) %>% unique(),
    protein_up = de_compare %>% dplyr::filter(.data$protein_logFC > 0) %>% dplyr::pull(gene_symbol) %>% unique(),
    protein_down = de_compare %>% dplyr::filter(.data$protein_logFC < 0) %>% dplyr::pull(gene_symbol) %>% unique()
  )
  enrich_results <- list()
  if (!is.null(string_net$string_db) && length(string_net$string_ids) >= 3L) {
    for (nm in names(enrich_dirs)) {
      genes <- enrich_dirs[[nm]]
      mapped <- string_net$mapped %>%
        dplyr::filter(.data$gene_symbol %in% genes)
      ids <- unique(mapped$STRING_id)
      er <- run_string_enrichment(string_net$string_db, ids)
      if (nrow(er) > 0L) {
        er$set <- nm
        readr::write_csv(er, file.path(dirs$enrichment, paste0(nm, "_string_enrichment.csv")))
        enrich_results[[nm]] <- er
      }
    }
  }

  top_enrich_terms <- function(er, n = 20L) {
    if (is.null(er) || nrow(er) == 0L || !"term" %in% names(er)) {
      return(character(0))
    }
    ord <- order(er$fdr)
    head(er$term[ord], n)
  }
  pathway_overlap <- tibble::tibble(
    comparison = c("up_in_high", "down_in_high"),
    jaccard_top20_terms = c(
      jaccard_top_terms(
        top_enrich_terms(enrich_results$rna_up),
        top_enrich_terms(enrich_results$protein_up)
      ),
      jaccard_top_terms(
        top_enrich_terms(enrich_results$rna_down),
        top_enrich_terms(enrich_results$protein_down)
      )
    ),
    n_rna_genes_in_set = c(
      length(enrich_dirs$rna_up),
      length(enrich_dirs$rna_down)
    ),
    n_protein_genes_in_set = c(
      length(enrich_dirs$protein_up),
      length(enrich_dirs$protein_down)
    ),
    n_panel_rna_fdr05 = sum(de_compare$rna_adj.P.Val < fdr_report, na.rm = TRUE),
    n_panel_protein_fdr05 = sum(de_compare$protein_adj.P.Val < fdr_report, na.rm = TRUE)
  )
  readr::write_csv(pathway_overlap, file.path(dirs$comparison, "pathway_overlap.csv"))

  export_de_scatter(
    panel_df,
    rna_de$results,
    prot_de$results,
    file.path(dirs$plots, "panel_de_rna_vs_protein_scatter.png"),
    title = "TMA panel: bulk RNA vs epithelial protein (High vs Low)"
  )

  if (!is.null(string_net$string_db) && length(string_net$string_ids) >= 2L) {
    png_path <- file.path(dirs$plots, "predicted_string_network.png")
    grDevices::png(png_path, width = 1400, height = 1000, res = 120)
    tryCatch(
      {
        string_net$string_db$plot_network(string_net$string_ids, add_link = FALSE)
        grDevices::dev.off()
        message("Wrote ", png_path)
      },
      error = function(e) {
        if (grDevices::dev.cur() > 1L) grDevices::dev.off()
        message("STRING plot failed: ", conditionMessage(e))
      }
    )
  }

  log_lines <- c(
    log_lines,
    "--- Networks ---",
    paste0("Predicted STRING gene edges: ", nrow(gene_edges)),
    paste0("Predicted STRING marker edges: ", nrow(pred_marker_edges)),
    paste0("Observed Spearman edges (|r|>=", corr_threshold, "): ", nrow(obs_edges)),
    paste0(
      "Edge concordance: supported=",
      sum(edge_disc$concordance == "supported", na.rm = TRUE),
      ", predicted_only=", sum(edge_disc$concordance == "predicted_only", na.rm = TRUE),
      ", observed_only=", sum(edge_disc$concordance == "observed_only", na.rm = TRUE)
    ),
    paste0(
      "Hub overlap (top ", hub_top_n, "): ",
      sum(hub_concordance$in_predicted_top & hub_concordance$in_observed_top, na.rm = TRUE),
      " markers in both"
    ),
    "",
    paste0("Output root: ", output_dir)
  )
  writeLines(log_lines, file.path(output_dir, "run_summary.txt"))
  cat(paste(log_lines, collapse = "\n"), "\n")
  invisible(list(
    panel = panel_df,
    de_compare = de_compare,
    predicted = list(nodes = pred_nodes, edges = pred_marker_edges),
    observed = list(nodes = obs_nodes, edges = obs_edges),
    edge_discordance = edge_disc,
    hub_concordance = hub_concordance,
    pathway_overlap = pathway_overlap
  ))
}
