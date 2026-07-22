# MC_13 (low risk) vs MC_8 (high risk) JAK-STAT investigation:
#   Option B — STRING PPI from spatial marker gene sets
#   Option C — Hallmark JAK-STAT regulon scores vs metacluster abundance
source("analysis_helpers.R")
source("progeny_helpers.R")
source("gsea_helpers.R")
source("metacluster_omics_correlation_helpers.R")

DEFAULT_MC_TARGETS <- c("MC_13", "MC_8")
DEFAULT_STRING_SCORE <- 400L
DEFAULT_TOP_MARKERS <- 12L

MARKER_MANUAL_GENE_MAP <- c(
  CD3 = "CD3D",
  CD8 = "CD8A",
  CD4 = "CD4",
  CD31 = "PECAM1",
  CD34 = "CD34",
  CD68 = "CD68",
  PD1 = "PDCD1",
  COLIV = "COL4A1",
  BCATENIN = "CTNNB1",
  CMYC = "MYC",
  KI67 = "MKI67",
  GRANZYMEB = "GZMB",
  HLA1 = "HLA-A",
  BCLXL = "BCL2L1",
  FLIP_CS = "CFLAR",
  CIAP1 = "BIRC2",
  S6 = "RPS6",
  SMA = "ACTA2",
  CASP3CLEAVED = "CASP3",
  PCK26 = "KRT26",
  SR2B = "HTR2B",
  MUC5 = "MUC5AC",
  AE1 = "SLC4A1",
  PCAD = "CDH17",
  NAK = "NAGK",
  GRP78 = "HSPA5",
  GLUT1 = "SLC2A1",
  BIM = "BCL2L11"
)

JAKSTAT_CORE_GENES <- c(
  "JAK1", "JAK2", "JAK3", "TYK2",
  "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
  "IL6", "IL6ST", "IL10", "IL10RB", "IFNGR1", "IFNGR2",
  "SOCS1", "SOCS2", "SOCS3", "PIAS1", "PIAS3",
  "MCL1", "BCL2", "BCL2L1", "CFLAR", "BIRC2"
)

HALLMARK_JAKSTAT_REGULONS <- c(
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_HYPOXIA",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_APOPTOSIS"
)

marker_to_gene <- function(marker) {
  marker <- as.character(marker)
  if (marker %in% names(MARKER_MANUAL_GENE_MAP)) {
    return(unname(MARKER_MANUAL_GENE_MAP[[marker]]))
  }
  toupper(marker)
}

load_mc_enriched_markers <- function(
  mc,
  mc_corr_dir,
  top_n = DEFAULT_TOP_MARKERS
) {
  path <- file.path(mc_corr_dir, paste0(mc, "_top_markers_enriched_vs_other.csv"))
  if (!file.exists(path)) {
    stop("Missing marker file: ", path)
  }
  readr::read_csv(path, show_col_types = FALSE) %>%
    dplyr::arrange(dplyr::desc(.data$log2fc)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      metacluster = mc,
      gene_symbol = vapply(.data$marker, marker_to_gene, character(1))
    )
}

build_marker_gene_table <- function(marker_df) {
  marker_df %>%
    dplyr::distinct(.data$marker, .data$gene_symbol, .data$metacluster) %>%
    dplyr::filter(!is.na(.data$gene_symbol), nzchar(.data$gene_symbol))
}

init_string_db_mc <- function(
  score_threshold = DEFAULT_STRING_SCORE,
  string_version = "12.0"
) {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("Package 'STRINGdb' is required.")
  }
  message("Initialising STRINGdb (first run may download the database)...")
  STRINGdb::STRINGdb$new(
    version = string_version,
    species = 9606L,
    score_threshold = score_threshold
  )
}

map_genes_to_string <- function(string_db, gene_tbl) {
  if (nrow(gene_tbl) == 0L) {
    return(data.frame())
  }
  string_db$map(
    as.data.frame(gene_tbl %>% dplyr::distinct(.data$gene_symbol)),
    "gene_symbol",
    removeUnmappedRows = TRUE
  )
}

string_interactions_to_edges <- function(interactions, mapped) {
  if (is.null(interactions) || nrow(interactions) == 0L || nrow(mapped) == 0L) {
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
  out %>%
    dplyr::mutate(
      gene_a = pmin(.data$from_gene, .data$to_gene),
      gene_b = pmax(.data$from_gene, .data$to_gene)
    ) %>%
    dplyr::group_by(.data$gene_a, .data$gene_b) %>%
    dplyr::slice_max(.data$combined_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      from_gene = .data$gene_a,
      to_gene = .data$gene_b,
      combined_score = .data$combined_score
    )
}

build_mc_string_network <- function(
  marker_df,
  string_db,
  jakstat_genes = JAKSTAT_CORE_GENES
) {
  gene_tbl <- build_marker_gene_table(marker_df)
  mapped <- map_genes_to_string(string_db, gene_tbl)
  if (nrow(mapped) < 2L) {
    return(list(
      marker_df = marker_df,
      gene_tbl = gene_tbl,
      mapped = mapped,
      edges = data.frame(),
      hubs = data.frame(),
      enrichment = data.frame(),
      jakstat_edges = data.frame()
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
  edges <- string_interactions_to_edges(interactions, mapped)
  if (nrow(edges) > 0L) {
    edges <- edges %>%
      dplyr::mutate(
        touches_jakstat = .data$from_gene %in% jakstat_genes |
          .data$to_gene %in% jakstat_genes
      )
  }
  jakstat_edges <- if (nrow(edges) > 0L) {
    edges %>% dplyr::filter(.data$touches_jakstat)
  } else {
    data.frame()
  }
  hubs <- compute_gene_hubs(edges, unique(gene_tbl$gene_symbol))
  enrichment <- tryCatch(
    string_db$get_enrichment(string_ids),
    error = function(e) {
      message("STRING enrichment failed: ", conditionMessage(e))
      data.frame()
    }
  )
  list(
    marker_df = marker_df,
    gene_tbl = gene_tbl,
    mapped = mapped,
    edges = edges,
    hubs = hubs,
    enrichment = enrichment,
    jakstat_edges = jakstat_edges
  )
}

compute_gene_hubs <- function(edges, nodes) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.")
  }
  if (length(nodes) == 0L) {
    return(data.frame())
  }
  node_df <- data.frame(gene_symbol = nodes, stringsAsFactors = FALSE)
  if (is.null(edges) || nrow(edges) == 0L) {
    node_df$degree <- 0L
    node_df$betweenness <- 0
    return(node_df %>% dplyr::arrange(dplyr::desc(.data$degree)))
  }
  g <- igraph::graph_from_data_frame(
    edges[, c("from_gene", "to_gene"), drop = FALSE],
    directed = FALSE,
    vertices = nodes
  )
  node_df$degree <- igraph::degree(g)
  node_df$betweenness <- igraph::betweenness(g, directed = FALSE)
  node_df %>%
    dplyr::arrange(dplyr::desc(.data$degree), dplyr::desc(.data$betweenness))
}

edge_key <- function(from_gene, to_gene) {
  paste(pmin(from_gene, to_gene), pmax(from_gene, to_gene), sep = "|")
}

compare_mc_string_edges <- function(edges_a, edges_b, label_a, label_b) {
  key_a <- if (nrow(edges_a) > 0L) edge_key(edges_a$from_gene, edges_a$to_gene) else character(0)
  key_b <- if (nrow(edges_b) > 0L) edge_key(edges_b$from_gene, edges_b$to_gene) else character(0)
  all_keys <- unique(c(key_a, key_b))
  if (length(all_keys) == 0L) {
    return(data.frame())
  }
  split_keys <- strsplit(all_keys, "|", fixed = TRUE)
  out <- data.frame(
    from_gene = vapply(split_keys, `[`, character(1), 1),
    to_gene = vapply(split_keys, `[`, character(1), 2),
    in_a = all_keys %in% key_a,
    in_b = all_keys %in% key_b,
    stringsAsFactors = FALSE
  )
  score_a <- if (nrow(edges_a) > 0L) {
    setNames(edges_a$combined_score, key_a)
  } else {
    setNames(numeric(0), character(0))
  }
  score_b <- if (nrow(edges_b) > 0L) {
    setNames(edges_b$combined_score, key_b)
  } else {
    setNames(numeric(0), character(0))
  }
  out[[paste0("score_", label_a)]] <- score_a[all_keys]
  out[[paste0("score_", label_b)]] <- score_b[all_keys]
  out$category <- dplyr::case_when(
    out$in_a & out$in_b ~ "shared",
    out$in_a & !out$in_b ~ paste0(label_a, "_only"),
    !out$in_a & out$in_b ~ paste0(label_b, "_only"),
    TRUE ~ "neither"
  )
  jak <- JAKSTAT_CORE_GENES
  out$touches_jakstat <- out$from_gene %in% jak | out$to_gene %in% jak
  out[order(out$category, out$touches_jakstat, decreasing = TRUE), , drop = FALSE]
}

export_string_network_png <- function(
  edges,
  mapped,
  out_path,
  title = NULL
) {
  if (nrow(edges) < 1L || nrow(mapped) < 2L) {
    message("Skipping network plot (<2 nodes or no edges): ", out_path)
    return(invisible(NULL))
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    warning("igraph not installed; skipping network plot.")
    return(invisible(NULL))
  }
  nodes <- unique(c(edges$from_gene, edges$to_gene))
  g <- igraph::graph_from_data_frame(
    edges[, c("from_gene", "to_gene"), drop = FALSE],
    directed = FALSE,
    vertices = nodes
  )
  jak <- JAKSTAT_CORE_GENES
  node_cols <- ifelse(nodes %in% jak, "#E41A1C", "#377EB8")
  edge_w <- if ("combined_score" %in% names(edges)) {
    scales::rescale(edges$combined_score, to = c(0.5, 3))
  } else {
    1
  }
  grDevices::png(out_path, width = 1400, height = 1000, res = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(
    g,
    vertex.color = node_cols,
    vertex.label.cex = 0.7,
    vertex.size = 12,
    edge.width = edge_w,
    main = title
  )
  legend(
    "bottomleft",
    legend = c("JAK-STAT core", "MC marker gene"),
    fill = c("#E41A1C", "#377EB8"),
    bty = "n",
    cex = 0.8
  )
  message("Saved: ", out_path)
  invisible(out_path)
}

hallmark_net_from_pathways <- function(pathway_names) {
  pathways <- hallmark_pathways()
  missing <- setdiff(pathway_names, names(pathways))
  if (length(missing) > 0L) {
    stop("Hallmark pathways not found: ", paste(missing, collapse = ", "))
  }
  dplyr::bind_rows(lapply(pathway_names, function(pw) {
    data.frame(
      source = pw,
      target = unique(pathways[[pw]]),
      weight = 1,
      stringsAsFactors = FALSE
    )
  })) %>%
    dplyr::distinct(.data$source, .data$target, .keep_all = TRUE)
}

run_hallmark_mlm <- function(expr_mat, pathway_names = HALLMARK_JAKSTAT_REGULONS) {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Package 'decoupleR' is required.")
  }
  net <- hallmark_net_from_pathways(pathway_names)
  genes <- intersect(unique(net$target), rownames(expr_mat))
  net <- net %>% dplyr::filter(.data$target %in% genes)
  decoupleR::run_mlm(
    mat = expr_mat,
    net = net,
    .source = "source",
    .target = "target",
    .mor = "weight",
    minsize = 5L
  )
}

hallmark_activities_to_long <- function(activities, sample_meta) {
  activities %>%
    dplyr::rename(
      sample_id = condition,
      regulon = source,
      score = score
    ) %>%
    dplyr::inner_join(
      sample_meta %>%
        dplyr::select(
          dplyr::any_of(c("sample_id", "patient_id", "risk_grp", "cohort"))
        ),
      by = "sample_id"
    )
}

run_mc_hallmark_regulon_correlation <- function(
  output_dir,
  target_mcs = DEFAULT_MC_TARGETS,
  regulons = HALLMARK_JAKSTAT_REGULONS,
  min_n = 10L,
  fisher_fdr = 0.05
) {
  p <- resolve_progeny_paths()
  survival_path <- file.path(p$proteomics, "survival_df_with_risk.rds")
  survival_df <- readRDS(survival_path)
  mc_df <- load_patient_metaclusters(survival_df)
  mc_cols <- grep(paste0("^", MC_PREFIX), names(mc_df), value = TRUE)
  target_cols <- mc_cols[vapply(mc_cols, function(col) mc_short_name(col) %in% target_mcs, logical(1))]

  cohort_loaders <- list(
    Stage2 = function() load_stage2_expression_meta(p, survival_df),
    Retrospective = function() load_retrospective_expression_meta(p, survival_df),
    Colossus = function() load_colossus_expression_meta(p, survival_df),
    Taxonomy = function() load_taxonomy_expression_meta(p, survival_path)
  )

  cohort_results <- list()
  for (cn in names(cohort_loaders)) {
    message("Hallmark regulon + MC correlations: ", cn)
    dat <- cohort_loaders[[cn]]()
    acts <- run_hallmark_mlm(dat$expr, regulons)
    acts_long <- hallmark_activities_to_long(acts, dat$meta)
    omics_wide <- aggregate_omics_to_patient(
      acts_long,
      patient_col = "patient_id",
      feature_col = "regulon",
      value_col = "score",
      sample_col = "sample_id"
    )
    if (ncol(omics_wide) == 0L) {
      warning("Skipping ", cn, " (no regulon scores).")
      next
    }
    combined <- join_mc_with_omics(mc_df, omics_wide)
    cor_df <- spearman_mc_omics(
      combined,
      mc_cols = intersect(target_cols, colnames(combined)),
      feature_cols = colnames(omics_wide),
      cohort = cn,
      omics_type = "Hallmark_regulon",
      min_n = min_n
    )
    if (nrow(cor_df) > 0L) {
      readr::write_csv(
        cor_df,
        file.path(output_dir, paste0(cn, "_mc_hallmark_regulon_spearman.csv"))
      )
      cohort_results[[cn]] <- cor_df
    }
  }

  if (length(cohort_results) < 2L) {
    stop("Need at least 2 cohorts with regulon correlations.")
  }
  meta <- fisher_meta_mc_omics(cohort_results)
  meta <- meta %>%
    dplyr::filter(.data$feature != "n_samples", is.finite(.data$fisher_p))
  meta$regulon_short <- sub("^HALLMARK_", "", meta$feature)
  readr::write_csv(meta, file.path(output_dir, "MC_hallmark_regulon_fisher_meta.csv"))

  sig <- meta %>%
    dplyr::filter(.data$fisher_padj <= fisher_fdr) %>%
    dplyr::filter(.data$sign_consistent %in% TRUE | .data$n_cohorts < 2L)
  readr::write_csv(sig, file.path(output_dir, "MC_hallmark_regulon_fisher_meta_sig.csv"))

  list(
    cohort = cohort_results,
    meta = meta,
    sig = sig
  )
}

export_regulon_heatmap <- function(meta_df, out_path, title) {
  if (nrow(meta_df) == 0L) {
    message("Skipping regulon heatmap (no rows): ", out_path)
    return(invisible(NULL))
  }
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    warning("pheatmap not installed; skipping heatmap.")
    return(invisible(NULL))
  }
  hm <- meta_df %>%
    dplyr::filter(.data$feature != "n_samples") %>%
    dplyr::transmute(
      metacluster = .data$metacluster,
      regulon = .data$regulon_short,
      mean_rho = .data$mean_rho
    ) %>%
    tidyr::pivot_wider(
      names_from = regulon,
      values_from = mean_rho
    ) %>%
    tibble::column_to_rownames("metacluster") %>%
    as.matrix()
  grDevices::png(out_path, width = 1400, height = 500, res = 120)
  pheatmap::pheatmap(
    hm,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    main = title,
    border_color = NA
  )
  grDevices::dev.off()
  message("Saved: ", out_path)
  invisible(out_path)
}

run_mc_jakstat_string_and_regulon <- function(
  output_dir = NULL,
  target_mcs = DEFAULT_MC_TARGETS,
  top_markers = DEFAULT_TOP_MARKERS,
  string_score = DEFAULT_STRING_SCORE
) {
  p <- resolve_progeny_paths()
  if (is.null(output_dir)) {
    output_dir <- file.path(
      p$base,
      "analysis_output",
      "metacluster_correlation",
      "jakstat_mc8_mc13"
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  string_dir <- file.path(output_dir, "string_networks")
  regulon_dir <- file.path(output_dir, "hallmark_regulons")
  dir.create(string_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(regulon_dir, recursive = TRUE, showWarnings = FALSE)

  mc_corr_dir <- file.path(p$base, "analysis_output", "metacluster_correlation")
  marker_tables <- stats::setNames(
    lapply(target_mcs, load_mc_enriched_markers, mc_corr_dir = mc_corr_dir, top_n = top_markers),
    target_mcs
  )
  combined_markers <- dplyr::bind_rows(marker_tables)
  readr::write_csv(combined_markers, file.path(string_dir, "mc_enriched_markers_input.csv"))

  string_db <- init_string_db_mc(score_threshold = string_score)
  networks <- lapply(target_mcs, function(mc) {
    build_mc_string_network(marker_tables[[mc]], string_db)
  })
  names(networks) <- target_mcs

  for (mc in target_mcs) {
    net <- networks[[mc]]
    prefix <- file.path(string_dir, mc)
    readr::write_csv(net$gene_tbl, paste0(prefix, "_marker_gene_map.csv"))
    readr::write_csv(net$mapped, paste0(prefix, "_string_mapped.csv"))
    readr::write_csv(net$edges, paste0(prefix, "_string_edges.csv"))
    readr::write_csv(net$jakstat_edges, paste0(prefix, "_jakstat_touching_edges.csv"))
    readr::write_csv(net$hubs, paste0(prefix, "_network_hubs.csv"))
    if (is.data.frame(net$enrichment) && nrow(net$enrichment) > 0L) {
      readr::write_csv(net$enrichment, paste0(prefix, "_string_enrichment.csv"))
    }
    export_string_network_png(
      net$edges,
      net$mapped,
      paste0(prefix, "_network.png"),
      title = paste0(mc, " marker STRING network")
    )
  }

  if (length(target_mcs) == 2L) {
    cmp <- compare_mc_string_edges(
      networks[[target_mcs[1]]]$edges,
      networks[[target_mcs[2]]]$edges,
      label_a = target_mcs[1],
      label_b = target_mcs[2]
    )
    readr::write_csv(cmp, file.path(string_dir, "MC_string_edge_comparison.csv"))
    jak_cmp <- cmp %>%
      dplyr::filter(.data$touches_jakstat)
    readr::write_csv(
      jak_cmp,
      file.path(string_dir, "MC_string_jakstat_edge_comparison.csv")
    )
  }

  regulon_out <- run_mc_hallmark_regulon_correlation(
    output_dir = regulon_dir,
    target_mcs = target_mcs
  )
  export_regulon_heatmap(
    regulon_out$meta,
    file.path(regulon_dir, "MC_hallmark_regulon_heatmap.png"),
    title = "MC vs Hallmark JAK-STAT regulons (Fisher meta rho)"
  )

  summary_lines <- c(
    "JAK-STAT MC_13 (low risk) vs MC_8 (high risk) — Options B + C",
    paste0("Output: ", output_dir),
    "",
    "=== Option B: STRING networks from enriched spatial markers ===",
    paste0("STRING score threshold: ", string_score),
    paste0("Top markers per MC: ", top_markers),
    ""
  )
  for (mc in target_mcs) {
    net <- networks[[mc]]
    top_hubs <- utils::head(net$hubs, 5)
    jak_n <- nrow(net$jakstat_edges)
    summary_lines <- c(
      summary_lines,
      paste0(mc, ":"),
      paste0("  markers: ", paste(marker_tables[[mc]]$marker, collapse = ", ")),
      paste0("  genes mapped to STRING: ", nrow(net$mapped)),
      paste0("  PPI edges: ", nrow(net$edges)),
      paste0("  JAK-STAT-touching edges: ", jak_n),
      if (nrow(top_hubs) > 0L) {
        paste0(
          "  top hubs: ",
          paste0(top_hubs$gene_symbol, " (deg=", top_hubs$degree, ")", collapse = ", ")
        )
      } else {
        "  top hubs: none"
      },
      ""
    )
  }
  if (length(target_mcs) == 2L && exists("cmp")) {
    shared_n <- sum(cmp$category == "shared", na.rm = TRUE)
    only_a <- sum(cmp$category == paste0(target_mcs[1], "_only"), na.rm = TRUE)
    only_b <- sum(cmp$category == paste0(target_mcs[2], "_only"), na.rm = TRUE)
    jak_shared <- sum(cmp$category == "shared" & cmp$touches_jakstat, na.rm = TRUE)
    summary_lines <- c(
      summary_lines,
      "Edge comparison:",
      paste0("  shared: ", shared_n, " (JAK-STAT-touching: ", jak_shared, ")"),
      paste0("  ", target_mcs[1], "-only: ", only_a),
      paste0("  ", target_mcs[2], "-only: ", only_b),
      ""
    )
  }

  summary_lines <- c(
    summary_lines,
    "=== Option C: Hallmark regulon vs MC abundance (Fisher meta) ===",
    paste0("Regulons: ", paste(HALLMARK_JAKSTAT_REGULONS, collapse = ", ")),
    ""
  )
  meta <- regulon_out$meta %>%
    dplyr::arrange(.data$metacluster, .data$fisher_padj)
  for (mc in target_mcs) {
    sub <- meta %>% dplyr::filter(.data$metacluster == mc)
    summary_lines <- c(summary_lines, paste0(mc, " top regulon links:"))
    if (nrow(sub) == 0L) {
      summary_lines <- c(summary_lines, "  (none tested)")
    } else {
      top <- sub %>% dplyr::arrange(.data$fisher_padj) %>% utils::head(5)
      summary_lines <- c(
        summary_lines,
        paste0(
          "  ",
          top$regulon_short,
          ": rho=",
          signif(top$mean_rho, 3),
          ", padj=",
          signif(top$fisher_padj, 3)
        )
      )
    }
    summary_lines <- c(summary_lines, "")
  }

  writeLines(summary_lines, file.path(output_dir, "summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(
    networks = networks,
    edge_comparison = if (exists("cmp")) cmp else NULL,
    regulon = regulon_out,
    output_dir = output_dir
  ))
}

MC_JAKSTAT_PANEL_COLORS <- list(
  MC_13 = "#2166AC",
  MC_8 = "#B2182B",
  jakstat = "#E41A1C",
  edge_other = "#BDBDBD",
  edge_jakstat = "#D95F02"
)

resolve_jakstat_panel_dir <- function() {
  p <- resolve_progeny_paths()
  file.path(
    p$base,
    "analysis_output",
    "metacluster_correlation",
    "jakstat_mc8_mc13"
  )
}

load_jakstat_panel_data <- function(panel_dir = NULL) {
  if (is.null(panel_dir)) {
    panel_dir <- resolve_jakstat_panel_dir()
  }
  string_dir <- file.path(panel_dir, "string_networks")
  regulon_dir <- file.path(panel_dir, "hallmark_regulons")
  list(
    panel_dir = panel_dir,
    regulon_meta = readr::read_csv(
      file.path(regulon_dir, "MC_hallmark_regulon_fisher_meta.csv"),
      show_col_types = FALSE
    ),
    edge_cmp = readr::read_csv(
      file.path(string_dir, "MC_string_edge_comparison.csv"),
      show_col_types = FALSE
    ),
    MC_13 = list(
      edges = readr::read_csv(
        file.path(string_dir, "MC_13_string_edges.csv"),
        show_col_types = FALSE
      ),
      hubs = readr::read_csv(
        file.path(string_dir, "MC_13_network_hubs.csv"),
        show_col_types = FALSE
      )
    ),
    MC_8 = list(
      edges = readr::read_csv(
        file.path(string_dir, "MC_8_string_edges.csv"),
        show_col_types = FALSE
      ),
      hubs = readr::read_csv(
        file.path(string_dir, "MC_8_network_hubs.csv"),
        show_col_types = FALSE
      )
    )
  )
}

regulon_label <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x, fixed = TRUE)
  tools::toTitleCase(tolower(x))
}

sig_stars <- function(padj) {
  dplyr::case_when(
    !is.finite(padj) ~ "",
    padj < 0.001 ~ "***",
    padj < 0.01 ~ "**",
    padj < 0.05 ~ "*",
    TRUE ~ ""
  )
}

plot_regulon_mc_comparison <- function(regulon_meta) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }
  plot_df <- regulon_meta %>%
    dplyr::filter(.data$feature != "n_samples", is.finite(.data$mean_rho)) %>%
    dplyr::mutate(
      regulon_label = regulon_label(.data$feature),
      stars = sig_stars(.data$fisher_padj),
      metacluster = factor(
        .data$metacluster,
        levels = c("MC_13", "MC_8"),
        labels = c("MC_13 (Low risk)", "MC_8 (High risk)")
      )
    )
  order_levels <- plot_df %>%
    dplyr::group_by(.data$regulon_label) %>%
    dplyr::summarise(avg = mean(.data$mean_rho), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$avg)) %>%
    dplyr::pull(.data$regulon_label)
  plot_df$regulon_label <- factor(plot_df$regulon_label, levels = rev(order_levels))

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = mean_rho,
      y = regulon_label,
      colour = metacluster,
      shape = metacluster
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4) +
    ggplot2::geom_point(size = 3.2, alpha = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = stars),
      hjust = -0.4,
      size = 3.5,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "MC_13 (Low risk)" = MC_JAKSTAT_PANEL_COLORS$MC_13,
        "MC_8 (High risk)" = MC_JAKSTAT_PANEL_COLORS$MC_8
      ),
      name = NULL
    ) +
    ggplot2::scale_shape_manual(
      values = c("MC_13 (Low risk)" = 16, "MC_8 (High risk)" = 17),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Hallmark regulon linkage",
      subtitle = "Fisher meta Spearman rho vs metacluster abundance",
      x = "Mean Spearman rho (4 cohorts)",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 12)
    )
}

plot_edge_category_summary <- function(edge_cmp) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }
  plot_df <- edge_cmp %>%
    dplyr::mutate(
      category = dplyr::case_when(
        .data$category == "MC_13_only" & .data$touches_jakstat ~ "MC_13 only (JAK-STAT)",
        .data$category == "MC_13_only" ~ "MC_13 only (other)",
        .data$category == "MC_8_only" ~ "MC_8 only",
        TRUE ~ "Shared"
      ),
      category = factor(
        .data$category,
        levels = c(
          "MC_13 only (JAK-STAT)",
          "MC_13 only (other)",
          "MC_8 only",
          "Shared"
        )
      )
    ) %>%
    dplyr::count(.data$category, name = "n_edges")

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = category, y = n_edges, fill = category)
  ) +
    ggplot2::geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = n_edges),
      vjust = -0.3,
      size = 3.5
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "MC_13 only (JAK-STAT)" = MC_JAKSTAT_PANEL_COLORS$edge_jakstat,
        "MC_13 only (other)" = MC_JAKSTAT_PANEL_COLORS$MC_13,
        "MC_8 only" = MC_JAKSTAT_PANEL_COLORS$MC_8,
        "Shared" = "grey70"
      ),
      guide = "none"
    ) +
    ggplot2::labs(
      title = "STRING PPI edge overlap",
      subtitle = "Predicted interactions among enriched spatial markers",
      x = NULL,
      y = "Edge count"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 12)
    )
}

plot_mc_string_network_gg <- function(
  edges,
  hubs,
  title,
  accent_color,
  jakstat_genes = JAKSTAT_CORE_GENES
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph is required.")
  }
  if (nrow(edges) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No STRING edges",
          size = 4
        ) +
        ggplot2::labs(title = title) +
        ggplot2::theme_void()
    )
  }

  nodes <- unique(c(edges$from_gene, edges$to_gene))
  g <- igraph::graph_from_data_frame(
    edges[, c("from_gene", "to_gene"), drop = FALSE],
    directed = FALSE,
    vertices = nodes
  )
  lay <- igraph::layout_with_fr(g)
  node_df <- data.frame(
    gene = nodes,
    x = lay[, 1],
    y = lay[, 2],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(hubs, by = c("gene" = "gene_symbol")) %>%
    dplyr::mutate(
      degree = dplyr::coalesce(.data$degree, 0L),
      node_type = dplyr::case_when(
        .data$gene %in% jakstat_genes ~ "JAK-STAT core",
        .data$degree >= max(.data$degree, na.rm = TRUE) - 1L & .data$degree > 0L ~ "Hub",
        TRUE ~ "Other"
      )
    )

  edge_df <- edges %>%
    dplyr::left_join(node_df, by = c("from_gene" = "gene")) %>%
    dplyr::rename(x1 = x, y1 = y) %>%
    dplyr::left_join(node_df, by = c("to_gene" = "gene")) %>%
    dplyr::rename(x2 = x, y2 = y) %>%
    dplyr::mutate(
      edge_type = if ("touches_jakstat" %in% names(.)) {
        dplyr::if_else(.data$touches_jakstat, "JAK-STAT", "Other")
      } else {
        "Other"
      }
    )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(
        x = x1,
        y = y1,
        xend = x2,
        yend = y2,
        colour = edge_type,
        linewidth = combined_score
      ),
      alpha = 0.75
    ) +
    ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(x = x, y = y, fill = node_type, size = degree + 1),
      shape = 21,
      colour = "grey20",
      stroke = 0.35
    ) +
    ggplot2::geom_text(
      data = node_df,
      ggplot2::aes(x = x, y = y, label = gene),
      size = 2.6,
      vjust = -1.1
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "JAK-STAT" = MC_JAKSTAT_PANEL_COLORS$edge_jakstat,
        "Other" = MC_JAKSTAT_PANEL_COLORS$edge_other
      ),
      name = "Edge"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "JAK-STAT core" = MC_JAKSTAT_PANEL_COLORS$jakstat,
        "Hub" = accent_color,
        "Other" = "grey92"
      ),
      name = "Node"
    ) +
    ggplot2::scale_linewidth_continuous(range = c(0.3, 1.4), guide = "none") +
    ggplot2::scale_size_continuous(range = c(3, 8), guide = "none") +
    ggplot2::labs(title = title) +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
}

export_mc_jakstat_figure_panel <- function(
  panel_dir = NULL,
  out_prefix = "MC_jakstat_panel"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required for the figure panel.")
  }

  data <- load_jakstat_panel_data(panel_dir)
  fig_dir <- file.path(data$panel_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  p_regulon <- plot_regulon_mc_comparison(data$regulon_meta)
  p_edges <- plot_edge_category_summary(data$edge_cmp)
  p_net13 <- plot_mc_string_network_gg(
    data$MC_13$edges,
    data$MC_13$hubs,
    title = "MC_13 - Low risk (immune)",
    accent_color = MC_JAKSTAT_PANEL_COLORS$MC_13
  )
  p_net8 <- plot_mc_string_network_gg(
    data$MC_8$edges,
    data$MC_8$hubs,
    title = "MC_8 - High risk (apoptosis)",
    accent_color = MC_JAKSTAT_PANEL_COLORS$MC_8
  )

  panel_top <- cowplot::plot_grid(
    p_regulon + ggplot2::theme(plot.margin = ggplot2::margin(5, 10, 5, 5)),
    p_edges + ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 10)),
    nrow = 1,
    rel_widths = c(1.35, 0.85),
    labels = c("A", "B")
  )
  panel_bottom <- cowplot::plot_grid(
    p_net13,
    p_net8,
    nrow = 1,
    rel_widths = c(1, 1),
    labels = c("C", "D")
  )
  panel <- cowplot::plot_grid(
    panel_top,
    panel_bottom,
    nrow = 2,
    rel_heights = c(1, 1.1)
  )

  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      "JAK-STAT links distinct spatial niches in low- vs high-risk metaclusters",
      fontface = "bold",
      size = 14,
      x = 0.01,
      hjust = 0
    )
  panel_final <- cowplot::plot_grid(title, panel, ncol = 1, rel_heights = c(0.05, 1))

  out_png <- file.path(fig_dir, paste0(out_prefix, ".png"))
  out_pdf <- file.path(fig_dir, paste0(out_prefix, ".pdf"))
  ggplot2::ggsave(out_png, panel_final, width = 14, height = 10, dpi = 300, bg = "white")
  ggplot2::ggsave(out_pdf, panel_final, width = 14, height = 10, bg = "white")

  index_lines <- c(
    "MC_13 vs MC_8 JAK-STAT figure panel",
    paste0("Generated: ", Sys.time()),
    "",
    paste0(out_prefix, ".png"),
    paste0(out_prefix, ".pdf"),
    "",
    "Panels:",
    "  A — Hallmark regulon Fisher meta rho (MC_13 vs MC_8)",
    "  B — STRING edge overlap (marker-derived PPIs)",
    "  C — MC_13 STRING network (immune / BCL2 hub)",
    "  D — MC_8 STRING network (apoptosis / CASP3 hub)"
  )
  writeLines(index_lines, file.path(fig_dir, "figures_index.txt"))
  message("Saved: ", out_png)
  message("Saved: ", out_pdf)
  invisible(list(png = out_png, pdf = out_pdf, fig_dir = fig_dir))
}
