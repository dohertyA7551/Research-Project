# STRING PPI networks from differentially expressed proteins (all compartments).
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

source("epi_protein_cluster_helpers.R")

DEFAULT_STRING_SCORE <- 400L
DEFAULT_P_CUTOFF <- 0.1
DEFAULT_FC_CUTOFF <- 0
DEFAULT_MAX_DE_MARKERS <- 60L

ALL_PROTEIN_COMPARTMENTS <- data.frame(
  suffix = c(
    "NonImmuneEpithelium",
    "NonImmuneStroma",
    "helperT",
    "cytT",
    "regT"
  ),
  label = c(
    "Non-immune epithelium (tumour)",
    "Stromal cells",
    "Helper T cells",
    "Cytotoxic T cells",
    "Regulatory T cells"
  ),
  stringsAsFactors = FALSE
)

compartment_display_label <- function(suffix) {
  idx <- match(suffix, ALL_PROTEIN_COMPARTMENTS$suffix)
  if (is.na(idx)) suffix else ALL_PROTEIN_COMPARTMENTS$label[idx]
}

DE_NETWORK_COLORS <- list(
  up = "#B2182B",
  down = "#2166AC",
  edge = "grey55"
)

init_de_string_db <- function(
  score_threshold = DEFAULT_STRING_SCORE,
  string_version = "12.0"
) {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("Package 'STRINGdb' is required.")
  }
  message("Initialising STRINGdb...")
  STRINGdb::STRINGdb$new(
    version = string_version,
    species = 9606L,
    score_threshold = score_threshold
  )
}

select_de_markers_for_string <- function(
  de_df,
  p_cutoff = DEFAULT_P_CUTOFF,
  fc_cutoff = DEFAULT_FC_CUTOFF,
  max_markers = DEFAULT_MAX_DE_MARKERS
) {
  id_col <- if ("marker" %in% names(de_df)) "marker" else "feature"

  df <- de_df %>%
    dplyr::filter(
      is.finite(.data$logFC),
      is.finite(.data$P.Value),
      !is.na(.data[[id_col]]),
      nzchar(.data[[id_col]]),
      .data$P.Value <= p_cutoff
    )
  if (fc_cutoff > 0) {
    df <- df %>% dplyr::filter(abs(.data$logFC) >= fc_cutoff)
  }
  df <- df %>%
    dplyr::arrange(.data$P.Value, dplyr::desc(abs(.data$logFC))) %>%
    dplyr::slice_head(n = max_markers) %>%
    dplyr::mutate(gene_symbol = vapply(.data[[id_col]], marker_primary_gene, character(1)))

  list(
    markers = df,
    n_input = nrow(df),
    p_mode = paste0("P<=", p_cutoff, if (fc_cutoff > 0) paste0(", |logFC|>=", fc_cutoff) else ""),
    id_col = id_col
  )
}

string_interactions_to_gene_edges <- function(interactions, mapped) {
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

compute_gene_hubs <- function(edges, nodes) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.")
  }
  node_df <- data.frame(gene_symbol = nodes, stringsAsFactors = FALSE)
  if (length(nodes) == 0L) {
    return(node_df)
  }
  if (is.null(edges) || nrow(edges) == 0L) {
    node_df$degree <- 0L
    return(node_df)
  }
  g <- igraph::graph_from_data_frame(
    edges[, c("from_gene", "to_gene"), drop = FALSE],
    directed = FALSE,
    vertices = nodes
  )
  node_df$degree <- igraph::degree(g)
  node_df
}

build_de_string_network <- function(de_markers, string_db) {
  gene_tbl <- de_markers %>%
    dplyr::distinct(.data$marker, .data$gene_symbol, .data$logFC, .data$P.Value, .data$adj.P.Val)
  mapped <- if (nrow(gene_tbl) > 0L) {
    string_db$map(
      as.data.frame(gene_tbl %>% dplyr::distinct(.data$gene_symbol)),
      "gene_symbol",
      removeUnmappedRows = TRUE
    )
  } else {
    data.frame()
  }
  if (nrow(mapped) < 1L) {
    return(list(de_markers = de_markers, mapped = mapped, edges = data.frame(), hubs = data.frame()))
  }
  string_ids <- unique(mapped$STRING_id)
  interactions <- if (length(string_ids) >= 2L) {
    tryCatch(
      string_db$get_interactions(string_ids),
      error = function(e) {
        message("STRING get_interactions failed: ", conditionMessage(e))
        data.frame()
      }
    )
  } else {
    data.frame()
  }
  edges <- string_interactions_to_gene_edges(interactions, mapped)
  nodes <- unique(gene_tbl$gene_symbol[gene_tbl$gene_symbol %in% mapped$gene_symbol])
  list(
    de_markers = de_markers,
    mapped = mapped,
    edges = edges,
    hubs = compute_gene_hubs(edges, nodes)
  )
}

run_compartment_de_high_vs_low <- function(compartment_suffix) {
  prot <- load_epi_protein_matrix(compartment_suffix = compartment_suffix)
  aligned <- subset_protein_high_low(prot)
  de <- run_limma_high_vs_low_matrix(
    aligned$matrix,
    aligned$meta,
    feature_col = "marker"
  )
  de$results %>%
    dplyr::mutate(compartment = compartment_suffix)
}

run_compartment_de_g1_vs_g4 <- function(
  compartment_suffix,
  groups_rds = DEFAULT_TMA_GROUPS_RDS,
  group_low = G1_GROUP,
  group_high = G4_GROUP
) {
  prot <- load_epi_protein_with_tma_groups(
    groups_rds = groups_rds,
    compartment_suffix = compartment_suffix
  )
  aligned <- subset_protein_two_groups(
    prot,
    group_col = "tma_comp_grp",
    group_low = group_low,
    group_high = group_high,
    min_patients = 10L
  )
  de <- run_limma_two_group_matrix(
    aligned$matrix,
    aligned$meta,
    group_col = "tma_comp_grp",
    group_low = group_low,
    group_high = group_high,
    feature_col = "marker"
  )
  de$results %>%
    dplyr::mutate(compartment = compartment_suffix)
}

plot_de_string_network_gg <- function(
  network,
  title,
  subtitle = NULL,
  logfc_limits = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph is required.")
  }

  de_markers <- network$de_markers
  edges <- network$edges
  gene_logfc <- de_markers %>%
    dplyr::group_by(.data$gene_symbol) %>%
    dplyr::summarise(
      logFC = mean(.data$logFC, na.rm = TRUE),
      markers = paste(unique(.data$marker), collapse = "/"),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$gene_symbol %in% network$mapped$gene_symbol)

  if (nrow(gene_logfc) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No mapped DE proteins", size = 4) +
        ggplot2::labs(title = title, subtitle = subtitle) +
        ggplot2::theme_void()
    )
  }

  nodes <- gene_logfc$gene_symbol
  if (nrow(edges) == 0L) {
    lay <- cbind(
      seq(0, 1, length.out = length(nodes)),
      rep(0.5, length(nodes))
    )
  } else {
    g <- igraph::graph_from_data_frame(
      edges[, c("from_gene", "to_gene"), drop = FALSE],
      directed = FALSE,
      vertices = nodes
    )
    lay <- igraph::layout_with_fr(g)
  }

  node_df <- data.frame(gene = nodes, x = lay[, 1], y = lay[, 2], stringsAsFactors = FALSE) %>%
    dplyr::left_join(gene_logfc, by = c("gene" = "gene_symbol")) %>%
    dplyr::left_join(network$hubs, by = c("gene" = "gene_symbol")) %>%
    dplyr::mutate(
      degree = dplyr::coalesce(.data$degree, 0L),
      label = dplyr::if_else(
        !is.na(.data$markers) & .data$markers != .data$gene,
        paste0(.data$gene, "\n(", .data$markers, ")"),
        .data$gene
      )
    )

  if (is.null(logfc_limits)) {
    lim <- max(abs(node_df$logFC), na.rm = TRUE)
    if (!is.finite(lim) || lim < 0.01) lim <- 0.5
    logfc_limits <- c(-lim, lim)
  }

  if (is.null(subtitle)) {
    subtitle <- paste0(
      nrow(de_markers), " DE protein(s), ",
      length(nodes), " mapped gene(s), ",
      nrow(edges), " STRING edge(s)"
    )
  }

  p <- ggplot2::ggplot()
  if (nrow(edges) > 0L) {
    edge_df <- edges %>%
      dplyr::left_join(node_df, by = c("from_gene" = "gene")) %>%
      dplyr::rename(x1 = x, y1 = y) %>%
      dplyr::left_join(node_df, by = c("to_gene" = "gene")) %>%
      dplyr::rename(x2 = x, y2 = y)
    p <- p +
      ggplot2::geom_segment(
        data = edge_df,
        ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = combined_score),
        colour = DE_NETWORK_COLORS$edge,
        alpha = 0.8
      )
  }

  p +
    ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(x = x, y = y, fill = logFC, size = degree + 2),
      shape = 21,
      colour = "grey20",
      stroke = 0.35
    ) +
    ggplot2::geom_text(
      data = node_df,
      ggplot2::aes(x = x, y = y, label = label),
      size = 2.4,
      lineheight = 0.85,
      vjust = -1.1
    ) +
    ggplot2::scale_fill_gradient2(
      low = DE_NETWORK_COLORS$down,
      mid = "white",
      high = DE_NETWORK_COLORS$up,
      midpoint = 0,
      limits = logfc_limits,
      name = "logFC"
    ) +
    ggplot2::scale_size_continuous(range = c(3.5, 9), guide = "none") +
    ggplot2::scale_linewidth_continuous(range = c(0.35, 1.5), guide = "none") +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 10),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 7.5, colour = "grey35"),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.8, "cm")
    )
}

export_compartment_de_string_networks <- function(
  compartment_suffix,
  output_dir = NULL,
  string_db = NULL,
  p_cutoff = DEFAULT_P_CUTOFF,
  fc_cutoff = DEFAULT_FC_CUTOFF,
  string_score = DEFAULT_STRING_SCORE
) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required.")
  }
  label <- compartment_display_label(compartment_suffix)
  if (is.null(output_dir)) {
    output_dir <- file.path(
      INTEGRATED_DIR,
      "results",
      "protein_de_string_networks",
      compartment_suffix
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  de_hl <- run_compartment_de_high_vs_low(compartment_suffix)
  de_g4 <- run_compartment_de_g1_vs_g4(compartment_suffix)
  readr::write_csv(de_hl, file.path(output_dir, "limma_high_vs_low.csv"))
  readr::write_csv(de_g4, file.path(output_dir, "limma_G1_vs_G4.csv"))

  sel_hl <- select_de_markers_for_string(de_hl, p_cutoff = p_cutoff, fc_cutoff = fc_cutoff)
  sel_g4 <- select_de_markers_for_string(de_g4, p_cutoff = p_cutoff, fc_cutoff = fc_cutoff)
  readr::write_csv(
    sel_hl$markers,
    file.path(output_dir, "de_markers_high_vs_low_for_string.csv")
  )
  readr::write_csv(
    sel_g4$markers,
    file.path(output_dir, "de_markers_G1_vs_G4_for_string.csv")
  )

  if (is.null(string_db)) {
    string_db <- init_de_string_db(score_threshold = string_score)
  }
  net_hl <- build_de_string_network(sel_hl$markers, string_db)
  net_g4 <- build_de_string_network(sel_g4$markers, string_db)
  readr::write_csv(net_hl$edges, file.path(output_dir, "string_edges_high_vs_low.csv"))
  readr::write_csv(net_g4$edges, file.path(output_dir, "string_edges_G1_vs_G4.csv"))

  p_hl <- plot_de_string_network_gg(
    net_hl,
    title = "High vs Low risk",
    subtitle = paste0(
      sel_hl$n_input, " DE proteins (", sel_hl$p_mode, "); ",
      nrow(net_hl$edges), " STRING edges"
    )
  )
  p_g4 <- plot_de_string_network_gg(
    net_g4,
    title = "G1_all_low vs G4_high_3plus",
    subtitle = paste0(
      sel_g4$n_input, " DE proteins (", sel_g4$p_mode, "); ",
      nrow(net_g4$edges), " STRING edges"
    )
  )

  panel <- cowplot::plot_grid(
    p_hl,
    p_g4,
    nrow = 1,
    labels = c("A", "B"),
    label_size = 11
  )
  title_gg <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0("DE protein STRING network — ", label),
      fontface = "bold",
      size = 13,
      x = 0.01,
      hjust = 0
    )
  fig <- cowplot::plot_grid(title_gg, panel, ncol = 1, rel_heights = c(0.06, 1))

  out_png <- file.path(output_dir, "de_string_network_panel.png")
  out_pdf <- file.path(output_dir, "de_string_network_panel.pdf")
  ggplot2::ggsave(out_png, fig, width = 13, height = 6.5, dpi = 200, bg = "white")
  ggplot2::ggsave(out_pdf, fig, width = 13, height = 6.5, bg = "white")
  message("Saved: ", out_png)

  summary_lines <- c(
    paste0("DE STRING network: ", label, " (", compartment_suffix, ")"),
    paste0(
      "High vs Low: ", sel_hl$n_input, " DE proteins (", sel_hl$p_mode, "), ",
      nrow(net_hl$edges), " edges | ",
      paste(sel_hl$markers$marker, collapse = ", ")
    ),
    paste0(
      "G1 vs G4: ", sel_g4$n_input, " DE proteins (", sel_g4$p_mode, "), ",
      nrow(net_g4$edges), " edges | ",
      paste(sel_g4$markers$marker, collapse = ", ")
    ),
    paste0("Figure: ", out_png)
  )
  writeLines(summary_lines, file.path(output_dir, "summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(
    de_hl = de_hl,
    de_g4 = de_g4,
    sel_hl = sel_hl,
    sel_g4 = sel_g4,
    net_hl = net_hl,
    net_g4 = net_g4,
    panel = fig
  ))
}

export_all_compartment_de_string_networks <- function(
  output_root = NULL,
  compartments = ALL_PROTEIN_COMPARTMENTS,
  p_cutoff = DEFAULT_P_CUTOFF,
  fc_cutoff = DEFAULT_FC_CUTOFF,
  string_score = DEFAULT_STRING_SCORE
) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required.")
  }
  if (is.null(output_root)) {
    output_root <- file.path(INTEGRATED_DIR, "results", "protein_de_string_networks")
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  string_db <- init_de_string_db(score_threshold = string_score)

  results <- list()
  panels_hl <- list()
  panels_g4 <- list()
  summary_rows <- list()

  for (i in seq_len(nrow(compartments))) {
    suffix <- compartments$suffix[i]
    label <- compartments$label[i]
    message("\n=== ", label, " ===")
    res <- export_compartment_de_string_networks(
      compartment_suffix = suffix,
      output_dir = file.path(output_root, suffix),
      string_db = string_db,
      p_cutoff = p_cutoff,
      fc_cutoff = fc_cutoff,
      string_score = string_score
    )
    results[[suffix]] <- res
    panels_hl[[suffix]] <- plot_de_string_network_gg(
      res$net_hl,
      title = label,
      subtitle = paste0("High vs Low (n=", res$sel_hl$n_input, ", edges=", nrow(res$net_hl$edges), ")")
    )
    panels_g4[[suffix]] <- plot_de_string_network_gg(
      res$net_g4,
      title = label,
      subtitle = paste0("G1 vs G4 (n=", res$sel_g4$n_input, ", edges=", nrow(res$net_g4$edges), ")")
    )
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      compartment = suffix,
      compartment_label = label,
      contrast = "high_vs_low",
      n_de = res$sel_hl$n_input,
      p_mode = res$sel_hl$p_mode,
      n_edges = nrow(res$net_hl$edges),
      de_markers = paste(res$sel_hl$markers$marker, collapse = ";"),
      stringsAsFactors = FALSE
    )
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      compartment = suffix,
      compartment_label = label,
      contrast = "G1_vs_G4",
      n_de = res$sel_g4$n_input,
      p_mode = res$sel_g4$p_mode,
      n_edges = nrow(res$net_g4$edges),
      de_markers = paste(res$sel_g4$markers$marker, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- dplyr::bind_rows(summary_rows)
  readr::write_csv(summary_df, file.path(output_root, "all_compartments_summary.csv"))

  col_hl <- cowplot::plot_grid(plotlist = panels_hl, ncol = 1)
  col_g4 <- cowplot::plot_grid(plotlist = panels_g4, ncol = 1)
  master <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_label("High vs Low risk", fontface = "bold", size = 14),
    col_hl,
    cowplot::ggdraw() +
      cowplot::draw_label("G1_all_low vs G4_high_3plus", fontface = "bold", size = 14),
    col_g4,
    nrow = 1,
    rel_widths = c(0.02, 1, 0.02, 1)
  )
  master_title <- cowplot::plot_grid(
    cowplot::ggdraw() +
      cowplot::draw_label(
        "DE protein STRING networks — all cell compartments",
        fontface = "bold",
        size = 15,
        x = 0.01,
        hjust = 0
      ),
    master,
    ncol = 1,
    rel_heights = c(0.04, 1)
  )
  master_png <- file.path(output_root, "de_string_network_all_compartments.png")
  master_pdf <- file.path(output_root, "de_string_network_all_compartments.pdf")
  ggplot2::ggsave(master_png, master_title, width = 16, height = 22, dpi = 180, bg = "white")
  ggplot2::ggsave(master_pdf, master_title, width = 16, height = 22, bg = "white")
  message("Saved master figure: ", master_png)

  epi_dir <- file.path(output_root, "NonImmuneEpithelium")
  epi_panel <- file.path(epi_dir, "de_string_network_panel.png")
  if (file.exists(epi_panel)) {
    file.copy(
      epi_panel,
      file.path(output_root, "cancer_epithelium_de_string_network.png"),
      overwrite = TRUE
    )
    message("Saved tumour epithelium highlight: ", file.path(output_root, "cancer_epithelium_de_string_network.png"))
  }

  master_lines <- c(
    "DE protein STRING networks — all compartments",
    paste0(
      "DE selection: nominal P<=", p_cutoff,
      " (no FDR; ", DEFAULT_MAX_DE_MARKERS, "-marker CODEX panel)"
    ),
    "Non-immune epithelium = tumour/cancer compartment (listed first).",
    ""
  )
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    master_lines <- c(
      master_lines,
      paste0(
        row$compartment_label, " [", row$contrast, "]: ",
        row$n_de, " DE (", row$p_mode, "), ", row$n_edges, " edges"
      )
    )
  }
  master_lines <- c(master_lines, "", paste0("Master figure: ", master_png))
  writeLines(master_lines, file.path(output_root, "master_summary.txt"))
  cat("\n", paste(master_lines, collapse = "\n"), "\n", sep = "")

  invisible(list(by_compartment = results, summary = summary_df, master = master_title))
}
