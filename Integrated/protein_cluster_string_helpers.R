# STRING PPI networks for top up/down protein DE clusters (combined figure panel).
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

source("panel_predicted_observed_network_helpers.R")

DEFAULT_STRING_SCORE <- 400L
CLUSTER_PANEL_COLORS <- list(
  up = "#B2182B",
  down = "#2166AC",
  edge = "grey55",
  node_neutral = "grey88"
)

init_cluster_string_db <- function(
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

map_markers_to_string <- function(string_db, marker_df) {
  gene_tbl <- marker_df %>%
    dplyr::mutate(
      gene_symbol = vapply(.data$marker, marker_primary_gene, character(1))
    ) %>%
    dplyr::filter(!is.na(.data$gene_symbol), nzchar(.data$gene_symbol)) %>%
    dplyr::distinct(.data$marker, .data$gene_symbol)
  if (nrow(gene_tbl) == 0L) {
    return(list(gene_tbl = gene_tbl, mapped = data.frame()))
  }
  mapped <- string_db$map(
    as.data.frame(gene_tbl %>% dplyr::distinct(.data$gene_symbol)),
    "gene_symbol",
    removeUnmappedRows = TRUE
  )
  list(gene_tbl = gene_tbl, mapped = mapped)
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
    node_df$betweenness <- 0
    return(node_df)
  }
  g <- igraph::graph_from_data_frame(
    edges[, c("from_gene", "to_gene"), drop = FALSE],
    directed = FALSE,
    vertices = nodes
  )
  node_df$degree <- igraph::degree(g)
  node_df$betweenness <- igraph::betweenness(g, directed = FALSE)
  node_df
}

build_cluster_string_network <- function(marker_df, string_db) {
  maps <- map_markers_to_string(string_db, marker_df)
  gene_tbl <- maps$gene_tbl
  mapped <- maps$mapped
  if (nrow(mapped) < 1L) {
    return(list(
      marker_df = marker_df,
      gene_tbl = gene_tbl,
      mapped = mapped,
      edges = data.frame(),
      hubs = data.frame()
    ))
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
  hubs <- compute_gene_hubs(edges, nodes)
  list(
    marker_df = marker_df,
    gene_tbl = gene_tbl,
    mapped = mapped,
    edges = edges,
    hubs = hubs
  )
}

load_top_cluster_markers <- function(contrast_dir) {
  up_path <- file.path(contrast_dir, "top_upregulated_cluster_markers.csv")
  down_path <- file.path(contrast_dir, "top_downregulated_cluster_markers.csv")
  if (!file.exists(up_path) || !file.exists(down_path)) {
    stop("Missing cluster marker files under: ", contrast_dir)
  }
  list(
    up = readr::read_csv(up_path, show_col_types = FALSE),
    down = readr::read_csv(down_path, show_col_types = FALSE)
  )
}

trim_cluster_for_network <- function(marker_df, max_markers = 12L) {
  if (nrow(marker_df) <= max_markers) {
    return(list(df = marker_df, trimmed = FALSE, n_total = nrow(marker_df)))
  }
  out <- marker_df %>%
    dplyr::arrange(dplyr::desc(abs(.data$logFC))) %>%
    dplyr::slice_head(n = max_markers)
  list(df = out, trimmed = TRUE, n_total = nrow(marker_df))
}

plot_cluster_string_network_gg <- function(
  network,
  title,
  direction = c("up", "down"),
  logfc_limits = NULL,
  n_total_markers = NULL
) {
  direction <- match.arg(direction)
  accent <- if (direction == "up") {
    CLUSTER_PANEL_COLORS$up
  } else {
    CLUSTER_PANEL_COLORS$down
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph is required.")
  }

  marker_df <- network$marker_df
  gene_tbl <- network$gene_tbl
  edges <- network$edges
  gene_logfc <- gene_tbl %>%
    dplyr::left_join(
      marker_df %>% dplyr::select(.data$marker, .data$logFC),
      by = "marker"
    ) %>%
    dplyr::group_by(.data$gene_symbol) %>%
    dplyr::summarise(
      logFC = mean(.data$logFC, na.rm = TRUE),
      markers = paste(unique(.data$marker), collapse = "/"),
      .groups = "drop"
    )

  if (nrow(gene_logfc) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No mapped genes",
          size = 4
        ) +
        ggplot2::labs(title = title) +
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

  node_df <- data.frame(
    gene = nodes,
    x = lay[, 1],
    y = lay[, 2],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(gene_logfc, by = c("gene" = "gene_symbol")) %>%
    dplyr::left_join(network$hubs, by = c("gene" = "gene_symbol")) %>%
    dplyr::mutate(
      degree = dplyr::coalesce(.data$degree, 0L),
      label = dplyr::if_else(
        !is.na(.data$markers) & .data$markers != .data$gene,
        paste0(.data$gene, "\n(", .data$markers, ")"),
        .data$gene
      ),
      node_type = dplyr::case_when(
        .data$degree >= max(.data$degree, na.rm = TRUE) - 1L & .data$degree > 0L ~ "Hub",
        TRUE ~ "Member"
      )
    )

  if (is.null(logfc_limits)) {
    lim <- max(abs(node_df$logFC), na.rm = TRUE)
    if (!is.finite(lim) || lim < 0.01) lim <- 0.5
    logfc_limits <- c(-lim, lim)
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
        ggplot2::aes(
          x = x1,
          y = y1,
          xend = x2,
          yend = y2,
          linewidth = combined_score
        ),
        colour = CLUSTER_PANEL_COLORS$edge,
        alpha = 0.8
      )
  }

  p +
    ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(
        x = x,
        y = y,
        fill = logFC,
        size = degree + 2
      ),
      shape = 21,
      colour = "grey20",
      stroke = 0.4
    ) +
    ggplot2::geom_text(
      data = node_df,
      ggplot2::aes(x = x, y = y, label = label),
      size = 2.5,
      lineheight = 0.85,
      vjust = -1.15
    ) +
    ggplot2::scale_fill_gradient2(
      low = CLUSTER_PANEL_COLORS$down,
      mid = "white",
      high = CLUSTER_PANEL_COLORS$up,
      midpoint = 0,
      limits = logfc_limits,
      name = "logFC"
    ) +
    ggplot2::scale_size_continuous(range = c(4, 10), guide = "none") +
    ggplot2::scale_linewidth_continuous(range = c(0.4, 1.6), guide = "none") +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        if (!is.null(n_total_markers) && n_total_markers > length(nodes)) {
          paste0("Top ", length(nodes), " of ", n_total_markers, " cluster markers by |logFC|. ")
        } else {
          ""
        },
        if (nrow(edges) == 0L) {
          paste0(length(nodes), " gene(s); no STRING edges above threshold")
        } else {
          paste0(nrow(edges), " STRING edge(s), ", length(nodes), " gene(s)")
        }
      )
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 10.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8.5, colour = "grey35"),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.9, "cm")
    )
}

export_protein_cluster_string_panel <- function(
  cluster_root,
  compartment_label = NULL,
  output_dir = NULL,
  out_prefix = "protein_cluster_string_panel",
  string_score = DEFAULT_STRING_SCORE,
  max_network_markers = 12L,
  save_edge_tables = TRUE,
  string_db = NULL
) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required for the combined figure.")
  }
  high_low_dir <- file.path(cluster_root, "high_vs_low")
  g1_g4_dir <- file.path(cluster_root, "G1_vs_G4")
  if (!dir.exists(high_low_dir) || !dir.exists(g1_g4_dir)) {
    stop("Expected high_vs_low/ and G1_vs_G4/ under: ", cluster_root)
  }
  if (is.null(output_dir)) {
    output_dir <- file.path(cluster_root, "string_panels")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  hl <- load_top_cluster_markers(high_low_dir)
  g4 <- load_top_cluster_markers(g1_g4_dir)
  if (is.null(string_db)) {
    string_db <- init_cluster_string_db(score_threshold = string_score)
  }

  hl_up_trim <- trim_cluster_for_network(hl$up, max_network_markers)
  hl_down_trim <- trim_cluster_for_network(hl$down, max_network_markers)
  g4_up_trim <- trim_cluster_for_network(g4$up, max_network_markers)
  g4_down_trim <- trim_cluster_for_network(g4$down, max_network_markers)

  nets <- list(
    hl_up = build_cluster_string_network(hl_up_trim$df, string_db),
    hl_down = build_cluster_string_network(hl_down_trim$df, string_db),
    g4_up = build_cluster_string_network(g4_up_trim$df, string_db),
    g4_down = build_cluster_string_network(g4_down_trim$df, string_db)
  )

  if (save_edge_tables) {
    for (nm in names(nets)) {
      readr::write_csv(
        nets[[nm]]$edges,
        file.path(output_dir, paste0(nm, "_string_edges.csv"))
      )
    }
  }

  hl_up_id <- hl$up$cluster[1]
  hl_down_id <- hl$down$cluster[1]
  g4_up_id <- g4$up$cluster[1]
  g4_down_id <- g4$down$cluster[1]

  comp_tag <- if (!is.null(compartment_label)) paste0(compartment_label, " — ") else ""

  p_hl_up <- plot_cluster_string_network_gg(
    nets$hl_up,
    title = paste0(comp_tag, "High vs Low: UP cluster ", hl_up_id),
    direction = "up",
    n_total_markers = hl_up_trim$n_total
  )
  p_hl_down <- plot_cluster_string_network_gg(
    nets$hl_down,
    title = paste0(comp_tag, "High vs Low: DOWN cluster ", hl_down_id),
    direction = "down",
    n_total_markers = hl_down_trim$n_total
  )
  p_g4_up <- plot_cluster_string_network_gg(
    nets$g4_up,
    title = paste0(comp_tag, "G1 vs G4: UP cluster ", g4_up_id),
    direction = "up",
    n_total_markers = g4_up_trim$n_total
  )
  p_g4_down <- plot_cluster_string_network_gg(
    nets$g4_down,
    title = paste0(comp_tag, "G1 vs G4: DOWN cluster ", g4_down_id),
    direction = "down",
    n_total_markers = g4_down_trim$n_total
  )

  panel <- cowplot::plot_grid(
    p_hl_up,
    p_hl_down,
    p_g4_up,
    p_g4_down,
    nrow = 2,
    labels = c("A", "B", "C", "D"),
    label_size = 12
  )

  title_gg <- cowplot::ggdraw() +
    cowplot::draw_label(
      if (!is.null(compartment_label)) {
        paste0("STRING networks — top DE protein clusters (", compartment_label, ")")
      } else {
        "STRING networks — top DE protein clusters"
      },
      fontface = "bold",
      size = 14,
      x = 0.01,
      hjust = 0
    )

  fig <- cowplot::plot_grid(
    title_gg,
    panel,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )

  out_png <- file.path(output_dir, paste0(out_prefix, ".png"))
  out_pdf <- file.path(output_dir, paste0(out_prefix, ".pdf"))
  ggplot2::ggsave(out_png, fig, width = 14, height = 11, dpi = 200, bg = "white")
  ggplot2::ggsave(out_pdf, fig, width = 14, height = 11, bg = "white")
  message("Saved: ", out_png)

  summary_lines <- c(
    paste0("STRING panel: ", cluster_root),
    if (!is.null(compartment_label)) paste0("Compartment: ", compartment_label),
    paste0("High vs Low UP cluster ", hl_up_id, ": ", paste(hl$up$marker, collapse = ", "),
           " | edges=", nrow(nets$hl_up$edges)),
    paste0("High vs Low DOWN cluster ", hl_down_id, ": ", paste(hl$down$marker, collapse = ", "),
           " | edges=", nrow(nets$hl_down$edges)),
    paste0("G1 vs G4 UP cluster ", g4_up_id, ": ", paste(g4$up$marker, collapse = ", "),
           " | edges=", nrow(nets$g4_up$edges)),
    paste0("G1 vs G4 DOWN cluster ", g4_down_id, ": ", paste(g4$down$marker, collapse = ", "),
           " | edges=", nrow(nets$g4_down$edges)),
    paste0("Figure: ", out_png)
  )
  writeLines(summary_lines, file.path(output_dir, paste0(out_prefix, "_summary.txt")))
  cat(paste(summary_lines, collapse = "\n"), "\n")

  invisible(list(networks = nets, panel = fig, output_dir = output_dir))
}

export_all_compartment_cluster_string_panels <- function(
  epithelium_root = file.path(INTEGRATED_DIR, "results", "protein_de_clusters"),
  immune_root = file.path(INTEGRATED_DIR, "results", "protein_de_clusters_immune"),
  string_score = DEFAULT_STRING_SCORE,
  max_network_markers = 12L
) {
  results <- list()
  string_db <- init_cluster_string_db(score_threshold = string_score)

  results$epithelium <- export_protein_cluster_string_panel(
    cluster_root = epithelium_root,
    compartment_label = "Non-immune epithelium",
    out_prefix = "protein_cluster_string_panel_epi",
    max_network_markers = max_network_markers,
    string_db = string_db
  )

  if (dir.exists(immune_root)) {
    for (suffix in IMMUNE_COMPARTMENTS$suffix) {
      comp_root <- file.path(immune_root, suffix)
      if (!dir.exists(file.path(comp_root, "high_vs_low"))) next
      label <- compartment_label(suffix)
      results[[suffix]] <- export_protein_cluster_string_panel(
        cluster_root = comp_root,
        compartment_label = label,
        out_prefix = paste0("protein_cluster_string_panel_", suffix),
        max_network_markers = max_network_markers,
        string_db = string_db
      )
    }
  }

  invisible(results)
}
