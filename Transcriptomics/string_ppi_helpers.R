# STRINGdb protein-protein interaction analysis from limma High vs Low DE genes.
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}

DEFAULT_STRING_COHORTS <- c("Stage2", "Retrospective", "Colossus", "Taxonomy")

resolve_string_base_dir <- function() {
  if (file.exists("analysis_output/stage2/Stage2_limma_high_vs_low_all_genes.csv")) {
    "."
  } else if (file.exists("Transcriptomics/analysis_output/stage2/Stage2_limma_high_vs_low_all_genes.csv")) {
    "Transcriptomics"
  } else {
    stop("Run from Transcriptomics/ or its parent directory.")
  }
}

string_cohort_configs <- function(base_dir) {
  list(
    Stage2 = list(
      prefix = "Stage2",
      dir = "stage2",
      limma = "Stage2_limma_high_vs_low_all_genes.csv"
    ),
    Retrospective = list(
      prefix = "Retrospective",
      dir = "retrospective",
      limma = "Retrospective_limma_high_vs_low_all_genes.csv"
    ),
    Colossus = list(
      prefix = "Colossus",
      dir = "colossus",
      limma = "Colossus_limma_high_vs_low_all_genes.csv"
    ),
    Taxonomy = list(
      prefix = "Taxonomy",
      dir = "taxonomy",
      limma = "Taxonomy_limma_high_vs_low_all_genes.csv"
    )
  )
}

is_protein_coding_symbol <- function(genes) {
  genes <- as.character(genes)
  ok <- grepl("^[A-Za-z][A-Za-z0-9-]*$", genes)
  bad_prefix <- grepl(
    "^(RP|RN|RNU|RPL|RPS|MT-|MTRNR|AC[0-9]|AL[0-9]|AP[0-9]|LINC|LOC|MIR|SNOR|U[0-9]|Y_RNA|TR)",
    genes,
    ignore.case = TRUE
  )
  bad_suffix <- grepl("P[0-9]+$|-AS[0-9]+$", genes, ignore.case = TRUE)
  ok & !bad_prefix & !bad_suffix & nchar(genes) >= 2L & nchar(genes) <= 15L
}

select_limma_genes_for_string <- function(
  limma_df,
  fdr_cutoff = 0.05,
  fc_cutoff = 0.5,
  direction = c("all", "high_up", "high_down"),
  max_genes = 200L,
  min_genes = 15L,
  p_use = c("fdr", "nominal", "auto")
) {
  direction <- match.arg(direction)
  p_use <- match.arg(p_use)
  p_col <- if ("adj.P.Val" %in% names(limma_df)) "adj.P.Val" else "padj"
  if (!p_col %in% names(limma_df)) {
    stop("limma table needs adj.P.Val column.")
  }
  if (!"P.Value" %in% names(limma_df)) {
    stop("limma table needs P.Value column for nominal fallback.")
  }

  filter_genes <- function(p_column, p_cutoff) {
    df <- limma_df %>%
      dplyr::filter(
        is.finite(.data$logFC),
        is.finite(.data[[p_column]]),
        !is.na(.data$gene),
        nzchar(.data$gene),
        .data[[p_column]] <= p_cutoff,
        abs(.data$logFC) >= fc_cutoff,
        is_protein_coding_symbol(.data$gene)
      )

    if (direction == "high_up") {
      df <- df %>% dplyr::filter(.data$logFC > 0)
      df <- df %>% dplyr::arrange(dplyr::desc(.data$logFC))
    } else if (direction == "high_down") {
      df <- df %>% dplyr::filter(.data$logFC < 0)
      df <- df %>% dplyr::arrange(.data$logFC)
    } else {
      df <- df %>% dplyr::arrange(.data[[p_column]], dplyr::desc(abs(.data$logFC)))
    }

    if (nrow(df) > max_genes) {
      df <- df %>% dplyr::slice_head(n = max_genes)
    }
    df
  }

  p_mode <- p_use
  df <- filter_genes(p_col, fdr_cutoff)
  if (p_use == "auto" && nrow(df) < min_genes) {
    df <- filter_genes("P.Value", fdr_cutoff)
    p_mode <- "nominal"
  } else if (p_use == "nominal") {
    df <- filter_genes("P.Value", fdr_cutoff)
    p_mode <- "nominal"
  } else if (p_use == "auto" && nrow(df) >= min_genes) {
    p_mode <- "fdr"
  }

  if (nrow(df) < min_genes) {
    return(list(
      genes = df,
      n_input = nrow(df),
      sufficient = FALSE,
      p_mode = p_mode
    ))
  }

  list(genes = df, n_input = nrow(df), sufficient = TRUE, p_mode = p_mode)
}

init_string_db <- function(
  score_threshold = 400L,
  string_version = "12.0",
  species = 9606L
) {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop(
      "Package 'STRINGdb' is required. Install with:\n",
      "  BiocManager::install('STRINGdb')"
    )
  }
  message("Initialising STRINGdb (first run may download the database)...")
  STRINGdb::STRINGdb$new(
    version = string_version,
    species = species,
    score_threshold = score_threshold
  )
}

logfc_to_color <- function(logfc) {
  pos <- logfc > 0
  cols <- rep("#CCCCCC", length(logfc))
  cols[pos] <- "#B2182B"
  cols[!pos] <- "#2166AC"
  cols
}

export_string_network_png <- function(
  string_db,
  string_ids,
  mapped_df,
  out_path,
  title = NULL
) {
  if (length(string_ids) < 2L) {
    message("Skipping network plot (<2 mapped proteins): ", out_path)
    return(invisible(NULL))
  }

  payload_id <- NULL
  if ("logFC" %in% names(mapped_df) && all(c("STRING_id", "logFC") %in% names(mapped_df))) {
    halo <- mapped_df %>%
      dplyr::filter(.data$STRING_id %in% string_ids) %>%
      dplyr::distinct(.data$STRING_id, .data$logFC) %>%
      dplyr::mutate(color = logfc_to_color(.data$logFC))
    if (nrow(halo) >= 2L) {
      payload_id <- tryCatch(
        string_db$post_payload(halo$STRING_id, colors = halo$color),
        error = function(e) {
          message("post_payload failed: ", conditionMessage(e))
          NULL
        }
      )
    }
  }

  grDevices::png(out_path, width = 1400, height = 1000, res = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (!is.null(title)) {
    graphics::plot.new()
    graphics::title(main = title, line = -1, adj = 0)
  }
  if (is.null(payload_id)) {
    string_db$plot_network(string_ids, add_link = FALSE)
  } else {
    string_db$plot_network(string_ids, payload_id = payload_id, add_link = FALSE)
  }
  message("Wrote ", out_path)
  invisible(out_path)
}

run_string_ppi_for_genes <- function(
  gene_df,
  output_dir,
  file_prefix,
  direction_label = "all_sig",
  score_threshold = 400L,
  string_version = "12.0"
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!isTRUE(gene_df$sufficient)) {
    msg <- paste0(
      "Skipping ", file_prefix, " (", direction_label, "): only ",
      gene_df$n_input, " genes after filtering (need >= 15)."
    )
    message(msg)
    return(list(skipped = TRUE, reason = msg, n_genes = gene_df$n_input))
  }

  de <- gene_df$genes %>%
    dplyr::distinct(.data$gene, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(abs(.data$logFC)))
  p_mode <- if (is.null(gene_df$p_mode)) "fdr" else gene_df$p_mode
  readr::write_csv(
    de %>%
      dplyr::mutate(string_p_mode = p_mode) %>%
      dplyr::select(dplyr::any_of(c("gene", "logFC", "adj.P.Val", "P.Value", "t", "string_p_mode"))),
    file.path(output_dir, paste0(file_prefix, "_string_", direction_label, "_input_genes.csv"))
  )

  string_db <- init_string_db(
    score_threshold = score_threshold,
    string_version = string_version
  )

  # STRINGdb$map() requires a base data.frame (tibbles trigger dimension errors).
  mapped <- string_db$map(as.data.frame(de), "gene", removeUnmappedRows = TRUE)
  if (nrow(mapped) < 15L) {
    msg <- paste0(
      "Skipping ", file_prefix, " (", direction_label, "): only ",
      nrow(mapped), " genes mapped to STRING."
    )
    message(msg)
    return(list(skipped = TRUE, reason = msg, n_mapped = nrow(mapped)))
  }

  readr::write_csv(
    mapped,
    file.path(output_dir, paste0(file_prefix, "_string_", direction_label, "_mapped.csv"))
  )

  string_ids <- unique(mapped$STRING_id)
  interactions <- tryCatch(
    string_db$get_interactions(string_ids),
    error = function(e) {
      message("get_interactions failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(interactions) && nrow(interactions) > 0L) {
    readr::write_csv(
      interactions,
      file.path(output_dir, paste0(file_prefix, "_string_", direction_label, "_interactions.csv"))
    )
  }

  enrichment <- tryCatch(
    string_db$get_enrichment(string_ids),
    error = function(e) {
      message("get_enrichment failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(enrichment) && nrow(enrichment) > 0L) {
    readr::write_csv(
      enrichment,
      file.path(output_dir, paste0(file_prefix, "_string_", direction_label, "_enrichment.csv"))
    )
  }

  ppi_enrich <- tryCatch(
    string_db$get_ppi_enrichment(string_ids),
    error = function(e) NULL
  )
  ppi_p <- NA_real_
  if (!is.null(ppi_enrich)) {
    if (is.data.frame(ppi_enrich) && "p_value" %in% names(ppi_enrich)) {
      ppi_p <- ppi_enrich$p_value[1]
    } else if (is.numeric(ppi_enrich)) {
      ppi_p <- as.numeric(ppi_enrich)[1]
    }
  }

  network_path <- file.path(
    output_dir,
    paste0(file_prefix, "_string_", direction_label, "_network.png")
  )
  export_string_network_png(
    string_db,
    string_ids,
    mapped,
    network_path,
    title = paste0(file_prefix, " High vs Low — ", direction_label)
  )

  list(
    skipped = FALSE,
    n_input = nrow(de),
    n_mapped = nrow(mapped),
    n_interactions = if (is.null(interactions)) 0L else nrow(interactions),
    n_enrichment_terms = if (is.null(enrichment)) 0L else nrow(enrichment),
    ppi_enrichment_p = ppi_p,
    p_mode = p_mode,
    network_path = network_path
  )
}

run_string_high_vs_low_cohort <- function(
  limma_csv,
  output_dir,
  file_prefix,
  fdr_cutoff = 0.05,
  fc_cutoff = 0.5,
  max_genes = 200L,
  score_threshold = 400L,
  string_version = "12.0",
  p_use = c("auto", "fdr", "nominal")
) {
  p_use <- match.arg(p_use)
  if (!file.exists(limma_csv)) {
    stop("Missing limma CSV: ", limma_csv)
  }
  limma_df <- readr::read_csv(limma_csv, show_col_types = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  directions <- c("all_sig", "high_up", "high_down")
  results <- list()
  for (dirn in directions) {
    direction_arg <- switch(
      dirn,
      all_sig = "all",
      high_up = "high_up",
      high_down = "high_down"
    )
    sel <- select_limma_genes_for_string(
      limma_df,
      fdr_cutoff = fdr_cutoff,
      fc_cutoff = fc_cutoff,
      direction = direction_arg,
      max_genes = max_genes,
      p_use = p_use
    )
    results[[dirn]] <- run_string_ppi_for_genes(
      sel,
      output_dir = output_dir,
      file_prefix = file_prefix,
      direction_label = dirn,
      score_threshold = score_threshold,
      string_version = string_version
    )
  }
  results
}

run_string_high_vs_low_all <- function(
  output_dir = NULL,
  cohorts = DEFAULT_STRING_COHORTS,
  fdr_cutoff = 0.05,
  fc_cutoff = 0.5,
  max_genes = 200L,
  score_threshold = 400L,
  string_version = "12.0",
  p_use = c("auto", "fdr", "nominal")
) {
  p_use <- match.arg(p_use)
  if (!requireNamespace("readr", quietly = TRUE)) library(readr)
  if (!requireNamespace("dplyr", quietly = TRUE)) library(dplyr)

  base_dir <- resolve_string_base_dir()
  if (is.null(output_dir)) {
    output_dir <- file.path(base_dir, "analysis_output", "string_ppi")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cfgs <- string_cohort_configs(base_dir)
  all_results <- list()

  for (cn in cohorts) {
    cfg <- cfgs[[cn]]
    limma_csv <- file.path(base_dir, "analysis_output", cfg$dir, cfg$limma)
    cohort_out <- file.path(output_dir, cfg$dir)
    message("\n=== STRING PPI: ", cfg$prefix, " (High vs Low) ===")
    all_results[[cn]] <- run_string_high_vs_low_cohort(
      limma_csv = limma_csv,
      output_dir = cohort_out,
      file_prefix = cfg$prefix,
      fdr_cutoff = fdr_cutoff,
      fc_cutoff = fc_cutoff,
      max_genes = max_genes,
      score_threshold = score_threshold,
      string_version = string_version,
      p_use = p_use
    )
  }

  summary_lines <- c(
    "STRINGdb PPI — High vs Low RSF risk (limma DE genes)",
    "",
    paste0("Contrast: High vs Low (positive logFC = higher in High risk)"),
    paste0("Gene filters: adj.P.Val <= ", fdr_cutoff,
           " (auto-fallback to P.Value if <15 genes), |logFC| >= ", fc_cutoff,
           ", protein-coding-like symbols, max ", max_genes, " genes"),
    paste0("STRING score threshold: ", score_threshold, " (medium confidence)"),
    paste0("Output root: ", output_dir),
    "",
    "Per cohort / direction:",
    "  {prefix}_string_{all_sig|high_up|high_down}_input_genes.csv",
    "  {prefix}_string_*_mapped.csv",
    "  {prefix}_string_*_interactions.csv",
    "  {prefix}_string_*_enrichment.csv",
    "  {prefix}_string_*_network.png  (red = up in High, blue = down in High)",
    ""
  )

  for (cn in cohorts) {
    summary_lines <- c(summary_lines, paste0("--- ", cn, " ---"))
    res <- all_results[[cn]]
    for (dirn in names(res)) {
      r <- res[[dirn]]
      if (isTRUE(r$skipped)) {
        summary_lines <- c(summary_lines, paste0("  ", dirn, ": SKIPPED — ", r$reason))
      } else {
        summary_lines <- c(
          summary_lines,
          paste0(
            "  ", dirn, ": ", r$n_mapped, " mapped, ",
            r$n_interactions, " interactions, ",
            r$n_enrichment_terms, " enrichment terms",
            if (!is.null(r$p_mode)) paste0(" [", r$p_mode, " p]") else "",
            if (is.finite(r$ppi_enrichment_p)) {
              paste0(", PPI enrichment p=", signif(r$ppi_enrichment_p, 3))
            } else {
              ""
            }
          )
        )
      }
    }
    summary_lines <- c(summary_lines, "")
  }

  writeLines(summary_lines, file.path(output_dir, "string_ppi_summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")
  invisible(all_results)
}
