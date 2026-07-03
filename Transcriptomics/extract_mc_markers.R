#!/usr/bin/env Rscript
# Extract top CODEX markers for selected metaclusters.
suppressPackageStartupMessages(library(dplyr))

path <- "/work_space/files/Features/Run2_all/Manual_classes/act4_metacluster_markers.csv"
out_dir <- "/work_space/files/Transcriptomics/analysis_output/metacluster_correlation"
df <- read.csv(path, check.names = FALSE)

parse_marker <- function(col) {
  sub("_MC_[0-9]+_mean$", "", sub("^mc_marker_Mean\\.Cell\\.", "", col))
}

parse_mc <- function(col) {
  sub(".*_(MC_[0-9]+)_mean$", "\\1", col)
}

all_mc_cols <- grep("^mc_marker_Mean\\.Cell\\..*_MC_[0-9]+_mean$", names(df), value = TRUE)

long <- bind_rows(lapply(all_mc_cols, function(c) {
  data.frame(
    marker = parse_marker(c),
    mc = parse_mc(c),
    val = mean(df[[c]], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

rank_mc_markers <- function(mc, top_n = 15L) {
  cols <- grep(
    paste0("^mc_marker_Mean\\.Cell\\..*", mc, "_mean$"),
    names(df),
    value = TRUE
  )
  tibble(
    marker = vapply(cols, parse_marker, character(1)),
    mean_intensity = colMeans(df[, cols, drop = FALSE], na.rm = TRUE)
  ) %>%
    arrange(desc(.data$mean_intensity)) %>%
    head(top_n)
}

enriched_vs_other <- function(mc_id, top_n = 15L) {
  sub <- long %>%
    dplyr::filter(.data$mc == mc_id) %>%
    dplyr::select(marker, mc_mean = val)
  other <- long %>%
    dplyr::filter(.data$mc != mc_id) %>%
    dplyr::group_by(.data$marker) %>%
    dplyr::summarise(other_mean = mean(.data$val), .groups = "drop")
  sub %>%
    dplyr::left_join(other, by = "marker") %>%
    dplyr::mutate(
      log2fc = log2((.data$mc_mean + 1e-3) / (.data$other_mean + 1e-3)),
      diff = .data$mc_mean - .data$other_mean
    ) %>%
    dplyr::arrange(dplyr::desc(.data$log2fc)) %>%
    head(top_n)
}

targets <- c("MC_8", "MC_13", "MC_17")
summary_lines <- c(
  "Top spatial (CODEX) markers for PROGENy-linked metaclusters",
  paste0("Source: ", path),
  paste0("TMA spots averaged: n=", nrow(df)),
  ""
)

all_enriched <- list()
for (mc in targets) {
  top_mean <- rank_mc_markers(mc, 12L)
  top_enr <- enriched_vs_other(mc, 12L)
  all_enriched[[mc]] <- top_enr %>%
    mutate(metacluster = mc)

  readr::write_csv(
    top_mean,
    file.path(out_dir, paste0(mc, "_top_markers_by_intensity.csv"))
  )
  readr::write_csv(
    top_enr,
    file.path(out_dir, paste0(mc, "_top_markers_enriched_vs_other.csv"))
  )

  summary_lines <- c(
    summary_lines,
    paste0("=== ", mc, " (top by mean intensity) ==="),
    paste0("  ", top_mean$marker, " (mean=", signif(top_mean$mean_intensity, 4), ")", collapse = "\n"),
    "",
    paste0("=== ", mc, " (most enriched vs other metaclusters) ==="),
    paste0(
      "  ",
      top_enr$marker,
      " (log2FC=",
      signif(top_enr$log2fc, 3),
      ", mc=",
      signif(top_enr$mc_mean, 3),
      " vs other=",
      signif(top_enr$other_mean, 3),
      ")",
      collapse = "\n"
    ),
    ""
  )
}

readr::write_csv(
  bind_rows(all_enriched),
  file.path(out_dir, "MC_8_13_17_markers_enriched_combined.csv")
)
writeLines(summary_lines, file.path(out_dir, "MC_8_13_17_marker_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
