# Separate PROGENy heatmaps per metacluster (MC_13, MC_8).
# Full colour by rho; weak |rho| stays near-white; black rings = significant.
#   cd /work_space/files/Transcriptomics
#   /work_space/envs/transcriptomics2/bin/Rscript run_mc_progeny_sig_heatmap.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

base_dir <- if (file.exists("analysis_output/metacluster_correlation/MC_progeny_fisher_meta.csv")) {
  "."
} else if (file.exists("Transcriptomics/analysis_output/metacluster_correlation/MC_progeny_fisher_meta.csv")) {
  "Transcriptomics"
} else {
  stop("Run from Transcriptomics/ or its parent directory.")
}

out_dir <- file.path(base_dir, "analysis_output", "metacluster_correlation", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta_path <- file.path(base_dir, "analysis_output", "metacluster_correlation", "MC_progeny_fisher_meta.csv")
meta <- readr::read_csv(meta_path, show_col_types = FALSE)

COHORTS <- c("Stage2", "Retrospective", "Colossus", "Taxonomy")
RHO_LIMIT <- 0.5
RHO_DEADZONE <- 0.15

mc_specs <- list(
  MC_13 = list(
    label = "MC_13 (Low risk)",
    subtitle = "Immune-enriched metacluster",
    file_tag = "MC13"
  ),
  MC_8 = list(
    label = "MC_8 (High risk)",
    subtitle = "Stress/apoptosis-enriched metacluster",
    file_tag = "MC8"
  )
)

pathway_order <- meta %>%
  filter(.data$metacluster %in% names(mc_specs), .data$feature != "n_samples") %>%
  group_by(.data$feature) %>%
  summarise(avg = mean(abs(.data$mean_rho), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(.data$avg)) %>%
  pull(.data$feature)

prepare_mc_heatmap_df <- function(meta_row_tbl) {
  cohort_long <- meta_row_tbl %>%
    select(
      feature,
      dplyr::starts_with("rho_"),
      dplyr::starts_with("pval_")
    ) %>%
    pivot_longer(
      cols = dplyr::starts_with("rho_"),
      names_to = "cohort",
      names_prefix = "rho_",
      values_to = "rho"
    ) %>%
    left_join(
      meta_row_tbl %>%
        select(feature, dplyr::starts_with("pval_")) %>%
        pivot_longer(
          cols = dplyr::starts_with("pval_"),
          names_to = "cohort",
          names_prefix = "pval_",
          values_to = "pval"
        ),
      by = c("feature", "cohort")
    ) %>%
    mutate(
      panel = cohort,
      sig = is.finite(.data$pval) & .data$pval <= 0.05
    )

  meta_col <- meta_row_tbl %>%
    transmute(
      feature = .data$feature,
      panel = "Meta",
      rho = .data$mean_rho,
      pval = .data$fisher_padj,
      sig = is.finite(.data$fisher_padj) &
        .data$fisher_padj <= 0.05 &
        .data$sign_consistent %in% TRUE
    )

  bind_rows(cohort_long, meta_col) %>%
    filter(is.finite(.data$rho)) %>%
    mutate(
      panel = factor(.data$panel, levels = c(COHORTS, "Meta")),
      feature = factor(.data$feature, levels = rev(pathway_order)),
      label = sprintf("%.2f", .data$rho),
      star = case_when(
        !.data$sig ~ "",
        .data$panel == "Meta" & .data$pval < 0.001 ~ "***",
        .data$panel == "Meta" & .data$pval < 0.01 ~ "**",
        .data$panel == "Meta" & .data$pval < 0.05 ~ "*",
        TRUE ~ "*"
      )
    )
}

rho_fill_scale <- function() {
  breaks <- c(-RHO_LIMIT, -0.3, -RHO_DEADZONE, 0, RHO_DEADZONE, 0.3, RHO_LIMIT)
  colours <- c(
    "#2166AC",
    "#67A9CF",
    "#F5F5F5",
    "#F5F5F5",
    "#FDBF6F",
    "#EF8A62",
    "#B2182B"
  )
  scale_fill_gradientn(
    colours = colours,
    values = scales::rescale(breaks, to = c(0, 1)),
    limits = c(-RHO_LIMIT, RHO_LIMIT),
    oob = scales::squish,
    name = "Spearman rho",
    breaks = c(-0.5, -0.25, 0, 0.25, 0.5)
  )
}

plot_mc_progeny_heatmap <- function(plot_df, title, subtitle) {
  plot_df <- plot_df %>%
    mutate(
      x_num = as.numeric(.data$panel),
      y_num = as.numeric(.data$feature),
      x_star = x_num + 0.5 - 0.16,
      y_star = y_num + 0.5 - 0.16
    )
  star_df <- plot_df %>% filter(nzchar(.data$star))

  ggplot(plot_df, aes(x = panel, y = feature)) +
    geom_tile(aes(fill = rho), colour = "grey92", linewidth = 0.45) +
    geom_text(aes(label = label), size = 2.7, colour = "grey15") +
    geom_label(
      data = star_df,
      aes(x = x_star, y = y_star, label = star),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 2.3,
      fontface = "bold",
      colour = "black",
      fill = "white",
      linewidth = 0,
      label.padding = grid::unit(0.08, "lines")
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    rho_fill_scale() +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "PROGENy pathway",
      caption = "* cohort p<=0.05; Meta: * p<0.05, ** p<0.01, *** p<0.001"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold", angle = 30, hjust = 1),
      plot.title = element_text(face = "bold", size = 12),
      plot.caption = element_text(size = 8, hjust = 0, colour = "grey30"),
      legend.position = "right",
      legend.key.height = grid::unit(0.45, "cm"),
      legend.key.width = grid::unit(0.25, "cm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      plot.margin = ggplot2::margin(5, 10, 8, 5)
    )
}

for (mc in names(mc_specs)) {
  spec <- mc_specs[[mc]]
  mc_df <- meta %>%
    filter(.data$metacluster == mc, .data$feature != "n_samples")

  plot_df <- prepare_mc_heatmap_df(mc_df)
  p <- plot_mc_progeny_heatmap(
    plot_df,
    title = paste0(spec$label, " vs PROGENy"),
    subtitle = paste0(
      spec$subtitle,
      " | pale = weak |rho|; stars = significant"
    )
  )

  out_png <- file.path(out_dir, paste0(spec$file_tag, "_progeny_heatmap.png"))
  out_pdf <- file.path(out_dir, paste0(spec$file_tag, "_progeny_heatmap.pdf"))
  ggsave(out_png, p, width = 7.8, height = 6.5, dpi = 300, bg = "white")
  ggsave(out_pdf, p, width = 7.8, height = 6.5, bg = "white")
  message("Saved: ", out_png)
  message("Saved: ", out_pdf)
}

index_lines <- c(
  "PROGENy heatmaps (one per metacluster)",
  paste0("Generated: ", Sys.time()),
  "",
  "MC13_progeny_heatmap.png — MC_13 (Low risk)",
  "MC8_progeny_heatmap.png — MC_8 (High risk)",
  "",
  "Colour: |rho| < 0.15 stays near-white",
  "Stars: top-right of cell (* cohort p<=0.05; **/*** Meta FDR)"
)
writeLines(index_lines, file.path(out_dir, "mc_progeny_heatmap_index.txt"))

cat("\nDone.\n")
