# Run preranked fgsea GSEA for all four transcriptomic cohorts.
# Requires limma DE CSVs (run_limma_de_all.R or cohort limma Rmds).
#   /work_space/envs/transcriptomics2/bin/Rscript run_gsea_all.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("analysis_helpers.R")
source("gsea_helpers.R")

resolve_paths <- function() {
  if (file.exists("analysis_output/stage2/Stage2_limma_high_vs_low_all_genes.csv")) {
    list(base = ".", proteomics = "../proteomics")
  } else if (file.exists("Transcriptomics/analysis_output/stage2/Stage2_limma_high_vs_low_all_genes.csv")) {
    list(base = "Transcriptomics", proteomics = "proteomics")
  } else {
    stop("Run from Transcriptomics/ or its parent directory.")
  }
}

p <- resolve_paths()
base_dir <- p$base

cohorts <- list(
  stage2 = list(
    prefix = "Stage2",
    title = "Stage 2 Transcriptomics",
    dir = "stage2",
    limma = "Stage2_limma_high_vs_low_all_genes.csv"
  ),
  retrospective = list(
    prefix = "Retrospective",
    title = "Retrospective Transcriptomics",
    dir = "retrospective",
    limma = "Retrospective_limma_high_vs_low_all_genes.csv"
  ),
  colossus = list(
    prefix = "Colossus",
    title = "Colossus Transcriptomics",
    dir = "colossus",
    limma = "Colossus_limma_high_vs_low_all_genes.csv"
  ),
  taxonomy = list(
    prefix = "Taxonomy",
    title = "Taxonomy RNA-Seq",
    dir = "taxonomy",
    limma = "Taxonomy_limma_high_vs_low_all_genes.csv"
  )
)

all_hallmark <- list()
all_kegg <- list()

for (cohort in names(cohorts)) {
  cfg <- cohorts[[cohort]]
  dirs <- setup_analysis_dirs(base_dir, cfg$dir)
  limma_csv <- file.path(dirs$output_dir, cfg$limma)
  if (!file.exists(limma_csv)) {
    stop(
      "Missing limma results for ", cfg$prefix, ": ", limma_csv,
      "\nRun run_limma_de_all.R or ", cfg$prefix, "_limma_de_by_risk_group.Rmd first."
    )
  }
  cat("\n=== ", cfg$prefix, " ===\n", sep = "")
  out <- run_gsea_from_limma_csv(
    limma_csv    = limma_csv,
    output_dir   = dirs$output_dir,
    file_prefix  = cfg$prefix,
    cohort_title = cfg$title,
    collections  = c("hallmark", "kegg")
  )
  all_hallmark[[cfg$prefix]] <- out$hallmark
  all_kegg[[cfg$prefix]] <- out$kegg
}

meta_dir <- file.path(base_dir, "analysis_output")
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n=== Fisher meta (Hallmark, 4 cohorts — fgsea preranked) ===\n")
run_fgsea_meta_fisher(
  fgsea_results_list = all_hallmark,
  output_dir = meta_dir,
  file_prefix = "Transcriptomics",
  collection = "hallmark"
)

cat("\n=== Fisher meta (KEGG, 4 cohorts — fgsea preranked) ===\n")
run_fgsea_meta_fisher(
  fgsea_results_list = all_kegg,
  output_dir = meta_dir,
  file_prefix = "Transcriptomics",
  collection = "kegg"
)

cat("\nGSEA complete (all cohorts: limma DE → fgsea preranked).\n")
