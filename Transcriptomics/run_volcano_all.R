suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(EnhancedVolcano)
})
source("analysis_helpers.R")

run_volcano <- function(differential_results_clean, title, subtitle, out_path) {
  diff <- differential_results_clean[
    is.finite(differential_results_clean$log2FC) & !is.na(differential_results_clean$pval),
  ]
  axes <- volcano_axis_limits(diff)
  p <- suppressWarnings(
    EnhancedVolcano(
      diff,
      lab       = diff$cell_type,
      x         = "log2FC",
      y         = "pval",
      pCutoff   = 0.05001,
      FCcutoff  = 0.5,
      title     = title,
      subtitle  = subtitle,
      xlab      = "Log2 Fold Change (High / Low)",
      ylab      = "-Log10(p-value)",
      pointSize = 2.0,
      labSize   = 3.5,
      drawConnectors        = TRUE,
      widthConnectors       = 0.4,
      maxoverlapsConnectors = Inf,
      xlim = axes$xlim,
      ylim = axes$ylim
    ) +
    theme(plot.margin = margin(20, 60, 20, 60), legend.position = "bottom")
  )
  ggsave(out_path, plot = p, width = 12, height = 8, dpi = 300, bg = "white")
  cat("Saved:", out_path, " xlim=", paste(round(axes$xlim, 2), collapse = ","),
      " ylim=", paste(round(axes$ylim, 2), collapse = ","), "\n", sep = "")
  diff[order(abs(diff$log2FC), decreasing = TRUE), c("cell_type", "log2FC", "pval")]
}

count_risk <- function(scores_long) {
  n_ct <- n_distinct(scores_long$cell_type)
  list(
    n = n_distinct(scores_long$sample),
    low = sum(scores_long$risk_grp == "Low") / n_ct,
    high = sum(scores_long$risk_grp == "High") / n_ct
  )
}

dirs_s2 <- setup_analysis_dirs(".", "stage2")
w_s2 <- ensure_working_copies(c(
  xcell = "Stage2_transcriptomics_xcell.csv",
  clinical = "stage2/clinical.csv",
  mapping = "../spatial/master_patient_mapping.csv",
  survival = "../proteomics/survival_df_with_risk.rds",
  ge = "stage2/ge_project_rcsi_cohort_clinical.csv"
), dirs_s2$work_dir)

stage2_res <- read_csv(w_s2["xcell"], show_col_types = FALSE)
clinical <- read_csv(w_s2["clinical"], show_col_types = FALSE)
master_patient_mapping <- read_csv(w_s2["mapping"], show_col_types = FALSE)
survival_df <- readRDS(w_s2["survival"])
stage2 <- stage2_res
r_codes <- clinical$r_code[match(colnames(stage2)[-1], clinical$patient_id)]
patient_id_g <- master_patient_mapping$patient_id_g[match(r_codes, master_patient_mapping$r_code)]
id_map <- data.frame(t_id = colnames(stage2)[-1], g_id = patient_id_g, stringsAsFactors = FALSE)
scores_s2 <- stage2_res %>%
  rename(cell_type = 1) %>%
  pivot_longer(-cell_type, names_to = "sample", values_to = "enrichment") %>%
  left_join(id_map %>% rename(sample = t_id, patient_id = g_id), by = "sample") %>%
  left_join(survival_df %>% select(patient_id, risk_grp, rsf_risk), by = "patient_id") %>%
  filter(risk_grp %in% c("Low", "High"))
rk <- count_risk(scores_s2)
cat("\n=== Stage 2 === n=", rk$n, " Low=", rk$low, " High=", rk$high, "\n", sep = "")
diff_s2 <- wilcox_enrichment_diff(scores_s2)
print(head(run_volcano(
  diff_s2,
  "Cell Type Enrichment: High vs Low Risk",
  paste0("Stage 2 Transcriptomics | xCell2 scores | n = ", rk$n,
         " (Low ", rk$low, ", High ", rk$high, ")"),
  file.path(dirs_s2$output_dir, "Stage2_risk_enrichment_volcano.png")
), 8))
print(export_top_log2fc_table(
  diff_s2,
  file.path(dirs_s2$output_dir, "Stage2_top5_log2fc_up_down.csv")
))

dirs_r <- setup_analysis_dirs(".", "retrospective")
w_r <- ensure_working_copies(c(
  xcell = "reterospective_transcriptomics_xcell_results.csv",
  survival = "../proteomics/survival_df_with_risk.rds"
), dirs_r$work_dir)
xcell_r <- read_csv(w_r["xcell"], show_col_types = FALSE)
survival_df <- readRDS(w_r["survival"])
scores_r <- xcell_r %>%
  rename(cell_type = 1) %>%
  pivot_longer(-cell_type, names_to = "sample", values_to = "enrichment") %>%
  mutate(patient_id = sample) %>%
  left_join(survival_df %>% select(patient_id, risk_grp, rsf_risk), by = "patient_id") %>%
  filter(risk_grp %in% c("Low", "High"))
rk <- count_risk(scores_r)
cat("\n=== Retrospective === n=", rk$n, " Low=", rk$low, " High=", rk$high, "\n", sep = "")
diff_r <- wilcox_enrichment_diff(scores_r)
print(head(run_volcano(
  diff_r,
  "Cell Type Enrichment: High vs Low Risk",
  paste0("Retrospective Transcriptomics | xCell2 scores | n = ", rk$n,
         " (Low ", rk$low, ", High ", rk$high, ")"),
  file.path(dirs_r$output_dir, "Retrospective_risk_enrichment_volcano.png")
), 8))
print(export_top_log2fc_table(
  diff_r,
  file.path(dirs_r$output_dir, "Retrospective_top5_log2fc_up_down.csv")
))

dirs_c <- setup_analysis_dirs(".", "colossus")
w_c <- ensure_working_copies(c(
  xcell = "collosus_transcriptomics_xcell_results.csv",
  survival = "../proteomics/survival_df_with_risk.rds",
  master = "../Integrated /All_patinet_IDS_concatonated_masterdoc.csv"
), dirs_c$work_dir)
xcell_c <- read_csv(w_c["xcell"], show_col_types = FALSE)
survival_df <- readRDS(w_c["survival"])
master_df <- read_csv(w_c["master"], show_col_types = FALSE)
strip_fp <- function(x) str_replace(x, "-FP.*$", "")
master_col_map <- bind_rows(
  master_df %>% filter(!is.na(patient_id_g), patient_id_g != "") %>%
    transmute(patient_id = patient_id_g, col_id = strip_fp(COLOSSUS_ID)),
  master_df %>% filter(!is.na(patient_id_g), patient_id_g != "") %>%
    transmute(patient_id = patient_id_g, col_id = strip_fp(Old_COLOSSUS_ID)),
  master_df %>% filter(!is.na(patient_id_g), patient_id_g != "") %>%
    transmute(patient_id = patient_id_g, col_id = strip_fp(Sample_ID)),
  master_df %>% filter(!is.na(patient_id_g), patient_id_g != "") %>%
    transmute(patient_id = patient_id_g, col_id = strip_fp(Alternative_ID))
) %>% filter(!is.na(col_id), col_id != "", col_id != "NA") %>% distinct(col_id, .keep_all = TRUE)
samples_c <- colnames(xcell_c)[-1]
id_map_c <- data.frame(sample = samples_c, stringsAsFactors = FALSE) %>%
  left_join(master_col_map, by = c("sample" = "col_id"))
scores_c <- xcell_c %>%
  rename(cell_type = 1) %>%
  pivot_longer(-cell_type, names_to = "sample", values_to = "enrichment") %>%
  left_join(id_map_c, by = "sample") %>%
  left_join(survival_df %>% select(patient_id, risk_grp, rsf_risk), by = "patient_id") %>%
  filter(risk_grp %in% c("Low", "High"))
rk <- count_risk(scores_c)
cat("\n=== Colossus === n=", rk$n, " Low=", rk$low, " High=", rk$high, "\n", sep = "")
diff_c <- wilcox_enrichment_diff(scores_c)
print(head(run_volcano(
  diff_c,
  "Cell Type Enrichment: High vs Low Risk",
  paste0("Colossus Transcriptomics | xCell2 scores | n = ", rk$n,
         " (Low ", rk$low, ", High ", rk$high, ")"),
  file.path(dirs_c$output_dir, "Colossus_risk_enrichment_volcano.png")
), 8))
print(export_top_log2fc_table(
  diff_c,
  file.path(dirs_c$output_dir, "Colossus_top5_log2fc_up_down.csv")
))

cat("\nDone.\n")
