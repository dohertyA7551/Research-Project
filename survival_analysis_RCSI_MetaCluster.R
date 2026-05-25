###########################################################
# ARCHIVED EXPLORATORY SCRIPT
#
# Contains:
# - early survival analyses
# - debugging attempts
# - TMA/patient mapping code
# - exploratory biology integration
#
# WARNING:
# Objects and analyses may rely on incomplete or inconsistent
# metacluster mappings and should not be treated as final.
#
# Clean restart began: 2026-05-24
###########################################################



















library(ggplot2)
library(dplyr)
library(tidyr)
library(survival)
library(coxphf)
library(purrr)
library(glmnet)
library(readr)
library(broom)

clin_data<- read.csv("/work_space/files/spatial/ColonRO1_2020_RCSI_clinical.csv")
p_id <-read.csv("/work_space/files/spatial/master_patient_mapping.csv")

meta_df <- read.csv("Features/Run2_all/Manual_classes/act3_metacluster_abundance.csv")
master_features <- read_csv("Features/Run2_all/Manual_classes/master_features.csv")

View(master_features)
names(meta_df)
names(p_id)
names(clin_data)
names


#####merging
#selecting one Leiden param setting to look at (1.0 - 12 clusters)
meta_leiden1 <- meta_df[, c(
  "tma_id",
  names(meta_df)[startsWith(
    names(meta_df),
    "abundance_metacluster_leiden_1.0_"
  )]
)]

str(meta_leiden1)
head(meta_leiden1)

####attatching patient_id_t
meta_leiden1 <- meta_leiden1 %>%
  mutate(tma_id = as.character(tma_id))

p_id <- p_id %>%
  mutate(
    SS_id        = as.character(SS_id),
    patient_id_g = as.character(patient_id_g)
  )

## Attach patient_id_t onto each TMA sample
meta_with_pid_g <- meta_leiden1 %>%
  left_join(
    p_id %>% select(SS_id, patient_id_g),
    by = c("tma_id" = "SS_id")
  ) %>%
  filter(!is.na(patient_id_g))

### see how many i lost
# Number of TMAs that successfully mapped to a patient_id_g
n_mapped <- sum(!is.na(meta_with_pid_g$patient_id_g))
n_total  <- nrow(meta_with_pid_g)

n_mapped
n_total
n_mapped / n_total  # proportion mapped
###343 mapped but we check again idk

meta_mapped_only <- meta_with_pid_g %>%
  filter(!is.na(patient_id_g))

nrow(meta_mapped_only)  # should be 343
head(meta_mapped_only[, c("tma_id", "patient_id_g")])

# Leiden 1.0 meta-cluster columns
abund_cols <- names(meta_mapped_only)[startsWith(
  names(meta_mapped_only),
  "abundance_metacluster_leiden_1.0_"
)]

meta_patient <- meta_mapped_only %>%
  group_by(patient_id_g) %>%
  summarise(
    across(all_of(abund_cols), mean, na.rm = TRUE),
    .groups = "drop"
  )

nrow(meta_patient)  # number of patients with mapped TMAs
head(meta_patient)


meta_mapped_only %>%
  count(patient_id_g) %>%
  arrange(desc(n)) %>%
  head(20)

names(meta_patient)
###only like 36 people.....

# Make IDs character
clin_data$patient_id      <- as.character(clin_data$Patient)
meta_patient$patient_id_g <- as.character(meta_patient$patient_id_g)

# Join clinical + meta-cluster data
surv_df <- clin_data %>%
  inner_join(
    meta_patient,
    by = c("Patient" = "patient_id_g")
  )

nrow(surv_df)
names(surv_df)[1:20]

surv_df %>%
  select(patient_id, os_months, os_event, dfs_months, dfs_event) %>%
  head()



###### survival analysis
# Recode "yes"/"no" to 1/0
table(surv_df$os_event, useNA = "ifany")

surv_df$os_event_num <- ifelse(surv_df$os_event == "yes", 1, 0)

table(surv_df$os_event, surv_df$os_event_num, useNA = "ifany")



# Ensure time is numeric
surv_df$os_months <- as.numeric(surv_df$os_months)

surv_os <- Surv(
  time  = surv_df$os_months,
  event = surv_df$os_event_num
)

# Meta-cluster columns
abund_cols <- names(surv_df)[startsWith(
  names(surv_df),
  "abundance_metacluster_leiden_1.0_"
)]

# Optional: scale covariates
surv_df[abund_cols] <- scale(surv_df[abund_cols])
######Run to here then hop down####







##probelm solvdeo... bulk version but its a bit intense
# identify abundance columns that are not all NA
valid_abund_cols <- abund_cols[colSums(!is.na(surv_df[abund_cols])) > 0]

length(abund_cols)        # original number of metaclusters
length(valid_abund_cols)  # number remaining after filtering

# keep only the valid ones
abund_cols <- valid_abund_cols

# scale only the valid abundance columns
surv_df[abund_cols] <- scale(surv_df[abund_cols])

# recommended: use Surv() directly in the formula
fit_os <- coxph(
  Surv(os_months, os_event_num) ~ .,
  data = surv_df[, c("os_months", "os_event_num", abund_cols)]
)

summary(fit_os)


# Cox model
formula_os <- as.formula(
  paste("surv_os ~", paste(abund_cols, collapse = " + "))
)

fit_os <- coxph(formula_os, data = surv_df)
summary(fit_os)


######## The above method is stat dubous.. lassso for cluster selection then coxph can follow#######
# keep rows with non-missing time and event
keep <- !is.na(surv_df$os_months) & !is.na(surv_df$os_event_num)

table(keep)  # just to see how many are dropped, if any

surv_df_clean <- surv_df[keep, ]

keep_time_pos <- surv_df_clean$os_months > 0

table(keep_time_pos)

surv_df_lasso <- surv_df_clean[keep_time_pos, ]

y <- Surv(
  time  = surv_df_lasso$os_months,
  event = surv_df_lasso$os_event_num
)

X <- as.matrix(surv_df_lasso[abund_cols])

set.seed(123)
cvfit <- cv.glmnet(
  x = X,
  y = y,
  family   = "cox",
  alpha    = 1,
  cox.ties = "breslow"
)

beta_lasso <- coef(cvfit, s = "lambda.min")

beta_vec <- as.numeric(beta_lasso)
names(beta_vec) <- rownames(beta_lasso)

selected_clusters <- names(beta_vec)[beta_vec != 0]
selected_clusters
length(selected_clusters)


###LASSO min returned 15... will try with 15.. if still to sparse.. will try with 1se
surv_selected <- surv_df_lasso[, c("os_months", "os_event_num", selected_clusters)]

fit_os_selected <- coxph(
  Surv(os_months, os_event_num) ~ .,
  data = surv_selected
)

summary(fit_os_selected)


####plot attempts####



tidy_fit <- broom::tidy(fit_os_selected, conf.int = TRUE, exponentiate = TRUE)

# Keep only metacluster terms (drop intercept if present)
hr_table <- tidy_fit %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value)

hr_table

ggplot(hr_table, aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  scale_x_log10() +
  labs(
    x = "Hazard ratio (log scale)",
    y = "Metacluster",
    title = "Hazard ratios for selected metaclusters"
  ) +
  theme_bw()

#########Sub setting##############################

##now lets try with top few
complete_columns <- c(
  "abundance_metacluster_leiden_1.0_MC_14",
  "abundance_metacluster_leiden_1.0_MC_2",
  "abundance_metacluster_leiden_1.0_MC_3"
)
 

# Define variables needed for this model
vars_model <- c("os_months", "os_event_num", complete_columns )

# Identify complete cases across these variables
cc_idx <- complete.cases(surv_df[, vars_model])
sum(cc_idx)  # number of patients used in the model

surv_cc <- surv_df[cc_idx, ]

# Optional: scale covariates within the complete-case set
surv_cc[complete_columns ] <- scale(surv_cc[complete_columns])

# Survival object
surv_os <- Surv(
  time  = surv_cc$os_months,
  event = surv_cc$os_event_num
)

# Cox formula and fit
formula_os <- as.formula(
  paste("surv_os ~", paste(complete_columns , collapse = " + "))
)

fit_os <- coxph(formula_os, data = surv_cc)
summary(fit_os)
 
#####Additional Stats
tidy_fit <- tidy(fit_os_selected, conf.int = TRUE, exponentiate = TRUE)

hr_table <- tidy_fit %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value)

# Drop clusters with problematic CIs (0 or Inf)
hr_table_clean <- hr_table %>%
  filter(is.finite(conf.low),
         is.finite(conf.high),
         conf.low > 0,
         conf.high > 0)

hr_table_clean

# Rank by HR
hr_table_clean <- hr_table_clean %>%
  arrange(estimate)

# Top 5 best prognosis (lowest HR)
best5 <- hr_table_clean %>% slice_head(n = 5)
best5

# Top 5 worst prognosis (highest HR)
worst5 <- hr_table_clean %>% slice_tail(n = 5)
worst5






#########plottig
s_os <- summary(fit_os)

multi_os_df <- data.frame(
  feature = rownames(s_os$coefficients),
  HR      = s_os$coefficients[,"exp(coef)"],
  lower95 = s_os$conf.int[,"lower .95"],
  upper95 = s_os$conf.int[,"upper .95"],
  pval    = s_os$coefficients[,"Pr(>|z|)"],
  row.names = NULL
)

multi_os_df
write.csv(multi_os_df, "multivariable_OS_meta_clusters.csv", row.names = FALSE)

multi_os_df <- multi_os_df %>%
  arrange(HR) %>%
  mutate(feature = factor(feature, levels = feature))

ggplot(multi_os_df, aes(x = HR, y = feature)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.2) +
  scale_x_log10() +
  xlab("Hazard ratio (log scale)") +
  ylab("Meta-cluster (Leiden 1.0)") +
  theme_bw()

#####Run this code if coming from  top and bottom 5 section

###ok cool now lets add some biologia

selected_clusters <- unique(c(best5$term, worst5$term))


meta_markers_mc <- meta_markers %>%
  left_join(
    meta_mapped_only %>%
      select(tma_id, all_of(selected_clusters)),
    by = "tma_id"
  )

selected_clusters
is.character(selected_clusters)

names(meta_markers_mc)[
  names(meta_markers_mc) %in% selected_clusters
]

#statistically messy to do all...
##search for best and worst (previously identified)
mc12_cols <- grep(
  "_MC_12_mean$",
  names(meta_markers),
  value = TRUE
)

mc5_cols <- grep(
  "_MC_5_mean$",
  names(meta_markers),
  value = TRUE
)

length(mc12_cols)
length(mc5_cols)

head(mc12_cols)
head(mc5_cols)

mc5_proteins <- mc5_cols %>%
  sub("^mc_marker_Mean\\.Cell\\.", "", .) %>%
  sub("_MC_5_mean$", "", .)

mc12_proteins <- mc12_cols %>%
  sub("^mc_marker_Mean\\.Cell\\.", "", .) %>%
  sub("_MC_12_mean$", "", .)

shared_markers <- intersect(mc5_proteins, mc12_proteins)

length(shared_markers)
shared_markers


comparison_df <- tibble(
  marker = shared_markers,
  
  MC_5_col = paste0(
    "mc_marker_Mean.Cell.",
    shared_markers,
    "_MC_5_mean"
  ),
  
  MC_12_col = paste0(
    "mc_marker_Mean.Cell.",
    shared_markers,
    "_MC_12_mean"
  )
)

comparison_df <- comparison_df %>%
  rowwise() %>%
  mutate(
    MC_5_mean = mean(meta_markers[[MC_5_col]], na.rm = TRUE),
    
    MC_12_mean = mean(meta_markers[[MC_12_col]], na.rm = TRUE),
    
    diff_MC12_minus_MC5 = MC_12_mean - MC_5_mean
  ) %>%
  ungroup() %>%
  arrange(desc(diff_MC12_minus_MC5))

head(comparison_df, 20)

tail(comparison_df, 20)


top_diff <- comparison_df %>%
  slice_max(abs(diff_MC12_minus_MC5), n = 20)

ggplot(
  top_diff,
  aes(
    x = reorder(marker, diff_MC12_minus_MC5),
    y = diff_MC12_minus_MC5
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Markers differing between MC_12 and MC_5",
    x = "Marker",
    y = "MC_12 - MC_5 mean expression"
  )
#####one at a time bio data merge


meta_markers <- read_csv("Features/Run2_all/Manual_classes/act4_metacluster_markers.csv")
names(meta_markers)


meta_markers_mc <- meta_markers %>%
  left_join(
    meta_mapped_only %>%
      select(tma_id, all_of(complete_ish_columns)),
    by = "tma_id"
  )

# Check MC4 is now present
summary(meta_markers_mc$abundance_metacluster_leiden_1.0_MC_4)

# Check IDs in each table
length(unique(meta_markers$tma_id))
length(unique(meta_mapped_only$tma_id))

# How many TMAs matched?
sum(!is.na(meta_markers_mc$abundance_metacluster_leiden_1.0_MC_4))


###too small people in clusters for high low split
mc4_col <- "abundance_metacluster_leiden_1.0_MC_4"

meta_markers_mc <- meta_markers_mc %>%
  mutate(
    MC4_present = ifelse(
      !is.na(.data[[mc4_col]]) & .data[[mc4_col]] > 0,
      "MC4_present", "MC4_absent"
    )
  )

table(meta_markers_mc$MC4_present, useNA = "ifany")

colnames(meta_markers_mc)
marker_cols <- names(meta_markers_mc)[
  grepl("^mc_marker_Mean\\.Cell\\..*_MC_4_mean$", names(meta_markers_mc))
]

length(marker_cols)
head(marker_cols)
# 0) Grouping already set up:
# meta_markers_mc$MC4_present ("MC4_present" / "MC4_absent")

# 1) Select the MC4 marker columns
marker_cols <- names(meta_markers_mc)[
  grepl("^mc_marker_Mean\\.Cell\\..*_MC_4_mean$", names(meta_markers_mc))
]

# 2) Differential markers: MC4_present vs MC4_absent
cmp_df <- meta_markers_mc %>%
  filter(!is.na(MC4_present))

marker_stats <- lapply(marker_cols, function(m) {
  x <- cmp_df[[m]]
  g <- cmp_df$MC4_present
  
  if (all(is.na(x))) {
    return(data.frame(
      marker = m,
      pval   = NA_real_,
      log2FC = NA_real_
    ))
  }
  
  mean_present <- mean(x[g == "MC4_present"], na.rm = TRUE)
  mean_absent  <- mean(x[g == "MC4_absent"],  na.rm = TRUE)
  
  log2FC <- log2(mean_present + 1e-6) - log2(mean_absent + 1e-6)
  
  wt <- try(wilcox.test(x ~ g), silent = TRUE)
  pval <- if (inherits(wt, "try-error")) NA_real_ else wt$p.value
  
  data.frame(
    marker = m,
    pval   = pval,
    log2FC = log2FC
  )
})

marker_stats_df <- do.call(rbind, marker_stats)
marker_stats_df$FDR <- p.adjust(marker_stats_df$pval, method = "BH")
marker_stats_df <- marker_stats_df[order(marker_stats_df$FDR), ]

head(marker_stats_df)


#######plotting


volcano_df <- marker_stats_df %>%
  mutate(
    neg_log10_FDR = -log10(FDR),
    sig = ifelse(FDR < 0.2 & log2FC > 0, "Enriched_in_MC4_present",
                 ifelse(FDR < 0.2 & log2FC < 0, "Enriched_in_MC4_absent", "NS"))
  )

ggplot(volcano_df, aes(x = log2FC, y = neg_log10_FDR, colour = sig)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_colour_manual(values = c(
    "Enriched_in_MC4_present" = "#d73027",
    "Enriched_in_MC4_absent"  = "#4575b4",
    "NS"                      = "grey70"
  )) +
  geom_hline(yintercept = -log10(0.2), linetype = "dashed", colour = "grey50") +
  labs(
    x = "log2FC (MC4_present vs MC4_absent)",
    y = "-log10(FDR)",
    colour = NULL
  ) +
  theme_minimal()

###nothing was significant after FDR cor

##lets see another plot
sel_markers <- c(
  "mc_marker_Mean.Cell.APAF1_MC_4_mean",
  "mc_marker_Mean.Cell.BAX_MC_4_mean",
  "mc_marker_Mean.Cell.BCATENIN_MC_4_mean",
  "mc_marker_Mean.Cell.PCAD_MC_4_mean",
  "mc_marker_Mean.Cell.MUC5_MC_4_mean",
  "mc_marker_Mean.Cell.CASP3CLEAVED_MC_4_mean"
)
plot_df <- meta_markers_mc %>%
  filter(!is.na(MC4_present)) %>%
  select(MC4_present, all_of(sel_markers)) %>%
  pivot_longer(
    cols = all_of(sel_markers),
    names_to = "marker",
    values_to = "expr"
  )

ggplot(plot_df, aes(x = MC4_present, y = expr, fill = MC4_present)) +
  geom_violin(trim = FALSE, colour = NA, alpha = 0.4) +
  geom_boxplot(width = 0.25, outlier.size = 0.6) +
  facet_wrap(~ marker, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("MC4_present" = "#d73027", "MC4_absent" = "#4575b4")) +
  labs(
    x = NULL,
    y = "Expression in MC4 cells (per TMA)",
    title = "Selected MC4 markers by MC4 presence at TMA level"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

############pATHWAY ENRICHMENT GOOD V BAD CLUSTER COHORT
write.csv(best5, "best5.csv", row.names = FALSE)
write.csv(worst5, "worst5.csv", row.names = FALSE)
