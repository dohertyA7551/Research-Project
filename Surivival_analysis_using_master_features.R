library(ggplot2)
library(dplyr)
library(tidyr)
library(survival)
library(coxphf)
library(purrr)
library(glmnet)
library(readr)
library(broom)
library(survminer)
library(stringr)
library(tidyr)
library(tibble)
library(umap)
library(limma)



meta_markers <- read_csv("Features/Run2_all/Manual_classes/act4_metacluster_markers.csv")
names(meta_markers)

master_features <- read_csv("/work_space/files/Features/Run2_all/Manual_classes/master_features.csv")

master_features

clin_data<- read.csv("/work_space/files/spatial/ColonRO1_2020_RCSI_clinical.csv")
p_id <-read.csv("/work_space/files/spatial/master_patient_mapping.csv")


cols <- colnames(master_features)

# --- 1. MEAN EXPRESSION MARKERS ---
# Pattern: MC_<number>_mean or similar _mean$ suffix
expression_cols <- cols[grepl("_mean$|_mean_|MC_.*_mean", cols, ignore.case = TRUE)]
cat("=== MEAN EXPRESSION MARKER COLUMNS ===\n")
cat("Total:", length(expression_cols), "\n")
print(expression_cols)

# --- 2. LEIDEN CLUSTER COLUMNS ---
leiden_cols <- cols[grepl("leiden|Leiden|LEIDEN", cols, ignore.case = TRUE)]
cat("\n=== LEIDEN CLUSTER COLUMNS ===\n")
cat("Total:", length(leiden_cols), "\n")
print(leiden_cols)

####Merge feature with clin data 


# --- STEP 1: Aggregate master_features to patient level ---
# Average all numeric columns across TMA cores per patient
master_patient <- master_features %>%
  group_by(tma_id) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE),
            .groups = "drop")

cat("master_features rows:", nrow(master_features), "\n")
cat("After aggregation (patient level):", nrow(master_patient), "\n")


# --- STEP 2: Check ID formats match before merging ---
cat("\nSample tma_id values:\n")
print(head(master_patient$tma_id, 10))

cat("\nSample SS values from clin_data:\n")
print(head(clin_data$SS, 10))

# --- STEP 3: Merge ---
merged_df <- left_join(master_patient, clin_data, by = c("tma_id" = "SS"))

cat("\nDimensions after merge:", dim(merged_df), "\n")
cat("Patients matched:", sum(!is.na(merged_df$tma_id)), "\n")

# Check for unmatched patients
unmatched <- master_patient$tma_id[!master_patient$tma_id %in% clin_data$SS]
cat("Unmatched TMA IDs:", length(unmatched), "\n")
if(length(unmatched) > 0) print(head(unmatched, 10))


# Then immediately convert yes no's to binary
merged_df$os_event  <- ifelse(tolower(trimws(merged_df$os_event))  == "yes", 1, 0)
merged_df$dss_event <- ifelse(tolower(trimws(merged_df$dss_event)) == "yes", 1, 0)
merged_df$dfs_event <- ifelse(tolower(trimws(merged_df$dfs_event)) == "yes", 1, 0)

# Also ensure time columns are numeric
merged_df$os_months  <- as.numeric(merged_df$os_months)
merged_df$dss_months <- as.numeric(merged_df$dss_months)
merged_df$dfs_months <- as.numeric(merged_df$dfs_months)

# Verify all looks good
cat("os_event:\n");  print(table(merged_df$os_event,  useNA = "always"))
cat("dss_event:\n"); print(table(merged_df$dss_event, useNA = "always"))
cat("dfs_event:\n"); print(table(merged_df$dfs_event, useNA = "always"))
###### Lasso for var selection


# -------------------------------------------------------
# SETUP
# -------------------------------------------------------

# Get all leiden abundance columns
leiden_cols <- colnames(merged_df)[grepl("abundance_metacluster_leiden", colnames(merged_df))]
cat("Number of metacluster abundance columns:", length(leiden_cols), "\n")

#LASSO COX FUNCTION

run_lasso_cox <- function(data, time_col, event_col, endpoint_name) {
  
  df_clean <- data[!is.na(data[[time_col]]) & !is.na(data[[event_col]]) & data[[time_col]] > 0, ]
  cat("\n========================================\n")
  cat("ENDPOINT:", endpoint_name, "| N =", nrow(df_clean), "\n")
  cat("========================================\n")
  
  X      <- as.matrix(df_clean[, leiden_cols])
  y      <- Surv(df_clean[[time_col]], df_clean[[event_col]])
  
  set.seed(42)
  cv_fit <- cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10, cox.ties = "efron")
  
  cat("Lambda min:", round(cv_fit$lambda.min, 5), "\n")
  cat("Lambda 1se:", round(cv_fit$lambda.1se, 5), "\n")
  
  coefs    <- coef(cv_fit, s = "lambda.min")
  selected <- data.frame(
    metacluster = rownames(coefs),
    coefficient = as.numeric(coefs)
  ) %>%
    filter(coefficient != 0) %>%
    mutate(
      HR        = round(exp(coefficient), 3),
      direction = ifelse(coefficient > 0, "POOR prognosis", "GOOD prognosis")
    ) %>%
    arrange(desc(abs(coefficient)))
  
  cat("\nSelected metaclusters at lambda.min (n =", nrow(selected), "):\n")
  print(selected)
  
  #write.csv(selected, paste0("lasso_cox_", endpoint_name, ".csv"), row.names = FALSE)
  
  plot(cv_fit, main = paste("LASSO CV -", endpoint_name))
  
  return(list(cv_fit = cv_fit, selected = selected))
}

# -------------------------------------------------------
# REMOVE ZERO-TIME PATIENTS
# -------------------------------------------------------

merged_df <- merged_df %>%
  filter(is.na(dss_months) | dss_months > 0) %>%
  filter(is.na(dfs_months) | dfs_months > 0)

cat("Rows remaining after removing zero-time patients:", nrow(merged_df), "\n")

# -------------------------------------------------------
# RUN DSS AND DFS ONLY
# -------------------------------------------------------

res_dss <- run_lasso_cox(merged_df, "dss_months", "dss_event", "DSS")
res_dfs <- run_lasso_cox(merged_df, "dfs_months", "dfs_event", "DFS")

# -------------------------------------------------------
# OVERLAP BETWEEN DSS AND DFS
# -------------------------------------------------------

cat("\n========================================\n")
cat("CONSISTENT ACROSS DSS AND DFS\n")
cat("========================================\n")

common_dss_dfs <- intersect(res_dss$selected$metacluster, res_dfs$selected$metacluster)

if (length(common_dss_dfs) > 0) {
  common_df <- res_dss$selected %>%
    filter(metacluster %in% common_dss_dfs) %>%
    select(metacluster, DSS_HR = HR, DSS_direction = direction) %>%
    left_join(res_dfs$selected %>% 
                select(metacluster, DFS_HR = HR, DFS_direction = direction), 
              by = "metacluster")
  # ---- FIX METACLUSTER NAMES HERE ----
  library(stringr)
  
  common_df <- common_df %>%
    mutate(
      metacluster = str_replace(
        metacluster,
        "^abundance_metacluster_leiden_1\\.0_abundance_metacluster_leiden_1\\.0_",
        "abundance_metacluster_leiden_1.0_"
      )
    )
  
  # Derive good / poor sets using cleaned names
  good_prog <- common_df %>%
    filter(DSS_direction == "GOOD prognosis",
           DFS_direction == "GOOD prognosis")
  
  poor_prog <- common_df %>%
    filter(DSS_direction == "POOR prognosis",
           DFS_direction == "POOR prognosis")
  
  print(common_df)

 # write.csv(common_df, "lasso_cox_consistent_DSS_DFS.csv", row.names = FALSE)

###filter out which ones are in agreement
# Metaclusters consistently GOOD prognosis across both DSS and DFS
good_prog <- common_df %>%
  filter(DSS_direction == "GOOD prognosis" & DFS_direction == "GOOD prognosis")

# Metaclusters consistently POOR prognosis across both DSS and DFS
poor_prog <- common_df %>%
  filter(DSS_direction == "POOR prognosis" & DFS_direction == "POOR prognosis")

# Discordant
discordant <- common_df %>%
  filter(DSS_direction != DFS_direction)

cat("=== CONSISTENTLY GOOD PROGNOSIS ===\n"); print(good_prog)
cat("\n=== CONSISTENTLY POOR PROGNOSIS ===\n"); print(poor_prog)
cat("\n=== DISCORDANT (exclude) ===\n"); print(discordant)





########seeing if cluster splits actually retian bio meaning
####### 
'
# Remove the _mean suffix to match merged_df column names
good_cols <- gsub("_mean$", "", good_cols)
poor_cols <- gsub("_mean$", "", poor_cols)

# Verify
print(good_cols %in% colnames(merged_df))
print(poor_cols %in% colnames(merged_df))'

######if re running this instead###
# These are already the correct column names — no gsub needed
good_cols <- good_prog$metacluster
poor_cols <- poor_prog$metacluster

# Verify they exist in merged_df
print(good_cols %in% colnames(merged_df))
print(poor_cols %in% colnames(merged_df))
nrow(good_prog)  # how many consistently GOOD?
nrow(poor_prog)  # how many consistently POOR?

# --- Composite score (unweighted for now) ---
merged_df$good_score <- rowMeans(merged_df[, good_cols], na.rm = TRUE)
merged_df$poor_score <- rowMeans(merged_df[, poor_cols], na.rm = TRUE)
merged_df$prog_score <- merged_df$good_score - merged_df$poor_score

# --- Tertile split ---
tertiles <- quantile(merged_df$prog_score, probs = c(1/3, 2/3), na.rm = TRUE)

merged_df$prognosis_group <- case_when(
  merged_df$prog_score >= tertiles[2] ~ "GOOD",
  merged_df$prog_score <= tertiles[1] ~ "POOR",
  TRUE                                ~ "INTERMEDIATE"
)
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate

print(table(merged_df$prognosis_group))

}

# Find survival-related columns
print(grep("DSS|DFS|OS|surv|time|event|status|death|follow", 
           colnames(merged_df), value = TRUE, ignore.case = TRUE))

# KM validation
# Subset to GOOD and POOR only
km_df <- merged_df %>% filter(prognosis_group %in% c("GOOD", "POOR"))

surv_obj <- Surv(time = km_df$dss_months, event = km_df$dss_event)

km_fit <- survfit(surv_obj ~ prognosis_group, data = km_df)

ggsurvplot(
  km_fit,
  data = km_df,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("#2ecc71", "#e74c3c"),
  legend.labs = c("GOOD", "POOR"),
  title = "KM: Prognostic Groups by Metacluster Score",
  xlab = "Time (months)"
)

# Check how much the prog_score actually varies
summary(merged_df$prog_score)
hist(merged_df$prog_score, breaks = 30, main = "Distribution of Prognostic Score")

# Find the valley between the two peaks
hist(merged_df$prog_score, breaks = 30, main = "Distribution of Prognostic Score")

# Get a density estimate to find the trough
d <- density(merged_df$prog_score, na.rm = TRUE)
plot(d, main = "Density of Prognostic Score")

# What value sits between the two peaks?
# Eyeball it from the plot, or find local minima:
# Compute density
d <- density(merged_df$prog_score, na.rm = TRUE)

# Find the valley (local minimum) between the two peaks
# Set the interval to cover the range between the two peaks
valley <- optimize(approxfun(d$x, d$y), 
                   interval = c(min(d$x), max(d$x)))$minimum

cat("Split point (valley):", valley, "\n")

# Plot to verify the split point makes sense
plot(d, main = "Density of Prognostic Score")
abline(v = valley, col = "blue", lty = 2)
merged_df$prognosis_group <- case_when(
  merged_df$prog_score >= valley ~ "GOOD",
  TRUE                           ~ "POOR"
)

print(table(merged_df$prognosis_group))


saveRDS(common_df, "common_dss_dfs_metaclusters.rds")
saveRDS(good_prog, "good_prog_metaclusters.rds")
saveRDS(poor_prog, "poor_prog_metaclusters.rds")
saveRDS(discordant, "discordant_metaclusters.rds")
#####saving objects for future use
saveRDS(merged_df, "merged_df_with_valley_groups.rds")
merged_df <- readRDS("merged_df_with_valley_groups.rds")

patient_groups <- merged_df %>%
  dplyr::select(tma_id, prognosis_group)

write.csv(
  patient_groups,
  "patient_group_assignments_valley.csv",
  row.names = FALSE
)
grep("group|prog|valley", colnames(merged_df), value = TRUE)


str(merged_df)
















#######

km_df <- merged_df %>% filter(prognosis_group %in% c("GOOD", "POOR"))

surv_obj <- Surv(time = km_df$dss_months, event = km_df$dss_event)
km_fit <- survfit(surv_obj ~ prognosis_group, data = km_df)

ggsurvplot(
  km_fit,
  data = km_df,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("#2ecc71", "#e74c3c"),
  legend.labs = c("GOOD", "POOR"),
  title = "KM: Natural Split by Prognostic Score",
  xlab = "Time (months)"
)


survdiff(
  Surv(dss_months, dss_event) ~ prognosis_group,
  data = km_df
)

#########Can skip above step#####'
###below code is if you need to equallise groups a bit
# 1) Get good metacluster columns directly
good_abund_cols <- good_prog$metacluster   # already proper column names
poor_abund_cols <- poor_prog$metacluster

# Sanity check
print(good_abund_cols)
print(good_abund_cols %in% colnames(merged_df))
print(poor_abund_cols %in% colnames(merged_df))
stopifnot(all(good_abund_cols %in% colnames(merged_df)))
stopifnot(all(poor_abund_cols %in% colnames(merged_df)))

##remove nas
lasso_df <- merged_df %>%
  filter(
    !is.na(dss_months),
    !is.na(dss_event)
  )

nrow(lasso_df)
# 2) Design matrix for LASSO: use lasso_df, not merged_df
X_good <- as.matrix(lasso_df[, good_abund_cols])

# 3) Survival response: also from lasso_df
y_surv <- Surv(time = lasso_df$dss_months, event = lasso_df$dss_event)

# 4) LASSO Cox
set.seed(123)
cvfit_good <- cv.glmnet(
  x      = X_good,
  y      = y_surv,
  family = "cox",
  alpha  = 1,
  nfolds = 10
)

warnings() #guess i will ignore
##then rank coef;
coef_good <- coef(cvfit_good, s = "lambda.min")  # or "lambda.min"
coef_good_df <- data.frame(
  feature = rownames(coef_good),
  coef    = as.numeric(coef_good)
) %>%
  dplyr::filter(coef != 0)

coef_good_df

coef_ranked <- coef_good_df %>%
  mutate(
    abs_coef = abs(coef)
  ) %>%
  arrange(desc(abs_coef)) %>%
  mutate(
    rank = row_number()
  )

coef_ranked

############# get stats values to report


# Good prognosis MCs from your second LASSO result
good_ids <- coef_ranked %>%
  slice(1:3) %>%
  pull(feature) %>%
  str_extract("MC_\\d+")

# Poor prognosis MCs from your poor_prog table
poor_ids <- poor_prog %>%
  pull(metacluster) %>%
  str_extract("MC_\\d+")

# Combine selected MCs
core_mc_ids <- c(good_ids, poor_ids) %>%
  unique()

core_mc_ids

### ph ratio and p vals

get_cox_for_mc <- function(mc, time_col, event_col, endpoint) {
  
  feature_col <- grep(
    paste0("abundance_metacluster_leiden.*", mc, "$"),
    colnames(merged_df),
    value = TRUE
  )
  
  if (length(feature_col) == 0) return(NULL)
  
  fit <- coxph(
    as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ `", feature_col[1], "`")),
    data = merged_df
  )
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(
      endpoint = endpoint,
      metacluster = mc,
      feature = feature_col[1]
    ) %>%
    select(endpoint, metacluster, feature, estimate, conf.low, conf.high, p.value)
}

cox_report_DSS <- map_dfr(
  core_mc_ids,
  get_cox_for_mc,
  time_col = "dss_months",
  event_col = "dss_event",
  endpoint = "DSS"
)

cox_report_DFS <- map_dfr(
  core_mc_ids,
  get_cox_for_mc,
  time_col = "dfs_months",
  event_col = "dfs_event",
  endpoint = "DFS"
)

cox_report <- bind_rows(cox_report_DSS, cox_report_DFS)

cox_report
cox_report <- cox_report %>%
  rename(
    HR = estimate,
    CI_lower = conf.low,
    CI_upper = conf.high,
    p_value = p.value
  )

write.csv(cox_report, "cox_report_selected_metaclusters.csv", row.names = FALSE)

####Finally getting onto de protein exp
# Top 3 good metaclusters from LASSO
top_k <- 3
good_top3 <- coef_ranked %>%
  slice(1:top_k) %>%
  mutate(metacluster = str_extract(feature, "MC_\\d+"))

good_top3_ids <- good_top3$metacluster

# Poor metaclusters from your earlier filter
poor_ids <- str_extract(poor_prog$metacluster, "MC_\\d+")

good_top3_ids
poor_ids

core_mc_ids <- c(good_top3_ids, poor_ids)
core_mc_ids

###subset out marker exp for relevant columns
# Build regex for the 6 core metaclusters
core_pattern <- paste(core_mc_ids, collapse = "|")  # e.g. "MC_55|MC_16|MC_14|MC_0|MC_1|MC_... "

core_marker_cols <- grep(
  paste0("^mc_marker_Mean\\.Cell\\..*_(", core_pattern, ")_mean$"),
  colnames(merged_df),
  value = TRUE
)

length(core_marker_cols)
head(core_marker_cols)



####pca with MC 38 instead of 47



# Extract clean MC IDs and prognosis group from top6 table
mc_group_table <- top6_metaclusters %>%
  mutate(
    mc_id = str_extract(metacluster, "MC_\\d+"),
    group = case_when(
      Prognosis == "Good" ~ "GOOD_MC",
      Prognosis == "Poor" ~ "POOR_MC",
      TRUE ~ NA_character_
    )
  ) %>%
  select(mc_id, group)

# PCA scores with group labels
pca_scores <- as.data.frame(pca_res$x) %>%
  rownames_to_column("metacluster") %>%
  left_join(mc_group_table, by = c("metacluster" = "mc_id"))

# Check it worked
pca_scores

# PCA plot
ggplot(pca_scores, aes(PC1, PC2, color = group, label = metacluster)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  theme_classic()

# Add group labels to marker expression table
cluster_level_grp <- cluster_level %>%
  left_join(mc_group_table, by = c("metacluster" = "mc_id"))

# Marker-level linear model: GOOD vs POOR
fit_marker_lm <- function(mk) {
  
  sub <- cluster_level_grp %>%
    filter(marker == mk) %>%
    mutate(group = factor(group, levels = c("POOR_MC", "GOOD_MC")))
  
  if (nrow(sub) < 3) return(NULL)
  
  fit <- lm(mean_expr ~ group, data = sub)
  
  coef_group <- broom::tidy(fit) %>%
    filter(term == "groupGOOD_MC") %>%
    transmute(
      marker = mk,
      estimate = estimate,
      std_error = std.error,
      stat = statistic,
      p_value = p.value
    )
  
  means <- sub %>%
    group_by(group) %>%
    summarise(mean_expr = mean(mean_expr), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = group,
      values_from = mean_expr,
      names_prefix = "mean_"
    )
  
  bind_cols(coef_group, means)
}

marker_lm <- map_dfr(unique(cluster_level_grp$marker), fit_marker_lm) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

head(marker_lm, 20)
###PCA attempt 2 ends here 


### prep for pca
core_pattern <- paste(core_mc_ids, collapse = "|")

core_marker_cols <- grep(
  paste0("^mc_marker_Mean\\.Cell\\..*_(", core_pattern, ")_mean$"),
  colnames(merged_df),
  value = TRUE
)

# 2) Build info table for each feature
core_info <- tibble(feature = core_marker_cols) %>%
  mutate(
    marker      = str_extract(feature, "(?<=Mean\\.Cell\\.).*(?=_MC_)"),
    metacluster = str_extract(feature, "MC_\\d+")
  )

# 3) Compute cluster-level mean per (marker, metacluster) across patients
core_mat <- merged_df[, core_marker_cols]

cluster_level <- core_info %>%
  mutate(mean_expr = colMeans(core_mat, na.rm = TRUE))

# 4) Wide matrix: rows = metaclusters, columns = markers
cluster_wide <- cluster_level %>%
  select(metacluster, marker, mean_expr) %>%
  pivot_wider(
    names_from  = marker,
    values_from = mean_expr
  ) %>%
  arrange(metacluster)

# Row names are metaclusters
cluster_mat <- cluster_wide %>%
  column_to_rownames("metacluster") %>%
  as.matrix()

cluster_mat[1:3, 1:5]  # quick sanity check



# Scale markers (columns) before PCA
cluster_scaled <- scale(cluster_mat)

pca_res <- prcomp(cluster_scaled, center = FALSE, scale. = FALSE)

# Extract PC scores for plotting
pca_scores <- as.data.frame(pca_res$x[, 1:2])
pca_scores$metacluster <- rownames(pca_scores)
pca_scores$group <- ifelse(
  pca_scores$metacluster %in% good_top3_ids,
  "GOOD_MC",
  "POOR_MC"
)

# Example ggplot
library(ggplot2)

ggplot(pca_scores, aes(PC1, PC2, color = group, label = metacluster)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  theme_classic()




fit_marker_lm <- function(mk) {
  sub <- cluster_level_grp %>%
    filter(marker == mk)
  
  # make sure it's a factor with POOR as reference (so coef is GOOD - POOR)
  sub <- sub %>%
    mutate(group = factor(group, levels = c("POOR_MC", "GOOD_MC")))
  
  if (nrow(sub) < 3) return(NULL)  # should be 6 rows, but just in case
  
  fit <- lm(mean_expr ~ group, data = sub)
  
  # tidy and keep the group coefficient
  coef_group <- broom::tidy(fit) %>%
    filter(term == "groupGOOD_MC") %>%
    transmute(
      marker   = mk,
      estimate = estimate,  # effect of GOOD vs POOR
      std_error = std.error,
      stat     = statistic,
      p_value  = p.value
    )
  
  # add group-wise means for context
  means <- sub %>%
    group_by(group) %>%
    summarise(mean_expr = mean(mean_expr), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = group,
      values_from = mean_expr,
      names_prefix = "mean_"
    )
  
  dplyr::bind_cols(coef_group, means)
}

marker_lm <- map_dfr(unique(cluster_level_grp$marker), fit_marker_lm) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

head(marker_lm, 20)



 ####plotting done### look at exp


# 2.1 select marker columns for the 6 metaclusters
core_pattern <- paste(core_mc_ids, collapse = "|")

core_marker_cols <- grep(
  paste0("^mc_marker_Mean\\.Cell\\..*_(", core_pattern, ")_mean$"),
  colnames(merged_df),
  value = TRUE
)

# 2.2 long format: patient × marker × metacluster
marker_long <- merged_df %>%
  select(prognosis_group, all_of(core_marker_cols)) %>%
  mutate(patient_id = row_number()) %>%  # replace with your real patient ID if you prefer
  pivot_longer(
    cols      = all_of(core_marker_cols),
    names_to  = "feature",
    values_to = "expr"
  ) %>%
  mutate(
    marker      = str_extract(feature, "(?<=Mean\\.Cell\\.).*(?=_MC_)"),
    metacluster = str_extract(feature, "MC_\\d+")
  )

# 2.3 aggregate across the 6 metaclusters per patient per marker (mean; sum also possible)
patient_marker <- marker_long %>%
  group_by(patient_id, prognosis_group, marker) %>%
  summarise(
    expr = mean(expr, na.rm = TRUE),
    .groups = "drop"
  )

# 2.4 wide: patients × markers (for modelling), and a matching group vector
patient_marker_mat <- patient_marker %>%
  select(patient_id, marker, expr) %>%
  pivot_wider(
    names_from  = marker,
    values_from = expr
  ) %>%
  arrange(patient_id)

marker_matrix <- as.matrix(patient_marker_mat[, -1])  # numeric matrix
group_vec     <- patient_marker %>%
  distinct(patient_id, prognosis_group) %>%
  arrange(patient_id) %>%
  pull(prognosis_group)

###simple DE...apparently#
library(limma)

group_factor <- factor(group_vec, levels = c("POOR", "GOOD"))  # POOR = reference
design <- model.matrix(~ group_factor)  # intercept + groupGOOD

fit <- lmFit(t(marker_matrix), design)  # genes/markers in rows
fit <- eBayes(fit)

marker_de <- topTable(
  fit,
  coef   = "group_factorGOOD",
  number = Inf,
  sort.by = "P"
)

head(marker_de) ##IMPORTANT TO YOUR LIFE
marker_de_out <- marker_de %>%
  as.data.frame() %>%
  tibble::rownames_to_column("marker")

write.csv(
  marker_de_out,
  file = "marker_de_patient_level_limma.csv",
  row.names = FALSE
)
library(dplyr)
library(readr)
library(tibble)


###fromcsv if coming back to it
# Read patient-level limma results from saved CSV


library(dplyr)
library(readr)
library(ggplot2)

# 1. Load patient-level limma proteomic results
marker_de_tbl <- readr::read_csv(
  "marker_de_patient_level_limma.csv",
  show_col_types = FALSE
)

# 2. Tidy metacluster-level proteomic results
marker_lm_tidy <- marker_lm %>%
  dplyr::select(
    marker,
    estimate,
    p_value,
    p_adj,
    mean_POOR_MC,
    mean_GOOD_MC
  )

# 3. Join patient-level and metacluster-level results
marker_cross <- marker_de_tbl %>%
  dplyr::inner_join(
    marker_lm_tidy,
    by = "marker"
  )

# 4. Add direction labels
marker_cross <- marker_cross %>%
  dplyr::mutate(
    dir_patient = dplyr::case_when(
      logFC_patient > 0 ~ "GOOD_high",
      logFC_patient < 0 ~ "POOR_high",
      TRUE ~ "no_change"
    ),
    dir_cluster = dplyr::case_when(
      estimate > 0 ~ "GOOD_high",
      estimate < 0 ~ "POOR_high",
      TRUE ~ "no_change"
    ),
    concordance = dplyr::case_when(
      dir_patient == dir_cluster ~ "concordant",
      dir_patient != dir_cluster ~ "discordant"
    )
  )

# 5. Check results
table(marker_cross$dir_patient, marker_cross$dir_cluster)
table(marker_cross$concordance)

# 6. Extract concordant good and poor markers
good_enriched <- marker_cross %>%
  dplyr::filter(
    dir_patient == "GOOD_high",
    dir_cluster == "GOOD_high"
  ) %>%
  dplyr::arrange(adj.P.Val_patient, p_adj)

poor_enriched <- marker_cross %>%
  dplyr::filter(
    dir_patient == "POOR_high",
    dir_cluster == "POOR_high"
  ) %>%
  dplyr::arrange(adj.P.Val_patient, p_adj)

discordant_markers <- marker_cross %>%
  dplyr::filter(concordance == "discordant") %>%
  dplyr::arrange(adj.P.Val_patient, p_adj)

# 7. View
head(good_enriched, 20)
head(poor_enriched, 20)
head(discordant_markers, 20)



































###cross val of patient level with meta cluster level
marker_de_tbl <- marker_de %>%
  as.data.frame() %>%
  tibble::rownames_to_column("marker") %>%
  select(marker, logFC, AveExpr, P.Value, adj.P.Val) %>%
  rename(
    logFC_patient   = logFC,
    AveExpr_patient = AveExpr,
    P.Value_patient = P.Value,
    adj.P.Val       = adj.P.Val   # keep this name if you like, or rename to adj.P.Val_patient
  )
#tidy cluster results#
marker_lm_tidy <- marker_lm %>%
  select(marker, estimate, p_value, p_adj, mean_POOR_MC, mean_GOOD_MC)

#join patient and cluster res
marker_cross <- marker_de_tbl %>%
  inner_join(
    marker_lm_tidy,
    by = "marker"
  )

marker_cross <- marker_cross %>%
  mutate(
    dir_patient = ifelse(logFC_patient > 0, "GOOD_high", "POOR_high"),
    dir_cluster = ifelse(estimate       > 0, "GOOD_high", "POOR_high")
  )
# quick look at top rows
head(marker_cross, 20)
##extract good enriched
good_enriched <- marker_cross %>%
  filter(
    dir_patient == "GOOD_high",
    dir_cluster == "GOOD_high"
  ) %>%
  arrange(adj.P.Val, p_adj)

head(good_enriched, 20)

##extract bad enriched
poor_enriched <- marker_cross %>%
  filter(
    dir_patient == "POOR_high",
    dir_cluster == "POOR_high"
  ) %>%
  arrange(adj.P.Val, p_adj)

head(poor_enriched, 20)


##save to csv
# assuming good_enriched and poor_enriched already defined and arranged

write.csv(
  good_enriched,
  file = "good_enriched_all.csv",
  row.names = FALSE
)

write.csv(
  poor_enriched,
  file = "poor_enriched_all.csv",
  row.names = FALSE
)

####pulling out patient groups
# Your real patient ID column appears to be tma_id
patient_groups <- merged_df %>%
  select(tma_id, prognosis_group) %>%
  filter(!is.na(prognosis_group))

print(patient_groups)
print(table(patient_groups$prognosis_group))





########Pulling out info post inital run 02/06/26


marker_de_tbl <- read_csv("marker_de_patient_level_limma.csv")

sig_markers <- marker_de_tbl %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(
    enriched_in = case_when(
      logFC > 0 ~ "GOOD prognosis",
      logFC < 0 ~ "POOR prognosis"
    )
  ) %>%
  arrange(adj.P.Val)

sig_good <- sig_markers %>%
  filter(enriched_in == "GOOD prognosis")

sig_poor <- sig_markers %>%
  filter(enriched_in == "POOR prognosis")

write.csv(sig_markers, "significant_limma_markers_all.csv", row.names = FALSE)
write.csv(sig_good, "significant_good_enriched_markers.csv", row.names = FALSE)
write.csv(sig_poor, "significant_poor_enriched_markers.csv", row.names = FALSE)

head(sig_markers, 30)
library(dplyr)
library(readr)

marker_de_tbl <- read_csv("marker_de_patient_level_limma.csv")

marker_de_tbl %>%
  filter(adj.P.Val < 0.05) %>%
  mutate(
    enriched_in = ifelse(
      logFC > 0,
      "GOOD prognosis",
      "POOR prognosis"
    )
  ) %>%
  select(
    marker,
    enriched_in,
    logFC,
    P.Value,
    adj.P.Val
  ) %>%
  arrange(adj.P.Val)

# GOOD prognosis markers
marker_de_tbl %>%
  filter(adj.P.Val < 0.05,
         logFC > 0) %>%
  select(marker, logFC, P.Value, adj.P.Val) %>%
  arrange(adj.P.Val)

# POOR prognosis markers
marker_de_tbl %>%
  filter(adj.P.Val < 0.05,
         logFC < 0) %>%
  select(marker, logFC, P.Value, adj.P.Val) %>%
  arrange(adj.P.Val)
print(sig_poor, n = Inf)
# Save full results tables
write_csv(
  sig_good,
  "significant_good_prognosis_markers.csv"
)

write_csv(
  sig_poor,
  "significant_poor_prognosis_markers.csv"
)

# Save marker-only lists for downstream analyses
write_csv(
  tibble(marker = sig_good$marker),
  "significant_good_prognosis_marker_list.csv"
)

write_csv(
  tibble(marker = sig_poor$marker),
  "significant_poor_prognosis_marker_list.csv"
)




#########Day 2 pulling out info 03/06/26
lasso_consistent <-read_csv("Features/lasso_cox_consistent_DSS_DFS.csv")

names(lasso_consistent)
head(lasso_consistent)
dim(lasso_consistent)
names(cox_report)




# Pull directly from your Cox results
forest_df <- cox_report %>%
  mutate(
    metacluster = gsub("MC_", "MC ", metacluster),
    endpoint = factor(endpoint, levels = c("DSS", "DFS"))
  )

# Optional: remove degenerate estimates that break log-scale plots
forest_df_clean <- forest_df %>%
  filter(
    HR > 0,
    CI_lower > 0,
    is.finite(CI_upper)
  )

# Forest plot
ggplot(
  forest_df_clean,
  aes(
    x = HR,
    y = reorder(metacluster, HR)
  )
) +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = CI_lower, xmax = CI_upper),
    height = 0.2
  ) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  facet_wrap(~ endpoint) +
  theme_classic() +
  labs(
    title = "Hazard ratios for selected metaclusters",
    x = "Hazard ratio (log scale)",
    y = "Metacluster"
  )


#########################
'comparing proteomic abundace between prognosis groups at a cell type level'



baseline_mean_cols <- grep(
  "^baseline.*mean",
  colnames(master_features),
  value = TRUE,
  ignore.case = TRUE
)

cat("Total baseline mean columns:", length(baseline_mean_cols), "\n\n")
print(baseline_mean_cols)

merged_df <- readRDS("merged_df_with_valley_groups.rds")

# Quick check
cat("Rows:", nrow(merged_df), "\n")
cat("Cols:", ncol(merged_df), "\n")
cat("Prognosis groups:\n")
print(table(merged_df$prognosis_group))

### get column names
# Check exact naming pattern for each cell type

# Print ALL matches - no truncation
nonie_cols <- grep("NonImmuneEpithelium", colnames(merged_df), 
                   value = TRUE, ignore.case = TRUE)

cat("Total columns:", length(nonie_cols), "\n\n")
writeLines(nonie_cols)  # prints every single one, no truncation

# Cancer cells
cancer_cols <- grep("Mean\\.Cell.*NonImmuneEpithelium|NonImmuneEpithelium.*Mean\\.Cell", 
                    colnames(merged_df), value = TRUE, ignore.case = TRUE)

# Stromal cells
stroma_cols <- grep("Mean\\.Cell.*NonImmuneStroma|NonImmuneStroma.*Mean\\.Cell", 
                    colnames(merged_df), value = TRUE, ignore.case = TRUE)

# Immune cells (all subtypes)
immune_cols <- grep("Mean\\.Cell.*(cytT|helperT|regT|OtherImmune|\\bT\\b|\\bB\\b)|
                    (cytT|helperT|regT|OtherImmune).*Mean\\.Cell", 
                    colnames(merged_df), value = TRUE, ignore.case = TRUE)

cat("Cancer cols:", length(cancer_cols), "\n"); writeLines(cancer_cols)
cat("\nStroma cols:", length(stroma_cols), "\n"); writeLines(stroma_cols)
cat("\nImmune cols:", length(immune_cols), "\n"); writeLines(immune_cols)



######now cells are sorted out into groups... look at DE for each cell type
# How many patients total have a prognosis group?
cat("Total patients with prognosis group:\n")
print(table(merged_df$prognosis_group, useNA = "always"))

# How many have ANY immune column data (not all NA)?
immune_data <- merged_df %>%
  filter(prognosis_group %in% c("GOOD", "POOR")) %>%
  dplyr::select(tma_id, prognosis_group, all_of(immune_cols))

cat("\nTotal patients with prognosis group:", nrow(immune_data), "\n")

# Count NAs per patient across immune cols
immune_data$n_missing <- rowSums(is.na(immune_data[, immune_cols]))
immune_data$n_present <- rowSums(!is.na(immune_data[, immune_cols]))

cat("\nNA summary per patient:\n")
print(summary(immune_data$n_missing))

cat("\nPatients with ALL immune cols missing:", 
    sum(immune_data$n_present == 0), "\n")

cat("\nPatients with at least 1 immune col present:", 
    sum(immune_data$n_present > 0), "\n")

cat("\nPatients with ALL immune cols present (complete.cases):", 
    sum(immune_data$n_missing == 0), "\n")

# Filter immune cols to only mean, high, low
immune_cols_filtered <- immune_cols[grepl("mean|high|low", 
                                          immune_cols, 
                                          ignore.case = TRUE)]

cat("Total immune cols (all):", length(immune_cols), "\n")
cat("After filtering to mean/high/low:", length(immune_cols_filtered), "\n")
writeLines(immune_cols_filtered)

cancer_cols_filtered <- cancer_cols[grepl("mean|high|low", 
                                          cancer_cols, 
                                          ignore.case = TRUE)]

stroma_cols_filtered <- stroma_cols[grepl("mean|high|low", 
                                          stroma_cols, 
                                          ignore.case = TRUE)]

cat("Cancer cols (mean/high/low):", length(cancer_cols_filtered), "\n")
cat("Stroma cols (mean/high/low):", length(stroma_cols_filtered), "\n")


# Patients with prognosis group
immune_check <- merged_df %>%
  filter(prognosis_group %in% c("GOOD", "POOR")) %>%
  dplyr::select(tma_id, prognosis_group, all_of(immune_cols_filtered))

total_immune_cols <- length(immune_cols_filtered)

# Count present (non-NA) columns per patient
immune_check$n_present <- rowSums(!is.na(immune_check[, immune_cols_filtered]))
immune_check$pct_present <- immune_check$n_present / total_immune_cols * 100

# Summary
cat("Total immune cols:", total_immune_cols, "\n")
cat("Total patients:", nrow(immune_check), "\n\n")

cat("Patients with >50% immune cols present:", 
    sum(immune_check$pct_present > 50), "\n")

cat("Patients with >75% immune cols present:", 
    sum(immune_check$pct_present > 75), "\n")

cat("Patients with 100% complete:", 
    sum(immune_check$pct_present == 100), "\n")

# Distribution
hist(immune_check$pct_present, 
     breaks = 20,
     main = "% Immune columns present per patient",
     xlab = "% columns present")



#####Immune cells will need to be done individually


#cancer and stroma first
library(limma)
library(tibble)
library(dplyr)

# Check what rownames look like BEFORE cleaning in the function
# Run this outside the function on your cancer cols
test_mat <- t(as.matrix(
  merged_df[1:5, cancer_cols_filtered]
))

cat("Raw rownames (first 10):\n")
writeLines(rownames(test_mat)[1:10])

# -------------------------------------------------------
# LIMMA DE: GOOD vs POOR — Cancer and Stroma cells
# -------------------------------------------------------

run_de <- function(feature_cols, celltype_name) {
  
  cat("\n========================================\n")
  cat("Cell type:", celltype_name, "\n")
  
  # Subset to GOOD/POOR patients with these columns
  de_df <- merged_df %>%
    filter(prognosis_group %in% c("GOOD", "POOR")) %>%
    dplyr::select(tma_id, prognosis_group, all_of(feature_cols)) %>%
    filter(complete.cases(.))
  
  cat("Patients:", nrow(de_df), "| Markers:", length(feature_cols), "\n")
  
  # Expression matrix: markers x patients
  expr_mat <- t(as.matrix(de_df[, feature_cols]))
  
  # Clean marker names
  rownames(expr_mat) <- feature_cols %>%
    str_remove(".*Mean\\.Cell\\.") %>%
    str_remove(paste0("_", celltype_name, ".*$")) %>%
    str_remove("_mean$|_high$|_low$")
  
  # Design: POOR = reference, GOOD = effect
  group_factor <- factor(de_df$prognosis_group, levels = c("POOR", "GOOD"))
  design <- model.matrix(~ group_factor)
  
  fit <- lmFit(expr_mat, design)
  fit <- eBayes(fit)
  
  topTable(fit, coef = "group_factorGOOD", number = Inf, sort.by = "P") %>%
    as.data.frame() %>%
    rownames_to_column("marker") %>%
    mutate(
      cell_type   = celltype_name,
      enriched_in = case_when(
        logFC > 0 ~ "GOOD prognosis",
        logFC < 0 ~ "POOR prognosis"
      )
    )
}

# -------------------------------------------------------
# RUN
# -------------------------------------------------------

de_cancer <- run_de(cancer_cols_filtered, "Cancer")
de_stroma <- run_de(stroma_cols_filtered, "Stroma")

# -------------------------------------------------------
# COMBINE AND SAVE
# -------------------------------------------------------

de_cancer_stroma <- bind_rows(de_cancer, de_stroma)

write.csv(de_cancer_stroma,
          "de_cancer_stroma_GOOD_vs_POOR.csv",
          row.names = FALSE)

# -------------------------------------------------------
# RESULTS SUMMARY
# -------------------------------------------------------

# Count significant hits per cell type
de_cancer_stroma %>%
  filter(adj.P.Val < 0.05) %>%
  count(cell_type, enriched_in) %>%
  tidyr::pivot_wider(names_from = enriched_in, 
                     values_from = n, 
                     values_fill = 0) %>%
  print()

# Top 10 per cell type
de_cancer_stroma %>%
  filter(adj.P.Val < 0.05) %>%
  group_by(cell_type) %>%
  slice_min(adj.P.Val, n = 10) %>%
  dplyr::select(cell_type, marker, logFC, adj.P.Val, enriched_in) %>%
  print(n = Inf)


# -------------------------------------------------------
# EXTRACT: Top proteins higher in GOOD prognosis
# -------------------------------------------------------

good_cancer <- de_cancer %>%
  filter(adj.P.Val < 0.05, logFC > 0) %>%
  arrange(desc(logFC)) %>%
  dplyr::select(marker, logFC, adj.P.Val, AveExpr)

good_stroma <- de_stroma %>%
  filter(adj.P.Val < 0.05, logFC > 0) %>%
  arrange(desc(logFC)) %>%
  dplyr::select(marker, logFC, adj.P.Val, AveExpr)

# -------------------------------------------------------
# EXTRACT: Top proteins higher in POOR prognosis
# -------------------------------------------------------

poor_cancer <- de_cancer %>%
  filter(adj.P.Val < 0.05, logFC < 0) %>%
  arrange(logFC) %>%  
  dplyr::select(marker, logFC, adj.P.Val, AveExpr)

poor_stroma <- de_stroma %>%
  filter(adj.P.Val < 0.05, logFC < 0) %>%
  arrange(logFC) %>%
  dplyr::select(marker, logFC, adj.P.Val, AveExpr)

# -------------------------------------------------------
# PRINT
# -------------------------------------------------------

cat("=== CANCER CELLS: Higher in GOOD prognosis ===\n"); print(good_cancer)
cat("\n=== CANCER CELLS: Higher in POOR prognosis ===\n"); print(poor_cancer)
cat("\n=== STROMA: Higher in GOOD prognosis ===\n");      print(good_stroma)
cat("\n=== STROMA: Higher in POOR prognosis ===\n");      print(poor_stroma)

# -------------------------------------------------------
# SAVE
# -------------------------------------------------------

#write.csv(good_cancer, "good_prognosis_proteins_cancer_cells.csv", row.names = FALSE)
#write.csv(poor_cancer, "poor_prognosis_proteins_cancer_cells.csv", row.names = FALSE)
#write.csv(good_stroma, "good_prognosis_proteins_stroma.csv",       row.names = FALSE)
#write.csv(poor_stroma, "poor_prognosis_proteins_stroma.csv",       row.names = FALSE)
