
library(xCell2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(readr)
library(GSVA)
library(tidyverse)

#if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
#}
#devtools::install_github("AlmogAngel/xCell2")

  


setwd("/work_space/files/Transcriptomics")
##exp matrix in correct orientation 
data_1<- read_tsv("retrospective/Leuven_rlog_values.txt")
mat <-data_1
#head(data)

data_c <- data
data_c1 <-data_c[,-1]
mat <-t(as.matrix(data_c1))

####different wiggling around for retrospective
mat <- as.data.frame(mat)        # drop tibble class

rownames(mat) <- mat$Geneid
mat$Geneid <- NULL

mat <- as.matrix(mat)            # if xCell2 expects a matrix



##Trying to keep patient id
# Keep patient IDs from the first column
patient_ids <- data[, 1]

# Drop the ID column and convert to matrix
data_c1 <- data[, -1]
mat_2 <- t(as.matrix(data_c1))

# Assign patient IDs as column names of the transposed matrix
colnames(mat_2) <- patient_ids

####gene set 



data("PanCancer.xCell2Ref", package = "xCell2")
data("TMECompendium.xCell2Ref", package = "xCell2")
data("ImmuneCompendium.xCell2Ref", package = "xCell2")
data("BlueprintEncode.xCell2Ref", package = "xCell2")
data("LM22.xCell2Ref", package = "xCell2")


head(rownames(mat), 20)


##cancer ref so more overlap
data("PanCancer.xCell2Ref", package = "xCell2")
head(res)
dim(res)

####need to find a better reference
ref_genes <- xCell2::getGenesUsed(ImmuneCompendium.xCell2Ref)
length(ref_genes)                      # how many genes in the reference

mix_genes <- rownames(data)
length(intersect(mix_genes, ref_genes))   # how many overlap?

length(intersect(mix_genes, ref_genes)) / length(ref_genes)


# Run the analysis
res_2 <- xCell2Analysis(
  mix = mat, 
  xcell2object = PanCancer.xCell2Ref
)

head (res_2)
##save results
# If your results are in a data frame called `results`
write.csv(res_2, "reterospective_transcriptomics_xcell_results.csv", row.names =T )
# Your expression matrix and reference this one not good ref
res <- xCell2Analysis(
  mix = mat,                              # your genes × samples matrix
  xcell2object = BlueprintEncode.xCell2Ref,  # the reference you downloaded
  minSharedGenes = 0.8,                   # minimum gene overlap (default)
  rawScores = FALSE,                      # return calibrated scores (default)
  spillover = TRUE,                       # apply spillover correction (default)
  spilloverAlpha = 0.5                    # spillover alpha parameter (default)
)

####making df of just patient ID and msi
clinical_msi <- clinical %>%
  select(patient_id, microsatellite_status)
head(clinical_msi)

library(dplyr)

# Transpose xCell2 results so patients become rows
xcell_df <- as.data.frame(t(res_2))

# The rownames of xcell_df are the patient IDs; make them an explicit column
xcell_df$patient_id <- rownames(xcell_df)
rownames(xcell_df) <- NULL
# Reorder so patient_id is first
xcell_df <- xcell_df %>%
  select(patient_id, everything())

# Check it
head(xcell_df)
##loooking perfeclty reasonable

clinical_msi <- clinical %>%
  select(patient_id, microsatellite_status)

merged <- xcell_df %>%
  inner_join(clinical_msi, by = "patient_id")

head(merged)


###now lets plot cd8 vs msi.. think we can guess what that might look like
merged_long <- merged %>%
  pivot_longer(
    cols = -c(patient_id, microsatellite_status),  # keep ID and MSI, pivot everything else
    names_to = "cell_type",
    values_to = "enrichment_score"
  )

head(merged_long)

# Filter for CD8+ alpha-beta T cell
merged_long %>%
  filter(cell_type == "CD8-positive, alpha-beta T cell") %>%
  ggplot(aes(x = microsatellite_status, y = enrichment_score, fill = microsatellite_status)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.3, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2) +
  theme_bw(base_size = 14) +
  labs(
    title = "CD8+ Alpha-Beta T Cell Enrichment by MSI Status",
    x = "Microsatellite Status",
    y = "xCell2 Enrichment Score"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )



#############Proteomic : Transcriptomic correlation starts here
saveRDS(res_2, "res_2.rds")
res_2 <- readRDS("res_2.rds")
feature_df <-read.csv("/work_space/files/spatial/patient_level_proteomics_frommaster_features.csv")
View(feature_df)
head(res_2)
cd8scores <-as.numeric(res_2["CD8-positive, alpha-beta T cell",])

cd8_df <-data.frame(
  patient_id =names(cd8scores),
  xcell_cd8 =as.numeric(cd8scores))
  

##merge whole df's
# Transpose transcriptomics so patients become rows
res_2_t <- as.data.frame(t(res_2))

# Add patient_id column from rownames
res_2_t$patient_id <- rownames(res_2_t)

# Clean patient IDs
res_2_t$patient_id <- trimws(as.character(res_2_t$patient_id))
feature_df$patient_id <- trimws(as.character(feature_df$patient_id))

# Merge patient-level transcriptomics + feature data
merged_df <- merge(
  res_2_t,
  feature_df,
  by = "patient_id",
  all = FALSE
)
wih


cor.test(
  merged_df$"CD8-positive, alpha-beta T cell",
  merged_df$density_Tcyto,
  method ="spearman"
)

cor.test(
  merged_df$'malignant cell',
  merged_df$density_Cancer,
  method ="spearman"
)

cor.test(
  merged_df$"regulatory T cell",
  merged_df$density_Treg,
  method ="spearman"
)
cor.test(
  merged_df$"macrophage",
  merged_df$density_Macrophage,
  method ="spearman"
)

colnames(merged_df)
colnames(res_2_t)
colnames(feature_df)
class(merged_df$xcell_cd8)
str(merged_df$xcell_cd8)




####plots


# Malignant cells
ggplot(
  merged_df,
  aes(
    x = `malignant cell`,
    y = density_Cancer
  )
) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Transcriptomic vs Proteomic Malignant Cells",
    x = "Transcriptomic Malignant Cell Score",
    y = "Proteomic Cancer Density"
  ) +
  theme_minimal()

# Regulatory T cells
ggplot(
  merged_df,
  aes(
    x = `regulatory T cell`,
    y = density_Treg
  )
) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Transcriptomic vs Proteomic Regulatory T Cells",
    x = "Transcriptomic Regulatory T Cell Score",
    y = "Proteomic Treg Density"
  ) +
  theme_minimal()



# Macrophages
ggplot(
  merged_df,
  aes(
    x = macrophage,
    y = density_Macrophage
  )
) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Transcriptomic vs Proteomic Macrophages",
    x = "Transcriptomic Macrophage Score",
    y = "Proteomic Macrophage Density"
  ) +
  theme_minimal()
