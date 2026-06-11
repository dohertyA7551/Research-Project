
# Install via BiocManager
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
BiocManager::install("decoupleR")
install.packages("tidyverse")
BiocManager::install("OmnipathR")
library(OmnipathR)
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(readr)

setwd("/RAIDDISK/home/alanna/Documents/Transcriptomics")

###loading data
data<- read.csv("stage2/normalized_genes.csv")
problems(data)
head(problems)

head(data_c)
names(data_c)

#make matrix in correct shape
#sum(is.na(data_c))
data_c1 <-data_c[-1]
mat <-as.matrix(data_c1)
counts <-t(mat)
head(counts)

##Get progeny
net <- decoupleR::get_progeny(organism = 'human', 
                              top = 500)

net

##activity inference
# Run mlm
sample_acts <- decoupleR::run_mlm(mat = counts, 
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'weight', 
                                  minsize = 5)
sample_acts


#visualisation
# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  tidyr::pivot_wider(id_cols = 'condition', 
                     names_from = 'source',
                     values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix()

# Scale per feature
sample_acts_mat <- scale(sample_acts_mat)

# Color scale
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05,2, length.out = floor(100 / 2)))

# Plot
p<-pheatmap::pheatmap(mat = sample_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 20,
                   cellheight = 20,
                   treeheight_row = 20,
                   treeheight_col = 20)

p
