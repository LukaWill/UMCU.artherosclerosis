setwd("~/Graduation/Rstudio")

library(dplyr)
library(tidyr)
library(edgeR)
library(DESeq2)
library(tibble)
library(limma)
library(tidyverse)
library(purrr)

#input data, rawcounts 
rawcounts <- read.table("AE_bulk_RNA_batch1.minRib.PC_07042023.txt", header = T, row.names = 1)

#Normalization
data <- DGEList(counts = rawcounts, genes = rownames(rawcounts))
data <- calcNormFactors(data, method = "TMM")
counts <- cpm(data, log = FALSE)
counts.log <- cpm(data, log = TRUE)

counts.log <- as.data.frame(counts.log)
#data <- rownames_to_column(data, var = "gene")
normalized.counts.log =round(limma::normalizeQuantiles(counts.log))

counts <- as.data.frame(counts)
normalized.counts =round(limma::normalizeQuantiles(counts))

write.table(counts, file = "normalized_counts.txt")
write.table(normalized.counts, file = "quantile_normalization.txt")


#Testplot
counts.plot = "~/Graduation/Rstudio/plot_normalization.pdf"
pdf(file = counts.plot)

boxplot(rawcounts[1:50])
boxplot(counts.log[1:50])
boxplot(normalized.counts.log[1:50])

#input
cluster <- read.table("Seurat_clusters.txt", header = T)
#transpose
counts <- t(counts)
#table to datafame
counts <- as.data.frame(counts)
#adding a header to the rownames(first column)
counts <- counts %>% rownames_to_column(var = "study_number")
#joining the counts with the cluster, creating an extra column with the cluster number
counts.cluster <- dplyr::full_join(cluster, counts, by = "study_number")
counts.cluster <- drop_na(counts.cluster)

write.csv(counts.cluster, file = "counts_cluster.csv", fileEncoding = "UTF-8")
counts.cluster <- read.csv("counts_cluster.csv", sep = ",")
counts.cluster <- counts.cluster[, -1]


dev.off()
