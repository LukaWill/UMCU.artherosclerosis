setwd("~/Graduation/Rstudio")

library(dplyr)
library(tidyr)
library(edgeR)
library(DESeq2)
library(tibble)
library(limma)
library(tidyverse)
library(purrr)
library(reshape2)
library(gplots)
library(RColorBrewer)
library(dendextend)


counts_log <- read.table("normalized_counts_log.txt", header = T, row.names = 1)
cluster <- read.table("Seurat_clusters.txt", header = T)

# Define the directory path where the files are located, the Genes per network are available in separate directories
directory_path <- "~/Graduation/Rstudio/AOR"

# Get the file names into one list with a for loop
file_names <- list.files(directory_path, pattern = "^Genes_.*\\.txt$", full.names = TRUE)
networks <- list()
for (i in 1:length(file_names)) {
  network_name <- gsub("^Genes_(.*)\\.txt$", "\\1", basename(file_names[i]))
  network_name <- paste0("network", network_name)
  network_object <- read.table(file_names[i], header = TRUE, sep = "\t")
  networks[[network_name]] <- network_object
}

# Separate the dataframes in the networks list into individual dataframes
list2env(networks, envir = .GlobalEnv)


#test for all AOR Networks, RENAME the pdf file and network_names, based on the networks. 
heatmap_cor = "~/Graduation/Rstudio/all_correlation.pdf"
pdf(file = heatmap_cor)

#create list with al the correlation matrices
network_names <- c("network39", "network195", "network177", "network24", "network35", "network74", "network109", "network172", "network147")

#if you want to seperate the clusters from the networks (this is only for network 195)
#network_names <- c("cluster_1", "cluster_2")

correlation_matrices <- list()
counts_log <- read.table("normalized_counts_log.txt", header = T, row.names = 1)
counts_log <- rownames_to_column(counts_log, "ensembl")

#loop to make a dataframe with only the genes that belong to the selected network and plaque type in case needed. Then performs a correlation matrix
for (name in network_names) {
  network <- get(name)
  # Perform the correlation matrix calculation per network
  overlap <- dplyr::semi_join(counts_log, network, by = "ensembl")
  overlap <- column_to_rownames(overlap, var = "ensembl")
  overlap <- t(overlap)
  overlap <- as.data.frame(overlap)
  overlap <- rownames_to_column(overlap, "study_number")
  overlap_cluster <- dplyr::full_join(cluster, overlap, by = "study_number")
  overlap_cluster <- arrange(overlap_cluster, cluster)
  overlap_cluster <- column_to_rownames(overlap_cluster, var = 'study_number')
  correlation <- cor(overlap_cluster[, -1])
 
  # Add the correlation matrix to the list
  correlation_matrices[[name]]<- correlation
}
# colorscheme and parameters 
color.scheme <- rev(brewer.pal(10,"RdBu")) 
par(mar = c(5, 4, 4, 6) + 0.1)

cluster_genes <- list()

#loop to make a dendrogram from the correlation and makes a heatmap based on the correlation including the dendogram
for(name in network_names) {
  network <- get(name)
  correlation <- correlation_matrices[[name]]
  correlation.dist <- as.dist(1- correlation)
  correlation.tree <- hclust(correlation.dist, method = "complete")
  correlation.dend <- as.dendrogram(correlation.tree)
  
  #heatmap correlation
  heatmap.2(correlation, 
            Rowv = ladderize(correlation.dend), 
            Colv = ladderize(correlation.dend), 
            dendrogram = "both", 
            revC = TRUE,  # rev column order of dendrogram so conforms to natural representation
            trace = "none", 
            density.info = "none",
            col = color.scheme, keysize = 1,
            labRow = FALSE, labCol = FALSE, 
            xlab = "Genes", ylab = "Genes", 
            main = paste("Heatmap for", name), 
            margins = c(3, 3))
  
  #generating a list with the genes per cluster per network
  num_clusters <- 2
  cluster_labels <- cutree(correlation.tree, k = num_clusters)
  cluster_genes1 <- data.frame(ensembl = c(rownames(correlation)[cluster_labels == 1]), 
                               cluster = "cluster_1")
  cluster_genes2 <- data.frame(ensembl = c(rownames(correlation)[cluster_labels == 2]), 
                               cluster = "cluster_2")
  cluster_column <- rbind(cluster_genes1, cluster_genes2)
  network <- get(name)
  cluster_column <- dplyr::left_join(cluster_column, network, by = "ensembl")
  cluster_genes[[name]] <- cluster_column
  
} 

dev.off()

#sinks the list with the genes per cluster per network into a txt file. 
sink("Test_cluster_genes1.txt")
print(cluster_genes)
sink()

#makes the txt file have better visibility and transforms it into a csv file. 
cluster_genes_total <- read.csv("Test_cluster_genes1.txt", header= F)
cluster_genes_total <- separate(cluster_genes_total, V1, into = c("number", "ensembl", "dendogram_cluster", "tissue", "gene_symbol", "network"), sep = "\\s+")
cluster_genes_total <- na.omit(cluster_genes_total)
cluster_genes_total <- cluster_genes_total[cluster_genes_total$ensembl != "ensembl", ]

write.table(cluster_genes_total, file = "cluster_genes_total_correlation.csv", sep = "\t")


#testing network 195 seperately and rerun rest of the code 

cluster_genes_195 <- filter(cluster_genes_total, network == "195")
write.table(cluster_genes_195, file = "cluster_genes_total_correlation_195.csv", sep = "\t")

cluster_1 <- filter(cluster_genes_195, dendogram_cluster == "cluster_1")
cluster_2 <- filter(cluster_genes_195, dendogram_cluster == "cluster_2")

sink("Test_cluster_195.txt")
print(cluster_genes)
sink()






 
  
 






