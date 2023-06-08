setwd("~/Graduation/Rstudio/test")

#load config
source("config.R")

library(dplyr)
library(tidyr)
library(edgeR)
library(DESeq2)
library(tibble)
library(limma)
library(tidyverse)
library(purrr)
library(UCell)
library(ggplot2)
library(gridExtra)
library(car)
library(carData)
library(stringr)
library(reshape2)
library(gplots)
library(RColorBrewer)
library(dendextend)
library(DescTools)


#INPUT DATA, RAWCOUNTS
rawcounts <- read.table(counts_file, header = T, row.names = 1)

#Normalization
data <- DGEList(counts = rawcounts, genes = rownames(rawcounts))
data <- calcNormFactors(data, method = "TMM")
counts <- cpm(data, log = FALSE)
counts.log <- cpm(data, log = TRUE)
counts.log <- as.data.frame(counts.log)

#quantile normalization
quantile.normalized.log =round(limma::normalizeQuantiles(counts.log))

#writing the different normalization steps to file
write.table(counts, file = normalized_counts)
write.table(counts.log, file = normalized_counts_log)
write.table(quantile.normalized.log, file = normalized_quantile)

#plot to quickly check the normalization from the first 50 genes
counts.plot = "~/Graduation/Rstudio/test/plot_normalization.pdf"
pdf(file = counts.plot)

boxplot(rawcounts[1:50])
boxplot(counts.log[1:50])
boxplot(quantile.normalized.log[1:50])

dev.off()
#MODULE SCORE

#Selecting files with the ensembl IDS. (check the read.me for the properties of these files)

#DEFINE the directory path where the files are located
directory_path <- "~/Graduation/Rstudio/Test/Networks"

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

#read the counts file
counts <- read.table("normalized_counts.txt", header = T, row.names = 1)

#list with all networks
AORnetworks <- list()
for (name in network_names) {
  network_object <- get(paste0("network", name))
  AORnetworks[[name]] <- network_object$ensembl
}

#run module score
module_score <- ScoreSignatures_UCell(counts, features = AORnetworks, maxRank = 7000) #maxrank is changed for network77, CRLF1
module_score <- as.data.frame(module_score)

#add plaque types
cluster <- read.table("Seurat_clusters.txt", header = T)
module_score <- module_score %>% rownames_to_column(var = "study_number")
module_score_cluster <- dplyr::full_join(cluster, module_score, by = "study_number")

#write to table
write.table(module_score_cluster, file = module_score_file, sep = "\t")

#plotting

module_score_cluster <- read.csv(module_score_file, sep = "\t", header = T)
module_score_cluster <- module_score_cluster[, -1]
module_score_cluster <- data.frame(module_score_cluster)

module_score_cluster$cluster <- factor(module_score_cluster$cluster, levels = c("0", "1", "2", "3", "4"))

generate_plot <- function(column_index, title) {
  ggplot(module_score_cluster, aes(x = cluster, y = module_score_cluster[, column_index], fill = cluster)) +
    geom_violin() +
    geom_jitter(color = "black", size = 0.1, alpha = 0.9) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      plot.title = element_text(size = 10),
      legend.text = element_text(size = 10), legend.title = element_text(size=10)
    ) +
    labs(title = title, y = 'modulescore', x = 'cluster')
}

plots <- list()

for (i in seq_along(network_names)) {
  column_index <- i + 1
  title <- network_names[i]
  
  plot <- list(column_index = column_index, title = title)
  
  plots[[i]] <- plot
}

plots <- lapply(plots, function(plot) {
  generate_plot(plot$column_index, plot$title)
})

ml <- marrangeGrob(plots, nrow = 2, ncol = 2)
ggsave(violin_plot_module, ml)

#WHICH GENES DO NOT OVERLAP
ensembl <- row.names(counts)
ensembl <- as.data.frame(ensembl)

no_overlap_counts <- list()
for (name in network_names_N) {
  network_object <- get(name)
  no_overlap <- anti_join(network_object, ensembl, by = "ensembl") %>% select(ensembl)
  no_overlap_counts[[name]] <- no_overlap
}

# Initialize empty vectors to store the information

total_genes <- vector("numeric", length(network_names_N))
no_overlap_genes <- vector("numeric", length(network_names_N))
percentage_no_overlap <- vector("numeric", length(network_names_N))

# for loop that uses different vectors for result table
for (i in seq_along(network_names_N)) {
  name <- network_names_N[i]
  network_object <- get(name)
  no_overlap <- no_overlap_counts[[name]] #
  total_genes[i] <- nrow(network_object) #total genes
  no_overlap_genes[i] <- nrow(no_overlap) # no overlap genes
  percentage_no_overlap[i] <- (no_overlap_genes[i] / total_genes[i]) * 100 #percentage no overlap
}

# Create the table
result_table <- data.frame(Network = network_names_N,
                           Total_Genes = total_genes,
                           No_Overlap_Genes = no_overlap_genes,
                           Percentage_No_Overlap = percentage_no_overlap)

write.table(result_table, file = result_no_overlap, sep = "\t")

#ANOVA --> separate script to run

source("ANOVA.R")

#CONNECTIVITY

#input counts
counts_log <- read.table(normalized_counts_log, header = T, row.names = 1)
counts_log <- rownames_to_column(counts_log, "ensembl")

#list per plaquetype
correlation_matrices <- list()
correlation_matrix_0 <- list()
correlation_matrix_1 <- list()
correlation_matrix_2 <- list()
correlation_matrix_3 <- list()
correlation_matrix_4 <- list()

for (name in network_names_N) {
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
  
  #separate network based on plaque types
  cluster0 <- filter(overlap_cluster, cluster == '0')
  cluster1 <- filter(overlap_cluster, cluster == '1')
  cluster2 <- filter(overlap_cluster, cluster == '2')
  cluster3 <- filter(overlap_cluster, cluster == '3')
  cluster4 <- filter(overlap_cluster, cluster == '4')
  
  #perform correlation on all plaque types
  correlation0 <- cor(cluster0[, -1])
  correlation1 <- cor(cluster1[, -1])
  correlation2 <- cor(cluster2[, -1])
  correlation3 <- cor(cluster3[, -1])
  correlation4 <- cor(cluster4[, -1])
  
  #add the correlation matrices to a list per plaque types
  correlation_matrix_0[[name]] <- correlation0
  correlation_matrix_1[[name]] <- correlation1
  correlation_matrix_2[[name]] <- correlation2
  correlation_matrix_3[[name]] <- correlation3
  correlation_matrix_4[[name]] <- correlation4
}

#correlation heatmap

pdf(file = heatmap_cor)

# colorscheme and parameters 
color.scheme <- rev(brewer.pal(10,"RdBu")) 
par(mar = c(3, 3, 2, 2) + 0.1)

cluster_genes <- list()

#loop to make a dendrogram from the correlation and makes a heatmap based on the correlation including the dendogram
for(name in network_names_N) {
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

cluster_genes_total <- do.call(rbind, cluster_genes)
cluster_genes_total <- rownames_to_column(cluster_genes_total, var = "network" )
cluster_genes_total <- cluster_genes_total[, -1]
colnames(cluster_genes_total)<- c("ensembl", "dendogram_cluster", "tissue", "gene_symbol", "network")

#dendogram_cluster to file
write.table(cluster_genes_total, file = dendogram_cluster_genes, sep = "\t")


#connectivity based on average correlation

# function to calculate the average correlation 
calculate_average_correlation_r <- function(correlation_matrix) {
  correlation_r <- as.data.frame(correlation_matrix)
  correlation_r_long  <- gather(correlation_r)
  correlation_r_long <- correlation_r_long %>% filter
  correlation_r_long$value <- abs(correlation_r_long$value)
  average_correlation_r <- correlation_r_long %>% group_by(key) %>% summarize(average_correlation_R = mean(value))
  average_correlation_r <- as.data.frame(average_correlation_r)
  colnames(average_correlation_r) <- c("genes", "average_correlation")
  return(average_correlation_r)
}

#make the lists

average_correlation_all <- list()
average_correlation_total <- list()

for (name in network_names_N) {
  network <- get(name)
  
  # Calculate average correlation for correlation_matrices
  average_correlationTotal <- calculate_average_correlation_r(correlation_matrices[[name]])
  average_correlation0  <- calculate_average_correlation_r(correlation_matrix_0[[name]])
  average_correlation1  <- calculate_average_correlation_r(correlation_matrix_1[[name]])
  average_correlation2  <- calculate_average_correlation_r(correlation_matrix_2[[name]])
  average_correlation3  <- calculate_average_correlation_r(correlation_matrix_3[[name]])
  average_correlation4  <- calculate_average_correlation_r(correlation_matrix_4[[name]])
  
  average_correlation_r <- cbind(average_correlation0, average_correlation1$average_correlation,
                                 average_correlation2$average_correlation, average_correlation3$average_correlation, average_correlation4$average_correlation)
  average_correlation_r <- as.data.frame(average_correlation_r)
  colnames(average_correlation_r) <- c("genes", "0", "1", "2", "3", "4") 
  average_long <- gather(average_correlation_r, key = "cluster", value = "Correlation", -genes)
  
  average_correlation_Total <- cbind(average_correlationTotal, average_correlation0$average_correlation, average_correlation1$average_correlation,
                                     average_correlation2$average_correlation, average_correlation3$average_correlation, average_correlation4$average_correlation)
  average_correlation_Total <- as.data.frame(average_correlation_Total)
  colnames(average_correlation_Total) <- c("genes", "all", "0", "1", "2", "3", "4")
  average_long_all <- gather(average_correlation_Total, key = "cluster", value = "Correlation", -genes)
  
  average_correlation_total[[name]] <- average_long_all[, -1]
  average_correlation_all[[name]] <- average_long[, -1]
  
}

plot_list <- list()

for (name in network_names_N) {
  network <- get(name)
  average_correlations <- average_correlation_all[[name]]
  
  p <- ggplot(average_correlations, aes(x = cluster, y = Correlation, fill = cluster)) +
    geom_violin() +
    geom_jitter(color = "black", size = 0.1, alpha = 0.9) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      plot.title = element_text(size = 10),
      legend.text = element_text(size = 10), legend.title = element_text(size=10)
    ) +
    labs(title = name, y = "Average correlation", x = 'cluster')
  
  plot_list[[name]] <- p
  
  #print(p)
  
} 

plot_list_all <- list()

for (name in network_names_N) {
  network <- get(name)
  average_correlations <- average_correlation_total[[name]]
  
  p <- ggplot(average_correlations, aes(x = cluster, y = Correlation, fill = cluster)) +
    geom_violin() +
    geom_jitter(color = "black", size = 0.1, alpha = 0.9) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      plot.title = element_text(size = 10),
      legend.text = element_text(size = 10), legend.title = element_text(size=10)
    ) +
    labs(title = name, y = "Average correlation", x = 'cluster')
  
  plot_list_all[[name]] <- p
  
}

ml <- marrangeGrob(plot_list, nrow = 2, ncol = 2)
ggsave(connectivity_plot, ml)
pl <- marrangeGrob(plot_list_all, nrow = 2, ncol = 2)
ggsave(connectivity_plot_all, pl)

dev.off()


#If the networks need to be separated into multiple clusters 

ensembl_ids_a <- list()
ensembl_ids_b <- list()

dir.create(directory)

for (name in network_names) {
  filtered_cluster_genes <- filter(cluster_genes_total, network == name)
  cluster_1 <- filter(filtered_cluster_genes, dendogram_cluster == "cluster_1")
  cluster_2 <- filter(filtered_cluster_genes, dendogram_cluster == "cluster_2")
  
  cluster_1 <- cluster_1$ensembl
  cluster_2 <- cluster_2$ensembl
  extra_character <- "ensembl"
  
  ensembl_ids_a[[name]] <- c(extra_character, cluster_1)
  ensembl_ids_b[[name]] <- c(extra_character, cluster_2)
  
  for (i in seq_along(ensembl_ids_a)) {
    file_name_a <- file.path(directory, paste0("Genes_", name, "a", ".txt"))
    writeLines(ensembl_ids_a[[i]], file_name_a)
  }
  for (i in seq_along(ensembl_ids_b)) { 
    file_name_b <- file.path(directory, paste0("Genes_", name, "b", ".txt"))
    writeLines(ensembl_ids_b[[i]], file_name_b)
  }
}


