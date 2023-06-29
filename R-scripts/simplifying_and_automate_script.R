setwd("~/Graduation/Rstudio/test")

library(dplyr)
library(tidyr)
library(edgeR)
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

#load config
source("config.R")

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
pdf(file = counts.plot)

boxplot(rawcounts[1:50])
boxplot(counts.log[1:50])
boxplot(quantile.normalized.log[1:50])

dev.off()
#MODULE SCORE

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

#when working with clusters of a selected network (Networka.b)

selected_names <- character()
for (name in network_names_N) {
  # Check if data frame with the given name exists in the environment
  if (exists(paste0(name, "a"))) {
    # Add the "a" name to the vector of selected names
    selected_names <- c(selected_names, paste0(name, "a"))
  }
  if (exists(paste0(name, "b"))) {
    # Add the "b" name to the vector of selected names
    selected_names <- c(selected_names, paste0(name, "b"))
  }
}

#selecting network_names_N and network_names based on the clusters
if (cluster_network) {
  network_names_N <- selected_names
  network_names <- gsub("^network(\\d+)([ab])$", "\\1\\2", selected_names)
}


#read the counts file
counts <- read.table(normalized_counts, header = T, row.names = 1)


#WHICH GENES DO AND DO NOT OVERLAP
ensembl <- row.names(counts)
ensembl <- as.data.frame(ensembl)

no_overlap_counts <- list()
overlap_counts <- list()
for (name in network_names_N) {
  network_object <- get(name)
  no_overlap <- anti_join(network_object, ensembl, by = "ensembl") %>% select(ensembl)
  overlap <- semi_join(network_object, ensembl, by = "ensembl") %>% select(ensembl)
  no_overlap_counts[[name]] <- no_overlap
  overlap_counts[[name]] <- overlap
}

# Initialize empty vectors to store the information

total_genes <- vector("numeric", length(network_names_N))
no_overlap_genes <- vector("numeric", length(network_names_N))
overlap_genes <- vector("numeric", length(network_names_N))
percentage_no_overlap <- vector("numeric", length(network_names_N))
percentage_overlap <- vector("numeric", length(network_names_N))

# for loop that uses different vectors for result table
for (i in seq_along(network_names_N)) {
  name <- network_names_N[i]
  network_object <- get(name)
  no_overlap <- no_overlap_counts[[name]] 
  overlap <- overlap_counts[[name]]
  total_genes[i] <- nrow(network_object) #total genes
  no_overlap_genes[i] <- nrow(no_overlap) # no overlap genes
  overlap_genes[i] <- nrow(overlap) #overlap genes
  percentage_no_overlap[i] <- (no_overlap_genes[i] / total_genes[i]) * 100 #percentage no overlap
  percentage_overlap[i] <- (overlap_genes[i] / total_genes[i]) * 100 #percentage overlap
}

# Create the table
result_table <- data.frame(Network = network_names_N,
                           Total_Genes = total_genes,
                           No_Overlap_Genes = no_overlap_genes,
                           Overlap_genes = overlap_genes,
                           Percentage_No_Overlap = percentage_no_overlap, 
                           Percentage_Overlap = percentage_overlap)

write.table(result_table, file = result_no_overlap, sep = "\t")


result_table <- filter(result_table, Total_Genes <1000)
result_table <- filter(result_table, Percentage_Overlap > 80)
column_name <- result_table$Network

#here the network_names are selected based on the result_table, The network_names_N in the config.R file is used for the beginning and is now not used anymore. 
network_names_N <- result_table %>% pull(Network)
network_names <- str_extract(column_name, "\\d+")

#list with all networks
Listnetworks <- list()
for (name in network_names) {
  network_object <- get(paste0("network", name))
  Listnetworks[[name]] <- network_object$ensembl
}


#run module score
module_score <- ScoreSignatures_UCell(counts, features = Listnetworks, maxRank = 7000) #maxrank is changed for network77, CRLF1
module_score <- as.data.frame(module_score)
module_score <- module_score %>% rownames_to_column(var = "study_number")

#adding the features to the module score 
module_score_feature <- dplyr::full_join(features, module_score, by = "study_number")
module_score_feature <- na.omit(module_score_feature)
module_score_feature <- module_score_feature[, -1]
module_score_feature <- data.frame(module_score_feature)
names(module_score_feature)[1] <- "Feature"
#write to table
write.table(module_score_feature, file = module_score_file, sep = "\t")

#plotting
module_score_feature$Feature <- factor(module_score_feature$Feature, levels = unique(module_score_feature$Feature)) 
generate_plot <- function(column_index, title) {
  ggplot(module_score_feature, aes(x = Feature, y = module_score_feature[, column_index], fill = Feature)) +
    geom_violin() +
    geom_jitter(color = "black", size = 0.1, alpha = 0.9) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      plot.title = element_text(size = 10),
      legend.text = element_text(size = 10), legend.title = element_text(size=10)
    ) +
    labs(title = title, y = 'modulescore', x = 'Feature')
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

#plot 
pl1 <- marrangeGrob(plots, nrow = 2, ncol = 2)
ggsave(violin_plot_module, pl1, width = 5, height = 5)



#CONNECTIVITY

#input counts
counts_log <- read.table(normalized_counts_log, header = T, row.names = 1)
counts_log <- rownames_to_column(counts_log, "ensembl")

correlation_matrices <- list()
unique_features <- unique(features[, 2])
unique_features <- na.omit(unique_features)
unique_features <- as.character(unique_features)
correlations_list <- vector("list", length(unique_features))
correlations_list <- list()

for (name in network_names_N) {
  network <- get(name)
  overlap <- dplyr::semi_join(counts_log, network, by = 'ensembl')
  overlap <- column_to_rownames(overlap, var = "ensembl")
  overlap <- t(overlap)
  overlap <- as.data.frame(overlap)
  overlap <- rownames_to_column(overlap, "study_number")
  overlap_feature <- dplyr::full_join(features, overlap, by = "study_number")
  overlap_feature <- arrange(overlap_feature, var = "study_number")
  overlap_feature <- na.omit(overlap_feature)
  names(overlap_feature)[2] <- "Feature"
  correlation <- cor(overlap_feature[, -(1:2)])
  correlation_matrices[[name]]<- correlation
  
  filtered_features <- list()
  correlations <- list()
  
  for ( feature in unique_features) {
    filtered_features[[feature]] <- filter(overlap_feature, Feature == feature)
    correlations[[feature]] <- cor(filtered_features[[feature]][, -(1:2)])
  }
  
  correlations_list[[name]] <- correlations
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
num_columns <- ncol(cluster_genes_total)

if(num_columns == 2) {
  colnames(cluster_genes_total) <- c("ensembl", "dendogram_cluster")
} else if  (num_columns == 5) {
  colnames(cluster_genes_total)<- c("ensembl", "dendogram_cluster", "tissue", "gene_symbol", "network")
} else {
  colnames(cluster_genes_total)[1:2] <- c("ensembl", "dendogram_cluster")
}

#dendogram_cluster to file
write.table(cluster_genes_total, file = dendogram_cluster_genes, sep = "\t")


#connectivity based on average correlation

# function to calculate the average correlation  #NEEDS TO BE CHECKED FOR THE PLAQUES
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

average_correlation_total <- list()

for (name in network_names_N) {
  network <- get(name)
  
  # Calculate average correlation for correlation_matrices
  
  for (feature in unique_features) {
    average_correlation <- calculate_average_correlation_r(correlations_list[[name]])
    average_correlations <- as.character(average_correlation$genes)
    split_values <- strsplit(average_correlations, "\\.")
    
    # Extract the parts before and after the dot
    before_dot <- sapply(split_values, `[`, 1)
    after_dot <- sapply(split_values, `[`, 2)
    
    average_correlation <- data.frame(feature = before_dot, genes= after_dot, average_correlation = average_correlation[, 2])
    colnames(average_correlation) <- c("feature", "genes", "Correlation")
    if (any(grepl("X", average_correlation$feature))) {
      # Remove 'X' if it appears in column 1
      average_correlation$feature <- gsub("X", "", average_correlation$feature)
    }
  }
  average_correlation_total[[name]] <- average_correlation[, -2]
}

plot_list <- list()
  
for (name in network_names_N) {
  network <- get(name)
  average_correlations <- average_correlation_total[[name]]
  
  p <- ggplot(average_correlations, aes(x = feature, y = Correlation, fill = feature)) +
    geom_violin() +
    geom_jitter(color = "black", size = 0.1, alpha = 0.9) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      plot.title = element_text(size = 10),
      legend.text = element_text(size = 10), legend.title = element_text(size=10)
    ) +
    labs(title = name, y = "Average correlation", x = 'feature')
  
  plot_list[[name]] <- p
  
  #print(p)
  
}   



ml <- marrangeGrob(plot_list, nrow = 2, ncol = 2)
ggsave(connectivity_plot, ml, width = 5, height = 5)

#save average correlation

save(average_correlation_total, file = average_correlation_file)



#Statistics 
source("ANOVA_automatic.R")


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


