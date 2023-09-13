setwd("~/luka_graduation/final_results")

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
library(cluster)    
library(factoextra) 

counts_log <- read.table(normalized_counts_log, header = T, row.names = 1) #reads counts in log
random_gene_lists <- list()
network_names <- 1:1000
network_names <- as.character(network_names)

for (i in 1:1000) {
  random_genes <- sample(rownames(counts_log), 1000, replace = FALSE)
  random_genes <- as.data.frame(random_genes)
  colnames(random_genes) <- "ensembl"
  random_genes_name <- paste0('network', i)
  network_names_N[i] <- random_genes_name
  random_gene_lists[[random_genes_name]] <- as.data.frame(random_genes)
}

list2env(random_gene_lists, envir = .GlobalEnv)

counts_log <- read.table(normalized_counts_log, header = T, row.names = 1) #reads counts in log
counts_log <- rownames_to_column(counts_log, "ensembl") #set rownames to column

correlation_matrices <- list()#empty list
unique_features <- unique(features[, 2])#selects the unique features
unique_features <- na.omit(unique_features) #removes NAs
unique_features <- as.character(unique_features) #set unique features as character
correlations_list <- vector("list", length(unique_features))
correlations_list <- list()

for (name in network_names_N) {
  network <- get(name)
  overlap <- dplyr::semi_join(counts_log, network, by = 'ensembl') #join counts with ensembl IDs network
  overlap <- column_to_rownames(overlap, var = "ensembl") #column to rownames
  overlap <- t(overlap) #transpose
  overlap <- as.data.frame(overlap) #set dataframe
  overlap <- rownames_to_column(overlap, "study_number") #rownames to column
  overlap_feature <- dplyr::full_join(features, overlap, by = "study_number") #join features with ensembl IDS 
  overlap_feature <- arrange(overlap_feature, var = "study_number") #arrange on study_number
  overlap_feature <- na.omit(overlap_feature) #remove NAs
  names(overlap_feature)[2] <- "Feature" #rename column 1
  correlation <- cor(overlap_feature[, -(1:2)]) #calculate correlation for all ensembl IDs
  correlation_matrices[[name]]<- correlation #adds the correlation per network to the empty list
  
  #empty lists
  filtered_features <- list()
  correlations <- list()
  #loops over the unique features and filteres for each unique filter and then calculates the correlation
  for ( feature in unique_features) {
    filtered_features[[feature]] <- filter(overlap_feature, Feature == feature)
    correlations[[feature]] <- cor(filtered_features[[feature]][, -(1:2)])
  }
  
  correlations_list[[name]] <- correlations #adds all the calculated correlations to a list per network
}

save(correlations_list, file = "Output_files/correlation_list_random_genes.Rdata")

#loop to make a dendrogram from the correlation and makes a heatmap based on the correlation including the dendogram

#connectivity based on average correlation

# function to calculate the average correlation  #NEEDS TO BE CHECKED FOR THE PLAQUES
calculate_average_correlation_r <- function(correlation_matrix) {
  correlation_r <- as.data.frame(correlation_matrix) #set dataframe
  correlation_r_long  <- gather(correlation_r) #matrix to long dataframe
  correlation_r_long <- correlation_r_long %>% filter
  correlation_r_long$value <- abs(correlation_r_long$value) #absolute value
  average_correlation_r <- correlation_r_long %>% group_by(key) %>% summarize(average_correlation_R = mean(value)) #calculates average of the absolute value
  average_correlation_r <- as.data.frame(average_correlation_r) #sets data frame
  colnames(average_correlation_r) <- c("genes", "average_correlation") #add column names
  return(average_correlation_r)
}

average_correlation_total <- list() #empty list

#loops over all selected networks
for (name in network_names_N) {
  network <- get(name)
  # Calculate average correlation for correlation_matrices
  for (feature in unique_features) {
    average_correlation <- calculate_average_correlation_r(correlations_list[[name]]) #calculate average correlation for all networks for all features per cluster
    average_correlations <- as.character(average_correlation$genes) #sets genes as character
    #splits the first column on the dot and gives 3 new columns
    split_values <- strsplit(average_correlations, "\\.")
    before_dot <- sapply(split_values, `[`, 1)
    after_dot <- sapply(split_values, `[`, 2)
    average_correlation <- data.frame(feature = before_dot, genes= after_dot, average_correlation = average_correlation[, 2])
    colnames(average_correlation) <- c("feature", "genes", "Correlation")
    #removes an X that is generated within the clustername
    if (any(grepl("X", average_correlation$feature))) {
      # Remove 'X' if it appears in column 1
      average_correlation$feature <- gsub("X", "", average_correlation$feature)
    }
  }
  average_correlation_total[[name]] <- average_correlation[, -2]#adds all the average correlation (connectivity) per network to one list
}

median_correlation_per_network <- list()
average_correlation_per_network<- list()
for (name in network_names_N) {
  median_correlation <- average_correlation_total[[name]] %>% group_by(feature) %>% summarize(median_correlation = median(Correlation))
  average_correlation_network <- average_correlation_total[[name]] %>% group_by(feature) %>% summarize(average_correlation_network = mean(Correlation))
  median_correlation_per_network[[name]] <- median_correlation
  average_correlation_per_network[[name]] <- average_correlation_network
}

median_correlation_dataframe <- bind_rows(median_correlation_per_network, .id = "Network")
median_correlation_df <- median_correlation_dataframe %>% spread(feature, median_correlation)
median_correlation_df <- column_to_rownames(median_correlation_df, var = "Network")
average_median_random_genes <- colMeans(median_correlation_df)
average_median_random_genes <- as.data.frame(average_median_random_genes)
write.table(average_median_random_genes, file = "Output_files/average_median_random_genes.txt", sep = "\t")
