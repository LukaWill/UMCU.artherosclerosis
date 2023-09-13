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


#load config
source("Scripts/config.R")

load(average_correlation_file) #the average_correlation_file from the base script

#Base of the simplifying_and_automate_script.R to filter the names, because the network_names_N in config file contains all names and the config file is reloaded at the start of this script
# Get the file names into one list with a for loop
file_names <- list.files(directory_path, pattern = "^Genes_.*\\.txt$", full.names = TRUE) #looks for file names in selected directory_path
networks <- list() #empty list
for (i in 1:length(file_names)) { #checks all file_names and changes the name to network + number and adds it to the empty list.
  network_name <- gsub("^Genes_(.*)\\.txt$", "\\1", basename(file_names[i]))
  network_name <- paste0("network", network_name)
  network_object <- read.table(file_names[i], header = TRUE, sep = "\t")
  networks[[network_name]] <- network_object
}

# Separate the dataframes in the networks list into individual dataframes
list2env(networks, envir = .GlobalEnv)

#makes a list of all the keydrivers per network
if (exists("directory_path_keydrivers", envir = .GlobalEnv)) {
  file_names_keydrivers <- list.files(directory_path_keydrivers, pattern = "^Keydrivers_.*\\.txt$", full.names = TRUE) #looks for file names in selected directory_path
  network_keydrivers <- list()
  for (i in 1:length(file_names_keydrivers)) { #checks all file_names and changes the name to network + number and adds it to the empty list.
    network_name <- gsub("^Keydrivers_(.*)\\.txt$", "\\1", basename(file_names_keydrivers[i]))
    network_name <- paste0("network", network_name, "_keydrivers")
    network_object <- read.table(file_names_keydrivers[i], header = TRUE, sep = "\t")
    network_keydrivers[[network_name]] <- network_object
  }
  # Separate the dataframes in the networks_keydrivers list into individual dataframes
  list2env(network_keydrivers, envir = .GlobalEnv)
  
  #select names for keydrivers dataframes
  network_names_keydrivers <- paste0(network_names_N, "_keydrivers")
} 

#read the counts file
counts <- read.table(normalized_counts, header = T, row.names = 1)

#select all genes in countsfile
ensembl <- row.names(counts) 
ensembl <- as.data.frame(ensembl)
no_overlap_counts <- list()  #empty list
overlap_counts <- list() #empty list
overlap_counts_keydrivers <- list() #empty list

#for every network that is given in network_names_N (in config.R), the overlap with the countsfile is checked. 
for (name in network_names_N) {
  network_object <- get(name)
  no_overlap <- anti_join(network_object, ensembl, by = "ensembl") %>% select(ensembl)
  overlap <- semi_join(network_object, ensembl, by = "ensembl") %>% select(ensembl)
  no_overlap_counts[[name]] <- no_overlap
  overlap_counts[[name]] <- overlap
}

#for every network that is given in network_names_keydrivers, based on network_names_N, the overlap with the countsfile is checked
if (exists("network_names_keydrivers", envir = .GlobalEnv)) {
  for (name in network_names_keydrivers) {
    if (exists(name, envir = .GlobalEnv)) {
      network_object <- get(name)
      overlap_keydrivers <- semi_join(network_object, ensembl, by = "ensembl") %>% select(ensembl)
      overlap_counts_keydrivers[[name]] <- overlap_keydrivers
    } else {
      warning(paste("Network object", name, "not found. Skipping..."))
    }
  }
}

# Initialize empty vectors to store the information
total_genes <- vector("numeric", length(network_names_N))
no_overlap_genes <- vector("numeric", length(network_names_N))
overlap_genes <- vector("numeric", length(network_names_N))
percentage_no_overlap <- vector("numeric", length(network_names_N))
percentage_overlap <- vector("numeric", length(network_names_N))

if (exists("network_names_keydrivers", envir = .GlobalEnv)) {
  total_genes_keydrivers <- vector("numeric", length(network_names_keydrivers))
  overlap_genes_keydrivers <- vector("numeric", length(network_names_keydrivers))
  percentage_overlap_keydrivers <- vector("numeric", length(network_names_keydrivers))
}


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
# calculates the total genes of the keydrivers and the percentage overlap of the counts with the keydrivers per network
if (exists("network_names_keydrivers", envir = .GlobalEnv)) {
  for (i in seq_along(network_names_keydrivers)) {
    name <- network_names_keydrivers[i]
    if (exists(name, envir = .GlobalEnv)) {
      network_object <- get(name)
      overlap_keydrivers <- overlap_counts_keydrivers[[name]]
      total_genes_keydrivers[i] <- nrow(network_object) #total genes
      overlap_genes_keydrivers[i] <- nrow(overlap_keydrivers) #overlap genes
      percentage_overlap_keydrivers[i] <- (overlap_genes_keydrivers[i] / total_genes_keydrivers[i]) * 100 #percentage overlap
      
    } else {
      warning(paste("Network object", name, "not found. Skipping..."))
    }
  }
} else {
  warning(paste("No keydrivers available. Skipping ..."))
}


# Create the table
result_table <- data.frame(Network = network_names_N,
                           Total_Genes = total_genes,
                           No_Overlap_Genes = no_overlap_genes,
                           Overlap_genes = overlap_genes,
                           Percentage_No_Overlap = percentage_no_overlap, 
                           Percentage_Overlap = percentage_overlap)
#adds percentage overlap with keydrivers to the result table
if (exists("percentage_overlap_keydrivers", envir = .GlobalEnv)) {
  percentage_overlap_keydrivers <- get("percentage_overlap_keydrivers", envir = .GlobalEnv)
  result_table$Percentage_Overlap_keydrivers <- percentage_overlap_keydrivers
}

#selects based on minimum of 25 genes and a maximum of 1000 genes and an overlap of at least 80% with the countsfile.
result_table <- filter(result_table, Total_Genes >25)
result_table <- filter(result_table, Total_Genes <1000)
result_table <- filter(result_table, Percentage_Overlap > 80)

#filters on 0 overlap of keydrivers = meaning there are no keydrivers available, and on at least 80% of the keydrivers need to be present
if ("Percentage_Overlap_keydrivers" %in% colnames(result_table)) {
  result_table <- filter(result_table, Percentage_Overlap_keydrivers == 0 | Percentage_Overlap_keydrivers > 80)
} else {
  # Handle the case when the column doesn't exist
  # For example, you can display a message or perform an alternative action
  message("percentage_overlap_keydrivers column not found in result_table.")
}

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

#Start script median_network.R

median_correlation_per_network <- list() #list with feature (cluster) and median correlation for every network

for (name in network_names_N) {
  median_correlation <- average_correlation_total[[name]] %>% group_by(feature) %>% summarize(median_correlation = median(Correlation))
  median_correlation_per_network[[name]] <- median_correlation
}

median_correlation_dataframe <- bind_rows(median_correlation_per_network, .id = "Network") #dataframe with the same info as the list

random_genes_average <- read.table("Output_files/average_median_random_genes.txt", header = T) #the pre-calculated average median correlation from the random genes

#substract  (median correlation - median correlation random genes)
median_correlation_df <- median_correlation_dataframe %>% spread(feature, median_correlation) #one row for each network with 5 columns corresponding to the plaques
median_correlation_df <- column_to_rownames(median_correlation_df, var = "Network") #network name = rownames
median_correlation_random_genes <- data.frame(matrix(0, nrow = nrow(median_correlation_df), ncol = ncol(median_correlation_df))) #empty dataframe
#colnames and rownames for empty dataframe
colnames(median_correlation_random_genes) <- colnames(median_correlation_df[, -7])
rownames(median_correlation_random_genes) <- rownames(median_correlation_df[, -7])

#calculate the median correlation - average median correlation random genes
for (network in rownames(median_correlation_df)) {
  median_correlation_random_genes[network, "0"] <- median_correlation_df[network, "0"] - random_genes_average[1, ]
  median_correlation_random_genes[network, "1"] <- median_correlation_df[network, "1"] - random_genes_average[2, ]
  median_correlation_random_genes[network, "2"] <- median_correlation_df[network, "2"] - random_genes_average[3, ]
  median_correlation_random_genes[network, "3"] <- median_correlation_df[network, "3"] - random_genes_average[4, ]
  median_correlation_random_genes[network, "4"] <- median_correlation_df[network, "4"] - random_genes_average[5, ]
}

#kmean clustering
fviz_nbclust(median_correlation_random_genes, kmeans, method = "silhouette")#add k.max when group size is lower than 10
k <- kmeans(median_correlation_random_genes, centers = 12, nstart =25) #for future set seed --> here 12 clusters are used
kcluster <- k$cluster
fviz_cluster(k, data = median_correlation_random_genes) #plot the clusters --> for future set seed 

#making dataframe ready
median_correlation_random_genes <- cbind(median_correlation_random_genes, Cluster = k$cluster) #bind cluster to the networks
median_correlation_random_genes <- arrange(median_correlation_random_genes, by_group = Cluster) #and groups them by cluster
median_correlation_random_genes <- rownames_to_column(median_correlation_random_genes, var = "Network") #rownames to column
network_order <- median_correlation_random_genes$Network #set order based on the clusters
#median_correlation_long$Network <- factor(median_correlation_long$Network, levels = network_order)
median_correlation_long_random_genes <- pivot_longer(median_correlation_random_genes, cols = -c(Network, Cluster), names_to = "feature", values_to = "median_correlation") #makes long pivot for plot
median_correlation_long_random_genes$Network <- factor(median_correlation_long_random_genes$Network, levels = network_order) #set factor levels

#plotting heatmap
ggplot(median_correlation_long_random_genes, aes(x = feature, y = Network, fill = median_correlation)) +
  geom_tile() +
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(11, "Spectral")) +
  labs(x = "Feature", y = "Network", title = "Median correlation per feature", fill = 'Median Correlation') +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 4, hjust = 0))  # Adjust the font size and alignment of y-axis labels

ggsave("Output_files/heatmap_median_correlation_random_genes_All_plaque.pdf",width = 5, height = 7) #save the plot


write.table(median_correlation_long_random_genes, file = "Output_files/median_correlation_plaques_random_genes.txt", sep = "\t") #write to table  

#pathway

#new packages needed for pathway analysis
library("clusterProfiler") 
require("org.Hs.eg.db")
library("AnnotationDbi")


network_cluster <- median_correlation_long_random_genes[, 1:2] #selecting networks per cluster
network_cluster <- unique(network_cluster) #only every network 1 time

cluster_list <- split(network_cluster, network_cluster$Cluster) #list with all networks per cluster. 
cluster_list <- lapply(cluster_list, function(x) x[, - which(names(x) == "Cluster")]) #removes the cluster column from the dataframe for each cluster

combined_list <- list() #empty list

#combines the ensembl_ids from each network of the cluster together
for(cluster_name in names(cluster_list)) {
  network_name <- cluster_list[[cluster_name]]
  combined_ensembl_ids <- character(0)
  network_name_clean <- gsub("network", "", network_name$Network)
  for(name in network_name_clean) {
    if(name %in% names(Listnetworks)) {
      network_ensembl_ids <- Listnetworks[[name]]
      combined_ensembl_ids <- c(combined_ensembl_ids, network_ensembl_ids)
    }
  }
  combined_list[[cluster_name]] <- unique(combined_ensembl_ids)
  
}

#comparing pathways between clusters

Entrez_per_cluster <- list() #empty list

#generates the entrezids from the ensemblids 
for(cluster_name in names(combined_list)) {
  ensembl_ids <- combined_list[[cluster_name]]
  Entrez_ids <- mapIds(org.Hs.eg.db, keys = ensembl_ids, keytype = "ENSEMBL", column = "ENTREZID")
  Entrez_per_cluster[[cluster_name]] <- Entrez_ids
  
}

#calculated the enrichedGO terms based on the entrezIDs and compares the different clusters
enrich_result_list_compare <- compareCluster(geneCluster = Entrez_per_cluster, fun = "enrichGO", OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, ont = "BP", pAdjustMethod = "BH")

#plotting dotplot for pathway analysis. 
pdf(file = "Output_files/combined_pathway.pdf", width = 10, height = 10)
dotplot(enrich_result_list_compare, font.size = 10, showCategory = 3) + scale_y_discrete(labels=function(x) str_wrap(x, width=90))

dev.off()

#results from enriched pathway analyisis in table format
result_comparison <- enrich_result_list_compare@compareClusterResult


#pathway analysis separate. 
enrich_result_list <- list()

for (cluster_name in names(combined_list)) {
  ensembl_ids <- combined_list[[cluster_name]]
  enrich_result <- enrichGO(gene = ensembl_ids, 
                            OrgDb    = org.Hs.eg.db,  # Specify your organism database
                            keyType  = "ENSEMBL",    # Specify the gene identifier type
                            ont      = "BP",          # Choose the ontology (e.g., "BP" for biological process)
                            pAdjustMethod = "BH",     # Multiple testing correction method
                            pvalueCutoff  = 0.05      # Significance cutoff
  )
  enrich_result_list[[cluster_name]] <- enrich_result
}

#plotting the enriched pathway analysis separate 
plot_list <- list()

for(cluster_name in names(enrich_result_list)) {
  enrich_result_total <- enrich_result_list[[cluster_name]]
  p <- dotplot(enrich_result_total, showCategory = 20)
  plot_list [[cluster_name]] <- p
}

ml <- marrangeGrob(plot_list, nrow = 1, ncol = 1)
ggsave("output_files/pathway_plot.pdf", ml, width = 10, height = 10)

#PERCENTAGES
#list with all tissues per network
Listnetworks_tissue <- list()
for (name in network_names_N) {
  network <- get(name)
  Listnetworks_tissue[[name]] <- network$tissue
}

#calculated the percentages of tissue per network
network_percentages <- lapply(Listnetworks_tissue, function(network) {
  tissue_counts <- table(network)
  tissue_percentages <- (tissue_counts / length(network)) * 100 # Use the total of the first network as the denominator
  return(tissue_percentages)
})

#puts the calculated percentages in a dataframe
network_percentage_df <- do.call(rbind, lapply(names(network_percentages), function(network) {
  df <- as.data.frame(network_percentages[[network]])
  df$Network <- network
  return(df)
}))

colnames(network_percentage_df) <- c("Tissue", "Percentage", "Network")

#adds NA values to the tissues that are not present in the network
heatmap_data <- network_percentage_df %>%
  pivot_wider(names_from = Tissue, values_from = Percentage)

#combines the data in three columns, network, tissue and percentages
heatmap_data_long <- heatmap_data %>%
  pivot_longer(cols = -Network, names_to = "Tissue", values_to = "Percentage") %>%
  mutate(Percentage = replace_na(Percentage, 0))

names_cluster <- unique(network_cluster$Network) #gives unique values new name

heatmap_data_long$Network <- factor(heatmap_data_long$Network, levels = names_cluster) #sets Network as factor

heatmap_data_long <- heatmap_data_long %>%
  left_join(network_cluster, by = "Network") #adds the cluster categories to the data

#group the data on the clusters ad calculate the average percentage tissue per cluster
cluster_tissue_averages <- heatmap_data_long %>%
  group_by(Cluster, Tissue) %>%
  summarize(Average_Percentage = mean(Percentage, na.rm = TRUE))

#plotting the average tissue percentages
unique_clusters <- unique(cluster_tissue_averages$Cluster)
tissue_colors <- c("AOR" = "red", "MAM" = "darkorange", "VAF" = "purple", "SF" = "magenta", "BLOOD" = "blue", "LIV" = "brown", "SKLM" = "darkgreen" )

plot_list <- list()
for (i in unique_clusters) {
  cluster_data <- cluster_tissue_averages %>%
    filter(Cluster == i)
  
  cluster_data$legend_label <- paste0(cluster_data$Tissue, " ", "(",
                                      round(cluster_data$Average_Percentage, 2), "%", ")")
  
  plot <- ggplot(cluster_data, aes(x = "", y = Average_Percentage, fill = Tissue)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    labs(title = paste("Cluster", i)) +
    theme_void() +
    scale_fill_manual(values = tissue_colors, labels = cluster_data$legend_label ) 
  
  plot_list[[i]] <- plot
}
ml <- marrangeGrob(plot_list, nrow = 2, ncol = 1) # Set parameters
ggsave("Output_files/piechart_average_percentages.pdf", ml, width = 5, height = 6)




#CALCULATING THE PERCENTAGES PER NETWORK 

#making a list with tissue per cluster
tissue_per_cluster <- list()

for(i in 1:length(cluster_list)) {
  network_selection <- cluster_list[[i]]
  
  filtered_data <- data.frame()
  
  for (name in network_selection$Network) {
    if (name %in% names(Listnetworks_tissue)) {
      network_data <- data.frame(Network = name, Tissue = unlist(Listnetworks_tissue[[name]]))
      filtered_data <- rbind(filtered_data, network_data)
    }
  }
  tissue_per_cluster[[i]] <- filtered_data
}
#list with percentage tissue per cluster
cluster_tissue_percentages <- list()

for (i in 1:length(tissue_per_cluster)) {
  cluster_data <- tissue_per_cluster[[i]]
  
  # Calculate tissue percentages for the current cluster
  tissue_counts <- table(cluster_data$Tissue)
  tissue_percentages <- (tissue_counts / nrow(cluster_data)) * 100
  
  # Store tissue percentages in cluster_tissue_percentages
  cluster_tissue_percentages[[i]] <- tissue_percentages
}
 
save(cluster_tissue_percentages, file = "Output_files/cluster_tissue_percentages.RData") #save the cluster_tissue_percentages as .RData


#plotting the percentages in a piechart
plot_list <- list()
for (i in 1:length(cluster_tissue_percentages)) {
  cluster_tissue_percentages_df <- data.frame(
    Tissue = names(cluster_tissue_percentages[[i]]),
    Percentage = cluster_tissue_percentages[[i]]
  )
  
  cluster_tissue_percentages_df$legend_label <- paste0(cluster_tissue_percentages_df$Tissue, " ", "(",
                                                      round(cluster_tissue_percentages_df$Percentage.Freq, 2), "%", ")")#adding the percentages to the legend
  
  
  # Create a pie chart
  p <- ggplot(cluster_tissue_percentages_df, aes(x = "", y = Percentage.Freq, fill = Tissue)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    labs(title = paste("Cluster", i)) +
    theme_void() +
    scale_fill_manual(values = tissue_colors, labels = cluster_tissue_percentages_df$legend_label) 
    
  
  
  # Print the pie chart
  plot_list[[i]] <- p
}


ml <- marrangeGrob(plot_list, nrow = 2, ncol = 1) # Set parameters

ggsave("Output_files/piechart_percentages.pdf", ml, width = 5, height = 6)


