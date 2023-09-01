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
rawcounts <- read.table(counts_file, header = T, row.names = 1) #read rawcounts file

#Normalization
data <- DGEList(counts = rawcounts, genes = rownames(rawcounts)) #creates a DGEList object fro a table of counts
data <- calcNormFactors(data, method = "TMM") #calculates normalization factors to scale of the raq library sizes
counts <- cpm(data, log = FALSE) #counts per million
counts.log <- cpm(data, log = TRUE) #counts per million, log2
counts.log <- as.data.frame(counts.log) #set counts.log as data frame
quantile.normalized.log =round(limma::normalizeQuantiles(counts.log)) #quantile normalization

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

write.table(result_table, file = result_no_overlap, sep = "\t") #writes total to file

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


#run module score
module_score <- ScoreSignatures_UCell(counts, features = Listnetworks, maxRank = 7000) #calculates module score
module_score <- as.data.frame(module_score) #set data frame
module_score <- module_score %>% rownames_to_column(var = "study_number") #rownames to first columns 


#adding the features to the module score 
module_score_feature <- dplyr::full_join(features, module_score, by = "study_number") #join module score with feature
names(module_score_feature)[2] <- "Feature" #rename first column
module_score_feature <- arrange(module_score_feature, module_score_feature[, 2])
module_score_feature <- na.omit(module_score_feature) #remove all NAs
module_score_feature <- module_score_feature[, -1] #deselect first column
module_score_feature <- data.frame(module_score_feature) #set data.frame
#names(module_score_feature)[2] <- "Feature" #rename first column

#write to table
write.table(module_score_feature, file = module_score_file, sep = "\t") #write to table

#plotting
module_score_feature$Feature <- factor(module_score_feature$Feature, levels = unique(module_score_feature$Feature)) #set factor levels
#set function for violin plot for module score
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


plots <- list() #empty list
for (i in seq_along(network_names)) { #loops over all selected networks
  column_index <- i + 1
  title <- network_names[i]
  
  plot <- list(column_index = column_index, title = title)
  
  plots[[i]] <- plot
}

plots <- lapply(plots, function(plot) { #apply the function to all networks
  generate_plot(plot$column_index, plot$title)
})

#plot 
pl1 <- marrangeGrob(plots, nrow = 2, ncol = 2) #set parameters
ggsave(violin_plot_module, pl1, width = 5, height = 5) #save plots and adds width / height parameters



#CONNECTIVITY

#input counts
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

#correlation heatmap
pdf(file = heatmap_cor)
# colorscheme and parameters 
color.scheme <- rev(brewer.pal(10,"RdBu")) 
par(mar = c(3, 3, 2, 2) + 0.1)

cluster_genes <- list() #empty list

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

cluster_genes_total <- do.call(rbind, cluster_genes) #add all results in the list with dendogram to a dataframe
cluster_genes_total <- rownames_to_column(cluster_genes_total, var = "network" )#rownames to column
cluster_genes_total <- cluster_genes_total[, -1] #removes column 1
num_columns <- ncol(cluster_genes_total)#checks amount of columns

#gives header to data frame
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

plot_list <- list() #empty lists
#loops over all selected networks and generates a plot  
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
}   

ml <- marrangeGrob(plot_list, nrow = 2, ncol = 2) #set parameters
ggsave(connectivity_plot, ml, width = 5, height = 5) #save plots and add extra parameters

save(average_correlation_total, file = average_correlation_file) #saves average correlation list to RData file

#Statistics 
source("ANOVA_automatic.R") #runs ANOVA for both module score and connectivity









#separates the networks in different clusters based on the dendogram
#to rerun the script with these run them separate for a and b, removing the a or b from the name or rename all names with numbers only

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


