setwd("~/Graduation/Rstudio")

library(UCell)
library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

# Define the directory path where the files are located
directory_path <- "~/Graduation/Rstudio"

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

#making list with all networks
network_names <- c("39", "195", "177", "24", "35", "74", "109", "172", "147")
AORnetworks <- list()
for (name in network_names) {
  network_object <- get(paste0("network", name))
  AORnetworks[[name]] <- network_object$ensembl
}

#run module score
module_score <- ScoreSignatures_UCell(counts, features = AORnetworks)
module_score <- as.data.frame(module_score)

#add plaque types
cluster <- read.table("Seurat_clusters.txt", header = T)
module_score <- module_score %>% rownames_to_column(var = "study_number")
module_score_cluster <- dplyr::full_join(cluster, module_score, by = "study_number")

#write to table
write.table(module_score_cluster, file = "module_score_networks.csv", sep = "\t")

#all genes in count files
ensembl <- row.names(counts)
ensembl <- as.data.frame(ensembl)

#search for genes with no overlap in network
network_names <- c("network39", "network195", "network177", "network24", "network35", "network74", "network109", "network172", "network147")
no_overlap_counts <- list()
for (name in network_names) {
  network_object <- get(name)
  no_overlap <- anti_join(network_object, ensembl, by = "ensembl") %>% select(ensembl)
  no_overlap_counts[[name]] <- no_overlap
}

# Initialize empty vectors to store the information
#network_names <- c("network39", "network195", "network177", "network24", "network35", "network74", "network109", "network172", "network147")
total_genes <- vector("numeric", length(network_names))
no_overlap_genes <- vector("numeric", length(network_names))
percentage_no_overlap <- vector("numeric", length(network_names))

# for loop that uses different vectors for result table
for (i in seq_along(network_names)) {
  name <- network_names[i]
  network_object <- get(name)
  no_overlap <- no_overlap_counts[[name]] #
  total_genes[i] <- nrow(network_object) #total genes
  no_overlap_genes[i] <- nrow(no_overlap) # no overlap genes
  percentage_no_overlap[i] <- (no_overlap_genes[i] / total_genes[i]) * 100 #percentage no overlap
}

# Create the table
result_table <- data.frame(Network = network_names,
                           Total_Genes = total_genes,
                           No_Overlap_Genes = no_overlap_genes,
                           Percentage_No_Overlap = percentage_no_overlap)





