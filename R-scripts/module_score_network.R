setwd("~/Graduation/Rstudio")

library(UCell)
library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)


counts <- read.table("normalized_counts.txt", header = T, row.names = 1)

network39<- read.table("Genes_39.txt", header = T, "\t")
network195 <- read.table("Genes_195.txt", header = T, "\t")
network177 <- read.table("Genes_177.txt", header = T, "\t")
network24 <- read.table("Genes_24.txt", header = T, "\t")
network35 <- read.table("Genes_35.txt", header = T, "\t")
network74 <- read.table("Genes_74.txt", header = T, "\t")
network109 <- read.table("Genes_109.txt", header = T, "\t")
network172 <- read.table("Genes_172.txt", header = T, "\t")
network147 <- read.table("Genes_147.txt", header = T, "\t")

#making list with all networks
AORnetworks <- list(N39 = network39$ensembl, N195 = network195$ensembl, N177 = network177$ensembl, 
                    N24 = network24$ensembl, N35 = network35$ensembl, N74 = network74$ensembl, 
                    N109 = network109$ensembl, N172 = network172$ensembl, N147 = network147$ensembl)

module_score <- ScoreSignatures_UCell(counts, features = AORnetworks)
module_score <- as.data.frame(module_score)
cluster <- read.table("Seurat_clusters.txt", header = T)
module_score <- module_score %>% rownames_to_column(var = "study_number")
module_score_cluster <- dplyr::full_join(cluster, module_score, by = "study_number")

write.table(module_score_cluster, file = "module_score_networks.csv", sep = "\t")

#which genes are present in network but not present in counts
ensembl <- row.names(counts)
ensembl <- as.data.frame(ensembl)

no_overlap <- map(list(network39 = network39, network195 = network195, network177 = network177, 
                       network24 = network24, network35 = network35, network74 = network74, network109 = network109, 
                       network172 = network172, network147 = network147), 
                  ~ anti_join(.x, ensembl, by = "ensembl") %>% select(ensembl))

overlap <- map(list(network39 = network39, network195 = network195, network177 = network177, 
                       network24 = network24, network35 = network35, network74 = network74, network109 = network109, 
                       network172 = network172, network147 = network147), 
                  ~ semi_join(.x, ensembl, by = "ensembl") %>% select(ensembl))


#which genes are not present in network but are present in counts

no_overlap_counts <- map(list(network39 = network39, network195 = network195, network177 = network177, 
                       network24 = network24, network35 = network35, network74 = network74, network109 = network109, 
                       network172 = network172, network147 = network147), 
                  ~ anti_join(ensembl, .x, by = "ensembl") %>% select(ensembl))

overlap_counts <- map(list(network39 = network39, network195 = network195, network177 = network177, 
                       network24 = network24, network35 = network35, network74 = network74, network109 = network109, 
                       network172 = network172, network147 = network147), 
                  ~ semi_join(ensembl, .x, by = "ensembl") %>% select(ensembl))

