setwd("~/Graduation/Rstudio")

library(UCell)
library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
#library(Seurat)


Network<- read.table("Genes_39.txt", header = T, "\t")
counts <- read.table("normalized_counts.txt", header = T, row.names = 1)

#Genes <- read.csv("counts_cluster.csv", header = T, ",")
#Genes <- Genes[, -1]
counts <- counts %>% rownames_to_column(var = "gene")
#joining both both input data by Ensembl ID
counts_network <- dplyr::full_join(Network, counts, by = "gene")
counts_network <- drop_na(counts_network)
counts_network = subset(counts_network, select = -c(tissue, clust, gene_symbol ))
counts_network <- t(counts_network)
counts_network <- as.data.frame(counts_network)
colnames(counts_network) <- counts_network[1,]
counts_network <- counts_network[-1,]
counts_network<- counts_network %>% rownames_to_column(var = "study_number")

#joining the counts in network with the cluster(different plaque type)
cluster <- read.table("Seurat_clusters.txt", header = T)
counts.cluster.network <- dplyr::full_join(cluster, counts_network, by = "study_number")
write.table(counts.cluster.network, file = "counts_cluster_network.csv", sep = "\t")

#try 1 studynumber = rownames to columnnames. genes is column to row
rownames(counts.cluster.network) <- counts.cluster.network$study_number
counts.cluster.network <- counts.cluster.network[, -1]
counts.cluster.network.transpose <- t(counts.cluster.network)
counts.cluster.network.transpose <- as.data.frame(counts.cluster.network.transpose)
counts.cluster.network.transpose <- counts.cluster.network.transpose[-1, ]

#matrix <- as.matrix(counts.cluster.network.transpose)

#try 2 (rerun counts.cluster.network <- dplyr::.......)
counts.cluster.network <- counts.cluster.network[, -2]
rownames(counts.cluster.network) <- counts.cluster.network$study_number
counts.cluster.network <- counts.cluster.network[, -1]
matrix.count <- as.matrix(counts.cluster.network)

type0 <- cluster[cluster$cluster == "0",]
type0 <- type0[, -2]


type1 <- cluster[cluster$cluster == "1",]
type1 <- type1[, -2]

type2 <- cluster[cluster$cluster == "2",]
type2 <- type2[, -2]


type3 <- cluster[cluster$cluster == "3",]
type3 <- type3[, -2]


type4 <- cluster[cluster$cluster == "4",]
type4 <- type4[, -2]


signature <- list(type0 = type0, type1 = type1, type2 = type2, type3 = type3, type4 = type4)





#not correct
genes <- list(colnames(counts.cluster.network[, !names(counts.cluster.network) %in% c("study_number", "cluster")]))

module_score <- ScoreSignatures_UCell(counts.cluster.network.transpose, features = genes)

module_score <- ScoreSignatures_UCell(counts.cluster.network, features = signature)

module_score <- ScoreSignatures_UCell(matrix.count, features = signature)

#try with all the counts

counts <- read.table("normalized_counts.txt", header = T, row.names = 1)

Network<- read.table("Genes_39.txt", header = T, "\t")
network39 <- list(Network$gene)

module_score <- ScoreSignatures_UCell(counts, features = network39)

#probably correct solution

counts <- read.table("normalized_counts.txt", header = T, row.names = 1)

Network<- read.table("Genes_39.txt", header = T, "\t")

network39 <- Network$gene
network39 <- as.data.frame(network39)

# group are all the genes that are present in network39 separated in a separate group
group <- split(network39$network39, network39$network39)

module_score <- ScoreSignatures_UCell(counts, features = group)

cluster <- read.table("Seurat_clusters.txt", header = T)
module_score <-as.data.frame(module_score)
module_score <- module_score %>% rownames_to_column(var = "study_number")
module_score_cluster <- dplyr::full_join(cluster, module_score, by = "study_number")
write.table(module_score_cluster, file = "module_score_network39.csv", sep = "\t")


   
