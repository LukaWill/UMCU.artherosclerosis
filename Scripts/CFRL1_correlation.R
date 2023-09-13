setwd("~/luka_graduation/final_results")

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
library(DescTools)

negative_correlation <- read.table("Scripts/negative_correlation.txt", header = T, row.names = 1)

negative_correlation <- rownames_to_column(negative_correlation)
negative_correlation <- negative_correlation %>% separate(rowname, c("study_number", "ensembl"))
negative_correlation_CRLF1_symbol <- as.data.frame(negative_correlation$study_number)
networkNegative <- as.data.frame(negative_correlation$ensembl)
colnames(networkNegative) <- c("ensembl")
write.table(networkNegative, file = "negative_correlation_CRLF1_symbol.csv", sep = "\t") #Symbols for STARNET
write.table(networkNegative, file = "negative_correlation_CRLF1_ensembl.csv", sep = "\t") #ensembl ids for network



positive_correlation <- read.table("Scripts/positive_correlation.txt", header = T, row.names = 1)

positive_correlation <- rownames_to_column(positive_correlation)
positive_correlation <- positive_correlation %>% separate(rowname, c("study_number", "ensembl"))
positive_correlation_CRLF1_symbol <- as.data.frame(positive_correlation$study_number)
networkPositive <- as.data.frame(positive_correlation$ensembl)
colnames(networkPositive) <- c("ensembl")

write.table(positive_correlation_CRLF1_symbol, file = "Output_files/positive_correlation_CRLF1_symbol.csv", sep = "\t") #Symbols for STARNET
write.table(networkPositive, file = "Output_files/positive_correlation_CRLF1_ensembl.csv", sep = "\t") #ensembl ids for network


network1 <- read.table("Output_files/positive_correlation_CRLF1_ensembl.csv", header = T) #positive
network2 <- read.table("Output_files/negative_correlation_CRLF1_ensembl.csv", header = T) #negative


