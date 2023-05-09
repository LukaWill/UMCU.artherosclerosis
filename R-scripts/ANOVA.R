setwd("~/Graduation/Rstudio")

module_score_cluster <- read.csv("module_score_network39.csv", sep = "\t", header=T)
#module_score_cluster <- data.frame(module_score_cluster)
module_score_cluster$cluster <- factor(module_score_cluster$cluster, levels = c("0", "1", "2", "3", "4"))

library(car)
library(carData)
library(stringr)
library(tidyr)
library(dplyr)

#all Anova & Tukey data together
2:ncol(module_score_cluster)
module_score_clusterz <- rep(NA, ncol(module_score_cluster)) #creates a table, with the same number of columns as counts.cluster
sink("ANOVA-Tukeys-network39.csv")
for (i in 3:ncol(module_score_cluster)) {
  column <- names(module_score_cluster[i])
  module_score_clusterz <- summary(aov(module_score_cluster[, i] ~cluster, data = module_score_cluster))
  tk<- TukeyHSD((aov(module_score_cluster[,i]~cluster, data= module_score_cluster)))
print(column)
print(module_score_clusterz)
print(tk)
}
sink()

#only Tukey data
2:ncol(module_score_cluster)
module_score_clusterz <- rep(NA, ncol(module_score_cluster)) #creates a table, with the same number of columns as counts.cluster
sink("Tukeys-network39.csv")
for (i in 3:ncol(module_score_cluster)) {
  column <- names(module_score_cluster[i])
  module_score_clusterz <- summary(aov(module_score_cluster[, i] ~cluster, data = module_score_cluster))
  tk<- TukeyHSD((aov(module_score_cluster[,i]~cluster, data= module_score_cluster)))
  print(column)
  print(tk)
}
sink()

#Anova p- values 
2:ncol(module_score_cluster)
module_score_clusterz <- rep(NA, ncol(module_score_cluster)) #creates a table, with the same number of columns as counts.cluster
sink("ANOVA-P-values-network39.doc")
for (i in 3:ncol(module_score_cluster)) {
  column <- names(module_score_cluster[i])
  module_score_clusterz <- summary(aov(module_score_cluster[, i] ~cluster, data = module_score_cluster))[[1]][["Pr(>F)"]]
print(column)
print(module_score_clusterz)
}
sink()

Anova.module_score_cluster = read.csv("ANOVA-P-values-network39.doc", header= F)
gene <- Anova.module_score_cluster[seq(from = 1, to = nrow(Anova.module_score_cluster), by = 2), 1]
pvalue <- Anova.module_score_cluster[seq(from =2, to =nrow(Anova.module_score_cluster), by = 2), 1]
Anova.pvalue<- data.frame(gene, pvalue)
Anova.pvalue$gene <- gsub("^\\[1\\]\\s*|\\s*\\[1\\]$", "", Anova.pvalue$gene) #removes [1] before ensembl ID
Anova.pvalue$pvalue <- gsub("^\\[1\\]\\s*|\\s*\\[1\\]$", "", Anova.pvalue$pvalue) #removes [1] before ensembl ID
Anova.pvalue[c("pvalue", "N")] <- str_split_fixed(Anova.pvalue$pvalue, ' ', 2)
Anova.pvalue <- Anova.pvalue [ 1: ncol(Anova.pvalue)-1]
write.csv(Anova.pvalue, file = "p-value_table-network39.csv")
Anova.pvalue.0.05 <- Anova.pvalue[Anova.pvalue$pvalue <.05, ]
write.csv(Anova.pvalue.0.05, file = "p-value_0.05_table-network39.csv")


Tukey_all = read.csv("Tukeys-network39.csv", header = F) # file needs to be changed in a table separately to be able to use 

Tukey_all <- Tukey_all[!grepl("Tukey multiple comparisons of means", Tukey_all$V1),] # removes line "Tukey multiple comparisons of means"
Tukey_all <- Tukey_all[!grepl("95% family-wise confidence level", Tukey_all$V1),] # removes line "95% family-wise confidence level"
removed <- c("Fit:") # select one part of the row that needs to be removed
Tukey_all <- Tukey_all[!grepl(paste(removed, collapse= "|"), Tukey_all$V1),] #removes row that has "fit:" in it
Tukey_all <- Tukey_all[1] #select only first column
Tukey_all$V1 <- gsub("^\\[1\\]\\s*|\\s*\\[1\\]$", "", Tukey_all$V1) #removes [1] before ensembl ID


#Tukey <- data.frame(matrix(ncol= 6, nrow = 10))
#for (i in 1:length(Tukey_all)) {
  gene_name <- Tukey_all[[i]][1]
  comp_data <- Tukey_all[[i]][[4]]
  comp_name <- rownames(comp_data)
  #for (j in 1: (length(comp_data)-1)) {
    comp_name <- comp_name[j]
    diff_val <- comp_data[j + length (comp_name) ]
    lwr_val <- comp_data[j]
    upr_val <- comp_data[j]
    p_val <- comp_data[j]
    new_row <- data.frame(gene_name, comp_name, diff_val, lwr_val, upr_val, p_val)
    Tukey <- rbind(Tukey, new_row)
  #}
#}

colnames(Tukey) <- c("gene", "cluster", "Diff", "Lwr", "Upr", "Pval")







