source("config.R")

library(car)
library(carData)
library(stringr)
library(tidyr)
library(dplyr)

module_score_cluster_ANOVA <- read.csv(module_score_file, sep = "\t", header=T)
module_score_cluster_ANOVA$cluster <- factor(module_score_cluster_ANOVA$cluster, levels = c("0", "1", "2", "3", "4"))

#all Anova & Tukey data together
2:ncol(module_score_cluster_ANOVA)
module_score_clusterz <- rep(NA, ncol(module_score_cluster_ANOVA)) #creates a table, with the same number of columns as counts.cluster
sink(Anova_Tukey)
for (i in 3:ncol(module_score_cluster_ANOVA)) {
  column <- names(module_score_cluster_ANOVA[i])
  module_score_clusterz <- summary(aov(module_score_cluster_ANOVA[, i] ~cluster, data = module_score_cluster_ANOVA))
  tk<- TukeyHSD((aov(module_score_cluster_ANOVA[,i]~cluster, data= module_score_cluster_ANOVA)))
  print(column)
  print(module_score_clusterz)
  print(tk)
}
sink()

#only Tukey data
2:ncol(module_score_cluster_ANOVA)
module_score_clusterz <- rep(NA, ncol(module_score_cluster_ANOVA)) #creates a table, with the same number of columns as counts.cluster
sink(Tukey)
for (i in 3:ncol(module_score_cluster_ANOVA)) {
  column <- names(module_score_cluster_ANOVA[i])
  module_score_clusterz <- summary(aov(module_score_cluster_ANOVA[, i] ~cluster, data = module_score_cluster_ANOVA))
  tk<- TukeyHSD((aov(module_score_cluster_ANOVA[,i]~cluster, data= module_score_cluster_ANOVA)))
  print(column)
  print(tk)
}
sink()

#Anova p- values 
2:ncol(module_score_cluster_ANOVA)
module_score_clusterz <- rep(NA, ncol(module_score_cluster_ANOVA)) #creates a table, with the same number of columns as counts.cluster
sink(Anova_pvalue_file)
for (i in 3:ncol(module_score_cluster_ANOVA)) {
  column <- names(module_score_cluster_ANOVA[i])
  module_score_clusterz <- summary(aov(module_score_cluster_ANOVA[, i] ~cluster, data = module_score_cluster_ANOVA))[[1]][["Pr(>F)"]]
  print(column)
  print(module_score_clusterz)
}
sink()

#seperate the p values 
Anova_pvalue = read.csv(Anova_pvalue_file, header= F)
gene <- Anova_pvalue[seq(from = 1, to = nrow(Anova_pvalue), by = 2), 1]
pvalue <- Anova_pvalue[seq(from =2, to =nrow(Anova_pvalue), by = 2), 1]
Anova.pvalue<- data.frame(gene, pvalue)
Anova.pvalue$gene <- gsub("^\\[1\\]\\s*|\\s*\\[1\\]$", "", Anova.pvalue$gene) #removes [1] before ensembl ID
Anova.pvalue$pvalue <- gsub("^\\[1\\]\\s*|\\s*\\[1\\]$", "", Anova.pvalue$pvalue) #removes [1] before ensembl ID
Anova.pvalue[c("pvalue", "N")] <- str_split_fixed(Anova.pvalue$pvalue, ' ', 2)
Anova.pvalue <- Anova.pvalue [ 1: ncol(Anova.pvalue)-1]

#write pvalues to table
write.table(Anova.pvalue, file = pvalue_Anova, sep = "\t")
Anova.pvalue.0.05 <- Anova.pvalue[as.numeric(Anova.pvalue$pvalue) < 0.05, ]
write.table(Anova.pvalue.0.05, file = pvalue_Anova_0.05, sep = "\t")


#separate Tukey results 

Tukey_all = read.csv(Tukey, header = F) 
Tukey_all <- Tukey_all[!grepl("Tukey multiple comparisons of means", Tukey_all$V1),] # removes line "Tukey multiple comparisons of means"
Tukey_all <- Tukey_all[!grepl("95% family-wise confidence level", Tukey_all$V1),] # removes line "95% family-wise confidence level"
removed <- c("Fit:") # select one part of the row that needs to be removed
Tukey_all <- Tukey_all[!grepl(paste(removed, collapse= "|"), Tukey_all$V1),] #removes row that has "fit:" in it
Tukey_all <- Tukey_all[1] #select only first column
Tukey_all$V1 <- gsub("^\\[1\\]\\s*|\\s*\\[1\\]$", "", Tukey_all$V1) #removes [1] before ensembl ID
Tukey_all <- separate(Tukey_all, V1, into = c("comparison", "diff", "lwr", "upr", "p_adj"), sep = "\\s+")
Tukey_all <- mutate_all(Tukey_all, trimws)
#removes the two lines with NA and removes the line that contains the header, adds the network_numbers and writes it to one file
Tukey_all <- Tukey_all[complete.cases(Tukey_all), ]
Tukey_all <- Tukey_all[Tukey_all$diff != "diff", ]
Tukey_all$network_number <- rep(network_names, each = 10)

#write to table
write.table(Tukey_all, file = Tukey_comparison, sep = "\t")
Tukey_all_0.05 <- Tukey_all[as.numeric(Tukey_all$p_adj) < 0.05, ]
write.table(Tukey_all_0.05, file = Tukey_comparison_0.05, sep = "\t")





