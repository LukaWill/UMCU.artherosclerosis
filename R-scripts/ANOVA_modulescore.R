setwd("~/Graduation/Rstudio")

module_score_cluster_AOR <- read.csv("module_score_networks.csv", sep = "\t", header=T)
#module_score_cluster <- data.frame(module_score_cluster)
module_score_cluster_AOR$cluster <- factor(module_score_cluster_AOR$cluster, levels = c("0", "1", "2", "3", "4"))

library(car)
library(carData)
library(stringr)
library(tidyr)
library(dplyr)

#all Anova & Tukey data together
2:ncol(module_score_cluster_AOR)
module_score_clusterz <- rep(NA, ncol(module_score_cluster_AOR)) #creates a table, with the same number of columns as counts.cluster
sink("ANOVA-Tukeys-networks_AOR.csv")
for (i in 3:ncol(module_score_cluster_AOR)) {
  column <- names(module_score_cluster_AOR[i])
  module_score_clusterz <- summary(aov(module_score_cluster_AOR[, i] ~cluster, data = module_score_cluster_AOR))
  tk<- TukeyHSD((aov(module_score_cluster_AOR[,i]~cluster, data= module_score_cluster_AOR)))
print(column)
print(module_score_clusterz)
print(tk)
}
sink()

#only Tukey data
2:ncol(module_score_cluster_AOR)
module_score_clusterz <- rep(NA, ncol(module_score_cluster_AOR)) #creates a table, with the same number of columns as counts.cluster
sink("Tukeys-network_AOR.csv")
for (i in 3:ncol(module_score_cluster_AOR)) {
  column <- names(module_score_cluster_AOR[i])
  module_score_clusterz <- summary(aov(module_score_cluster_AOR[, i] ~cluster, data = module_score_cluster_AOR))
  tk<- TukeyHSD((aov(module_score_cluster_AOR[,i]~cluster, data= module_score_cluster_AOR)))
  print(column)
  print(tk)
}
sink()

#Anova p- values 
2:ncol(module_score_cluster_AOR)
module_score_clusterz <- rep(NA, ncol(module_score_cluster_AOR)) #creates a table, with the same number of columns as counts.cluster
sink("ANOVA-P-values-network_AOR.doc")
for (i in 3:ncol(module_score_cluster_AOR)) {
  column <- names(module_score_cluster_AOR[i])
  module_score_clusterz <- summary(aov(module_score_cluster_AOR[, i] ~cluster, data = module_score_cluster_AOR))[[1]][["Pr(>F)"]]
print(column)
print(module_score_clusterz)
}
sink()

Anova.module_score_cluster = read.csv("ANOVA-P-values-network_AOR.doc", header= F)
gene <- Anova.module_score_cluster[seq(from = 1, to = nrow(Anova.module_score_cluster), by = 2), 1]
pvalue <- Anova.module_score_cluster[seq(from =2, to =nrow(Anova.module_score_cluster), by = 2), 1]
Anova.pvalue<- data.frame(gene, pvalue)
Anova.pvalue$gene <- gsub("^\\[1\\]\\s*|\\s*\\[1\\]$", "", Anova.pvalue$gene) #removes [1] before ensembl ID
Anova.pvalue$pvalue <- gsub("^\\[1\\]\\s*|\\s*\\[1\\]$", "", Anova.pvalue$pvalue) #removes [1] before ensembl ID
Anova.pvalue[c("pvalue", "N")] <- str_split_fixed(Anova.pvalue$pvalue, ' ', 2)
Anova.pvalue <- Anova.pvalue [ 1: ncol(Anova.pvalue)-1]
write.table(Anova.pvalue, file = "p-value_table-network_AOR.csv", sep = "\t")
Anova.pvalue.0.05 <- Anova.pvalue[as.numeric(Anova.pvalue$pvalue) < 0.05, ]
write.table(Anova.pvalue.0.05, file = "p-value_0.05_table-network_AOR.csv", sep = "\t")
 

Tukey_all = read.csv("Tukeys-network_AOR.csv", header = F) # file needs to be changed in a table separately to be able to use 

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
network_numbers <- c("N39", "N195", "N177", "N24", "N35", "N74", "N109", "N172", "N147")
Tukey_all$network_number <- rep(network_numbers, each = 10)
write.table(Tukey_all, file = "comparison_Tukey_AOR.csv", sep = "\t")

Tukey_all_0.05 <- Tukey_all[as.numeric(Tukey_all$p_adj) < 0.05, ]
write.table(Tukey_all_0.05, file = "comparison_Tukey_AOR_0.05.csv", sep = "\t")
write.table(Tukey_all_0.05, file = "comparison_Tukey_AOR_0.05.txt", sep = "\t")

library(gridExtra)
library(knitr)
library(kableExtra)
table <- kable(Tukey_all)
table <- kable(Tukey_all) %>% 
  kable_styling(full_width = F) %>%
  add_header_above(c("TukeyHSD ANOVA" = 7)) %>%
  column_spec(1, bold = TRUE) %>%
  column_spec(3, width = "4em") %>%
  row_spec(0, bold = TRUE, color = "white", background = "darkblue")

pdf("~/Graduation/Rstudio/Tukey.pdf")
cat(table)
dev.off()

options(min.print = 90)
table <- tableGrob(Tukey_all)
pdf("~/Graduation/Rstudio/Tukey.pdf")
grid.arrange(table)
dev.off()


library(officer)
library(magrittr)
table <- data.frame(Tukey_all)
doc <- read_docx()
doc <- doc %>%
  body_add_table(table)

output_file <- "Tukey.docx"
print(doc, target = output_file)

pdf_output_file <- "Tukey.pdf"
convert_to(output_file, target = pdf_output_file)







