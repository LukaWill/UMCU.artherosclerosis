
library(car)
library(carData)
library(stringr)
library(tidyr)
library(dplyr)

module_score_feature_ANOVA <- read.csv(module_score_file, sep = "\t", header=T)

module_score_feature_ANOVA$Feature <- as.character(module_score_feature_ANOVA$Feature)
module_score_feature_ANOVA$Feature <- factor(module_score_feature_ANOVA$Feature, levels = unique(module_score_feature_ANOVA$Feature))

network <- colnames(module_score_feature_ANOVA[, -1])
networknames <- gsub('[^0-9]+', "", network)
module_score_list <- list()

for (i in seq_along(network)){
  name <- networknames[i]
  network_dataframe <- data.frame(feature = module_score_feature_ANOVA[, 1], module_scores = module_score_feature_ANOVA[[network[i]]])
  module_score_list[[name]] <- network_dataframe
  
}

#CHANGE THESE NAMES BASED ON THE QUESTION, ONLY CHANGE THE PART BETWEEN THE ""
runAnovaTukey <- function(input_data, file_prefix) {
  Anova_Tukey <- paste0(file_prefix, "_Anova_Tukey_ALL_plaque_cluster.txt")
  Tukey <- paste0(file_prefix, "_Tukey_ALL_symptoms_plaque_cluster.txt")
  Anova_pvalue_file <- paste0(file_prefix, "_Anova_pvalue_ALL_plaque_cluster.txt")
  pvalue_Anova <- paste0(file_prefix, "_pvalue_Anova_ALL_plaque_cluster.txt")
  pvalue_Anova_0.05 <- paste0(file_prefix, "_pvalue_Anova_0.05_ALL_plaque_cluster.txt")
  Tukey_comparison <- paste0(file_prefix, "_Tukey_comparison_ALL_plaque_cluster.txt")
  Tukey_comparison_0.05 <- paste0(file_prefix, "_Tukey_comparison_0.05_ALL_plaque_cluster.txt")
  
  #all 
  sink(Anova_Tukey)
  
  for (network_name in names(input_data)) {
    network_df <- input_data[[network_name]]
    column <- network_name
    empty_list <- summary(aov(network_df[, 2] ~network_df[, 1]))
    tk <- TukeyHSD(aov(network_df[, 2] ~network_df[, 1]))
    print(column)
    print(empty_list)
    print(tk)
  }
  
  sink()
  
  #Tukey
  sink(Tukey)
  
  for (network_name in names(input_data)) {
    network_df <- input_data[[network_name]]
    column <- network_name
    empty_list <- summary(aov(network_df[, 2] ~ network_df[, 1]))
    tk <- TukeyHSD(aov(network_df[, 2] ~ network_df[, 1]))
    print(column)
    print(tk)
  }
  
  sink()
  
  sink(Anova_pvalue_file)
  
  for (network_name in names(input_data)) {
    network_df <- input_data[[network_name]]
    column <- network_name
    empty_list <- summary(aov(network_df[, 2] ~ network_df[, 1]))[[1]][["Pr(>F)"]]
    print(column)
    print(empty_list)
  }
  
  sink()
  
  #seperate the p values 
  Anova_pvalue = read.csv(Anova_pvalue_file, header= F)
  Network <- Anova_pvalue[seq(from = 1, to = nrow(Anova_pvalue), by = 2), 1]
  pvalue <- Anova_pvalue[seq(from =2, to =nrow(Anova_pvalue), by = 2), 1]
  Anova.pvalue<- data.frame(Network, pvalue)
  Anova.pvalue$Network <- gsub("^\\[1\\]\\s*|\\s*\\[1\\]$", "", Anova.pvalue$Network) #removes [1] before ensembl ID
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
  Tukey_all <- Tukey_all[Tukey_all$diff != "Df", ]
  Tukey_all <- Tukey_all[Tukey_all$diff != "codes:", ]
  Tukey_all <- Tukey_all[Tukey_all$comparison != "Residuals", ]
  
  amount_feature <- length(unique(unique_features))
  comparisons = (amount_feature *(amount_feature -1))/2
  comparisons = as.character(comparisons)
  Tukey_all$network_number <- rep(network_names, each = comparisons)
  
  #write to table
  write.table(Tukey_all, file = Tukey_comparison, sep = "\t")
  Tukey_all_0.05 <- Tukey_all[as.numeric(Tukey_all$p_adj) < 0.05, ]
  write.table(Tukey_all_0.05, file = Tukey_comparison_0.05, sep = "\t")
}

runAnovaTukey(average_correlation_total, "Connectivity")
runAnovaTukey(module_score_list, "module_score")


