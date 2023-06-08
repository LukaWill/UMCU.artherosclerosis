# config.R

#input counts_file
counts_file <- "AE_bulk_RNA_batch1.minRib.PC_07042023.txt"

#normalized counts
normalized_counts <- "normalized_counts.txt"
normalized_counts_log <- "normalized_counts_log.txt"
normalized_quantile <- "quantile_normalization.txt"

#normalization plot
counts.plot = "~/Graduation/Rstudio/test/plot_normalization.pdf"

#the network names that need to be selected. ("1")
network_names <- c("39", "195", "177", "24", "35", "74", "109", "172", "147")

#name of the module score file that includes the clusters 
module_score_file <- "module_score_networks_AOR.csv"

#The network names that need to be selected, but it includes the word network. ("network1")
network_names_N <- c("network39", "network195", "network177", "network24", "network35", "network74", "network109", "network172", "network147")

#table that gives the percentage no overlap of the genes compared to the genes in the network. The no overlap is the amount of genes in the network that are not present in the data. 
result_no_overlap <- "network_genes_no_overlap_percentage_AOR.txt"

#violin plot module score
violin_plot_module <- "~/Graduation/Rstudio/violin_plot_module_score_AOR.pdf"

#ANOVA
#results Anova and Tukey together 
Anova_Tukey <- "ANOVA-Tukeys-networks_AOR.csv"
#Tukey only
Tukey <- "Tukeys-network_AOR.csv"
#Anova P-values in between file
Anova_pvalue_file <- "ANOVA-P-values-network_AOR.doc"

#pvalue to file
pvalue_Anova <- "p-value_table-network_AOR.csv"
pvalue_Anova_0.05 <- "p-value_0.05_table-network_AOR.csv"

Tukey_comparison <- "comparison_Tukey_AOR.csv"
Tukey_comparison_0.05 <- "comparison_Tukey_AOR_0.05.txt"

#Connectivity

heatmap_cor = "~/Graduation/Rstudio/Test/all_correlation_AOR.pdf"
dendogram_cluster_genes <- "cluster_genes_correlation_AOR.csv"

directory <- "~/Graduation/Rstudio/Test/clusters"

#connectivity plots
connectivity_plot <- "~/Graduation/Rstudio/Test/connectivity_grid_AOR.pdf"
connectivity_plot_all <- "~/Graduation/Rstudio/Test/connectivity_grid_AOR_all.pdf"