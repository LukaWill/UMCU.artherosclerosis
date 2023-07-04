# config.R

#file that gives all input and output names, and variable changes if needed

#input counts_file
counts_file <- "AE_bulk_RNA_batch1.minRib.PC_07042023.txt"

#output normalized counts
normalized_counts <- "normalized_counts.txt" #normalized counts
normalized_counts_log <- "normalized_counts_log.txt"  #normalized counts log
normalized_quantile <- "quantile_normalization.txt"
counts.plot = "~/Graduation/Rstudio/test/plot_normalization.pdf"#normalization plot for first 50 genes

#Input plaque types / clinical info
plaque_types <- "Seurat_clusters.txt"
clinical_info_file <- "clinical_data_good_selection.txt"

#directory path where the input networks are (based on ensembl IDs) 
directory_path <- "~/Graduation/Rstudio/All_genes"

#INPUT Features: 
cluster <- read.table(plaque_types, header = T)
clinical_info <- read.table(clinical_info_file, header = T, sep = "\t") 
names(clinical_info)[1] <- "study_number" #renames column 1 of the clinical_info


#if the question is based on the plaques only than fill in "plaques only" fill in "plaques combined" when both plaque and feature are used.
question = "plaque combined" #"plaque combined" "plaque only"  I used F if neither of them is used,meaning another feature of the clinical_info only

#Select column name of feature and input file, either clinical_info or cluster == plaque types, when comparing the feature within the clusters,  select the feature of the clinical info
if (question == "plaque only") {
  features <- cluster 
} else { #features is based on a variable within the clinical_info, 
  features <- clinical_info %>% select('study_number','secondary_events') #CHANGE the second variable 
}

  
#Output: name of the module score file, average expression per network per feature  
module_score_file <- "module_score_networks_ALL_plaque_secondary_events.csv"

#The network names that need to be selected, includes the word network. ("network1")
#All STARNET networks (1:225 except 12)
network_names_N <- c("network1", "network2", "network3", "network4", "network5", "network6", "network7", "network8", "network9", "network10", "network11", "network13", "network14", "network15", "network16", "network17", "network18", "network19", "network20", "network21", "network22", "network23", "network24", "network25", "network26", "network27", "network28", "network29", "network30", "network31", "network32", "network33", "network34", "network35", "network36", "network37", "network38", "network39", "network40", "network41", "network42", "network43", "network44", "network45", "network46", "network47", "network48", "network49", "network50", "network51", "network52", "network53", "network54", "network55", "network56", "network57", "network58", "network59", "network60", "network61", "network62", "network63", "network64", "network65", "network66", "network67", "network68", "network69", "network70", "network71", "network72", "network73", "network74", "network75", "network76", "network77", "network78", "network79", "network80", "network81", "network82", "network83", "network84", "network85", "network86", "network87", "network88", "network89", "network90", "network91", "network92", "network93", "network94", "network95", "network96", "network97", "network98", "network99", "network100", "network101", "network102", "network103", "network104", "network105", "network106", "network107", "network108", "network109", "network110", "network111", "network112", "network113", "network114", "network115", "network116", "network117", "network118", "network119", "network120", "network121", "network122", "network123", "network124", "network125", "network126", "network127", "network128", "network129", "network130", "network131", "network132", "network133", "network134", "network135", "network136", "network137", "network138", "network139", "network140", "network141", "network142", "network143", "network144", "network145", "network146", "network147", "network148", "network149", "network150", "network151", "network152", "network153", "network154", "network155", "network156", "network157", "network158", "network159", "network160", "network161", "network162", "network163", "network164", "network165", "network166", "network167", "network168", "network169", "network170", "network171", "network172", "network173", "network174", "network175", "network176", "network177", "network178", "network179", "network180", "network181", "network182", "network183", "network184", "network185", "network186", "network187", "network188", "network189", "network190", "network191", "network192", "network193", "network194", "network195", "network196", "network197", "network198", "network199", "network200", "network201", "network202", "network203", "network204", "network205", "network206", "network207", "network208", "network209", "network210", "network211", "network212", "network213", "network214", "network215", "network216", "network217", "network218", "network219", "network220", "network221", "network222", "network223", "network224", "network225")

#Output table, percentage overlap and no overlap of the input genes(countsfile) compared to the genes in the network. 
result_no_overlap <- "network_genes_no_overlap_percentage_ALL_plaque_secondary_events.txt"

#violin plot module score
violin_plot_module <- "~/Graduation/Rstudio/test/violin_plot_module_score_ALL_plaque_secondary_events.pdf"

#Connectivity

heatmap_cor = "~/Graduation/Rstudio/Test/all_correlation_ALL_plaque_secondary_events.pdf" #correlation heatmap for the total correlation of the network
dendogram_cluster_genes <- "cluster_genes_correlation_ALL_plaque_secondary_events.csv" #table of the dendogram clusters shown in the heatmap

directory <- "~/Graduation/Rstudio/Test/Networks" #directory to save the dendogram clusters separately

average_correlation_file <- "average_correlation_total_ALL_plaque_secondary_events.RData" #Rdata file that contains all average correlations
connectivity_plot <- "~/Graduation/Rstudio/Test/connectivity_grid_ALL_plaque_secondary_events.pdf" #connectivity plot that is based on the average correlation. 
