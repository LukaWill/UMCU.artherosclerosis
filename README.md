# UMCU.artherosclerosis

To run the base analysis there are three scripts needed. 
1. config.R
Here all the input and output names are present and the feature will be selected.
feature is either cluster or one of the variables of the clinical_info (such as sex). 
With question = "plaque only" only the plaque differences are compared.
With question = "plaque combined" both the plaque differences are compared in combination with sex differences
The network names are selected based on a number. for STARNET this is number 1 to 224.

2. ANOVA_automatic.R
This is the script that performs the statistical analysis, which can be done for both plaque only or plaque combined. and this is used to select one of the two functions.

3. plaque_type_analyis.R or plaque_feature_combined.R
Here the data is normalized, the networks are selected and the module score and connectivity analysis are performed. The networks need to be saved as Genes_NUMBER.txt

