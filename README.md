# TraumaRDB
This folder contains the data from Trauma-Related Database (TRDB) and the code developed for its analysis.


All of the TRDB, from the “Inflammation and the Host Response to Injury” (IHRI) research program, are freely available in the folder “data”, as downloaded from www.gluegrant.org.


The code published were fully developed in R. The following “.R” files have the code described:
1. Load_data: 
Contains the code used to import the excel and csv files necessary for the analysis (microarrays and clinical data).
2. AnalysisT0:
The gene expression were analysed at the day of injury (day 0). Survival analysis and dimensionality reduction was accomplished. 
3. AnalysisT1:
The gene expression were analysed at the first day after injury (day 1). Survival analysis and dimensionality reduction was accomplished. 
4. AnalysisT4:
The gene expression were analysed at the fourth day after injury (day 4). Survival analysis and dimensionality reduction was accomplished. 
5. Genes_selected_trajectory
The time trajectory of the intersection of genes previously selected at the first 4 days after injury are computed. 
5. Imputation:
Automated code used for missing data imputation (junction of LOCF and NOCB with linear interpolation). 
5. MTSclustering:
Clustering of multivariate time series (MTS clustering) with R package "dtwclust” and survival analysis of the results. 


The .xlsx file “ProbeID_to_Accession” was used to substitute the probe ID into the accession number (genbank reference number). 


Please contact the corresponding author Cláudia S. Constantino for more details if needed (claudia.s.constantino@tecnico.ulisboa.pt). 


