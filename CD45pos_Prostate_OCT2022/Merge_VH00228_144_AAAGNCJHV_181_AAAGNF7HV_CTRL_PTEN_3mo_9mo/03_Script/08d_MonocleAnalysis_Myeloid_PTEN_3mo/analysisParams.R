###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "08d_MonocleAnalysis_Myeloid_PTEN_3mo"
ANALYSIS_STEP_LITTERAL_DESCRIPTION = "Analysis of cell dynamics using Monocle 3"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Seurat object from Heterogeniety analysis
SEURAT_OBJECT_HETEROGENEITY_ANALYSIS = file.path( PATH_EXPERIMENT_OUTPUT, "03d_AnalyzeCellTypes_Myeloid_NoConta_NoDoublet", "Resolution_0.8", paste0( outputFilesPrefix, "seuratObject_final.RDS"));

# Define the cluster resolution
CLUSTER_RESOLUTION = 0.003

# Define the condition and date to use for analysis
SELECTED_CONDITION = CONDITION_NAME_PTEN
SELECTED_DATE = DATE_NAME_3MO

# Groups of clusters to analyze separatly
CLUSTER_GROUP_LIST = list( ClusterGroup1 = c( "5","6","1","3","12","11")#,
                           #ClusterGroup2 = c( "9","8","0","4"),
                           #ClusterGroup3 = c( "5","1","3","12"),
                           #ClusterGroup4 = c( "0","4","8","9","3")
                         )

## Genes monitored as modules (tsv file, one column for each group of genes)
MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, "05_Module", "Modules.csv"),
                                     sep = ",",
                                     header = TRUE,
                                     stringsAsFactors = FALSE,
                                     row.names = NULL, fill = TRUE));
MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

## Genes groups to analyze (csv file, one column for each group of genes)
GROUP_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, 
                                              "05_Module",
                                              "Modules.csv"),
                                   sep = ",", quote = '"',
                                   header = TRUE,
                                   stringsAsFactors = FALSE,
                                   row.names = NULL, fill = TRUE));
GROUP_GENES = Map('[', GROUP_GENES, lapply(GROUP_GENES, function(x){ which( nchar( x)>0)})) # Remove empty strings

