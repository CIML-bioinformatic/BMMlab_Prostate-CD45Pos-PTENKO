###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "07c_RNAVelocityAnalysis_Neutrophil_PTEN_3mo_9mo"
ANALYSIS_STEP_LITTERAL_DESCRIPTION = "Analysis of cell dynamics using Velocyto - Neutrophil cells"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Location of the velocyto loom files
VELOCYTO_LOOM_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "07_RNAVelocityAnalysis", "VelocytoLoom", paste0( EXPERIMENT_NAME, "_sorted.loom"))

# Retrieve result of CellRanger from the analysis of individual samples
SEURAT_OBJECT_HETEROGENEITY_ANALYSIS = file.path( PATH_EXPERIMENT_OUTPUT, "03c_AnalyzeCellTypes_Neutrophil_NoConta_NoDoublet", "Resolution_0.8", paste0( outputFilesPrefix, "seuratObject_final.RDS"));

# Define the condition and date to use for analysis
SELECTED_CONDITION = CONDITION_NAME_PTEN

# Groups of clusters to analyze separately
CLUSTER_GROUP_LIST = list( ClusterGroup1 = c( "0","1","5","6"),
                           ClusterGroup2 = c( "3","2","5","6"),
                           ClusterGroup3 = c( "2","5","6"),
                           ClusterGroup4 = c( "3","2","4","7"),
                           ClusterGroup5 = c( "0","1","3"),
                           ClusterGroup6 = c( "0","1","2","3","5","6"),
                           ClusterGroup7 = c( "0","1","2","3","4","5","6","7")
)

## Genes monitored as modules (tsv file, one column for each group of genes)
MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, "05_Module", "Modules.csv"),
                                     sep = ",",
                                     header = TRUE,
                                     stringsAsFactors = FALSE,
                                     row.names = NULL, fill = TRUE));
MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

#MONITORED_GENES = sort( unique( unlist( MODULES_GENES)))
