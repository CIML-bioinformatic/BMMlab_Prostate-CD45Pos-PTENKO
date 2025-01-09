###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "07d_RNAVelocityAnalysis_Myeloid_PTEN_9mo"
ANALYSIS_STEP_LITTERAL_DESCRIPTION = "Analysis of cell dynamics using Velocyto - Myeloid cells"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Location of the velocyto loom files
VELOCYTO_LOOM_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "07_RNAVelocityAnalysis", "VelocytoLoom", paste0( EXPERIMENT_NAME, "_sorted.loom"))

# Retrieve result of CellRanger from the analysis of individual samples
SEURAT_OBJECT_HETEROGENEITY_ANALYSIS = file.path( PATH_EXPERIMENT_OUTPUT, "03d_AnalyzeCellTypes_Myeloid_NoConta_NoDoublet", "Resolution_0.8", paste0( outputFilesPrefix, "seuratObject_final.RDS"));

# Define the condition and date to use for analysis
SELECTED_CONDITION = CONDITION_NAME_PTEN
SELECTED_DATE = DATE_NAME_9MO

# Groups of clusters to analyze separatly
CLUSTER_GROUP_LIST = list( ClusterGroup1 = c( "5","6","1","3","12","11"),
                           ClusterGroup2 = c( "9","8","0","4"),
                           ClusterGroup3 = c( "5","1","3","12"),
                           ClusterGroup4 = c( "0","4","8","9","3")
                         )

## Genes monitored as modules (tsv file, one column for each group of genes)
MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, "05_Module", "Modules.csv"),
                                     sep = ",",
                                     header = TRUE,
                                     stringsAsFactors = FALSE,
                                     row.names = NULL, fill = TRUE));
MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

#MONITORED_GENES = sort( unique( unlist( MODULES_GENES)))
