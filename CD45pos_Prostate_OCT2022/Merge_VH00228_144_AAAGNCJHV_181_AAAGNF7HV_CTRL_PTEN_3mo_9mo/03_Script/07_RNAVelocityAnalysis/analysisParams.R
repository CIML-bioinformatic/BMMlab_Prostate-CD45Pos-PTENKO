###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME = "07_RNAVelocityAnalysis"
ANALYSIS_STEP_LITTERAL_DESCRIPTION = "Analysis of cell dynamics using Velocyto"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Location of the velocyto loom files
VELOCYTO_LOOM_FILE = file.path( PATH_ANALYSIS_OUTPUT, "VelocytoLoom", paste0( EXPERIMENT_NAME, "_sorted.loom"))

# Seurat object from Heterogeniety analysis
SEURAT_OBJECT_HETEROGENEITY_ANALYSIS = file.path( PATH_EXPERIMENT_OUTPUT, "01_DatasetIntegration", "resolution_0.4", paste0( outputFilesPrefix, "seuratObject_final.RDS"))

## Genes monitored as modules (tsv file, one column for each group of genes)
MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, "05_Module", "Modules.csv"),
                                     sep = ",",
                                     header = TRUE,
                                     stringsAsFactors = FALSE,
                                     row.names = NULL, fill = TRUE));
MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

#MONITORED_GENES = sort( unique( unlist( MODULES_GENES)))
