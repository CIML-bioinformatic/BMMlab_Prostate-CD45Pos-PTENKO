###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "04a_CompareConditions_Bcell_CompositionAnalysis_3mo"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Comparison of WT and PTEN KO conditions at 3 months - Bcell - Analysis of differential composition"

# =========================
# DATA TO LOAD VARIABLES
# =========================

# Retrieve result of CellRanger from the analysis of individual samples
PATH_SEURAT_OBJECT_RDS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "04a_CompareConditions_Bcell", paste0( outputFilesPrefix, "seuratObject_final.RDS"));


# =========================
# ANALYSIS VARIABLES
# =========================

# The seed of random generator
SEED = 42

# Nb of cores when use of multiple cores using 'future' library, as implemented in Seurat3 when possible
NBCORES = 4

