###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "04a_CompareConditions_Bcell_3mo"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Comparison of WT and PTEN KO conditions at 3mo : B cells"

# =========================
# DATA TO LOAD VARIABLES
# =========================

# Retrieve result of CellRanger from the analysis of individual samples
PATH_SEURAT_OBJECT_RDS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "03a_AnalyzeCellTypes_Bcell_NoConta_NoDoublet_3mo", "Resolution_0.2", paste0( outputFilesPrefix, "seuratObject_final.RDS"));

# =========================
# ANALYSIS VARIABLES
# =========================

# Define the selected cell type
SELECTED_CELL_TYPE = "B cell"

# The seed of random generator
SEED = 42

# Nb of cores when use of multiple cores using 'future' library, as implemented in Seurat3 when possible
NBCORES = 4

# ========== LOAD MODULES OF GENES ======================================================================

## Genes monitored as modules (tsv file, one column for each group of genes)
MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, "05_Module", "Modules.csv"),
                                     sep = ",",
                                     header = TRUE,
                                     stringsAsFactors = FALSE,
                                     row.names = NULL, fill = TRUE));
MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

# Path to the Cell cycle gene lists
CELL_CYCLE_SPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "01_CellCycle", "S_phase_genes.csv")))
CELL_CYCLE_G2MPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "01_CellCycle", "G2M_phase_genes.csv")))

MODULES_GENES[[ "S_PHASE"]] = CELL_CYCLE_SPHASE_GENELIST
MODULES_GENES[[ "G2M_PHASE"]] = CELL_CYCLE_G2MPHASE_GENELIST

# ========== DEFINE VARIABLES FOR DEG ANALYSIS==========================================================

# Alpha error for the FindMakers analysis
MARKER_GENES_ALPHA_THRESHOLD = 0.05
DEG_TEST_USE = "wilcox"

# Parameter for enrichment analysis in GO Terms
ENRICHMENT_GO_PVALUECUTOFF = 0.05
ENRICHMENT_GO_QVALUECUTOFF = 0.05

## KEGG Pathways
KEGG_PATHWAY_DF = read.table( file.path( PATH_EXPERIMENT_REFERENCE,
                                         "04_KEGGPathway",
                                         "KeggPathway.csv"),
                              sep = ",",
                              header = TRUE,
                              stringsAsFactors = FALSE,
                              row.names = NULL, fill = TRUE)
