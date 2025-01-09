###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "02_CompareConditions_9mo"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Comparison of WT and PTEN KO conditions at 9 months"

# =========================
# DATA TO LOAD VARIABLES
# =========================

CLUSTER_RESOLUTION = 0.1
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, paste0( "Resolution_", CLUSTER_RESOLUTION))

# Retrieve result of CellRanger from the analysis of individual samples
PATH_SEURAT_OBJECT_RDS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "01_DatasetIntegration", paste0( "Resolution_", CLUSTER_RESOLUTION), paste0( outputFilesPrefix, "seuratObject_final.RDS"));

# =========================
# ANALYSIS VARIABLES
# =========================

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

## KEGG Pathways
KEGG_PATHWAY_DF = read.table( file.path( PATH_EXPERIMENT_REFERENCE,
                                         "04_KEGGPathway",
                                         "KeggPathway.csv"),
                              sep = ",",
                              header = TRUE,
                              stringsAsFactors = FALSE,
                              row.names = NULL, fill = TRUE)

# ========== DEFINE VARIABLES FOR DEG ANALYSIS==========================================================

# Alpha error for the FindMakers analysis
MARKER_GENES_ALPHA_THRESHOLD = 0.05
DEG_TEST_USE = "wilcox"

# Differential expression analysis
DEG_ALPHA_THRESHOLD = 0.05

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE      # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1       # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25      # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001     # PValue threshold for identification of significant markers
FINDMARKERS_SHOWTOP_TABLE     = 100   # Number of marker genes to show in report tables (NULL for all)
FINDMARKERS_SHOWTOP_HEATMAP   = 5     # Number of marker genes to show in report heatmaps (NULL for all)

# Parameter for enrichment analysis in GO Terms
ENRICHMENT_GO_PVALUECUTOFF = 0.05
ENRICHMENT_GO_QVALUECUTOFF = 0.05