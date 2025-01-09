###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "03c_AnalyzeCellTypes_Neutrophil_NoConta_NoDoublet_CTRL_3mo"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Neutrophil analysis for CTRL 3 months cells"

# =========================
# DATA TO LOAD VARIABLES
# =========================

# Retrieve result of CellRanger from the analysis of individual samples
PATH_SEURAT_OBJECT_INTEGRATION = file.path( PATH_PROJECT, 
                         "Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo",
                         "05_Output", 
                         "03c_AnalyzeCellTypes_Neutrophil_NoConta_NoDoublet",
                         "Resolution_0.8", 
                         "Merge_VH00228_144_AAAGNCJHV_181_AAAGNF7HV_CTRL_PTEN_3mo_9mo_seuratObject_final.RDS")

# =========================
# ANALYSIS VARIABLES
# =========================

# The seed of random generator
SEED = 42

# Nb of cores when use of multiple cores using 'future' library, as implemented in Seurat3 when possible
NBCORES = 4


#### Filtering / Normalization

# Switch from QC exploration mode and QC filtering
# In exploration mode, the cells that would be filtered by QC threshold are not filtered but marked
# They are shown in the later analysis with a different color
QC_EXPLORATION_MODE = FALSE

# Filters for loading seurat object
LOAD_MIN_CELLS     = 3;    # Retain cells with at least this many features (annotations)
LOAD_MIN_FEATURES  = 200;  # Retain annotations appearing in at least this many cells

# Cells with number of UMIs outside the range will be excluded
FILTER_UMI_MIN     = 0;
FILTER_UMI_MAX     = 70000;

# Cells with number of genes outside the range will be excluded
FILTER_FEATURE_MIN = 0;
FILTER_FEATURE_MAX = 8000;

# Cells with percentage of mitocohondrial genes above threshold will be excluded
FILTER_MITOPCT_MAX = 10;

# Cells with percentage of ribosomal genes below threshold will be excluded
FILTER_RIBOPCT_MIN = 0;

# Normalization parameters (see Seurat::NormalizeData())
DATA_NORM_METHOD      = "LogNormalize";
DATA_NORM_SCALEFACTOR = 10000;

# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE;
DATA_SCALE        = FALSE;
DATA_VARS_REGRESS = NULL;  # c("nCount_RNA") for UMIs (NULL to ignore)


#### Variable features / Dimensionality reduction

# Maximum number of variable features to keep
VARIABLE_FEATURES_MAXNB   = 2000
VARIABLE_FEATURES_SHOWTOP = 200   # For table in report
VARIABLE_FEATURE_METHOD = "vst"

# PCA parameters
PCA_NPC              = 50;  # Default number of dimensions to use for PCA (see Seurat::RunPCA())
PCA_PLOTS_NBDIMS     = 3;   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15;  # Number of'top' features to show when plotting PCA loadings

# Dimensionality reduction parameters (UMAP)
OPTIMIZE_UMAP_PARAMETERS = FALSE
BEST_UMAP_N_NEIGHBORS = 10
BEST_UMAP_MIN_DIST = 0.15
DIMREDUC_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results# PCA parameters

# Nearest-neighbor graph construction
FINDNEIGHBORS_K = 30

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION       = 0.8;
FINDCLUSTERS_RESOLUTION_RANGE = c( 0.1, 0.2, 0.4, 0.6, 0.8)
FINDCLUSTERS_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results
FINDCLUSTERS_ALGORITHM      = 1;   # 1 = Louvain; 2 = Louvain with multilevel refinement; 3 = SLM; 4 = Leiden

# Differential expression analysis
DEG_ALPHA_THRESHOLD = 0.05

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE      # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1       # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25      # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001     # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP_TABLE     = 100   # Number of marker genes to show in report tables (NULL for all)
FINDMARKERS_SHOWTOP_HEATMAP   = 5     # Number of marker genes to show in repot heatmaps (NULL for all)

# Module score analysis parameters
MODULES_CONTROL_SIZE = 100;

## Genes monitored as modules (tsv file, one column for each group of genes)
MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE, "05_Module", "Modules.csv"),
                                     sep = ",",
                                     header = TRUE,
                                     stringsAsFactors = FALSE,
                                     row.names = NULL, fill = TRUE))
MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

# Path to the Cell cycle gene lists
CELL_CYCLE_SPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "01_CellCycle", "S_phase_genes.csv")))
CELL_CYCLE_G2MPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "01_CellCycle", "G2M_phase_genes.csv")))

MODULES_GENES[[ "S_PHASE"]] = CELL_CYCLE_SPHASE_GENELIST
MODULES_GENES[[ "G2M_PHASE"]] = CELL_CYCLE_G2MPHASE_GENELIST
