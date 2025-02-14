###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#



ANALYSIS_STEP_NAME = "03_GlobalHeterogeneity_ByCellGroup"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "Quality Control, normalization and clustering"

# Retrieve result of CellRanger from previous analysis step
PATH_SEURAT_RDS_FILE = file.path( PATH_EXPERIMENT_OUTPUT,
                                            "02_GlobalHeterogeneity",
                                            "Resolution_0.8",
                                            paste0( "Prostate_Project_BMSH_DP_CD45pos_Prostate_OCT2022_230127_VH00228_144_AAAGNCJHV_", SAMPLE_NAME, "_seuratObject_final.RDS"))

# Path to the Cell cycle gene lists
CELL_CYCLE_SPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "01_CellCycle", "S_phase_genes.csv")))
CELL_CYCLE_G2MPHASE_GENELIST = unlist( use.names = FALSE, read.table( quote = NULL, header = TRUE, file = file.path( PATH_EXPERIMENT_REFERENCE, "01_CellCycle", "G2M_phase_genes.csv")))

# Path to Heat Shock stress genes list
PATH_HS_STRESS_MARKER_GENES_TABLE_FILE = file.path( PATH_EXPERIMENT_REFERENCE, "02_HeatShock", "coregene_df-FALSE-v3.csv")

#### General

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (for Seurat3 using 'future')
NBCORES = 4;

# Number of cells above which use ggplot instead of interactive plotly
PLOT_RASTER_NBCELLS_THRESHOLD = 20000;


#### Analysis parameters

# Maximum number of variable features to keep
VARIABLE_FEATURES_MAXNB   = 2000;  # For analysis (PCA)
VARIABLE_FEATURES_SHOWTOP = 200;   # For table in report

# Nearest-neighbor graph construction
FINDNEIGHBORS_K = 30

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION     = 0.1;
FINDCLUSTERS_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results
FINDCLUSTERS_ALGORITHM      = 1;   # 1 = Louvain; 2 = Louvain with multilevel refinement; 3 = SLM; 4 = Leiden

# RENAME_CLUSTERS = list()
# RENAME_CLUSTERS[[ SAMPLE_ALL]] = list()
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "0"]] = 0
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "1"]] = 5
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "2"]] = 1
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "3"]] = 3
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "4"]] = 6
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "5"]] = 2
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "6"]] = 4
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "7"]] = 7
  

# PCA parameters
PCA_NPC              = 50;  # Default number of dimensions to use for PCA (see Seurat::RunPCA())
PCA_PLOTS_NBDIMS     = 3;   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15;  # Number of'top' features to show when plotting PCA loadings

# Dimensionality reduction parameters (TSNE/UMAP)
DIMREDUC_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP_TABLE     = 100;       # Number of marker genes to show in report tables (NULL for all)
FINDMARKERS_SHOWTOP_HEATMAP   = 5;       # Number of marker genes to show in repot heatmaps (NULL for all)

# Parameter for enrichment analysis in GO Terms
ENRICHMENT_GO_PVALUECUTOFF = 0.05
ENRICHMENT_GO_QVALUECUTOFF = 0.05



#### List of group of clusters of interest

CLUSTER_GROUP_LIST = list( T_NK_cell = c( 0, 4, 5, 10, 12, 13, 15, 16),
                           Neutrophil_cell = c( 20),
                           B_cell = c( 3, 8),
                           Myeloid_cell = c( 1, 2, 6, 7, 9, 11, 14, 17, 18, 19))


#### Lists of genes of interest

# Number of genes under which genes of a module (MODULES_GENES) are transfered to be analyzed individually (MONITORED_GENES)
MONITORED_GENES_SMALL_MODULE = 5;
MODULES_CONTROL_SIZE = 100;


## Contamination-related genes (txt file, one gene name by line)
#CONTAMINATION_GENES = readLines( file.path( PATH_PROJECT_EXTERNALDATA, "ContaminationGenes.txt"))
CONTAMINATION_GENES=NULL

## Genes monitored individually (tsv file, one column for each group of genes)
MONITORED_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE,
                                                  "03_MonitoredGenes",
                                                  "Monitored.csv"),
                                       sep = ",",
                                       header = TRUE,
                                       stringsAsFactors = FALSE,
                                       row.names = NULL, fill = TRUE));
MONITORED_GENES = Map('[', MONITORED_GENES, lapply(MONITORED_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings

## Genes monitored individually (tsv file, one column for each group of genes)
CHOSEN_MARKER_GENES = read.table( file.path( PATH_EXPERIMENT_REFERENCE,
                                                  "03_MonitoredGenes",
                                                  "ChosenMarkerGenes.csv"),
                                       sep = ",",
                                       header = TRUE,
                                       stringsAsFactors = FALSE,
                                       row.names = NULL, fill = TRUE);
CHOSEN_MARKER_GENES = CHOSEN_MARKER_GENES$Gene

## Genes groups to analyze (csv file, one column for each group of genes)
GROUP_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE,
                                              "04_GeneGroups",
                                              "GeneGroups.csv"),
                                   sep = ",", quote = '"',
                                   header = TRUE,
                                   stringsAsFactors = FALSE,
                                   row.names = NULL, fill = TRUE));
GROUP_GENES = Map('[', GROUP_GENES, lapply(GROUP_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings


## Genes monitored as modules (tsv file, one column for each group of genes)
MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE,
                                                "05_Modules",
                                                "Modules.csv"),
                                        sep = ",",
                                        header = TRUE,
                                        stringsAsFactors = FALSE,
                                        row.names = NULL, fill = TRUE));
MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings
MODULES_GENES[[ "S_PHASE"]] = CELL_CYCLE_SPHASE_GENELIST
MODULES_GENES[[ "G2M_PHASE"]] = CELL_CYCLE_G2MPHASE_GENELIST


