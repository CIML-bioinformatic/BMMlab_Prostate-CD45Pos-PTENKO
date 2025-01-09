# #########################################
# This script reads and filters sc10x  data
# #########################################




# READ DATA
# ---------

# .....................................................
## @knitr loadData
# .....................................................

# Create the seurat object (RNA only)

cat("#### ANALYZING CELL TYPE :", SELECTED_CELL_TYPE)
cat("\n")
cat("\nLoad data from:\n", PATH_SEURAT_OBJECT_RDS_FILE)

# Load the data from previous analysis
merge_seurat_object = readRDS( file = PATH_SEURAT_OBJECT_RDS_FILE)

# Select the cells at the chosen timepoint
Idents( merge_seurat_object) = "orig.date"
merge_seurat_object = subset( merge_seurat_object, idents= DATE_NAME_3MO)

cat("<BR>Dispersion of cells among timepoints and conditions after selection")
table( merge_seurat_object$orig.condition, merge_seurat_object$orig.date) %>% kable() %>% kable_styling( full_width = FALSE)

# Plot the UMAP split by condition and colored by clusters
DimPlot( merge_seurat_object, 
         group.by = "seurat_clusters", 
         split.by = "orig.condition",
         cols = CLUSTER_COLOR_PANEL[ 1:length( unique( merge_seurat_object$seurat_clusters))],
         label = TRUE)
