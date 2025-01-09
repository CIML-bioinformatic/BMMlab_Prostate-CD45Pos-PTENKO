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

# Plot the UMAP split by condition and colored by clusters
DimPlot(merge_seurat_object, group.by = "seurat_clusters", split.by = "orig.condition")
