# #########################################
# This script reads and filters sc10x  data
# #########################################




# READ DATA
# ---------

# .....................................................
## @knitr loadData
# .....................................................

# Create the seurat object (RNA only)
merge_seurat_object = readRDS( file = PATH_SEURAT_OBJECT_RDS_FILE)


# Plot the cell in UMAP colored by clusters and split by condition
DimPlot( merge_seurat_object, 
         group.by = "seurat_clusters", 
         split.by = "orig.condition",
         cols = CLUSTER_COLOR_PANEL[ 1:length( unique( merge_seurat_object$seurat_clusters))],
         label = TRUE)
