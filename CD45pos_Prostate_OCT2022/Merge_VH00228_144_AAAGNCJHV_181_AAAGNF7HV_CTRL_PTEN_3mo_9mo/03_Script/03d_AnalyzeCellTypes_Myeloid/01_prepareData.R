# #########################################
# This script reads and filters sc10x  data
# #########################################




# READ DATA
# ---------

## @knitr loadData

# Create the seurat object (RNA only)

cat("#### ANALYZING CELL TYPE :", SELECTED_CELL_TYPE)
cat("\n")
cat("\nLoad data from:\n", PATH_SEURAT_OBJECT_RDS_FILE)

# Load Seurat object from file
merge_seurat_object = readRDS( PATH_SEURAT_OBJECT_RDS_FILE)

# Subset the selected cells by clusters
Idents( merge_seurat_object) = "seurat_clusters"
merge_seurat_object = subset( merge_seurat_object, idents = SELECTED_CLUSTERS)

# Relevel the cluster meta-data factor
merge_seurat_object = AddMetaData( merge_seurat_object, metadata = droplevels( Idents( merge_seurat_object)), col.name = "seurat_clusters")

# Copy the clusters to an other meta-data to keep trace of the original values
merge_seurat_object = AddMetaData( merge_seurat_object, metadata = droplevels( Idents( merge_seurat_object)), col.name = "seurat_clusters_allcells")

