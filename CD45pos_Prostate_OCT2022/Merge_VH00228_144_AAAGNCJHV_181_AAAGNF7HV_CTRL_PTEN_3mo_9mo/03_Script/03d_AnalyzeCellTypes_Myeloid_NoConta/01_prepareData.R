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

# If some clusters has been detected as contaminant, remove the corresponding cells
if( ! is.null( CONTAMINATION_CLUSTERS)){
  # Determine the cells that are contaminants
  Idents( merge_seurat_object) = "seurat_clusters"
  contamination_cells = Cells( merge_seurat_object)[ which( Idents( merge_seurat_object) %in% CONTAMINATION_CLUSTERS)]
  write.table( contamination_cells, file = file.path( PATH_ANALYSIS_EXTRA_OUTPUT, "Contamination_cells.csv"))
  
  # Remove the contamination cells
  merge_seurat_object = subset( merge_seurat_object, cells = setdiff( Cells( merge_seurat_object), contamination_cells))
  
  # Relevel the cluster meta-data factor
  merge_seurat_object = AddMetaData( merge_seurat_object, metadata = droplevels( Idents( merge_seurat_object)), col.name = "seurat_clusters")
}
